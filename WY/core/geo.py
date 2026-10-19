"""GEO 数据获取模块:E-utilities 检索 + series matrix 下载与自解析。

第 0 步前置检查结论:GEOparse 2.0.4 的下载器跨境不稳定(FTP 拉全量 SOFT、
大小校验失败),filepath 解析 series matrix 取不到数据表 —— 故本模块不用
GEOparse,采用 requests 下载(重试 + gzip 完整性校验 + 永久缓存)+ 自行解析。

设计要点(实施方案 §4.1/§4.2):
- GEO 是归档库:series matrix 一经发布不可变,文件层永久缓存(存在即命中,
  仅 force_refresh 重下);检索结果走 HTTP 缓存会话(24h TTL);
- 失败语义:响亮报错并给出路 —— 下载重试耗尽后,本地有旧缓存则降级沿用并
  警告,否则报错退出;不产生空壳产物;
- 解析防御:每条 !Sample_ 行 = 字段名 + 全部样本取值(带引号 TSV,必须逐值
  拆开);多行 characteristics 每行一个注释维度;空矩阵数据集记入 issues。
"""

import csv
import gzip
import io
import logging
import os
import time
from dataclasses import dataclass, field

import pandas as pd
import requests

_NCBI_HEADERS = {"User-Agent": "NetworkPharmacologyPlatform/1.0"}
# 检索式只保留表达谱类条目(芯片/高通量测序),排除 eQTL、甲基化等无关类型
_EXPRESSION_FILTERS = (
    '("expression profiling by array"[Filter] OR '
    '"expression profiling by high throughput sequencing"[Filter])'
)


def _config():
    import config
    return config


def _geo_bucket(accession: str) -> str:
    """GEO FTP 千位分桶:去掉编号末三位再拼 nnn(不足四位不补零,如 GPL570 → GPLnnn)。"""
    digits = accession.strip().upper()[3:]
    prefix = digits[:-3] if len(digits) > 3 else ""
    return accession.strip().upper()[:3] + prefix + "nnn"


def series_matrix_url(gse: str) -> str:
    """按 GEO FTP 目录规律拼 series matrix 直链(如 GSE63067 → GSE63nnn 桶)。"""
    gse = gse.strip().upper()
    return (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{_geo_bucket(gse)}/{gse}/"
            f"matrix/{gse}_series_matrix.txt.gz")


@dataclass
class GeoDataset:
    """一个 GSE 数据集的解析产物。

    samples: 样本注释表(行=GSM 编号,列=各注释字段;多行 characteristics
    已按注释维度分列,并把「键: 值」对进一步解析为独立列);
    table: 表达矩阵(行=探针或基因符号,列=GSM 编号,数值型)。
    """
    gse: str
    meta: dict
    samples: pd.DataFrame
    table: pd.DataFrame
    gpl: str = ""
    issues: list = field(default_factory=list)


# ---------------- 检索(E-utilities,限流并入现有节流表) ----------------


def search_geo_datasets(disease: str, organism: str = "Homo sapiens",
                        limit: int = 20) -> pd.DataFrame:
    """按疾病名检索 GEO 表达谱数据集,返回候选表供人工选型(不替用户做主)。

    列:gse/title/n_samples/platform/possible_cellline(启发式标注)/link。
    「可能为细胞系」由标题/摘要关键词启发式判定:本平台湿实验为细胞模型,
    组织来源数据的表达证据迁移有折扣,检索阶段即需可见。
    """
    import config
    from .http_client import request_with_retry

    disease = disease.strip()
    term = (f'"{disease}" AND "{organism}"[Organism] AND gse[ETYP] '
            f"AND {_EXPRESSION_FILTERS}")
    empty = pd.DataFrame(columns=["gse", "title", "n_samples",
                                  "platform", "possible_cellline", "link"])
    uids: list[str] = []
    params_base = {"db": "gds", "retmode": "json"}
    if config.NCBI_API_KEY:
        params_base["api_key"] = config.NCBI_API_KEY

    params = dict(params_base, term=term, retmax=limit)
    resp = request_with_retry(
        f"{config.NCBI_EUTILS_URL}/esearch.fcgi", method="get", params=params,
        max_retries=3, timeout=60, label=f"NCBI esearch ({disease})",
        throttle_key="NCBI_EUTILS", headers=_NCBI_HEADERS,
        raise_on_exhausted=True, use_cache=True,
    )
    uids = list(resp.json().get("esearchresult", {}).get("idlist", []))
    if not uids:
        logging.info("GEO search: no datasets for '%s'", disease)
        return empty

    params = dict(params_base, id=",".join(uids))
    resp = request_with_retry(
        f"{config.NCBI_EUTILS_URL}/esummary.fcgi", method="get", params=params,
        max_retries=3, timeout=60, label=f"NCBI esummary ({disease})",
        throttle_key="NCBI_EUTILS", headers=_NCBI_HEADERS,
        raise_on_exhausted=True, use_cache=True,
    )
    result = resp.json().get("result", {})
    out = []
    for uid in result.get("uids", []):
        d = result.get(uid, {})
        text = " ".join(str(d.get(k, "")) for k in ("title", "summary"))
        cell_flag = any(k in text.lower() for k in (
            "cell line", "cell lines", "cells in vitro", "hepg2", "hela",
            "hek293", "a549", "sh-sy5y", "raw 264", "bv2"))
        gse = d.get("accession", "")
        out.append({
            "gse": gse,
            "title": d.get("title", ""),
            "n_samples": d.get("n_samples", ""),
            "platform": d.get("gpl", ""),
            "possible_cellline": "是" if cell_flag else "",
            "link": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
        })
    return pd.DataFrame(out) if out else empty


# ---------------- 下载(重试 + 校验 + 永久缓存) ----------------


def _geo_cache_path(gse: str) -> str:
    config = _config()
    d = os.path.join(config.CACHE_DIR, config.GEO_SERIES_CACHE_SUBDIR)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{gse.strip().upper()}_series_matrix.txt.gz")


def _gzip_intact(path: str) -> bool:
    """gzip 完整性全量校验(第 0 步 GEOparse 的翻车点:截断文件必须视为失败)。"""
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except OSError:
        return False


def _download_with_retry(url: str, dest: str) -> None:
    """下载到 .part 再原子替换,不留半截文件;Content-Length 与 gzip 双重校验。"""
    config = _config()
    last_exc: Exception | None = None
    for attempt in range(1, config.GEO_DOWNLOAD_RETRIES + 1):
        try:
            got = 0
            with requests.get(url, timeout=config.GEO_DOWNLOAD_TIMEOUT,
                              headers=_NCBI_HEADERS, stream=True) as r:
                r.raise_for_status()
                expected = int(r.headers.get("content-length", 0) or 0)
                part = dest + ".part"
                with open(part, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
                        got += len(chunk)
            if expected and got != expected:
                raise OSError(f"size mismatch: got {got}, expected {expected}")
            if not _gzip_intact(part):
                raise OSError("gzip integrity check failed (truncated or corrupt)")
            os.replace(part, dest)  # 原子替换:要么完整新文件,要么保留旧缓存
            logging.info("GEO download OK: %s (%d bytes, attempt %d)",
                         url, got, attempt)
            return
        except Exception as exc:
            last_exc = exc
            logging.exception("GEO download failed (attempt %d/%d): %s",
                              attempt, config.GEO_DOWNLOAD_RETRIES, url)
            if attempt < config.GEO_DOWNLOAD_RETRIES:
                time.sleep(config.GEO_DOWNLOAD_RETRY_DELAY)
    if os.path.isfile(dest) and _gzip_intact(dest):
        logging.warning("下载耗尽,降级使用本地旧缓存: %s", dest)
        return
    raise RuntimeError(
        f"GEO 文件下载失败({config.GEO_DOWNLOAD_RETRIES} 次重试后): {url}。"
        f"可选出路:① 稍后重试(跨境带宽波动大);② 用浏览器手动下载该文件后"
        f"放到 {dest};③ 换一个规模更小的数据集。最后错误:{last_exc}")


def fetch_geo_series_matrix(gse: str, force_refresh: bool = False) -> str:
    """下载并永久缓存 series matrix,返回本地文件路径(归档不可变,存在即命中)。"""
    gse = gse.strip().upper()
    path = _geo_cache_path(gse)
    if os.path.isfile(path) and not force_refresh and _gzip_intact(path):
        logging.info("GEO series matrix 缓存命中: %s", path)
        return path
    _download_with_retry(series_matrix_url(gse), path)
    return path


# ---------------- 自解析(第 0 步验证的备选路线,现为正式路线) ----------------


def _split_tsv(line: str) -> list[str]:
    """按带引号的 TSV 规范拆行(series matrix 全文件使用该格式)。"""
    return next(csv.reader(io.StringIO(line.rstrip("\r\n")), delimiter="\t",
                           quotechar='"'))


def parse_series_matrix(path: str) -> GeoDataset:
    """解析 series matrix:元数据键值 + 样本注释(逐样本取值)+ 数据表。

    空矩阵(元数据型提交/超级系列)不清真地返回空表,而是记入 issues 并由
    调用方决定报错 —— 防住第 0 步实际踩到的 GSE84465 型数据集。
    """
    gse = os.path.basename(path).split("_")[0]
    meta: dict = {}
    sample_lines: dict[str, list] = {}   # 字段名 -> [每行的逐样本取值列表]
    table_rows: list[list[str]] = []
    table_header: list[str] | None = None
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                table_header = None
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                parts = _split_tsv(line)
                if table_header is None:
                    table_header = parts
                else:
                    table_rows.append(parts)
            elif line.startswith("!Sample_"):
                parts = _split_tsv(line)
                sample_lines.setdefault(parts[0], []).append(parts[1:])
            elif line.startswith("!Series_") or line.startswith("!Platform_"):
                key, _, val = line.partition("\t")
                meta.setdefault(key.lstrip("!"), val.strip('"'))

    # 样本注释表:行 = GSM(!Sample_geo_accession 首行的逐样本取值),列 = 各字段
    gsm_key = "!Sample_geo_accession"
    if gsm_key not in sample_lines or not sample_lines[gsm_key][0]:
        raise RuntimeError(f"{gse}: series matrix 缺少样本编号行,无法解析")
    gsms = [str(x) for x in sample_lines[gsm_key][0]]
    samples = pd.DataFrame(index=gsms)
    issues: list = []
    for key, lines in sample_lines.items():
        if key == gsm_key:
            continue
        base = key.removeprefix("!Sample_")
        for i, values in enumerate(lines):
            if len(values) != len(gsms):
                issues.append(f"字段 {base} 第 {i + 1} 行样本数不匹配"
                              f"({len(values)}!={len(gsms)}),已跳过")
                continue
            existing = [c for c in samples.columns
                        if c == base or c.startswith(base + "_")]
            name = base if not existing else f"{base}_{len(existing) + 1}"
            samples[name] = values

    # characteristics 的「键: 值」进一步解析为独立列(分组字段的实际来源)
    char_cols = [c for c in samples.columns
                 if c == "characteristics_ch1" or c.startswith("characteristics_ch1_")]
    for c in char_cols:
        parsed: dict[str, dict] = {}
        for gsm, cell in samples[c].items():
            k, sep, v = str(cell).partition(":")
            if sep:
                parsed.setdefault(k.strip(), {})[gsm] = v.strip()
        for k, mapping in parsed.items():
            if k and len(mapping) >= max(2, int(len(samples) * 0.9)):
                if k and k not in samples.columns:
                    samples[k] = [mapping.get(g, "") for g in samples.index]

    # 数据表:ID_REF 列为行标识,GSM 列为样本;数值化
    if not table_header or not table_rows:
        issues.append("数据表为空(可能为元数据型提交或超级系列,无表达矩阵)")
        table = pd.DataFrame()
    else:
        ids = [r[0].strip('"') for r in table_rows]
        data = {table_header[j]: [r[j] if j < len(r) else "" for r in table_rows]
                for j in range(1, len(table_header))}
        table = pd.DataFrame(data, index=ids)
        table = table.apply(pd.to_numeric, errors="coerce")
        table = table.dropna(how="all")
    gpl = ""
    if "!Sample_platform_id" in sample_lines and sample_lines["!Sample_platform_id"]:
        gpl = str(sample_lines["!Sample_platform_id"][0][0])
    elif meta.get("Series_platform_id"):
        gpl = meta["Series_platform_id"]
    return GeoDataset(gse=gse, meta=meta, samples=samples, table=table,
                      gpl=gpl, issues=issues)


def fetch_geo_dataset(gse: str, force_refresh: bool = False) -> GeoDataset:
    """下载(或命中缓存)+ 解析,一步到位。"""
    path = fetch_geo_series_matrix(gse, force_refresh=force_refresh)
    return parse_series_matrix(path)


# ---------------- 分组辅助(「全模块最脏环节」的人工决策支持) ----------------


def group_fields_overview(ds: GeoDataset) -> pd.DataFrame:
    """列出样本注释中所有字段的取值分布,供人工确定 case/control 分组正则。

    每行:字段名 | 不同取值数 | 取值样例(频次)。--case-pattern 匹配的对象
    是「样本注释全文」(各注释字段按 | 拼接)。
    """
    rows = []
    for col in ds.samples.columns:
        vals = ds.samples[col].astype(str)
        dist = vals.value_counts()
        sample = "; ".join(f"{v}(x{n})" for v, n in dist.head(5).items())
        rows.append({"field": col, "n_unique": int(dist.size), "values": sample})
    return pd.DataFrame(rows)


def assign_groups(ds: GeoDataset, case_pattern: str, control_pattern: str) -> pd.Series:
    """按正则把样本分为 case(1)/control(0);匹配对象为样本注释全文。

    四道硬边界(实施方案 §2 + 评审补充,缺一不可):
    - case/control 各自匹配 0 个样本 → 报错(正则写错);
    - 任一 pattern 匹配全部样本 → 报错(没有对照/病例组);
    - 两 pattern 命中集合重叠 → 报错(标签方向歧义);
    - 任一组样本数 < ML_MIN_SAMPLES_PER_GROUP → 报错(小样本无意义)。
    未被任何 pattern 命中的样本丢弃(如中间分级),并记录数量。
    """
    config = _config()
    text = ds.samples.apply(lambda row: " | ".join(str(v) for v in row.values), axis=1)
    case_hit = text.str.contains(case_pattern, case=False, regex=True)
    ctrl_hit = text.str.contains(control_pattern, case=False, regex=True)
    n_case, n_ctrl = int(case_hit.sum()), int(ctrl_hit.sum())
    if n_case == 0 or n_ctrl == 0:
        raise RuntimeError(
            f"分组正则匹配失败:case 匹配 {n_case} 个、control 匹配 {n_ctrl} 个样本。"
            "请先运行不带分组参数的 deg 命令查看字段取值后修正正则。")
    if n_case == len(text) or n_ctrl == len(text):
        raise RuntimeError(
            "分组正则匹配了全部样本(case 或 control 必须留出对照组/病例组),"
            "请把正则收紧到具体分组字段值。")
    overlap = text[case_hit & ctrl_hit]
    if len(overlap):
        raise RuntimeError(
            f"case 与 control 正则命中重叠({len(overlap)} 个样本,如 "
            f"{list(overlap.index[:3])}):同一样本不能既是病又是对照,标签方向歧义。")
    for name, n in (("case", n_case), ("control", n_ctrl)):
        if n < config.ML_MIN_SAMPLES_PER_GROUP:
            raise RuntimeError(
                f"{name} 组仅 {n} 个样本(低于 {config.ML_MIN_SAMPLES_PER_GROUP}):"
                "小样本下统计与 ML 均无意义,请换样本量更大的数据集或放宽分组定义。")
    dropped = int((~case_hit & ~ctrl_hit).sum())
    if dropped:
        logging.info("分组:丢弃未命中任何 pattern 的样本 %d 个", dropped)
    # 只保留命中任一 pattern 的样本 —— 未命中的"丢弃"必须真正从标签中剔除,
    # 否则会被默认值 0 误标为对照,污染对照组(2026-10-19 M1 审查发现并修复)
    labels = pd.Series(0, index=text.index[case_hit | ctrl_hit], dtype=int)
    labels[text.index[case_hit]] = 1
    return labels
