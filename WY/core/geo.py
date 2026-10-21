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
import glob
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
# 样本来源启发式关键词(方案 §4.1「样本来源」标注):本平台湿实验全为细胞模型,
# 组织来源数据的表达证据迁移有折扣,检索表与 shortlist(evidence_scope)共用同一判定
_CELL_LINE_KEYWORDS = (
    "cell line", "cell lines", "cells in vitro", "hepg2", "hela",
    "hek293", "a549", "sh-sy5y", "raw 264", "bv2")
_TISSUE_KEYWORDS = (
    "tissue", "biopsy", "biopsies", "plasma", "serum", "whole blood",
    "pbmc", "peripheral blood", "patient", "patients")


def sample_source_label(*texts: str) -> str:
    """标题/摘要 → 样本来源三值标注:细胞系/组织/未标注(关键词启发式,非保证)。

    检索表(possible_cellline 列)与 screen 链(shortlist 的 evidence_scope 列,
    方案 §4.3)共用此判定,保证两处口径一致;同时命中时细胞系优先。
    """
    text = " ".join(str(t) for t in texts if t).lower()
    if any(k in text for k in _CELL_LINE_KEYWORDS):
        return "细胞系"
    if any(k in text for k in _TISSUE_KEYWORDS):
        return "组织"
    return "未标注"


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
    # 本次文件为「下载耗尽后降级沿用的旧缓存」时置 True(方案失效场景②):
    # 降级不能只进日志,CLI/页面据此向用户打旧数据警告
    cache_fallback: bool = False
    cache_date: str = ""   # series matrix 文件的最后修改时间(降级警告中展示)


# ---------------- 检索(E-utilities,限流并入现有节流表) ----------------


def search_geo_datasets(disease: str, organism: str = "Homo sapiens",
                        limit: int = 20) -> pd.DataFrame:
    """按疾病名检索 GEO 表达谱数据集,返回候选表供人工选型(不替用户做主)。

    列:gse/title/n_samples/platform/type(芯片/测序)/symbol_level(行标识粒度
    提示,启发式)/possible_cellline(启发式标注)/link。
    - type/symbol_level 来自 esummary 的 gdstype:探针级数据集首跑要多下 GPL
      注释表(+几分钟),检索阶段标出让用户建议优先选「多为基因级」(方案 §2);
    - possible_cellline 由标题/摘要关键词启发式判定(sample_source_label,
      与 shortlist 的 evidence_scope 同一判定):本平台湿实验为细胞模型,
      组织来源数据的表达证据迁移有折扣,检索阶段即需可见。
    """
    import config
    from .http_client import request_with_retry

    disease = disease.strip()
    term = (f'"{disease}" AND "{organism}"[Organism] AND gse[ETYP] '
            f"AND {_EXPRESSION_FILTERS}")
    empty = pd.DataFrame(columns=["gse", "title", "n_samples", "platform",
                                  "type", "symbol_level", "possible_cellline",
                                  "link"])
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
        cache_ttl_hours=config.GEO_SEARCH_CACHE_TTL_HOURS,
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
        cache_ttl_hours=config.GEO_SEARCH_CACHE_TTL_HOURS,
    )
    result = resp.json().get("result", {})
    out = []
    for uid in result.get("uids", []):
        d = result.get(uid, {})
        gse = d.get("accession", "")
        gdstype = str(d.get("gdstype", "") or "")
        low = gdstype.lower()
        if "sequencing" in low:
            type_short, sym_hint = "测序", "多为基因级"
        elif "array" in low:
            type_short, sym_hint = "芯片", "多为探针级"
        else:
            type_short, sym_hint = (gdstype[:24] or "未知"), "未知"
        out.append({
            "gse": gse,
            "title": d.get("title", ""),
            "n_samples": d.get("n_samples", ""),
            "platform": d.get("gpl", ""),
            "type": type_short,
            "symbol_level": sym_hint,
            "possible_cellline": "是" if sample_source_label(
                d.get("title", ""), d.get("summary", "")) == "细胞系" else "",
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


def _read_lock_pid(lock: str) -> int | None:
    """读锁文件里记录的持有者 PID;空锁(旧版/写 PID 前被杀)返回 None。"""
    try:
        with open(lock, "r", encoding="ascii", errors="ignore") as f:
            content = f.read().strip()
        return int(content) if content else None
    except (OSError, ValueError):
        return None


def _pid_alive(pid: int) -> bool:
    """进程是否存活(Windows 走 OpenProcess 探测,Linux 走 kill 0;判不了按存活算)。"""
    if pid <= 0:
        return False
    if os.name == "nt":
        import ctypes
        PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
        handle = ctypes.windll.kernel32.OpenProcess(
            PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
        if handle:
            ctypes.windll.kernel32.CloseHandle(handle)
            return True
        return False
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except OSError:
        return True  # PermissionError 等无权限情形:进程存在


def acquire_download_lock(dest: str) -> bool:
    """同文件并发下载互斥(O_EXCL 锁文件,内含持有者 PID),返回是否拿到下载权。

    页面运行中重复点击、或 CLI 与页面同时跑同一数据集时,只允许一个写入者
    (两个线程交错写同一 .part 会产出坏文件)。接管判定按序:
    ① 成品文件已完整 → 直接用缓存(返回 False),顺手清掉遗留锁;
    ② 持有者进程已死(运行被强杀/崩溃的典型遗留)→ 立即接管;
    ③ PID 不可读(旧版空锁)或进程仍在但锁龄超过全部重试窗口 → 兜底接管
      (兜底窗口必须覆盖最慢的正常下载,故为 超时×重试+余量)。
    """
    config = _config()
    lock = dest + ".lock"
    stale_after = config.GEO_DOWNLOAD_TIMEOUT * config.GEO_DOWNLOAD_RETRIES + 60
    waits = 0
    while True:
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, str(os.getpid()).encode("ascii"))
            os.close(fd)
            return True
        except FileExistsError:
            if os.path.isfile(dest) and _gzip_intact(dest):
                release_download_lock(dest)
                return False
            holder_pid = _read_lock_pid(lock)
            reason = None
            if holder_pid is not None and not _pid_alive(holder_pid):
                reason = f"持有者 PID {holder_pid} 已不存在(上次运行被强杀/崩溃)"
            elif time.time() - os.path.getmtime(lock) > stale_after:
                reason = f"锁龄超过兜底窗口 {int(stale_after)}s"
            if reason:
                logging.warning("下载锁接管(%s): %s", reason, lock)
                try:
                    os.remove(lock)
                except OSError:
                    pass
                continue
            waits += 1
            if waits % 15 == 1:  # 2s 轮询,30s 记一条,避免日志刷屏
                logging.info("同文件下载进行中(持有者 PID %s),等待锁释放: %s",
                             holder_pid if holder_pid is not None else "未知", lock)
            time.sleep(2)


def release_download_lock(dest: str) -> None:
    try:
        os.remove(dest + ".lock")
    except OSError:
        pass


def _probe_range_support(url: str) -> int | None:
    """探测目标是否支持 Range 分段下载:返回文件总字节数;不支持返回 None。

    先 HEAD 看 Accept-Ranges;HEAD 不给准信再用 1 字节 Range GET 复探
    (看 206 与 Content-Range 总长)。两个探测请求都极小,不占实质带宽。
    """
    try:
        r = requests.head(url, timeout=(30, 30), headers=_NCBI_HEADERS,
                          allow_redirects=True)
        total = int(r.headers.get("content-length", 0) or 0)
        if r.ok and r.headers.get("accept-ranges", "").lower() == "bytes" and total:
            return total
    except requests.RequestException:
        pass
    try:
        r = requests.get(url, timeout=(30, 30),
                         headers=dict(_NCBI_HEADERS, Range="bytes=0-0"),
                         stream=True)
        try:
            if r.status_code == 206:
                cr = r.headers.get("content-range", "")  # 形如 bytes 0-0/8471521
                total = int(cr.rsplit("/", 1)[-1]) if "/" in cr else 0
                if total:
                    return total
        finally:
            r.close()
    except (requests.RequestException, ValueError):
        pass
    return None


def _fetch_segment(url: str, start: int, end: int, seg_path: str,
                   connect_timeout: int, read_timeout: int) -> None:
    """下载单个分段到 seg_path(断点续传:从段文件当前大小接着要)。

    读超时即断流判定:连续 read_timeout 秒无任何字节到达则抛异常,
    该段从已写入处续传,不动其他段。
    """
    have = os.path.getsize(seg_path) if os.path.isfile(seg_path) else 0
    want = end - start + 1
    if have >= want:
        return
    headers = dict(_NCBI_HEADERS, Range=f"bytes={start + have}-{end}")
    with requests.get(url, headers=headers,
                      timeout=(connect_timeout, read_timeout), stream=True) as r:
        if r.status_code == 200 and have > 0:
            have = 0  # 服务端无视 Range 发来全量:本段从头重写
        elif r.status_code != 206:
            r.raise_for_status()
        with open(seg_path, "ab" if have else "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 18):
                f.write(chunk)


def _download_one_pass(url: str, dest: str, config) -> None:
    """一轮完整下载:支持 Range 的大文件分段并发,否则单流;产物落 dest。"""
    import shutil
    from concurrent.futures import ThreadPoolExecutor

    total = _probe_range_support(url)
    n_conn = max(1, int(getattr(config, "GEO_DOWNLOAD_CONCURRENCY", 4)))
    stall = int(getattr(config, "GEO_DOWNLOAD_STALL_TIMEOUT", 90))
    part = dest + ".part"

    if total is None or total < n_conn * (1 << 20):
        # 不支持分段或文件太小:单流直下
        got = 0
        with requests.get(url, timeout=(30, stall),
                          headers=_NCBI_HEADERS, stream=True) as r:
            r.raise_for_status()
            expected = int(r.headers.get("content-length", 0) or 0)
            with open(part, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    f.write(chunk)
                    got += len(chunk)
        if expected and got != expected:
            raise OSError(f"size mismatch: got {got}, expected {expected}")
    else:
        # 均匀切段;段文件名编码布局(连接数×总长×段号):改并发数或换文件后
        # 旧布局段文件自动作废 —— 若复用旧段,`have >= want` 会把字节范围错位
        # 的段误判为已完成,拼装总数对不上且重试永远修不好(审查 BUG-13)
        seg_paths = [f"{dest}.p{n_conn}x{total}.{i}" for i in range(n_conn)]
        bounds = []
        step = total // n_conn
        for i in range(n_conn):
            start = i * step
            end = (total - 1) if i == n_conn - 1 else (start + step - 1)
            bounds.append((start, end, seg_paths[i]))
        with ThreadPoolExecutor(max_workers=n_conn) as pool:
            futures = [pool.submit(_fetch_segment, url, s, e, sp, 30, stall)
                       for s, e, sp in bounds]
            for fut in futures:
                fut.result()  # 任一段异常上抛 → 外层重试,已完成的段保留
        got_total = sum(os.path.getsize(sp) for sp in seg_paths)
        if got_total != total:
            raise OSError(f"segments total mismatch: {got_total}, expected {total}")
        with open(part, "wb") as out:
            for sp in seg_paths:
                with open(sp, "rb") as f:
                    shutil.copyfileobj(f, out, 1 << 20)
        # 成品已替换,清掉本布局段文件与历史遗留的其他布局段/.part 残留
        for leftover in glob.glob(glob.escape(dest) + ".p*"):
            try:
                os.remove(leftover)
            except OSError:
                pass

    if not _gzip_intact(part):
        raise OSError("gzip integrity check failed (truncated or corrupt)")
    os.replace(part, dest)  # 原子替换:要么完整新文件,要么保留旧缓存


def download_geo_file(url: str, dest: str) -> bool:
    """GEO 文件下载统一入口(series matrix 与 GPL 注释共用):
    并发锁 + 分段并发 + 断点续传 + 双重校验。

    锁防同文件并发写坏 .part;分段并发把跨境对 NCBI 单连接的限速聚合起来
    (方案 §4.1 原定"不做下载并发",实测单连接被压到个位数 KB/s、8.5MB 爬
    16 分钟后修订;E-utilities 的 3 req/s 限流只约束检索接口,不约束文件服
    务器,单文件 4 连接在 NCBI 容忍范围内)。每段独立落盘,断流 90s 即判,
    重试只补缺段;耗尽后本地有旧缓存则降级沿用。

    返回 True=本次拿到新文件(或并发另一请求已完成下载);False=重试耗尽后
    降级沿用了本地旧缓存 —— 调用方必须向用户暴露该情况(方案失效场景②),
    不得静默。
    """
    config = _config()
    if not acquire_download_lock(dest):
        logging.info("并发的同文件下载已完成,直接使用缓存: %s", dest)
        return True
    last_exc: Exception | None = None
    try:
        for attempt in range(1, config.GEO_DOWNLOAD_RETRIES + 1):
            try:
                _download_one_pass(url, dest, config)
                logging.info("GEO download OK: %s", url)
                return True
            except Exception as exc:
                last_exc = exc
                logging.exception("GEO download failed (attempt %d/%d): %s",
                                  attempt, config.GEO_DOWNLOAD_RETRIES, url)
                if attempt < config.GEO_DOWNLOAD_RETRIES:
                    time.sleep(config.GEO_DOWNLOAD_RETRY_DELAY)
        if os.path.isfile(dest) and _gzip_intact(dest):
            logging.warning("下载耗尽,降级使用本地旧缓存: %s", dest)
            return False
        raise RuntimeError(
            f"GEO 文件下载失败({config.GEO_DOWNLOAD_RETRIES} 次重试后): {url}。"
            f"可选出路:① 稍后重试(跨境带宽波动大);② 用浏览器手动下载该文件后"
            f"放到 {dest};③ 换一个规模更小的数据集。最后错误:{last_exc}")
    finally:
        release_download_lock(dest)


def fetch_geo_series_matrix(gse: str, force_refresh: bool = False) -> tuple[str, bool]:
    """下载并永久缓存 series matrix,返回(本地路径, 是否降级沿用了旧缓存)。"""
    gse = gse.strip().upper()
    path = _geo_cache_path(gse)
    if os.path.isfile(path) and not force_refresh and _gzip_intact(path):
        logging.info("GEO series matrix 缓存命中: %s", path)
        return path, False
    fell_back = not download_geo_file(series_matrix_url(gse), path)
    return path, fell_back


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
    """下载(或命中缓存)+ 解析,一步到位。

    下载耗尽降级沿用旧缓存时,在产物上打 cache_fallback 标并记缓存日期
    (方案失效场景②:CLI 打印警告、页面显示黄色警告条 + 强制更新按钮)。
    """
    path, fell_back = fetch_geo_series_matrix(gse, force_refresh=force_refresh)
    ds = parse_series_matrix(path)
    if fell_back:
        ds.cache_fallback = True
        ds.cache_date = time.strftime("%Y-%m-%d %H:%M",
                                      time.localtime(os.path.getmtime(path)))
    return ds


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
