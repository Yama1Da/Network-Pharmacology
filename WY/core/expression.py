"""表达谱预处理与差异表达(DEG)模块。

- log2 启发式检测 + 手动覆盖(不同数据集提交的标准化状态不一);
- 探针级数据按 GPL 注释映射基因符号,多探针对一基因取 IQR 最大者;
  symbol 级数据直通;映射覆盖率 < 50% 报错,不做激进猜测;
- Welch t 检验 + 自实现 BH 校正(不为此引 statsmodels)。
"""

import gzip
import logging
import os

import numpy as np
import pandas as pd
import requests

_PROBE_SUFFIX = "_at"
_GPL_HEADERS = {"User-Agent": "NetworkPharmacologyPlatform/1.0"}


def _config():
    import config
    return config


# ---------------- 预处理:log2 启发式 ----------------


def looks_symbol_level(index_values) -> bool:
    """行标识是否已是基因符号级。

    排除四类常见探针/非符号标识(全面审核 BUG-1 修复):
    - Affy 探针: 以 _at 结尾(如 1007_s_at)
    - Ensembl ID: ENSG/ENSMUSG 等前缀
    - Illumina 探针: ILMN_ 前缀
    - Agilent 探针: A_XX_P 前缀(如 A_23_P123456)
    """
    strs = [str(x) for x in index_values[:200]]
    if not strs:
        return False
    # 探针/非符号模式匹配
    probe_patterns = (
        _PROBE_SUFFIX,           # Affy _at
        "ENSG", "ENSMUSG", "ENSRNOG",  # Ensembl
        "ILMN_",                 # Illumina
        "A_",                    # Agilent (A_23_P...)
    )
    probe_like = sum(1 for s in strs
                     if any(s.startswith(p) or s.endswith(p)
                            for p in probe_patterns))
    digit_like = sum(1 for s in strs if s.replace(".", "").isdigit())
    alpha_like = sum(1 for s in strs if any(c.isalpha() for c in s))
    return (probe_like / len(strs) < 0.05
            and digit_like / len(strs) < 0.5
            and alpha_like / len(strs) > 0.8)


def preprocess(matrix: pd.DataFrame, force_log: bool | None = None,
               ) -> tuple[pd.DataFrame, str]:
    """log2 变换决策,返回(处理后的矩阵, 实际采用的分支说明)。

    启发式(全面审核 BUG-3 修复:尺度判定优先于负值):
    1. 95 分位 > 100 → 线性尺度,做 log2(x+1)(即使含负值也先做,
       因为负值可能来自背景校正而非已 log);
    2. 含负值且分位不高 → 已归一化数据(如 z-score),不再变换;
    3. 其余(分位不高、无负值)→ 视为已 log。
    force_log True/False 可覆盖(summary 记录实际分支,保证可追溯)。
    """
    num = matrix.select_dtypes(include=[np.number])
    if force_log is True:
        out, branch = np.log2(num.clip(lower=0) + 1), "log2(手动指定)"
    elif force_log is False:
        out, branch = num, "不变换(手动指定)"
    elif float(num.quantile(0.95).quantile(0.95)) > 100:
        out, branch = np.log2(num.clip(lower=0) + 1), "log2(x+1)(启发式:线性尺度)"
    elif (num < 0).any().any():
        out, branch = num, "不变换(负值且非大尺度,判定为已归一化)"
    else:
        out, branch = num, "不变换(启发式:判定为已 log)"
    return out, branch


# ---------------- 探针 → 基因符号 ----------------


def _gpl_annot_path(gpl: str) -> str:
    config = _config()
    d = os.path.join(config.CACHE_DIR, config.GEO_SERIES_CACHE_SUBDIR)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{gpl}_annot.gz")


def _gpl_annot_url(gpl: str) -> str:
    from .geo import _geo_bucket
    return (f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{_geo_bucket(gpl)}/{gpl}/"
            f"annot/{gpl}.annot.gz")


def _fetch_gpl_annotation(gpl: str) -> str:
    """下载 GPL 注释表(GEO 平台基因注释,含 ID ↔ Gene symbol 列),永久缓存。"""
    path = _gpl_annot_path(gpl)
    if os.path.isfile(path):
        return path
    config = _config()
    url = _gpl_annot_url(gpl)
    last: Exception | None = None
    for attempt in range(1, config.GEO_DOWNLOAD_RETRIES + 1):
        try:
            r = requests.get(url, timeout=config.GEO_DOWNLOAD_TIMEOUT,
                             headers=_GPL_HEADERS)
            r.raise_for_status()
            if not r.content[:2] == b"\x1f\x8b":
                raise OSError("非 gzip 内容(该平台可能无 annot 文件)")
            part = path + ".part"
            with open(part, "wb") as f:
                f.write(r.content)
            os.replace(part, path)
            logging.info("GPL annotation downloaded: %s", url)
            return path
        except Exception as exc:
            last = exc
            logging.exception("GPL annotation download failed (attempt %d): %s",
                              attempt, url)
            if attempt < config.GEO_DOWNLOAD_RETRIES:
                import time
                time.sleep(config.GEO_DOWNLOAD_RETRY_DELAY)
    raise RuntimeError(
        f"GPL 注释下载失败: {url}。该平台可能不提供 annot 文件 —— "
        "请换用 symbol 级数据集,或在 GEO 页面确认平台注释后手动放置注释文件到 "
        f"{path}。最后错误:{last}")


def map_probe_to_symbol(matrix: pd.DataFrame, gpl: str) -> pd.DataFrame:
    """把探针级矩阵按 GPL 注释映射为基因符号矩阵。

    多探针对一基因取 IQR 最大者(代表性与区分度最好的那条探针);
    映射覆盖率 < 50% 直接报错(平台冷门/注释缺失),不猜测。
    """
    if looks_symbol_level(matrix.index):
        logging.info("行标识为基因符号级,跳过映射")
        return matrix

    if not gpl:
        raise RuntimeError("探针级数据缺少平台(GPL)信息,无法映射基因符号;请换 symbol 级数据集")
    annot_path = _fetch_gpl_annotation(gpl.strip())
    mapping: dict[str, str] = {}
    with gzip.open(annot_path, "rt", encoding="utf-8", errors="replace") as f:
        platform_table = False
        header: list[str] | None = None
        id_col = sym_col = None
        for line in f:
            if line.startswith("!platform_table_begin"):
                platform_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if platform_table:
                import csv
                import io
                parts = next(csv.reader(io.StringIO(line.rstrip("\r\n")),
                                        delimiter="\t", quotechar='"'))
                if header is None:
                    header = parts
                    for cand_id in ("ID", "PROBE_ID"):
                        if cand_id in header:
                            id_col = cand_id
                    for cand_sym in ("Gene symbol", "Gene Symbol", "GENE_SYMBOL",
                                     "Symbol", "gene_assignment"):
                        if cand_sym in header and cand_sym != "gene_assignment":
                            sym_col = cand_sym
                    if id_col is None or sym_col is None:
                        raise RuntimeError(
                            f"GPL 注释缺少 ID 或 Gene symbol 列(表头:{header[:8]}...);"
                            "请换 symbol 级数据集")
                    continue
                probe = parts[header.index(id_col)]
                symbol = parts[header.index(sym_col)].strip()
                if probe and symbol and "---" not in symbol:
                    mapping[probe.strip()] = symbol

    idx = pd.Series([str(x).strip() for x in matrix.index], index=matrix.index)
    mapped = idx.map(mapping)
    coverage = float(mapped.notna().mean())
    if coverage < 0.5:
        raise RuntimeError(
            f"探针→基因符号映射覆盖率仅 {coverage:.0%}(< 50%):平台注释不完整,"
            "结果不可信 —— 请换 symbol 级数据集或注释完善的平台")
    keep = matrix.loc[mapped.notna()].copy()          # 行保持探针级(索引唯一)
    symbols = mapped[mapped.notna()]
    # 多探针对一基因:IQR 最大者胜出(代表该基因表达差异能力最强的一条探针)
    iqr = keep.quantile(0.75, axis=1) - keep.quantile(0.25, axis=1)
    cand = (pd.DataFrame({"symbol": symbols.values, "iqr": iqr.values},
                         index=keep.index)
            .sort_values("iqr", ascending=False)
            .drop_duplicates("symbol"))
    out = keep.loc[cand.index].copy()
    out.index = cand["symbol"]
    out = out[~out.index.duplicated(keep="first")].sort_index()
    logging.info("探针映射: %d 探针 → %d 基因(覆盖率 %.0f%%)",
                 len(keep), len(out), coverage * 100)
    return out


# ---------------- 差异表达:Welch t + BH ----------------


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """BH-FDR 校正(Benjamini & Hochberg 1995 标准四步,自实现 ~10 行)。"""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]  # 单调性修正
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def differential_expression(matrix: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    """按标签做两组 Welch t 检验 + BH 校正,输出 gene|log2FC|t|p|FDR|direction。

    log2FC = case 均值 − control 均值(输入应为预处理后的矩阵);
    过滤两组中全为 NaN 的基因后计算,方向列便于与机制假设比对。
    """
    from scipy.stats import ttest_ind

    labels = labels.loc[labels.index.intersection(matrix.columns)]
    case_cols = labels[labels == 1].index
    ctrl_cols = labels[labels == 0].index
    mat = matrix.loc[:, list(case_cols) + list(ctrl_cols)]
    mat = mat.dropna(how="all")
    a = mat[list(case_cols)].to_numpy(dtype=float)
    b = mat[list(ctrl_cols)].to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        t_stat, pvals = ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
        log2fc = (np.nanmean(a, axis=1) - np.nanmean(b, axis=1))
    valid = np.isfinite(pvals)
    res = pd.DataFrame({
        "gene": mat.index,
        "log2FC": log2fc,
        "t": t_stat,
        "p": pvals,
    })[valid].reset_index(drop=True)
    res["FDR"] = benjamini_hochberg(res["p"].to_numpy())
    res["direction"] = np.where(res["log2FC"] > 0, "up",
                                np.where(res["log2FC"] < 0, "down", "flat"))
    return res.sort_values("p").reset_index(drop=True)
