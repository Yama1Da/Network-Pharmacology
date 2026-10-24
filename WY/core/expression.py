"""表达谱预处理与差异表达(DEG)模块。

- log2 启发式检测 + 手动覆盖(不同数据集提交的标准化状态不一);
- 探针级数据按 GPL 注释映射基因符号,多探针对一基因取 IQR 最大者;
  symbol 级数据直通;映射覆盖率 < 50% 报错,不做激进猜测;
- Welch t 检验 + 自实现 BH 校正(不为此引 statsmodels)。
"""

import gzip
import logging
import os
import time

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


def normalize_gpl(gpl: str) -> str:
    """平台编号规范化:检索元数据里的 gpl 可能不带 GPL 前缀(如 "26963"),
    统一补成 GPL26963 —— 否则 annot URL 与缓存文件名都会拼错、探测恒 404。"""
    gpl = gpl.strip().upper()
    return gpl if gpl.startswith("GPL") else f"GPL{gpl}"


def _gpl_annot_path(gpl: str) -> str:
    config = _config()
    d = os.path.join(config.CACHE_DIR, config.GEO_SERIES_CACHE_SUBDIR)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{normalize_gpl(gpl)}_annot.gz")


def _gpl_annot_url(gpl: str) -> str:
    from .geo import _geo_bucket
    gpl = normalize_gpl(gpl)
    return (f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{_geo_bucket(gpl)}/{gpl}/"
            f"annot/{gpl}.annot.gz")


def _fetch_gpl_annotation(gpl: str) -> str:
    """下载 GPL 注释表(GEO 平台基因注释,含 ID ↔ Gene symbol 列),永久缓存。

    GEO 只为部分平台生成 annot 文件:新式 Agilent/Illumina 平台常缺失
    (404 属服务器端"文件不存在",确定性问题),这类平台表往往连符号列都
    没有、只有转录本 ID。404 直接给"换数据集"的准确报错,不走 3 轮重试
    ——重试与手动下载对不存在的文件都无意义。
    """
    from .geo import download_geo_file

    path = _gpl_annot_path(gpl)
    if os.path.isfile(path):
        return path
    url = _gpl_annot_url(gpl.strip())
    try:
        head = requests.head(url, timeout=(15, 30), headers=_GPL_HEADERS)
        if head.status_code == 404:
            raise RuntimeError(
                f"GEO 未为平台 {gpl} 提供基因注释文件(服务器 404,常见于新式 "
                "Agilent/Illumina 平台)。这不是网络问题,重试与手动下载均无效"
                " —— 请换 symbol 级数据集,或选用平台注释齐全的数据集。")
    except requests.RequestException:
        pass  # 探测失败不阻断:交给 download_geo_file 的重试链路给出网络类报错
    download_geo_file(url, path)
    return path


# annot 探测结果短时记忆:10 分钟 TTL,只存确定性结论(见 probe_platform_annot)
_ANNOT_PROBE_TTL = 600.0
_ANNOT_PROBE_CACHE: dict[str, tuple[float, str]] = {}


def probe_platform_annot(gpl: str) -> str:
    """探测平台注释文件可用性(选型预检,单次 HEAD,不落盘):
    返回 "cached"(本地已有)| "ok"(服务器存在)| "missing"(404,GEO 未生成,
    确定性缺失)| "unknown"(网络原因,开始分析时会经重试链路再判)。

    结果短时记忆(10 分钟 TTL):总览每次点击都会走到这里,对探针级平台反复
    HEAD 没有意义。只缓存确定性结论,"unknown" 不缓存 —— 网络恢复后可自愈。
    本地文件存在优先于记忆判断:手动放置注释文件后立即翻转为 "cached"。
    """
    gpl = normalize_gpl(gpl)
    if os.path.isfile(_gpl_annot_path(gpl)):
        return "cached"
    cached = _ANNOT_PROBE_CACHE.get(gpl)
    if cached and time.time() - cached[0] < _ANNOT_PROBE_TTL:
        return cached[1]
    try:
        status = requests.head(_gpl_annot_url(gpl),
                               timeout=(15, 30), headers=_GPL_HEADERS).status_code
    except requests.RequestException:
        return "unknown"
    if status == 200:
        result = "ok"
    elif status == 404:
        result = "missing"
    else:
        return "unknown"
    _ANNOT_PROBE_CACHE[gpl] = (time.time(), result)
    return result


# 预检的分组字段噪声词:字段名命中即不算候选分组列(年龄/性别等样本属性,非病/对照)
_GROUP_FIELD_DENY = ("age", "sex", "gender", "date", "batch")
# 性别取值词表:候选列的取值规约到该词表内即整列剔除(值级判断,覆盖
# "Sex: male/female"、单字母 "M/F"、中文 "男/女" 等编码 —— 字段名不含
# sex 字样时也拦得住;性别是样本属性,混进候选会误导用户对它写分组正则)
_SEX_VALUE_TOKENS = {"m", "f", "male", "female", "man", "woman", "男", "女"}


def assess_expression_dataset(ds) -> dict:
    """数据集可用性预检(实施方案 §4.1/§11:平台注释缺失、无病例对照设计
    → 选型阶段就排除,而不是下载分组之后才失败)。

    对已解析的 GeoDataset 体检,返回页面可直接消费的字典:
    matrix_empty  表达矩阵为空(元数据型提交/超级系列)→ 不可用
    symbol_level  行标识是否基因符号级(True=免平台注释映射)
    shape / gpl   矩阵形状与平台编号(供展示)
    annot         平台注释状态(cached/ok/missing/unknown;基因级为 None)
    group_fields  候选分组字段(取值 2-8 类的注释列,剔除年龄/性别等属性列),
                  空列表 = 可能没有病例/对照设计(启发式,人工仍需核对标题)
    n_samples     样本数(页面对照两组下限给软提示)
    usable        综合判定:矩阵非空 且 (基因级 或 注释可用) 且 有候选分组字段
    """
    n = len(ds.samples)
    if ds.table.empty:
        # 空矩阵 ≠ 判死:RNA-seq 数据集的 counts 在补充文件里,运行时会自动回退
        # DESeq2 路径 —— 预检此时探测 counts 可用性,有则放行(否则锁定才正确)
        try:
            from .geo import has_counts_supplement
            counts_available = has_counts_supplement(ds.gse)
        except Exception:
            counts_available = False
        return {"matrix_empty": True, "symbol_level": None, "shape": None,
                "gpl": ds.gpl, "annot": None, "group_fields": [],
                "counts_available": counts_available,
                "n_samples": n, "usable": counts_available}

    symbol = looks_symbol_level(ds.table.index)
    if symbol:
        annot = None           # 基因级矩阵免映射,注释文件无关紧要
    elif ds.gpl:
        annot = probe_platform_annot(ds.gpl)
    else:
        annot = "missing"      # 探针级又无平台信息,映射无从谈起

    group_fields = []
    for col in ds.samples.columns:
        if not (col.startswith("characteristics") or col in ("title", "source_name_ch1")):
            continue
        if any(k in col.lower() for k in _GROUP_FIELD_DENY):
            continue
        uniq = ds.samples[col].astype(str).unique()
        if not 2 <= len(uniq) <= 8:
            continue
        lows = {v.strip().lower() for v in uniq}
        # 值级性别判断:每个取值取冒号后最后一段词,全部落进性别词表 → 整列剔除
        # (覆盖 "Sex: male/female"、单字母 "M/F"、中文 "男/女" 等编码)
        toks = {v.split(":")[-1].strip().split()[-1]
                for v in lows if v.split(":")[-1].strip().split()}
        if toks and toks <= _SEX_VALUE_TOKENS:
            continue
        group_fields.append(f"{col}: {'/'.join(sorted(map(str, uniq))[:4])}")

    usable = bool(group_fields) and (symbol or annot in ("cached", "ok"))
    return {"matrix_empty": False, "symbol_level": symbol, "shape": tuple(ds.table.shape),
            "gpl": ds.gpl, "annot": annot, "group_fields": group_fields,
            "counts_available": None,
            "n_samples": n, "usable": usable}


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


# ---------------- raw counts 判定与 DESeq2 差异表达(二期 RNA-seq 路径) ----------------


def looks_count_matrix(matrix: pd.DataFrame) -> bool:
    """行=基因 矩阵是否为 raw counts(全为非负整数)。

    TPM/FPKM/归一化矩阵含小数,log 化矩阵可能含负值 —— DESeq2 只接受
    非负整数 counts,喂错输入统计上无效,必须在此挡下。
    """
    num = matrix.select_dtypes(include=[np.number])
    if num.empty:
        return False
    v = num.to_numpy(dtype=float)
    if v.size == 0 or np.isnan(v).any() or (v < 0).any():
        return False
    return bool((v % 1 == 0).all())


def _deseq2_run(counts: pd.DataFrame, labels: pd.Series):
    """DESeq2 共享实现:校验 + 一次拟合,返回 (results_df, 样本×基因 标准化counts)。"""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise RuntimeError(
            "pydeseq2 未安装(RNA-seq counts 路径依赖,见 requirements.txt "
            f"的 pydeseq2==0.5.4): {exc}")

    if not looks_count_matrix(counts):
        raise RuntimeError(
            "输入矩阵不是 raw counts(存在小数/负值/缺失,疑似 TPM/FPKM 或已"
            "归一化矩阵)—— DESeq2 只接受非负整数 counts,请提供原始 counts。")
    labels = labels.loc[labels.index.intersection(counts.columns)]
    case_cols = labels[labels == 1].index
    ctrl_cols = labels[labels == 0].index
    # 守卫:交集为空或某组不足 2 个样本时,pyDESeq2 会抛深堆栈的底层数值错误
    # (或静默产出空表),在此给出可读的报错
    if len(case_cols) == 0 or len(ctrl_cols) == 0:
        raise RuntimeError(
            "分组标签与 counts 矩阵的样本列无交集(或某组为 0)—— "
            f"case={len(case_cols)}, control={len(ctrl_cols)};请检查样本命名"
            "是否与 counts 列名一致")
    if len(case_cols) < 2 or len(ctrl_cols) < 2:
        raise RuntimeError(
            f"每组至少需要 2 个样本才能估计离散度(case={len(case_cols)}, "
            f"control={len(ctrl_cols)});请补充样本或换数据集")
    X = counts.loc[:, list(case_cols) + list(ctrl_cols)].T.astype(int)  # 样本×基因
    meta = pd.DataFrame(
        {"group": ["case"] * len(case_cols) + ["control"] * len(ctrl_cols)},
        index=X.index)
    X = X.loc[:, X.sum(axis=0) > 0]  # 全零基因:DESeq2 必置 NaN,提前剔除

    try:
        dds = DeseqDataSet(counts=X, metadata=meta, design="~group")
    except TypeError:  # 旧版参数名(design_factors)
        dds = DeseqDataSet(counts=X, metadata=meta, design_factors="group")
    dds.deseq2()
    stats = DeseqStats(dds, contrast=("group", "case", "control"))
    stats.summary()
    norm = pd.DataFrame(np.asarray(dds.layers["normed_counts"]),
                        index=X.index, columns=X.columns)
    return stats.results_df, norm


def _deg_table_from_results(res: pd.DataFrame) -> pd.DataFrame:
    out = res.copy()
    out.insert(0, "gene", out.index)
    out = out[["gene", "log2FoldChange", "stat", "pvalue", "padj"]].dropna(
        subset=["pvalue"])
    out.columns = ["gene", "log2FC", "stat", "p", "FDR"]
    out["direction"] = np.where(out["log2FC"] > 0, "up",
                                np.where(out["log2FC"] < 0, "down", "flat"))
    return out.sort_values("p").reset_index(drop=True)


def differential_expression_counts(counts: pd.DataFrame,
                                   labels: pd.Series) -> pd.DataFrame:
    """raw counts 的差异表达(pyDESeq2 负二项 GLM + Wald 检验,二期 RNA-seq 路径)。

    输入与 differential_expression 同构: counts 为 基因×样本 整数矩阵,
    labels 样本→1(case)/0(control)。输出同构 DEG 表:
    gene|log2FC|stat|p|FDR|direction —— stat 为 Wald 统计量;FDR 取
    pyDESeq2 的 padj(BH + 独立过滤,优于朴素 BH)。全零/未通过检验的基因
    (p 为 NaN)按现路径 valid 口径过滤。行标识可为 Ensembl ID,符号映射由
    调用方决定(职责同 map_probe_to_symbol)。需要 pydeseq2(见 requirements)。
    """
    res, _ = _deseq2_run(counts, labels)
    return _deg_table_from_results(res)


def rnaseq_deg_and_matrix(counts: pd.DataFrame,
                          labels: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    """DESeq2 一次拟合,同时产出 DEG 表与 ML 用表达矩阵。

    ML 不吃 raw counts(方差随均值增长,树/线性模型会被高表达基因主导),
    也不直接用 log2 raw counts(未去文库深度)—— 用 log2(尺寸因子标准化
    counts + 1),是 VST 的廉价折中(VST 更优但耗时,一期取折中)。
    返回 (deg 同构表, ml_matrix: 样本×基因)。
    """
    res, norm = _deseq2_run(counts, labels)
    deg = _deg_table_from_results(res)
    ml_matrix = np.log2(norm.clip(lower=0) + 1)
    return deg, ml_matrix


# ---------------- Ensembl → 符号映射(路线 A:mygene 批量,懒映射) ----------------

_MYGENE_URL = "https://mygene.info/v3/query"


def map_ids_to_symbol(gene_ids, chunk: int = 1000) -> dict:
    """批量把基因 ID 映射为 HGNC 符号(mygene.info REST,方案二期路线 A)。

    设计要点:
    - 懒映射:调用方只传 DEG 预过滤池(几百个 ID),不映射全矩阵;
    - 非 Ensembl ID(已是大写符号等)不经网络,直接原样直通 —— 零成本;
    - Ensembl ID 分块批量 POST;未命中/多义取第一个 symbol,未命中记 None;
    - 网络失败响亮抛错(调用方负责降级保留原 ID),HTTP 层走 request_with_retry
      的统一重试;POST 响应不被 requests-cache 缓存(GET-only),重复调用即重查,
      结果由调用方落盘 symbol_map.csv 保证可复现。
    返回 {输入ID: 符号 或 None}(只含入参中出现过的 ID)。
    """
    from .http_client import request_with_retry

    ids = [str(g).strip() for g in gene_ids]
    ens = [g for g in ids if g.upper().startswith("ENS")]
    mapping: dict = {g: g for g in ids if not g.upper().startswith("ENS")}
    for i in range(0, len(ens), chunk):
        part = ens[i:i + chunk]
        resp = request_with_retry(
            _MYGENE_URL, method="post", max_retries=3, timeout=60,
            label=f"mygene 批量符号映射({len(part)} 个 ID)",
            data={"q": "\n".join(part), "scopes": "ensembl.gene",
                  "fields": "symbol", "species": "all"},
            headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"},
            raise_on_exhausted=True)
        payload = resp.json()
        # 批量 POST 返回平铺的 hit 列表(每项自带 query 字段);兼容字典形态
        hits = payload if isinstance(payload, list) else payload.get("hits", [])
        for hit in hits:
            q = hit.get("query")
            if q in mapping or q not in part:   # 同 ID 多 hit 取首个
                continue
            mapping[q] = hit.get("symbol")      # 未命中为 None
    return mapping


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
