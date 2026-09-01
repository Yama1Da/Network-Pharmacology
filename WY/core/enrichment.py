"""富集分析模块：g:Profiler 优先 + Enrichr 降级，含熔断器与结果绘图。"""

import hashlib
import logging
import os

import pandas as pd

from .caching import _is_cache_fresh, _read_cache_csv

# gseapy (Enrichr) 基因集库名与 g:Profiler sources 的映射关系：
#   g:Profiler source    →  Enrichr 库名
#   GO:BP                →  GO_Biological_Process_2021
#   GO:MF                →  GO_Molecular_Function_2021
#   GO:CC                →  GO_Cellular_Component_2021
#   KEGG                 →  KEGG_2021_Human
ENRICHR_LIB_MAP = {
    "GO:BP": "GO_Biological_Process_2021",
    "GO:MF": "GO_Molecular_Function_2021",
    "GO:CC": "GO_Cellular_Component_2021",
    "KEGG": "KEGG_2021_Human",
}

# g:Profiler 物种名 → Enrichr 物种名。
# 注意：Enrichr 只支持 human / mouse 等，大鼠（rnorvegicus）无对应库，
# 因此大鼠场景下降级 Enrichr 会被跳过，仅依赖 g:Profiler。
_ENRICHR_ORGANISM_MAP = {
    "hsapiens": "human",
    "mmusculus": "mouse",
}


def _enrichment_cache_key(genes: list[str], organism: str, sources, pvalue: float) -> str:
    """根据基因列表与查询参数生成富集缓存文件名（md5）：相同输入命中同一缓存文件。"""
    payload = (
        "\r".join(sorted(set(genes)))
        + f"|{organism}|{','.join(sources)}|{pvalue}"
    )
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


# ---------- g:Profiler 熔断（circuit breaker） ----------
# g:Profiler 服务器不可达时，单次连接可能空等 20s+ 才超时（socket.setdefaulttimeout
# 对库内部 requests 调用不生效），重试次数越多拖得越久。因此失败时写一个带时间戳的
# 标记文件，TTL 内直接跳过 g:Profiler 走 Enrichr；TTL 过期后自动重试，服务恢复即复用。


def _gprofiler_circuit_file() -> str:
    """g:Profiler 熔断标记文件路径（内容为上次失败时间戳）。"""
    import config

    return os.path.join(config.ENRICHMENT_CACHE_DIR, ".gprofiler_unreachable")


def _gprofiler_circuit_open() -> bool:
    """熔断是否开启：失败标记存在且未超过 TTL。"""
    import config
    import time

    path = _gprofiler_circuit_file()
    if not os.path.exists(path):
        return False
    try:
        with open(path, "r", encoding="utf-8") as f:
            failed_at = float(f.read().strip())
    except Exception:
        return False
    return (time.time() - failed_at) < config.GPROFILER_CIRCUIT_TTL_SECONDS


def _gprofiler_circuit_open_trip() -> None:
    """记录 g:Profiler 本次失败时间（开启熔断）。"""
    import time

    path = _gprofiler_circuit_file()
    try:
        with open(path, "w", encoding="utf-8") as f:
            f.write(str(time.time()))
    except Exception:
        logging.exception("Failed to write g:Profiler circuit breaker file")


def _gprofiler_circuit_close() -> None:
    """g:Profiler 请求成功时清除熔断标记。"""
    path = _gprofiler_circuit_file()
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception:
        logging.exception("Failed to remove g:Profiler circuit breaker file")


def fetch_enrichment_gprofiler(
    genes: list[str],
    organism: str = "hsapiens",
    sources: tuple = ("GO:BP", "GO:MF", "GO:CC", "KEGG"),
    pvalue: float = 0.05,
    label: str | None = None,
    status: dict | None = None,
) -> pd.DataFrame:
    """
    通过 g:Profiler 官方客户端（gprofiler-official）获取功能富集结果。

    结果缓存在 Static/data/cache/enrichment 下，文件名为 gProfiler_{label}_{md5}.csv
    （label 为来源标识，如 药物名_疾病名），缓存首行为 # 元数据注释头；
    命中缓存时直接读取本地文件，不访问远端 API。

    Parameters
    ----------
    genes : list[str]
        输入基因符号列表。
    organism : str
        g:Profiler 物种名，如 "hsapiens"。
    sources : tuple
        g:Profiler 数据源，如 ("GO:BP", "GO:MF", "GO:CC", "KEGG")。
    pvalue : float
        显著性阈值（g:Profiler 的 user_threshold，默认即 g_SCS 校正后阈值）。
    label : str | None
        来源标识（如 "aspirin_cancer" 或 "manual"），用于缓存文件名前缀，便于用户识别。
    status : dict | None
        可选出参。提供时写入本次调用状态，供调用方区分"真无显著项"与"服务失败"：
          ok: 布尔；reason 取值
            "success"        请求成功且有显著结果
            "no_significant" 请求成功但无显著条目
            "cache"          命中本地缓存
            "request_failed" 重试后仍失败（网络/服务不可用）
            "circuit_open"   熔断开启，本次跳过 g:Profiler（未发起请求）
            "empty_input"    输入基因列表为空

    Returns
    -------
    pd.DataFrame
        统一列富集结果（source / term_id / term_name / p_value /
        adjusted_p_value / intersection_size / query_size / genes）；
        请求失败或无显著结果时返回空 DataFrame。
    """
    import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes:
        if status is not None:
            status.update(ok=False, reason="empty_input")
        return pd.DataFrame()
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    md5 = _enrichment_cache_key(genes, organism, sources, pvalue)
    name = f"gProfiler_{label}_{md5}.csv" if label else f"gProfiler_{md5}.csv"
    cache_file = os.path.join(config.ENRICHMENT_CACHE_DIR, name)
    if _is_cache_fresh(cache_file):
        logging.info("Loading cached g:Profiler enrichment for %d genes", len(genes))
        if status is not None:
            status.update(ok=True, reason="cache")
        return _read_cache_csv(cache_file)

    # 熔断：gProfiler 近期失败过则直接跳过，不再空等连接超时（走 Enrichr 兜底）
    if _gprofiler_circuit_open():
        logging.info("g:Profiler circuit breaker open — skipping, Enrichr fallback")
        if status is not None:
            status.update(ok=False, reason="circuit_open")
        return pd.DataFrame()

    import time
    import socket

    # 请求保护：gprofiler-official 内部 requests.post 未传 timeout、也无重试，
    # 这里在调用层统一补上【请求间隔 + 失败重试 + socket 层超时兜底】。
    # socket.setdefaulttimeout 兜底库内部的 requests 调用（未显式传 timeout 时生效），
    # try/finally 确保调用结束后恢复，不影响其他请求。
    from gprofiler import GProfiler

    res = None
    for attempt in range(config.ENRICHMENT_MAX_RETRIES):
        time.sleep(config.API_REQUEST_DELAY)  # 每次请求前等待，避免快速连发被限流
        socket.setdefaulttimeout(config.API_TIMEOUT)
        try:
            gp = GProfiler(return_dataframe=True)
            res = gp.profile(
                query=genes,
                organism=organism,
                sources=list(sources),
                user_threshold=pvalue,
                no_evidences=False,  # 需要返回交集基因（intersections），供 genes 列使用
            )
            break  # 成功
        except Exception:
            logging.exception(
                "g:Profiler enrichment failed (attempt %d/%d)",
                attempt + 1, config.ENRICHMENT_MAX_RETRIES,
            )
            _gprofiler_circuit_open_trip()  # 记录失败时间，开启熔断
            if attempt < config.ENRICHMENT_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)  # 失败退避后重试
        finally:
            socket.setdefaulttimeout(None)

    if res is None:
        # 重试后仍失败（网络/服务不可用）
        logging.warning("g:Profiler unreachable after %d retries for %d genes",
                        config.ENRICHMENT_MAX_RETRIES, len(genes))
        if status is not None:
            status.update(ok=False, reason="request_failed")
        return pd.DataFrame()
    if len(res) == 0:
        # 请求成功但无显著条目（真"空结果"，非失败）
        logging.info("g:Profiler: no significant terms (p < %s) for %d genes", pvalue, len(genes))
        if status is not None:
            status.update(ok=True, reason="no_significant")
        return pd.DataFrame()

    _gprofiler_circuit_close()  # 请求成功，关闭熔断
    df = _normalize_enrichment_gprofiler(res)
    if df.empty:
        # 防御：返回行全部无法归一化时不写缓存，避免"异常空表"被后续 TTL 命中
        logging.warning("g:Profiler result normalized to empty — skip caching (%d genes)", len(genes))
        if status is not None:
            status.update(ok=False, reason="normalize_empty")
        return df
    if status is not None:
        status.update(ok=True, reason="success")
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    logging.info("g:Profiler enrichment for %d genes: %d terms cached", len(genes), len(df))
    return df


def _normalize_enrichment_gprofiler(df: pd.DataFrame) -> pd.DataFrame:
    """把 g:Profiler 返回的 DataFrame 整理为统一列。

    g:Profiler 默认已做 g_SCS 多重检验校正，p_value 即校正后 p 值，
    因此 adjusted_p_value 与 p_value 相同。

    兼容 gprofiler-official 0.3.x 与 1.x 的列名差异：
    0.3.x 返回 term_id / term_name / intersection；
    1.x 返回 native / name / intersections（需 no_evidences=False）。
    """
    df = df.reset_index()  # return_dataframe=True 时 term_id 常在 index 中

    def _pick(*names) -> pd.Series:
        """按优先级取第一个存在的列，都不存在时返回全 NA 序列。"""
        for n in names:
            if n in df.columns:
                return df[n]
        return pd.Series(pd.NA, index=df.index, dtype=object)

    out = pd.DataFrame({
        "source": df["source"].astype(str),
        "term_id": _pick("term_id", "native"),
        "term_name": _pick("term_name", "name"),
        "p_value": pd.to_numeric(df["p_value"], errors="coerce"),
        "adjusted_p_value": pd.to_numeric(df["p_value"], errors="coerce"),
        "intersection_size": pd.to_numeric(df["intersection_size"], errors="coerce"),
        "query_size": pd.to_numeric(df["query_size"], errors="coerce"),
        "genes": _pick("intersection", "intersections").fillna(""),
    })
    out = out.dropna(subset=["term_id", "p_value"]).sort_values("p_value")
    out["term_id"] = out["term_id"].astype(str)
    out["term_name"] = out["term_name"].astype(str)
    return out.reset_index(drop=True)


def fetch_enrichment_enrichr(
    genes: list[str],
    libs: list[str],
    organism: str = "human",
    pvalue: float = 0.05,
    label: str | None = None,
    status: dict | None = None,
) -> pd.DataFrame:
    """
    通过 gseapy（Enrichr 官方 wrapper）获取功能富集结果（作为 g:Profiler 的降级方案）。

    结果同样缓存在 Static/data/cache/enrichment 下，文件名为 Enrichr_{label}_{md5}.csv，
    缓存首行为 # 元数据注释头。

    Parameters
    ----------
    genes : list[str]
        输入基因符号列表。
    libs : list[str]
        Enrichr 基因集库名（见 ENRICHR_LIB_MAP 的映射）。
    organism : str
        Enrichr 物种名（"human" / "mouse"）。
    pvalue : float
        显著性阈值（cutoff）。
    label : str | None
        来源标识（如 "aspirin_cancer" 或 "manual"），用于缓存文件名前缀，便于用户识别。
    status : dict | None
        可选出参，语义同 fetch_enrichment_gprofiler 的 status
        （reason 取值：success / no_significant / cache / request_failed / empty_input）。

    Returns
    -------
    pd.DataFrame
        统一列富集结果；请求失败或无显著结果时返回空 DataFrame。
    """
    import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes or not libs:
        if status is not None:
            status.update(ok=False, reason="empty_input")
        return pd.DataFrame()
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    md5 = _enrichment_cache_key(genes, organism, libs, pvalue)
    name = f"Enrichr_{label}_{md5}.csv" if label else f"Enrichr_{md5}.csv"
    cache_file = os.path.join(config.ENRICHMENT_CACHE_DIR, name)
    if _is_cache_fresh(cache_file):
        logging.info("Loading cached Enrichr enrichment for %d genes", len(genes))
        if status is not None:
            status.update(ok=True, reason="cache")
        return _read_cache_csv(cache_file)

    import time
    import socket

    # 请求保护：gseapy 内部 session 虽自带 5xx 重试，但部分请求 timeout=None
    # 无超时、且无请求间隔；这里在调用层统一补上【请求间隔 + 失败重试 + socket 层超时兜底】。
    # socket.setdefaulttimeout 兜底库内部的 requests 调用（未显式传 timeout 时生效），
    # try/finally 确保调用结束后恢复，不影响其他请求。
    import gseapy as gp

    res = None
    for attempt in range(config.ENRICHMENT_MAX_RETRIES):
        time.sleep(config.API_REQUEST_DELAY)  # 每次请求前等待，避免快速连发被限流
        socket.setdefaulttimeout(config.API_TIMEOUT)
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=list(libs),
                organism=organism,
                outdir=None,  # 不写盘，结果由我们自行缓存
                cutoff=pvalue,
            )
            res = enr.results
            break  # 成功
        except Exception:
            logging.exception(
                "Enrichr enrichment failed (attempt %d/%d)",
                attempt + 1, config.ENRICHMENT_MAX_RETRIES,
            )
            if attempt < config.ENRICHMENT_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)  # 失败退避后重试
        finally:
            socket.setdefaulttimeout(None)

    if res is None:
        # 重试后仍失败（网络/服务不可用）
        logging.warning("Enrichr unreachable after %d retries for %d genes",
                        config.ENRICHMENT_MAX_RETRIES, len(genes))
        if status is not None:
            status.update(ok=False, reason="request_failed")
        return pd.DataFrame()
    if len(res) == 0:
        # 请求成功但无显著条目（真"空结果"，非失败）
        if status is not None:
            status.update(ok=True, reason="no_significant")
        return pd.DataFrame()

    df = _normalize_enrichment_enrichr(res)
    if df.empty:
        # 防御：返回行全部无法归一化时不写缓存，避免"异常空表"被后续 TTL 命中
        logging.warning("Enrichr result normalized to empty — skip caching (%d genes)", len(genes))
        if status is not None:
            status.update(ok=False, reason="normalize_empty")
        return df
    if status is not None:
        status.update(ok=True, reason="success")
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    return df


def _normalize_enrichment_enrichr(df: pd.DataFrame) -> pd.DataFrame:
    """把 Enrichr 返回的 DataFrame 整理为统一列（并把库名反查回 g:Profiler 的 source 命名）。"""
    import re

    def _parse_overlap(s) -> tuple[int, int]:
        try:
            a, b = str(s).split("/")
            return int(a), int(b)
        except (ValueError, TypeError):
            return 0, 0

    def _extract_term_id(term) -> str:
        # Enrichr 的 Term 形如 "peptidyl-tyrosine modification (GO:0018212)"，提取末尾括号内的 ID
        m = re.search(r"\(([^()]+)\)\s*$", str(term))
        return m.group(1) if m else str(term)

    def _strip_suffix(term) -> str:
        return re.sub(r"\s*\([^()]*\)\s*$", "", str(term)).strip() or str(term)

    parsed = df["Overlap"].map(_parse_overlap)
    source_rev = {v: k for k, v in ENRICHR_LIB_MAP.items()}  # 库名 → g:Profiler source 命名
    out = pd.DataFrame({
        "source": df["Gene_set"].map(source_rev).fillna(df["Gene_set"]),
        "term_id": df["Term"].map(_extract_term_id),
        "term_name": df["Term"].map(_strip_suffix),
        "p_value": pd.to_numeric(df["P-value"], errors="coerce"),
        "adjusted_p_value": pd.to_numeric(df["Adjusted P-value"], errors="coerce"),
        "intersection_size": parsed.map(lambda t: t[0]),
        "query_size": parsed.map(lambda t: t[1]),
        "genes": df["Genes"].fillna(""),
    })
    out = out.dropna(subset=["term_id", "p_value"]).sort_values("p_value")
    return out.reset_index(drop=True)


def run_enrichment(
    genes: list[str],
    organism: str = "hsapiens",
    pvalue: float = 0.05,
    sources: tuple | None = None,
    label: str | None = None,
) -> tuple[pd.DataFrame, str | None]:
    """
    富集分析主入口：g:Profiler 优先，失败降级 Enrichr。

    Parameters
    ----------
    genes : list[str]
        输入基因符号列表（通常为药物∩疾病交集靶点）。
    organism : str
        g:Profiler 物种名（hsapiens / mmusculus / rnorvegicus）。
    pvalue : float
        显著性阈值。
    sources : tuple | None
        g:Profiler 数据源，默认取 config.ENRICHMENT_SOURCES（GO三分类 + KEGG）。
    label : str | None
        来源标识（如 "aspirin_cancer" 或 "manual"），用于缓存文件名前缀，便于用户识别。

    Returns
    -------
    tuple[pd.DataFrame, str | None]
        (统一列富集结果, 实际使用的数据源 "gProfiler" / "Enrichr")；
        两者都失败时返回 (空表, None)。
    """
    import config

    sources = sources or config.ENRICHMENT_SOURCES

    # 主：g:Profiler
    df = fetch_enrichment_gprofiler(genes, organism=organism, sources=sources, pvalue=pvalue, label=label)
    if not df.empty:
        return df, "gProfiler"

    # 备：Enrichr（不支持大鼠，跳过）
    enrichr_organism = _ENRICHR_ORGANISM_MAP.get(organism)
    if enrichr_organism:
        libs = [ENRICHR_LIB_MAP[s] for s in sources if s in ENRICHR_LIB_MAP]
        df = fetch_enrichment_enrichr(genes, libs=libs, organism=enrichr_organism, pvalue=pvalue, label=label)
        if not df.empty:
            return df, "Enrichr"
    return pd.DataFrame(), None


def render_enrichment_plots(df: pd.DataFrame) -> tuple:
    """
    绘制富集结果图：柱状图（每个 source 取前 N 条，-log10 校正 p）+ 气泡图（全局前 N 条）。

    Returns
    -------
    tuple
        (柱状图 Figure, 气泡图 Figure)；df 无有效行时返回 (None, None)。
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import config

    df = df[df["adjusted_p_value"].notna() & (df["adjusted_p_value"] > 0)].copy()
    if df.empty:
        return None, None
    df["neg_log10_p"] = -np.log10(df["adjusted_p_value"])

    source_colors = {
        "GO:BP": "#4C72B0", "GO:MF": "#DD8452",
        "GO:CC": "#55A868", "KEGG": "#C44E52",
    }

    # ---- 柱状图：每个 source 取校正 p 最小的前 N 条 ----
    top = (
        df.sort_values("adjusted_p_value")
        .groupby("source", as_index=False)
        .head(config.ENRICHMENT_TOP_PER_SOURCE)
    )
    fig_bar, ax_bar = plt.subplots(figsize=(9, max(5, 0.28 * len(top))))
    colors = [source_colors.get(s, "#7f7f7f") for s in top["source"]]
    ax_bar.barh(top["term_name"], top["neg_log10_p"], color=colors, height=0.6)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("-log10(adjusted p-value)")
    ax_bar.set_title("Functional Enrichment (top terms per source)")
    ax_bar.tick_params(axis="y", labelsize=7)
    fig_bar.tight_layout()

    # ---- 气泡图：全局校正 p 最小的前 N 条 ----
    bubble = df.sort_values("adjusted_p_value").head(config.PATHWAY_DISPLAY_NUM)
    sizes = 40 + 120 * (bubble["intersection_size"] / max(bubble["intersection_size"].max(), 1))
    fig_dot, ax_dot = plt.subplots(figsize=(10, max(6, 0.3 * len(bubble))))
    sc = ax_dot.scatter(
        bubble["neg_log10_p"], bubble["term_name"],
        s=sizes, c=bubble["adjusted_p_value"], cmap="viridis_r",
        edgecolors="grey", alpha=0.85,
    )
    cb = fig_dot.colorbar(sc, ax=ax_dot)
    cb.set_label("adjusted p-value")
    ax_dot.set_xlabel("-log10(adjusted p-value)")
    ax_dot.set_title("Enrichment bubble plot")
    ax_dot.tick_params(axis="y", labelsize=7)
    fig_dot.tight_layout()

    return fig_bar, fig_dot


# ---------- 排序富集（GSEA prerank，仅智人） ----------
# 与 ORA 的区别：输入是"基因|分数"排序列表（如 pChEMBL 亲和力、Open Targets 关联分），
# 而非二元基因列表，能利用打分信息检验"高分靶点整体偏向哪些通路"。
# 基因集来自 Enrichr 库的本地 GMT 快照（首次使用时下载一次，之后永不重复联网）；
# 失败语义与 ORA 链路一致：失败返回空表 + status，且绝不把失败结果写进缓存。


def _gsea_gmt_path(lib: str) -> str | None:
    """返回该 Enrichr 库的本地 GMT 快照路径；无快照时下载一次并写盘，失败返回 None。

    快照命名 {lib}_Human_{YYYY-MM-DD}.gmt，取字典序最新（即最近下载）的一份；
    快照存在则永不重新下载，删除 GSEA_GMT_DIR 目录即可强制重建。
    """
    import glob
    import config

    os.makedirs(config.GSEA_GMT_DIR, exist_ok=True)
    snaps = sorted(glob.glob(os.path.join(config.GSEA_GMT_DIR, f"{lib}_Human_*.gmt")))
    if snaps:
        return snaps[-1]

    import gseapy

    try:
        lib_dict = gseapy.get_library(lib, organism="Human")
    except Exception:
        logging.exception("Failed to download Enrichr library %s for GSEA", lib)
        return None
    if not lib_dict:
        logging.warning("Enrichr library %s downloaded empty", lib)
        return None
    from datetime import date

    path = os.path.join(config.GSEA_GMT_DIR, f"{lib}_Human_{date.today().isoformat()}.gmt")
    try:
        with open(path, "w", encoding="utf-8") as f:
            for term, gene_set in lib_dict.items():
                f.write("\t".join([str(term), "na", *[str(g) for g in gene_set]]) + "\n")
    except Exception:
        logging.exception("Failed to write GMT snapshot for %s", lib)
        return None
    logging.info("GSEA GMT snapshot saved: %s (%d gene sets)",
                 os.path.basename(path), len(lib_dict))
    return path


def gsea_gene_sets() -> dict[str, str]:
    """返回可用的 {库名: GMT路径}；单库缺失自动降级跳过，全部缺失时返回空 dict。"""
    import config

    out: dict[str, str] = {}
    for lib in config.GSEA_ENRICHR_LIBS:
        path = _gsea_gmt_path(lib)
        if path:
            out[lib] = path
        else:
            logging.warning("GSEA gene set library unavailable, skipped: %s", lib)
    return out


def build_ranked_list(targets_df: pd.DataFrame, gene_col: str, score_col: str) -> tuple[pd.DataFrame, dict]:
    """
    把带打分的靶点表转为 gseapy.prerank 的输入（gene | score，按分数降序）。

    规则：基因统一大写去空白；同基因多行取最高分；打分缺失的行剔除（不补零，
    补零会制造假尾部信号）。drop_ratio = 被剔除行数 / 输入行数，供页面预警。

    Returns
    -------
    tuple[pd.DataFrame, dict]
        (rnk: 两列 DataFrame[gene, score]，按 score 降序,
         stats: n_input / n_kept / n_dropped / drop_ratio / missing_cols)
    """
    if targets_df is None or targets_df.empty:
        return pd.DataFrame(columns=["gene", "score"]), {
            "n_input": 0, "n_kept": 0, "n_dropped": 0, "drop_ratio": 0.0,
            "missing_cols": [c for c in (gene_col, score_col)],
        }
    missing_cols = [c for c in (gene_col, score_col) if c not in targets_df.columns]
    if missing_cols:
        return pd.DataFrame(columns=["gene", "score"]), {
            "n_input": len(targets_df), "n_kept": 0, "n_dropped": len(targets_df),
            "drop_ratio": 1.0, "missing_cols": missing_cols,
        }

    df = targets_df[[gene_col, score_col]].rename(columns={gene_col: "gene", score_col: "score"})
    n_input = len(df)
    df["gene"] = df["gene"].astype(str).str.strip().str.upper()
    df = df[df["gene"] != ""]
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])
    df = df.groupby("gene", as_index=False)["score"].max()  # 同基因取最高分
    df = df.sort_values("score", ascending=False).reset_index(drop=True)
    n_kept = len(df)
    stats = {
        "n_input": n_input,
        "n_kept": n_kept,
        "n_dropped": n_input - n_kept,
        "drop_ratio": (n_input - n_kept) / n_input if n_input else 0.0,
        "missing_cols": [],
    }
    logging.info("Ranked list built: %d inputs -> %d scored genes (drop_ratio=%.2f)",
                 n_input, n_kept, stats["drop_ratio"])
    return df, stats


def _normalize_enrichment_gsea(res2d: pd.DataFrame) -> pd.DataFrame:
    """
    把 gseapy prerank 的 res2d 整理为统一列。

    gseapy 1.x 列名：Term / ES / NES / NOM p-val / FDR q-val / Tag % / Gene % / Lead_genes；
    这里按"大小写不敏感精确匹配 → 包含子串"两级查找，兼容后续版本的更名。
    """
    if res2d is None or res2d.empty:
        return pd.DataFrame()
    cols = list(res2d.columns)

    def pick(*cands: str) -> str | None:
        for c in cands:  # 先精确（大小写不敏感）
            for col in cols:
                if col.lower() == c.lower():
                    return col
        for c in cands:  # 再子串
            for col in cols:
                if c.lower() in col.lower():
                    return col
        return None

    name_of = {
        "term_name": pick("Term"),
        "es": pick("ES"),
        "nes": pick("NES"),
        "nominal_p_value": pick("NOM p-val", "NOM p"),
        "fdr_q_value": pick("FDR q-val", "FDR q", "FDR"),
        "tag_pct": pick("Tag %"),
        "gene_pct": pick("Gene %"),
        "lead_genes": pick("Lead_genes", "Lead genes"),
    }
    if name_of["term_name"] is None or name_of["nes"] is None:
        logging.warning("GSEA res2d missing expected columns: %s", cols)
        return pd.DataFrame()

    out = pd.DataFrame({
        "term_name": res2d[name_of["term_name"]].astype(str),
        "es": pd.to_numeric(res2d[name_of["es"]], errors="coerce") if name_of["es"] else pd.NA,
        "nes": pd.to_numeric(res2d[name_of["nes"]], errors="coerce"),
        "nominal_p_value": (pd.to_numeric(res2d[name_of["nominal_p_value"]], errors="coerce")
                            if name_of["nominal_p_value"] else pd.NA),
        "fdr_q_value": (pd.to_numeric(res2d[name_of["fdr_q_value"]], errors="coerce")
                        if name_of["fdr_q_value"] else pd.NA),
        "tag_pct": res2d[name_of["tag_pct"]] if name_of["tag_pct"] else "",
        "gene_pct": res2d[name_of["gene_pct"]] if name_of["gene_pct"] else "",
        "lead_genes": res2d[name_of["lead_genes"]].fillna("") if name_of["lead_genes"] else "",
        "source": res2d["source"],
    })
    out = out.dropna(subset=["term_name", "nes"])
    # 统一按 |NES| 降序（页面再按 FDR 筛选展示）
    out = out.reindex(out["nes"].abs().sort_values(ascending=False).index)
    return out.reset_index(drop=True)


def _gsea_cache_file(rnk: pd.DataFrame, libs: dict[str, str], label: str | None) -> str:
    """GSEA 结果缓存路径：键 = 排序基因+分数 + 库快照文件名 + 全部计算参数。"""
    import hashlib
    import config

    gene_part = "\r".join(rnk["gene"].astype(str))
    score_part = "\r".join(f"{v:.6g}" for v in rnk["score"].astype(float))
    libs_part = "|".join(f"{k}:{os.path.basename(v)}" for k, v in sorted(libs.items()))
    payload = (gene_part + "||" + score_part + "||" + libs_part
               + f"|{config.GSEA_SET_MIN_SIZE}|{config.GSEA_SET_MAX_SIZE}"
               + f"|{config.GSEA_PERMUTATION_NUM}|{config.GSEA_SEED}")
    md5 = hashlib.md5(payload.encode("utf-8")).hexdigest()
    name = f"GSEA_{label}_{md5}.csv" if label else f"GSEA_{md5}.csv"
    return os.path.join(config.ENRICHMENT_CACHE_DIR, name)


def fetch_enrichment_gsea_prerank(
    rnk: pd.DataFrame,
    libs: dict[str, str] | None = None,
    label: str | None = None,
    status: dict | None = None,
) -> pd.DataFrame:
    """
    GSEA prerank 主入口（gseapy 本地计算，逐库运行以标记 source）。

    Parameters
    ----------
    rnk : pd.DataFrame
        build_ranked_list 的输出（gene | score 两列，降序）。
    libs : dict[str, str] | None
        {库名: GMT路径}；None 时自动取 gsea_gene_sets()。
    label : str | None
        缓存文件名前缀（如 aspirin_cancer_drug）。
    status : dict | None
        可选出参，语义与 fetch_enrichment_gprofiler 一致：
          success / no_significant / cache / request_failed（GMT 全缺或计算异常）/
          insufficient_input（打分基因少于 GSEA_MIN_RANKED_GENES）/ empty_input

    Returns
    -------
    pd.DataFrame
        统一列结果表（source / term_name / es / nes / nominal_p_value /
        fdr_q_value / tag_pct / gene_pct / lead_genes），按 |NES| 降序；
        失败或无显著结果时返回空 DataFrame（失败不写缓存）。
    """
    import config

    if libs is None:
        libs = gsea_gene_sets()
    if rnk is None or len(rnk) == 0:
        if status is not None:
            status.update(ok=False, reason="empty_input")
        return pd.DataFrame()
    if len(rnk) < config.GSEA_MIN_RANKED_GENES:
        logging.info("GSEA skipped: only %d scored genes (< %d)",
                     len(rnk), config.GSEA_MIN_RANKED_GENES)
        if status is not None:
            status.update(ok=False, reason="insufficient_input")
        return pd.DataFrame()
    if not libs:
        logging.warning("GSEA skipped: no gene set library available")
        if status is not None:
            status.update(ok=False, reason="request_failed")
        return pd.DataFrame()

    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    cache_file = _gsea_cache_file(rnk, libs, label)
    if _is_cache_fresh(cache_file):
        logging.info("Loading cached GSEA prerank result for %d genes", len(rnk))
        if status is not None:
            status.update(ok=True, reason="cache")
        return _read_cache_csv(cache_file)

    import gseapy

    frames: list[pd.DataFrame] = []
    try:
        for lib, gmt_path in libs.items():
            res = gseapy.prerank(
                rnk=rnk, gene_sets=gmt_path, outdir=None, verbose=False,
                min_size=config.GSEA_SET_MIN_SIZE, max_size=config.GSEA_SET_MAX_SIZE,
                permutation_num=config.GSEA_PERMUTATION_NUM, seed=config.GSEA_SEED,
            )
            r2d = res.res2d
            if r2d is not None and len(r2d) > 0:
                frames.append(r2d.assign(source=lib))
    except Exception:
        logging.exception("GSEA prerank failed (%d genes, libs=%s)", len(rnk), list(libs))
        if status is not None:
            status.update(ok=False, reason="request_failed")
        return pd.DataFrame()

    if not frames:
        # 各库都算完但没有（通过规模窗口的）显著基因集
        logging.info("GSEA prerank: no gene sets passed for %d genes", len(rnk))
        if status is not None:
            status.update(ok=True, reason="no_significant")
        return pd.DataFrame()

    raw = pd.concat(frames, ignore_index=True)
    df = _normalize_enrichment_gsea(raw)
    if df.empty:
        # 防御：归一化清零不写缓存（避免异常空表被 TTL 命中）
        logging.warning("GSEA result normalized to empty — skip caching (%d genes)", len(rnk))
        if status is not None:
            status.update(ok=False, reason="normalize_empty")
        return df
    if status is not None:
        status.update(ok=True, reason="success")
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    logging.info("GSEA prerank for %d genes: %d terms cached", len(rnk), len(df))
    return df


def render_gsea_plot(df: pd.DataFrame) -> "matplotlib.figure.Figure | None":
    """
    绘制 GSEA 结果条形图：按 |NES| 取前 N 条基因集，
    NES > 0 红色（高分靶点中富集/激活）、NES < 0 蓝色（抑制）。

    Returns
    -------
    matplotlib.figure.Figure | None
        df 无有效行时返回 None。
    """
    import matplotlib.pyplot as plt
    import config

    df = df[df["nes"].notna()].copy()
    if df.empty:
        return None
    top = df.reindex(df["nes"].abs().sort_values(ascending=False).index)
    top = top.head(config.GSEA_TOP_DISPLAY)
    colors = ["#C44E52" if v > 0 else "#4C72B0" for v in top["nes"]]
    fig, ax = plt.subplots(figsize=(9, max(4, 0.3 * len(top))))
    ax.barh(top["term_name"], top["nes"], color=colors, height=0.6)
    ax.axvline(0, color="#888888", linewidth=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("NES (normalized enrichment score)")
    ax.set_title("GSEA prerank (red = enriched in high-score targets)")
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    return fig
