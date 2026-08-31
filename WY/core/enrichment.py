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
        return pd.DataFrame()
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    md5 = _enrichment_cache_key(genes, organism, sources, pvalue)
    name = f"gProfiler_{label}_{md5}.csv" if label else f"gProfiler_{md5}.csv"
    cache_file = os.path.join(config.ENRICHMENT_CACHE_DIR, name)
    if _is_cache_fresh(cache_file):
        logging.info("Loading cached g:Profiler enrichment for %d genes", len(genes))
        return _read_cache_csv(cache_file)

    # 熔断：gProfiler 近期失败过则直接跳过，不再空等连接超时（走 Enrichr 兜底）
    if _gprofiler_circuit_open():
        logging.info("g:Profiler circuit breaker open — skipping, Enrichr fallback")
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

    if res is None or len(res) == 0:
        logging.info("g:Profiler: no significant terms (p < %s) for %d genes", pvalue, len(genes))
        return pd.DataFrame()

    _gprofiler_circuit_close()  # 请求成功，关闭熔断
    df = _normalize_enrichment_gprofiler(res)
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

    Returns
    -------
    pd.DataFrame
        统一列富集结果；请求失败或无显著结果时返回空 DataFrame。
    """
    import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes or not libs:
        return pd.DataFrame()
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    md5 = _enrichment_cache_key(genes, organism, libs, pvalue)
    name = f"Enrichr_{label}_{md5}.csv" if label else f"Enrichr_{md5}.csv"
    cache_file = os.path.join(config.ENRICHMENT_CACHE_DIR, name)
    if _is_cache_fresh(cache_file):
        logging.info("Loading cached Enrichr enrichment for %d genes", len(genes))
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

    if res is None or len(res) == 0:
        return pd.DataFrame()

    df = _normalize_enrichment_enrichr(res)
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
