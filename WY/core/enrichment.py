"""GO/KEGG 功能富集分析：主用 g:Profiler，连不上时自动改用 Enrichr，并绘制结果图。"""

import hashlib
import logging
import os

import pandas as pd

from .caching import _is_cache_fresh, _read_cache_csv

# 两个数据源的基因集库名对应关系（GO 生物过程/分子功能/细胞组分、KEGG 通路）。
ENRICHR_LIB_MAP = {
    "GO:BP": "GO_Biological_Process_2021",
    "GO:MF": "GO_Molecular_Function_2021",
    "GO:CC": "GO_Cellular_Component_2021",
    "KEGG": "KEGG_2021_Human",
}

# g:Profiler 物种名 → Enrichr 物种名。Enrichr 只支持人/小鼠等，大鼠没有对应库，
# 因此大鼠场景下 Enrichr 降级会被跳过，只能用 g:Profiler。
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


# ---------- g:Profiler 熔断 ----------
# g:Profiler 暂时连不上时，短期内自动改用备选数据源 Enrichr，避免反复苦等超时；
# 冷却期过后自动恢复尝试，服务恢复即回到正常流程。


def _gprofiler_circuit_file() -> str:
    """熔断标记文件路径（记录上次失败时间）。"""
    import config

    return os.path.join(config.ENRICHMENT_CACHE_DIR, ".gprofiler_unreachable")


def _gprofiler_circuit_open() -> bool:
    """熔断是否生效：上次失败距今仍在冷却期内则为 True。"""
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
    """记录本次失败时间，开启熔断。"""
    import time

    path = _gprofiler_circuit_file()
    try:
        with open(path, "w", encoding="utf-8") as f:
            f.write(str(time.time()))
    except Exception:
        logging.exception("Failed to write g:Profiler circuit breaker file")


def _gprofiler_circuit_close() -> None:
    """请求成功时清除熔断标记。"""
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
    用 g:Profiler 查询基因列表的功能富集（GO/KEGG）。结果整理为统一列后缓存到本地，
    相同基因列表与参数再次查询时直接读缓存。

    请求失败或没有显著条目时返回空表；status（若提供）标明本次属于哪种情况，
    便于区分"真无显著项"与"服务失败"。
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

    # 熔断：g:Profiler 近期失败过则直接跳过，走 Enrichr 兜底
    if _gprofiler_circuit_open():
        logging.info("g:Profiler circuit breaker open — skipping, Enrichr fallback")
        if status is not None:
            status.update(ok=False, reason="circuit_open")
        return pd.DataFrame()

    import time
    import socket

    # 请求保护：查询库自身不带超时和重试，这里统一补上请求间隔、失败重试和超时限制。
    from gprofiler import GProfiler

    res = None
    for attempt in range(config.ENRICHMENT_MAX_RETRIES):
        time.sleep(config.API_REQUEST_DELAY)  # 每次请求前稍等，避免请求过快被服务器拒绝
        socket.setdefaulttimeout(config.API_TIMEOUT)
        try:
            gp = GProfiler(return_dataframe=True)
            res = gp.profile(
                query=genes,
                organism=organism,
                sources=list(sources),
                user_threshold=pvalue,
                no_evidences=False,  # 需要返回每个通路命中的基因列表
            )
            break  # 成功
        except Exception:
            logging.exception(
                "g:Profiler enrichment failed (attempt %d/%d)",
                attempt + 1, config.ENRICHMENT_MAX_RETRIES,
            )
            _gprofiler_circuit_open_trip()  # 记录失败时间，开启熔断
            if attempt < config.ENRICHMENT_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)  # 失败后稍等再重试
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
        # 结果整理后为空时不写缓存，避免下次把异常空表当成有效结果
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
    """把 g:Profiler 返回的结果整理为统一列。

    g:Profiler 返回的 p 值已做多重检验校正，因此校正前后 p 值相同；
    同时兼容不同版本返回的列名差异。
    """
    df = df.reset_index()  # 通路 ID 有时在行索引里，先还原成列

    def _pick(*names) -> pd.Series:
        """按优先级取第一个存在的列，都没有时返回空值。"""
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
    用 Enrichr 查询功能富集，作为 g:Profiler 连不上时的备选方案。
    结果同样缓存到本地，相同基因列表与参数直接读缓存；
    请求失败或没有显著条目时返回空表，status 含义与 g:Profiler 版本一致。
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

    # 请求保护：查询库自身的超时和重试不完整，这里统一补上请求间隔、失败重试和超时限制。
    import gseapy as gp

    res = None
    for attempt in range(config.ENRICHMENT_MAX_RETRIES):
        time.sleep(config.API_REQUEST_DELAY)  # 每次请求前稍等，避免请求过快被服务器拒绝
        socket.setdefaulttimeout(config.API_TIMEOUT)
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=list(libs),
                organism=organism,
                outdir=None,  # 不让库自行写盘，结果由本程序统一缓存
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
                time.sleep(config.API_RETRY_DELAY)  # 失败后稍等再重试
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
        # 结果整理后为空时不写缓存，避免下次把异常空表当成有效结果
        logging.warning("Enrichr result normalized to empty — skip caching (%d genes)", len(genes))
        if status is not None:
            status.update(ok=False, reason="normalize_empty")
        return df
    if status is not None:
        status.update(ok=True, reason="success")
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    return df


def _normalize_enrichment_enrichr(df: pd.DataFrame) -> pd.DataFrame:
    """把 Enrichr 返回的结果整理为统一列，库名换回与 g:Profiler 一致的写法。"""
    import re

    def _parse_overlap(s) -> tuple[int, int]:
        try:
            a, b = str(s).split("/")
            return int(a), int(b)
        except (ValueError, TypeError):
            return 0, 0

    def _extract_term_id(term) -> str:
        # 通路名形如 "XXX (GO:0018212)"，提取末尾括号内的通路 ID
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
    富集分析主入口：输入基因列表（通常是药物与疾病的交集靶点），
    先请求 g:Profiler，失败自动改用 Enrichr。

    返回 (富集结果表, 实际使用的数据源)；两个数据源都失败时返回 (空表, None)。
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
    绘制富集结果图：柱状图（每类数据源各取前 N 条通路，横轴为 -log10 校正 p 值）
    和气泡图（全局前 N 条，点越大表示命中的基因越多）。
    没有有效结果时返回 (None, None)。
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
# 与普通富集的区别：输入是"基因|分数"排序列表（如药物亲和力、疾病关联分），
# 而非单纯的基因列表，用于检验高分靶点整体偏向哪些通路。
# 基因集来自 Enrichr 库的本地文件（首次使用时联网下载一次，之后不再重复下载）；
# 失败时返回空表，且失败结果不会写进缓存。


def _gsea_gmt_path(lib: str) -> str | None:
    """返回 Enrichr 库对应的本地基因集文件（GMT）路径；本地没有时联网下载一次，
    之后一直复用，下载失败返回 None。删除 GSEA_GMT_DIR 目录可强制重新下载。
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
    """返回可用的基因集库（库名 → 文件路径）；单个库缺失自动跳过，全部缺失返回空字典。"""
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
    把带打分的靶点表转为排序富集的输入（基因 | 分数，按分数降序）。

    处理规则：基因名统一大写；同一基因取最高分；没有分数的行直接剔除
    （不补零，补零会制造假信号），剔除比例记入 stats 供页面提示。
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
    """把排序富集的原始结果整理为统一列；兼容不同版本软件的列名差异。"""
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
    """GSEA 结果缓存路径：由基因排序与分数、基因集库、计算参数共同决定，相同输入命中同一缓存。"""
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
    排序富集（GSEA prerank）主入口：在本地计算，逐个基因集库运行后合并结果。

    输入为"基因|分数"表（build_ranked_list 的输出）；结果按 |NES| 降序，
    失败或无显著结果时返回空表（失败不写缓存），status 含义与 g:Profiler 版本一致。
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
        # 各基因集库都算完了，但没有显著富集的通路
        logging.info("GSEA prerank: no gene sets passed for %d genes", len(rnk))
        if status is not None:
            status.update(ok=True, reason="no_significant")
        return pd.DataFrame()

    raw = pd.concat(frames, ignore_index=True)
    df = _normalize_enrichment_gsea(raw)
    if df.empty:
        # 结果整理后为空时不写缓存
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
    绘制 GSEA 结果条形图：按 |NES| 取前 N 条通路，
    NES > 0 红色（在高分靶点中富集/激活），NES < 0 蓝色（抑制）；
    没有有效结果时返回 None。
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
