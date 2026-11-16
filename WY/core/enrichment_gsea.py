"""排序富集(GSEA prerank,仅智人)。

输入是"基因|分数"排序列表(如药物亲和力、疾病关联分),检验高分靶点整体偏向哪些通路。
基因集来自 Enrichr 库的本地 GMT 文件(首次使用联网下载一次,之后不再重复下载);
失败时返回空表,且失败结果不会写进缓存。
"""

import hashlib
import logging
import os

import pandas as pd

from .enrichment_common import (_clean_cache_dir_limited, _hit_csv_cache,
                                _set_status, enrichment_cache_filename,
                                save_enrichment_result)


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
    from WY import config

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
    from WY import config

    out: dict[str, str] = {}
    for lib in config.GSEA_ENRICHR_LIBS:
        path = _gsea_gmt_path(lib)
        if path:
            out[lib] = path
        else:
            logging.warning("GSEA gene set library unavailable, skipped: %s", lib)
    return out


def build_ranked_list(targets_df: pd.DataFrame, gene_col: str, score_col: str) -> tuple[pd.DataFrame, dict]:
    """把带打分的靶点表转为排序富集的输入（基因 | 分数，按分数降序）。
    基因名统一大写、同一基因取最高分;没有分数的行直接剔除不补零（补零会制造假信号）,剔除比例记入 stats。
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
    df = df.groupby("gene", as_index=False)["score"].max()
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


def _gsea_cache_md5(rnk: pd.DataFrame, libs: dict[str, str]) -> str:
    """GSEA 缓存键:由基因排序与分数、基因集库、计算参数共同决定,相同输入命中同一缓存。"""
    from WY import config

    gene_part = "\r".join(rnk["gene"].astype(str))
    score_part = "\r".join(f"{v:.6g}" for v in rnk["score"].astype(float))
    libs_part = "|".join(f"{k}:{os.path.basename(v)}" for k, v in sorted(libs.items()))
    payload = (gene_part + "||" + score_part + "||" + libs_part
               + f"|{config.GSEA_SET_MIN_SIZE}|{config.GSEA_SET_MAX_SIZE}"
               + f"|{config.GSEA_PERMUTATION_NUM}|{config.GSEA_SEED}")
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def _gsea_cache_file(rnk: pd.DataFrame, libs: dict[str, str], label: str | None,
                     md5: str | None = None) -> str:
    """GSEA 结果缓存路径:同名覆盖式命名 GSEA_侧_标识.csv(键记录在文件头部,见 _hit_csv_cache)。
    label 约定"侧别_分析标识",侧别放在最前,避免与疾病名里可能含的"disease"字样混淆。
    """
    from WY import config

    if md5 is None:
        md5 = _gsea_cache_md5(rnk, libs)
    return os.path.join(config.ENRICHMENT_CACHE_DIR,
                        enrichment_cache_filename("GSEA", label, md5))


def fetch_enrichment_gsea_prerank(
    rnk: pd.DataFrame,
    libs: dict[str, str] | None = None,
    label: str | None = None,
    status: dict | None = None,
) -> pd.DataFrame:
    """排序富集（GSEA prerank）主入口：在本地计算,逐个基因集库运行后合并,结果按 |NES| 降序。
    输入为 build_ranked_list 输出的"基因|分数"表;失败或无显著结果返回空表（失败不写缓存）,status 含义与 g:Profiler 版本一致。
    """
    from WY import config

    if libs is None:
        libs = gsea_gene_sets()
    if rnk is None or len(rnk) == 0:
        _set_status(status, False, "empty_input")
        return pd.DataFrame()
    if len(rnk) < config.GSEA_MIN_RANKED_GENES:
        logging.info("GSEA skipped: only %d scored genes (< %d)",
                     len(rnk), config.GSEA_MIN_RANKED_GENES)
        _set_status(status, False, "insufficient_input")
        return pd.DataFrame()
    if not libs:
        logging.warning("GSEA skipped: no gene set library available")
        _set_status(status, False, "request_failed")
        return pd.DataFrame()

    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    _clean_cache_dir_limited()
    gsea_md5 = _gsea_cache_md5(rnk, libs)
    cache_file = _gsea_cache_file(rnk, libs, label, md5=gsea_md5)
    cached = _hit_csv_cache(cache_file, gsea_md5, f"GSEA prerank result for {len(rnk)} genes")
    if cached is not None:
        _set_status(status, True, "cache")
        return cached

    import gseapy

    frames: list[pd.DataFrame] = []
    n_error = 0
    for lib, gmt_path in libs.items():
        try:
            res = gseapy.prerank(
                rnk=rnk, gene_sets=gmt_path, outdir=None, verbose=False,
                min_size=config.GSEA_SET_MIN_SIZE, max_size=config.GSEA_SET_MAX_SIZE,
                permutation_num=config.GSEA_PERMUTATION_NUM, seed=config.GSEA_SEED,
            )
        except LookupError:
            # 该库与本次基因列表交集过小，没有基因集通过 min/max 过滤：
            # 属数据特征而非服务失败，跳过该库，让其余库继续计算
            logging.info("GSEA prerank: no gene sets in %s matched %d genes, skipped",
                         lib, len(rnk))
            continue
        except Exception:
            logging.exception("GSEA prerank failed for %s (%d genes)", lib, len(rnk))
            n_error += 1
            continue
        r2d = res.res2d
        if r2d is not None and len(r2d) > 0:
            frames.append(r2d.assign(source=lib))

    if not frames:
        if n_error:
            _set_status(status, False, "request_failed")
            return pd.DataFrame()
        # 各基因集库都没有可用的基因集或没有显著富集的通路——真"无显著结果"，非失败
        logging.info("GSEA prerank: no gene sets passed for %d genes", len(rnk))
        _set_status(status, True, "no_significant")
        return pd.DataFrame()

    raw = pd.concat(frames, ignore_index=True)
    df = _normalize_enrichment_gsea(raw)
    return save_enrichment_result(df, cache_file, gsea_md5, status,
                                  log_prefix="GSEA prerank", n_genes=len(rnk))


def render_gsea_plot(df: pd.DataFrame) -> "matplotlib.figure.Figure | None":
    """绘制 GSEA 结果条形图：按 |NES| 取前 N 条通路,NES > 0 红色（在高分靶点中富集/激活）、NES < 0 蓝色（抑制）。
    没有有效结果时返回 None。
    """
    import matplotlib.pyplot as plt
    from WY import config

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
