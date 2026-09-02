"""靶点汇总：整合药物靶点（ChEMBL+DGIdb）与疾病靶点（Open Targets+DisGeNET），并取交集。"""

import logging
import os

import pandas as pd

from .caching import open_csv_cache, save_csv_cache
from .chembl import fetch_drug_targets_chembl
from .disgenet import fetch_disease_targets_disgenet
from .dgidb import fetch_drug_targets_dgidb
from .opentargets import get_disease_targets
from .text_utils import safe_name


def _max_score_map(df: pd.DataFrame, score_col: str) -> dict[str, float]:
    """按基因符号分组取分数最大值；符号统一大写，分数统一转数值（非数值置空）。"""
    if df.empty or score_col not in df.columns:
        return {}
    tmp = df.dropna(subset=["approvedSymbol"]).copy()
    tmp["approvedSymbol"] = tmp["approvedSymbol"].str.upper()
    tmp[score_col] = pd.to_numeric(tmp[score_col], errors="coerce")
    return tmp.groupby("approvedSymbol")[score_col].max().to_dict()


def _merge_source_pair(
    name: str,
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    *,
    cache_path: str,
    kind: str,
    src_a: str,
    src_b: str,
    score_a: tuple[str, str],
    score_b: tuple[str, str],
    sources_desc: tuple[str, str],
) -> pd.DataFrame:
    """两个数据源的靶点按基因符号取并集去重，并保留各源分数最大值。

    kind 为 "drug"/"disease"，用于日志；src_a/src_b 为来源标记列名；
    score_a/score_b 为 (输出列名, 源分数列名)；两源都没有数据时返回空表。
    """
    genes_a = set(df_a["approvedSymbol"].dropna().str.upper()) if not df_a.empty else set()
    genes_b = set(df_b["approvedSymbol"].dropna().str.upper()) if not df_b.empty else set()
    all_genes = sorted(genes_a | genes_b)
    if not all_genes:
        logging.warning("No %s targets from either source for '%s'", kind, name)
        return pd.DataFrame()

    col_a, raw_a = score_a
    col_b, raw_b = score_b
    map_a = _max_score_map(df_a, raw_a)
    map_b = _max_score_map(df_b, raw_b)

    records: list[dict] = []
    for g in all_genes:
        in_a = g in genes_a
        in_b = g in genes_b
        records.append({
            "approvedSymbol": g,
            src_a: int(in_a),
            src_b: int(in_b),
            "source_count": int(in_a) + int(in_b),
            col_a: map_a.get(g),
            col_b: map_b.get(g),
        })

    merged = pd.DataFrame(records)
    save_csv_cache(merged, cache_path, skip_empty=False)
    logging.info(
        "Merged %s '%s': %d targets (%s %d, %s %d)",
        kind, name, len(merged), sources_desc[0], len(genes_a), sources_desc[1], len(genes_b),
    )
    return merged


def get_drug_targets_merged(drug_name: str) -> pd.DataFrame:
    """整合 ChEMBL 与 DGIdb 两个数据库的药物靶点，按基因符号取并集去重。

    为每个基因标记来自哪些数据库（source_chmbl / source_dgidb），
    并保留各库的活性/交互分数。两个库都没有数据时返回空表；
    结果存为本地 CSV，已有缓存则直接复用。
    """
    import config

    cache_path = os.path.join(config.MERGED_CACHE_DIR, f"merged_drug_{safe_name(drug_name)}.csv")
    cached = open_csv_cache(cache_path, label=f"merged drug targets for '{drug_name}'")
    if cached is not None:
        return cached

    chmbl_df, _ = fetch_drug_targets_chembl(drug_name)
    dgidb_df = fetch_drug_targets_dgidb(drug_name)
    return _merge_source_pair(
        drug_name, chmbl_df, dgidb_df,
        cache_path=cache_path, kind="drug",
        src_a="source_chmbl", src_b="source_dgidb",
        score_a=("score_pchembl", "pchembl_value"),
        score_b=("score_dgidb", "interaction_score"),
        sources_desc=("ChEMBL", "DGIdb"),
    )


def get_disease_targets_merged(disease_name: str) -> pd.DataFrame:
    """整合 Open Targets 与 DisGeNET 两个数据库的疾病靶点，按基因符号取并集去重。

    为每个基因标记来自哪些数据库（source_ot / source_disgenet），
    并保留各库的关联分数。两个库都没有数据时返回空表；
    结果存为本地 CSV，已有缓存则直接复用。
    """
    import config

    cache_path = os.path.join(config.MERGED_CACHE_DIR, f"merged_disease_{safe_name(disease_name)}.csv")
    cached = open_csv_cache(cache_path, label=f"merged disease targets for '{disease_name}'")
    if cached is not None:
        return cached

    ot_df = get_disease_targets(disease_name)
    disgenet_df = fetch_disease_targets_disgenet(disease_name)
    return _merge_source_pair(
        disease_name, ot_df, disgenet_df,
        cache_path=cache_path, kind="disease",
        src_a="source_ot", src_b="source_disgenet",
        score_a=("score_ot", "association_score"),
        score_b=("score_disgenet", "disgenet_score"),
        sources_desc=("OpenTargets", "DisGeNET"),
    )


def get_common_targets(drug_name: str, disease_name: str, top_n: int = 1000) -> list[str]:
    """
    取药物靶点与疾病靶点的交集基因，作为后续 PPI 网络分析的输入。

    按字母序排序，最多返回 top_n 个（默认 1000，对齐 STRING 查询上限）；
    任一侧无数据时返回空列表。
    """
    drug_df = get_drug_targets_merged(drug_name)
    disease_df = get_disease_targets_merged(disease_name)
    # 提取药物和疾病两侧的基因符号
    drug_set = set(drug_df["approvedSymbol"].dropna()) if not drug_df.empty else set()
    disease_set = set(disease_df["approvedSymbol"].dropna()) if not disease_df.empty else set()
    # 取交集，按字母序排序并截断到 top_n
    return sorted(drug_set & disease_set)[:top_n]
