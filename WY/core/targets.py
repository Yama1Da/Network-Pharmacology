"""多来源靶点整合：药物（ChEMBL+DGIdb）、疾病（Open Targets+DisGeNET）并集与交集。"""

import logging
import os

import pandas as pd

from .chembl import fetch_drug_targets_chembl
from .disgenet import fetch_disease_targets_disgenet
from .dgidb import fetch_drug_targets_dgidb
from .opentargets import get_disease_targets
from .text_utils import safe_name


def get_drug_targets_merged(drug_name: str) -> pd.DataFrame:
    """整合 ChEMBL 与 DGIdb 两个来源的药物靶点（并集去重，保留来源证据）。

    两个来源分别抓取（各自独立缓存），按基因符号（approvedSymbol，统一大写）取并集，
    去重后为每个基因标记来源（source_chmbl / source_dgidb）与来源数量，
    并保留各来源的分数（score_pchembl / score_dgidb）。

    Parameters
    ----------
    drug_name : str
        药物英文名。

    Returns
    -------
    pd.DataFrame
        列: approvedSymbol、source_chmbl、source_dgidb、source_count、
        score_pchembl、score_dgidb；两来源均无数据时返回空表。

    缓存
    ----
    结果持久化为 CSV（Static/data/cache/merged/merged_drug_{drug_name}.csv）。
    """
    import config

    os.makedirs(config.MERGED_CACHE_DIR, exist_ok=True)
    merged_path = os.path.join(config.MERGED_CACHE_DIR, f"merged_drug_{safe_name(drug_name)}.csv")
    if os.path.exists(merged_path):
        logging.info("Loading cached merged drug targets for '%s'", drug_name)
        return pd.read_csv(merged_path)

    chmbl_df, _ = fetch_drug_targets_chembl(drug_name)
    dgidb_df = fetch_drug_targets_dgidb(drug_name)

    chmbl_genes = set(chmbl_df["approvedSymbol"].dropna().str.upper()) if not chmbl_df.empty else set()
    dgidb_genes = set(dgidb_df["approvedSymbol"].dropna().str.upper()) if not dgidb_df.empty else set()

    all_genes = sorted(chmbl_genes | dgidb_genes)
    if not all_genes:
        logging.warning("No drug targets from either source for '%s'", drug_name)
        return pd.DataFrame()

    # 各来源分数聚合（同一基因多条记录取最大值）
    chmbl_score: dict[str, float] = {}
    if not chmbl_df.empty and "pchembl_value" in chmbl_df.columns:
        tmp = chmbl_df.dropna(subset=["approvedSymbol"]).copy()
        tmp["approvedSymbol"] = tmp["approvedSymbol"].str.upper()
        tmp["pchembl_value"] = pd.to_numeric(tmp["pchembl_value"], errors="coerce")
        chmbl_score = tmp.groupby("approvedSymbol")["pchembl_value"].max().to_dict()

    dgidb_score: dict[str, float] = {}
    if not dgidb_df.empty and "interaction_score" in dgidb_df.columns:
        tmp = dgidb_df.dropna(subset=["approvedSymbol"]).copy()
        tmp["approvedSymbol"] = tmp["approvedSymbol"].str.upper()
        tmp["interaction_score"] = pd.to_numeric(tmp["interaction_score"], errors="coerce")
        dgidb_score = tmp.groupby("approvedSymbol")["interaction_score"].max().to_dict()

    records: list[dict] = []
    for g in all_genes:
        in_chmbl = g in chmbl_genes
        in_dgidb = g in dgidb_genes
        records.append({
            "approvedSymbol": g,
            "source_chmbl": int(in_chmbl),
            "source_dgidb": int(in_dgidb),
            "source_count": int(in_chmbl) + int(in_dgidb),
            "score_pchembl": chmbl_score.get(g),
            "score_dgidb": dgidb_score.get(g),
        })

    merged = pd.DataFrame(records)
    merged.to_csv(merged_path, index=False, encoding="utf-8-sig")
    logging.info(
        "Merged drug '%s': %d targets (ChEMBL %d, DGIdb %d)",
        drug_name, len(merged), len(chmbl_genes), len(dgidb_genes),
    )
    return merged


def get_disease_targets_merged(disease_name: str) -> pd.DataFrame:
    """整合 Open Targets 与 DisGeNET 两个来源的疾病靶点（并集去重，保留来源证据）。

    两个来源分别抓取（各自独立缓存），按基因符号（approvedSymbol，统一大写）取并集，
    去重后为每个基因标记来源（source_ot / source_disgenet）与来源数量，
    并保留各来源的分数（score_ot / score_disgenet）。

    Parameters
    ----------
    disease_name : str
        疾病英文名。

    Returns
    -------
    pd.DataFrame
        列: approvedSymbol、source_ot、source_disgenet、source_count、
        score_ot、score_disgenet；两来源均无数据时返回空表。

    缓存
    ----
    结果持久化为 CSV（Static/data/cache/merged/merged_disease_{disease_name}.csv）。
    """
    import config

    os.makedirs(config.MERGED_CACHE_DIR, exist_ok=True)
    merged_path = os.path.join(config.MERGED_CACHE_DIR, f"merged_disease_{safe_name(disease_name)}.csv")
    if os.path.exists(merged_path):
        logging.info("Loading cached merged disease targets for '%s'", disease_name)
        return pd.read_csv(merged_path)

    ot_df = get_disease_targets(disease_name)
    disgenet_df = fetch_disease_targets_disgenet(disease_name)

    ot_genes = set(ot_df["approvedSymbol"].dropna().str.upper()) if not ot_df.empty else set()
    disgenet_genes = set(disgenet_df["approvedSymbol"].dropna().str.upper()) if not disgenet_df.empty else set()

    all_genes = sorted(ot_genes | disgenet_genes)
    if not all_genes:
        logging.warning("No disease targets from either source for '%s'", disease_name)
        return pd.DataFrame()

    ot_score: dict[str, float] = {}
    if not ot_df.empty and "association_score" in ot_df.columns:
        tmp = ot_df.dropna(subset=["approvedSymbol"]).copy()
        tmp["approvedSymbol"] = tmp["approvedSymbol"].str.upper()
        ot_score = tmp.groupby("approvedSymbol")["association_score"].max().to_dict()

    disgenet_score: dict[str, float] = {}
    if not disgenet_df.empty and "disgenet_score" in disgenet_df.columns:
        tmp = disgenet_df.dropna(subset=["approvedSymbol"]).copy()
        tmp["approvedSymbol"] = tmp["approvedSymbol"].str.upper()
        tmp["disgenet_score"] = pd.to_numeric(tmp["disgenet_score"], errors="coerce")
        disgenet_score = tmp.groupby("approvedSymbol")["disgenet_score"].max().to_dict()

    records: list[dict] = []
    for g in all_genes:
        in_ot = g in ot_genes
        in_dis = g in disgenet_genes
        records.append({
            "approvedSymbol": g,
            "source_ot": int(in_ot),
            "source_disgenet": int(in_dis),
            "source_count": int(in_ot) + int(in_dis),
            "score_ot": ot_score.get(g),
            "score_disgenet": disgenet_score.get(g),
        })

    merged = pd.DataFrame(records)
    merged.to_csv(merged_path, index=False, encoding="utf-8-sig")
    logging.info(
        "Merged disease '%s': %d targets (OpenTargets %d, DisGeNET %d)",
        disease_name, len(merged), len(ot_genes), len(disgenet_genes),
    )
    return merged


def get_common_targets(drug_name: str, disease_name: str, top_n: int = 1000) -> list[str]:
    """
    计算药物靶点与疾病靶点的交集基因列表（PPI 网络输入）。

    药物靶点整合 ChEMBL + DGIdb（get_drug_targets_merged），
    疾病靶点整合 Open Targets + DisGeNET（get_disease_targets_merged），
    再取两表 approvedSymbol 的交集，按字母序排序并截断到 top_n（默认对齐 STRING 查询上限 1000）。

    Parameters
    ----------
    drug_name : str
        药物英文名。
    disease_name : str
        疾病英文名。
    top_n : int
        返回的交集基因数量上限，默认 1000。

    Returns
    -------
    list[str]
        交集基因符号列表；任一来源无数据时返回空列表。
    """
    drug_df = get_drug_targets_merged(drug_name)
    disease_df = get_disease_targets_merged(disease_name)
    # 从药物和疾病靶点表中提取基因符号（approvedSymbol），并转换为集合
    drug_set = set(drug_df["approvedSymbol"].dropna()) if not drug_df.empty else set()
    disease_set = set(disease_df["approvedSymbol"].dropna()) if not disease_df.empty else set()
    # 取药物和疾病靶点的交集，按字母序排序并截断到 top_n
    return sorted(drug_set & disease_set)[:top_n]
