"""ChEMBL 数据源：查询药物的全部人类活性靶点，并提供缓存有效期调整。"""

import logging
import os

import pandas as pd
import streamlit as st

from .caching import open_csv_cache, save_csv_cache
from .text_utils import safe_name


def _sort_drug_df_by_pchembl(df: pd.DataFrame) -> pd.DataFrame:
    """按 pChEMBL 活性值从高到低排列药物靶点（无活性值的排最后）。"""
    if df.empty or "pchembl_value" not in df.columns:
        return df
    df = df.copy()
    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    return (
        df.sort_values("pchembl_value", ascending=False, na_position="last")
        .reset_index(drop=True)
    )


def fetch_drug_targets_chembl(drug_name: str) -> tuple[pd.DataFrame, str | None]:
    """
    从 ChEMBL 获取药物的全部活性靶点（全靶点模式）：不按活性值过滤，
    靶点更全，但会包含 ALB、CYP450 等弱关联/非特异性结合的蛋白；
    仅保留人类靶点，按 pChEMBL 活性值降序排列。
    已有本地缓存时直接复用；请求失败时返回空表。
    """
    import config

    csv_path = os.path.join(config.DRUG_CACHE_DIR, f"ChEMBL_{safe_name(drug_name)}.csv")

    # 已有缓存则直接读取本地文件
    cached = open_csv_cache(csv_path, label=f"ChEMBL targets for '{drug_name}'")
    if cached is not None:
        return _sort_drug_df_by_pchembl(cached), csv_path

    # 无缓存时联网查询
    try:
        # 延迟导入：该库一导入就会联网，放在函数内可避免启动时
        # ChEMBL 服务器不可用导致整个程序崩溃。
        from chembl_webresource_client.new_client import new_client

        # 第一步：按药物名查找化合物
        mols = list(new_client.molecule.filter(pref_name__iexact=drug_name.strip()))
        if not mols:
            # 找不到时改用不区分大小写的模糊搜索
            mols = list(new_client.molecule.filter(pref_name__icontains=drug_name.strip()))
        if not mols:
            logging.warning("No ChEMBL molecule found for '%s'", drug_name)
            return pd.DataFrame(), None

        chembl_id = mols[0]["molecule_chembl_id"]

        # 第二步：取该药物在人类中的全部活性记录（不按 pChEMBL 过滤，靶点更全）。
        # 如需只保留有 pChEMBL 定量值的高置信度靶点，可在 filter 中
        # 加回 pchembl_value__isnull=False 参数。
        acts = list(new_client.activity.filter(
            molecule_chembl_id=chembl_id,
            organism="Homo sapiens",
        ))

        # 第三步：逐个获取靶点详情（同一靶点只查一次）
        import time
        import config
        target_cache: dict[str, dict | None] = {}
        for act in acts:
            tid = act.get("target_chembl_id")
            if tid and tid not in target_cache:
                # 每次请求前稍作等待，避免请求过快被 ChEMBL 限流
                time.sleep(config.API_REQUEST_DELAY)
                try:
                    t = new_client.target.get(tid)
                except Exception:
                    logging.exception("Failed to fetch target details for '%s'", tid)
                    st.toast("靶点获取错误，停止信息获取", icon="⚠️")
                    return pd.DataFrame(), None  # 失败即停止，且不保存缓存
                if not t:
                    logging.error("Target '%s' not found, aborting fetch", tid)
                    st.toast("靶点获取错误，停止信息获取", icon="⚠️")
                    return pd.DataFrame(), None  # 失败即停止，且不保存缓存
                target_cache[tid] = t

        def _extract_gene_name(t: dict) -> str | None:
            """从靶点信息中提取基因符号，没有则退而用 HGNC 名称，都没有时返回 None"""
            if "target_components" in t:
                for comp in t["target_components"]:
                    for syn in comp.get("target_component_synonyms", []):
                        if syn.get("syn_type") == "GENE_SYMBOL":
                            return syn.get("component_synonym")
                    for xref in comp.get("target_component_xrefs", []):
                        if xref.get("xref_src_db") == "HGNC":
                            return xref.get("xref_name")
            return None

        # 第四步：给每条活性记录补上靶点信息，整理成表格行
        rows: list[dict] = []
        seen: set[tuple] = set()
        for act in acts:
            tid = act.get("target_chembl_id")
            if not tid:
                continue

            t = target_cache.get(tid)
            if not t:
                # 正常不会走到这里；若发生说明数据异常，停止并提示
                logging.error("Target details missing for '%s', aborting fetch", tid)
                st.toast("靶点获取错误，停止信息获取", icon="⚠️")
                return pd.DataFrame(), None

            gene_name = _extract_gene_name(t)
            if not gene_name:
                continue

            pref_name = t.get("pref_name", "") or ""
            organism = t.get("organism", "") or ""
            if organism != "Homo sapiens":
                # 双重保险：只保留人类靶点
                continue
            assay_type = act.get("assay_type", "") or ""
            pchembl_value = act.get("pchembl_value")

            key = (gene_name, assay_type, pchembl_value)
            if key in seen:
                continue
            seen.add(key)

            rows.append({
                "approvedSymbol": gene_name,
                "target_id": tid,
                "approvedName": pref_name,
                "organism": organism,
                "assay_type": assay_type,
                "pchembl_value": pchembl_value,
            })

        df = pd.DataFrame(rows).reset_index(drop=True)
        # 不设数量上限，按 pChEMBL 活性值降序排列（无值的排最后）
        df = _sort_drug_df_by_pchembl(df)

        # 保存到本地缓存（空结果也落盘，避免反复查询同一查不到的药物）
        save_csv_cache(df, csv_path, skip_empty=False)

        return df, csv_path

    except Exception:
        logging.exception("fetch_drug_targets_chembl failed for '%s'", drug_name)
        return pd.DataFrame(), None


def set_chembl_cache_ttl(hours: int) -> None:
    """调整 ChEMBL 客户端本地缓存的有效期（默认 24 小时）。

    hours 为缓存小时数：设得越小，数据越新、联网请求越频繁；
    超过有效期的旧数据会自动重新联网获取。
    """
    try:
        from chembl_webresource_client.settings import Settings
        Settings.Instance().CACHE_EXPIRE = hours * 3600
        from chembl_webresource_client.new_client import new_client
        for name in ("molecule", "activity", "target"):
            qs = getattr(new_client, name, None)
            query = getattr(qs, "query", None) if qs is not None else None
            session = getattr(query, "session", None) if query is not None else None
            if session is not None:
                session.expire_after = hours * 3600
    except Exception:
        logging.exception("Failed to set ChEMBL cache TTL to %d hours", hours)
