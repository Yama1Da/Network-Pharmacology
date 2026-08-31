"""ChEMBL 数据源：药物全部活性靶点（全靶点模式，仅人类），含本地 HTTP 缓存 TTL 调整。"""

import logging
import os

import pandas as pd
import streamlit as st

from .paths import BASE_DIR
from .text_utils import safe_name


def _sort_drug_df_by_pchembl(df: pd.DataFrame) -> pd.DataFrame:
    """按 pchembl_value 降序排列药物靶点表格（空值排最后）。"""
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
    Fetch drug targets from ChEMBL via the webresource client.

    If a cached CSV already exists for this drug, loads it directly
    instead of calling the remote API.

    当前为【全靶点模式】：获取药物在 ChEMBL 中的全部活性记录
    （不限数量、不过滤 pChEMBL 值），再提取基因符号，并按 pChEMBL 降序排列。
    靶点更全，但会包含 ALB、CYP450 等弱关联/非特异性结合靶点。
    【物种限定】：仅保留人类（organism = Homo sapiens）的靶点记录。

    若后续需要【高置信度模式】，只需在下方 Step 2 的 activity.filter 中
    加回 pchembl_value__isnull=False 参数，即可只保留有 pChEMBL 定量值的靶点
    （数量更少，但均为经过活性定量的高置信度靶点）。

    Parameters
    ----------
    drug_name : str
        Drug name in English (e.g. \"Aspirin\").

    Returns
    -------
    tuple[pd.DataFrame, str | None]
        (DataFrame with columns [approvedSymbol, target_id, approvedName,
         organism, assay_type, pchembl_value] sorted by pchembl_value desc,
         path to the cached CSV file, or None if the request failed).
    """
    cache_dir = os.path.join(BASE_DIR, "Static", "data", "cache", "drug")
    os.makedirs(cache_dir, exist_ok=True)
    csv_path = os.path.join(cache_dir, f"ChEMBL_{safe_name(drug_name)}.csv")

    # ---- cache hit: load local ----
    if os.path.exists(csv_path):
        logging.info("Loading cached ChEMBL targets for '%s'", drug_name)
        df = _sort_drug_df_by_pchembl(pd.read_csv(csv_path))
        return df, csv_path

    # ---- cache miss: call API ----
    try:
        # 延迟导入：chembl_webresource_client 在 import 的瞬间就会向远端
        # (https://www.ebi.ac.uk/chembl/api/data/spore) 请求 schema。
        # 如果放在模块顶层，ChEMBL 服务器不可用会导致整个应用启动崩溃；
        # 放在函数内 + 外层 try-except 兜底，请求失败时只返回空数据、不影响页面。
        from chembl_webresource_client.new_client import new_client

        # ---- Step 1: find molecule by name ----
        mols = list(new_client.molecule.filter(pref_name__iexact=drug_name.strip()))
        if not mols:
            # fallback: case-insensitive contains search
            mols = list(new_client.molecule.filter(pref_name__icontains=drug_name.strip()))
        if not mols:
            logging.warning("No ChEMBL molecule found for '%s'", drug_name)
            return pd.DataFrame(), None

        chembl_id = mols[0]["molecule_chembl_id"]

        # ---- Step 2: 获取全部活性数据（全靶点模式，仅人类） ----
        # 【全靶点模式】不筛选 pChEMBL 值，取药物全部活性记录，靶点更全。
        # 【物种限定】organism="Homo sapiens" 在服务端只返回人类实验物种的活性记录，
        # 减少下载量与后续靶点详情请求数。
        # 后续如需切换为【高置信度模式】，在此 filter 中加回 pchembl_value__isnull=False，
        # 即可只保留有 pChEMBL 定量值的靶点（更少但更可靠）。
        acts = list(new_client.activity.filter(
            molecule_chembl_id=chembl_id,
            organism="Homo sapiens",
        ))

        # ---- Step 3: get target details (cache by target_id) ----
        import time
        import config
        target_cache: dict[str, dict | None] = {}
        for act in acts:
            tid = act.get("target_chembl_id")
            if tid and tid not in target_cache:
                # 请求间隔：每个靶点请求前先等待，避免连续快速请求被 ChEMBL 限流/封禁
                time.sleep(config.API_REQUEST_DELAY)
                try:
                    t = new_client.target.get(tid)
                except Exception:
                    logging.exception("Failed to fetch target details for '%s'", tid)
                    st.toast("靶点获取错误，停止信息获取", icon="⚠️")
                    return pd.DataFrame(), None  # 失败即停止，且不保存 CSV
                if not t:
                    logging.error("Target '%s' not found, aborting fetch", tid)
                    st.toast("靶点获取错误，停止信息获取", icon="⚠️")
                    return pd.DataFrame(), None  # 失败即停止，且不保存 CSV
                target_cache[tid] = t

        def _extract_gene_name(t: dict) -> str | None:
            """从靶点字典中提取基因符号或HGNC名称，返回第一个找到的（如果找到），否则返回None"""
            if "target_components" in t:
                for comp in t["target_components"]:
                    for syn in comp.get("target_component_synonyms", []):
                        if syn.get("syn_type") == "GENE_SYMBOL":
                            return syn.get("component_synonym")
                    for xref in comp.get("target_component_xrefs", []):
                        if xref.get("xref_src_db") == "HGNC":
                            return xref.get("xref_name")
            return None

        # ---- Step 4: build one row per activity with target info ----
        rows: list[dict] = []
        seen: set[tuple] = set()
        for act in acts:
            tid = act.get("target_chembl_id")
            if not tid:
                continue

            t = target_cache.get(tid)
            if not t:
                # 正常情况下不会走到这里（Step 3 已保证所有 target 有效）；
                # 若发生则视为获取失败，停止并提示，保证数据完整性。
                logging.error("Target details missing for '%s', aborting fetch", tid)
                st.toast("靶点获取错误，停止信息获取", icon="⚠️")
                return pd.DataFrame(), None

            gene_name = _extract_gene_name(t)
            if not gene_name:
                continue

            pref_name = t.get("pref_name", "") or ""
            organism = t.get("organism", "") or ""
            if organism != "Homo sapiens":
                # 兜底过滤：即使 activity 层 organism 过滤有遗漏，也只保留人类靶点
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
        # 全量展示：不设数量上限，按 pchembl_value 降序排列（空值排最后）
        df = _sort_drug_df_by_pchembl(df)

        # ---- cache to CSV ----
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")

        return df, csv_path

    except Exception:
        logging.exception("fetch_drug_targets_chembl failed for '%s'", drug_name)
        return pd.DataFrame(), None


def set_chembl_cache_ttl(hours: int) -> None:
    """调整 chembl_webresource_client 的 HTTP 响应缓存有效期（TTL）。

    chembl 客户端默认把每个 API 响应缓存在本地 sqlite（CACHE_EXPIRE=24h）。
    本函数把 TTL 改为指定小时数，并同步更新当前进程内已创建的会话，
    使其对之后保存的响应立即生效（超过 TTL 的响应会自动重新联网请求）。

    Parameters
    ----------
    hours : int
        缓存小时数。数值越小数据越新、真实请求越频繁。
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
