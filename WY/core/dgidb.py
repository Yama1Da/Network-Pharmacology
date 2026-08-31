"""DGIdb v5 数据源：药物-基因交互（GraphQL API，原生 requests 调用）。"""

import logging
import os

import pandas as pd
import requests

from .paths import BASE_DIR
from .text_utils import safe_name


def fetch_drug_targets_dgidb(drug_name: str) -> pd.DataFrame:
    """从 DGIdb v5 GraphQL API 获取药物的交互基因（drug-gene interaction）。

    直接通过 requests 调用 DGIdb 的 GraphQL 端点（dgipy 0.1.0 与 gql 4.x
    存在 execute 签名不兼容的 bug，故弃用第三方库、改用原生 requests）。
    返回药物与基因的交互记录，含交互分数与来源证据。

    Parameters
    ----------
    drug_name : str
        药物英文名，如 "aspirin"。

    Returns
    -------
    pd.DataFrame
        列: approvedSymbol（基因符号）、drug_name、interaction_score、
        interaction_sources、interaction_pmids；无数据或失败时返回空表。

    缓存
    ----
    结果持久化为 CSV（Static/data/cache/drug/DGIdb_{drug_name}.csv）。
    """
    import time
    import config

    cache_dir = os.path.join(BASE_DIR, "Static", "data", "cache", "drug")
    os.makedirs(cache_dir, exist_ok=True)
    csv_path = os.path.join(cache_dir, f"DGIdb_{safe_name(drug_name)}.csv")

    # ---- cache hit: load local ----
    if os.path.exists(csv_path):
        logging.info("Loading cached DGIdb targets for '%s'", drug_name)
        return pd.read_csv(csv_path)

    query = """
    query getInteractionsByDrug($names: [String!]) {
      drugs(names: $names) {
        nodes {
          interactions {
            gene { name longName conceptId }
            drug { name conceptId approved }
            interactionScore
            interactionClaims {
              publications { pmid }
              source { sourceDbName }
            }
          }
        }
      }
    }
    """

    try:
        # 请求带重试：网络抖动 / 临时 5xx 时按退避重试，避免一次失败就返回空表。
        # 单次分析仅 1 次请求，用轻量重试次数（2），避免网络不佳时拖慢更新。
        payload = None
        for attempt in range(config.LIGHT_API_MAX_RETRIES):
            try:
                resp = requests.post(
                    config.DGIDB_GRAPHQL_URL,
                    json={"query": query, "variables": {"names": [drug_name.strip()]}},
                    timeout=30,
                )
                resp.raise_for_status()
                payload = resp.json()
                break
            except Exception:
                logging.exception(
                    "DGIdb request failed (attempt %d/%d)",
                    attempt + 1, config.LIGHT_API_MAX_RETRIES,
                )
                if attempt < config.LIGHT_API_MAX_RETRIES - 1:
                    time.sleep(config.API_RETRY_DELAY)

        if payload is None:
            return pd.DataFrame()

        if "errors" in payload or "data" not in payload:
            logging.error("DGIdb query failed: %s", payload.get("errors"))
            return pd.DataFrame()

        records: list[dict] = []
        for node in payload.get("data", {}).get("drugs", {}).get("nodes", []):
            for inter in node.get("interactions", []):
                gene = inter.get("gene") or {}
                drug = inter.get("drug") or {}
                sources = [
                    c.get("source", {}).get("sourceDbName")
                    for c in inter.get("interactionClaims", [])
                    if c.get("source", {}).get("sourceDbName")
                ]
                pmids = [
                    p.get("pmid")
                    for c in inter.get("interactionClaims", [])
                    for p in c.get("publications", [])
                    if p.get("pmid")
                ]
                records.append({
                    "approvedSymbol": gene.get("name"),
                    "drug_name": drug.get("name"),
                    "interaction_score": inter.get("interactionScore"),
                    "interaction_sources": ";".join(dict.fromkeys(sources)),
                    "interaction_pmids": ";".join(str(p) for p in dict.fromkeys(pmids)),
                })

        if not records:
            logging.warning("DGIdb returned no interactions for '%s'", drug_name)
            return pd.DataFrame()

        df = (
            pd.DataFrame(records)
            .drop_duplicates(subset=["approvedSymbol", "drug_name"])
            .reset_index(drop=True)
        )
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")
        logging.info("DGIdb '%s': %d targets fetched", drug_name, len(df))
        return df

    except Exception:
        logging.exception("fetch_drug_targets_dgidb failed for '%s'", drug_name)
        return pd.DataFrame()
