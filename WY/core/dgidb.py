"""DGIdb 数据源：查询药物的相互作用基因（drug-gene interaction）。"""

import logging

import pandas as pd

from .caching import open_named_csv_cache, save_csv_cache
from .http_client import request_with_retry
from .text_utils import remember_name_mapping


def fetch_drug_targets_dgidb(drug_name: str) -> pd.DataFrame:
    """从 DGIdb 获取药物的相互作用基因，含交互分数与证据来源（文献、数据库）。

    失败或无数据时返回空表；结果存为本地 CSV，已有缓存则直接复用。
    """
    import config

    # 已有缓存则直接读取本地文件
    csv_path, cached = open_named_csv_cache(
        config.DRUG_CACHE_DIR, "DGIdb_", drug_name,
        label=f"DGIdb targets for '{drug_name}'",
        remember=remember_name_mapping,  # 记录安全文件名→原始名,供管理页还原
    )
    if cached is not None:
        return cached

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
        # 网络偶发失败时做少量重试，避免一次失败就返回空表
        resp = request_with_retry(
            config.DGIDB_GRAPHQL_URL,
            json={"query": query, "variables": {"names": [drug_name.strip()]}},
            max_retries=config.LIGHT_API_MAX_RETRIES, timeout=30,
            label="DGIdb request", use_cache=True,
        )
        payload = resp.json() if resp is not None else None

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
        save_csv_cache(df, csv_path)
        logging.info("DGIdb '%s': %d targets fetched", drug_name, len(df))
        return df

    except Exception:
        logging.exception("fetch_drug_targets_dgidb failed for '%s'", drug_name)
        return pd.DataFrame()
