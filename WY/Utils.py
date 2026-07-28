import os
import logging

import pandas as pd
import requests
import streamlit as st
from chembl_webresource_client.new_client import new_client

# 读取疾病靶点列表（测试假数据）
@st.cache_data
def get_disease_targets(disease_name: str) -> list[str]:
    """读取 static/data/{disease_name}.csv，返回 'Gene' 列的靶点列表。"""
    df = pd.read_csv(f"Static/data/{disease_name}.csv")
    return df["Gene"].tolist()


OPENTARGETS_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def _search_drug(drug_name: str) -> str | None:
    """Search Open Targets for a drug and return its ChEMBL ID, or None."""
    query = """
    query SearchDrug($queryString: String!) {
        search(queryString: $queryString) {
            hits {
                id
                entity
            }
        }
    }
    """
    resp = requests.post(
        OPENTARGETS_API_URL,
        json={"query": query, "variables": {"queryString": drug_name.strip()}},
        timeout=30,
    )
    resp.raise_for_status()
    payload = resp.json()
    for hit in payload.get("data", {}).get("search", {}).get("hits", []):
        if hit.get("entity") == "drug":
            return hit["id"]
    return None


def fetch_drug_targets_opentargets(drug_name: str) -> tuple[pd.DataFrame, str | None]:
    """
    Fetch drug targets from the Open Targets Platform GraphQL API.

    Uses a two-step approach:
      1. Search for the drug by name to get its ChEMBL ID.
      2. Query mechanismsOfAction to retrieve target genes.

    Parameters
    ----------
    drug_name : str
        Drug name in English (e.g. \"Aspirin\").

    Returns
    -------
    tuple[pd.DataFrame, str | None]
        (DataFrame with columns [approvedSymbol, ensembl_id, approvedName],
         path to the cached CSV file, or None if the request failed).
    """
    try:
        # ---- Step 1: search ----
        chembl_id = _search_drug(drug_name)
        if chembl_id is None:
            logging.warning("No drug found for '%s'", drug_name)
            return pd.DataFrame(), None

        # ---- Step 2: get targets via mechanismsOfAction ----
        query = """
        query DrugTargets($chemblId: String!) {
            drug(chemblId: $chemblId) {
                id
                name
                mechanismsOfAction {
                    rows {
                        targets {
                            id
                            approvedSymbol
                            approvedName
                        }
                    }
                }
            }
        }
        """
        resp = requests.post(
            OPENTARGETS_API_URL,
            json={"query": query, "variables": {"chemblId": chembl_id}},
            timeout=30,
        )
        resp.raise_for_status()
        payload = resp.json()

        rows: list[dict] = []
        for row in (
            payload.get("data", {})
            .get("drug", {})
            .get("mechanismsOfAction", {})
            .get("rows", [])
        ):
            for tgt in row.get("targets", []):
                rows.append({
                    "approvedSymbol": tgt.get("approvedSymbol"),
                    "ensembl_id": tgt.get("id"),
                    "approvedName": tgt.get("approvedName"),
                })

        df = pd.DataFrame(rows).drop_duplicates(subset=["ensembl_id"]).reset_index(drop=True)

        # ---- cache to CSV ----
        cache_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache"
        )
        os.makedirs(cache_dir, exist_ok=True)
        csv_path = os.path.join(cache_dir, f"OpenTargets_{drug_name}.csv")
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")

        return df, csv_path

    except Exception:
        logging.exception("fetch_drug_targets_opentargets failed for '%s'", drug_name)
        return pd.DataFrame(), None


def fetch_drug_targets_chembl(drug_name: str) -> tuple[pd.DataFrame, str | None]:
    """
    Fetch drug targets from ChEMBL via the webresource client.

    If a cached CSV already exists for this drug, loads it directly
    instead of calling the remote API.

    Uses activity data with pChEMBL values to identify target proteins,
    then extracts gene symbols from target component synonyms.

    Parameters
    ----------
    drug_name : str
        Drug name in English (e.g. \"Aspirin\").

    Returns
    -------
    tuple[pd.DataFrame, str | None]
        (DataFrame with columns [approvedSymbol, target_id, approvedName],
         path to the cached CSV file, or None if the request failed).
    """
    cache_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache"
    )
    os.makedirs(cache_dir, exist_ok=True)
    csv_path = os.path.join(cache_dir, f"ChEMBL_{drug_name}.csv")

    # ---- cache hit: load local ----
    if os.path.exists(csv_path):
        logging.info("Loading cached ChEMBL targets for '%s'", drug_name)
        df = pd.read_csv(csv_path)
        return df, csv_path

    # ---- cache miss: call API ----
    try:
        # ---- Step 1: find molecule by name ----
        mols = list(new_client.molecule.filter(pref_name__iexact=drug_name.strip()))
        if not mols:
            # fallback: case-insensitive contains search
            mols = list(new_client.molecule.filter(pref_name__icontains=drug_name.strip()))
        if not mols:
            logging.warning("No ChEMBL molecule found for '%s'", drug_name)
            return pd.DataFrame(), None

        chembl_id = mols[0]["molecule_chembl_id"]

        # ---- Step 2: get activities with pChEMBL ----
        acts = list(new_client.activity.filter(
            molecule_chembl_id=chembl_id,
            pchembl_value__isnull=False,
        ))

        # ---- Step 3: collect unique target IDs ----
        target_ids: set[str] = set()
        for act in acts:
            tid = act.get("target_chembl_id")
            if tid:
                target_ids.add(tid)

        # ---- Step 4: get target details & gene symbols ----
        rows: list[dict] = []
        for tid in sorted(target_ids):
            try:
                t = new_client.target.get(tid)
            except Exception:
                continue
            if not t:
                continue

            pref_name = t.get("pref_name", "") or ""

            # extract gene symbol from target_component synonyms
            gene_name = None
            if "target_components" in t:
                for comp in t["target_components"]:
                    for syn in comp.get("target_component_synonyms", []):
                        if syn.get("syn_type") == "GENE_SYMBOL":
                            gene_name = syn.get("component_synonym")
                            break
                    if gene_name:
                        break

            # fallback: HGNC xref
            if not gene_name and "target_components" in t:
                for comp in t["target_components"]:
                    for xref in comp.get("target_component_xrefs", []):
                        if xref.get("xref_src_db") == "HGNC":
                            gene_name = xref.get("xref_name")
                            break
                    if gene_name:
                        break

            if gene_name:
                rows.append({
                    "approvedSymbol": gene_name,
                    "target_id": tid,
                    "approvedName": pref_name,
                })

        df = pd.DataFrame(rows).drop_duplicates(subset=["approvedSymbol"]).reset_index(drop=True)

        # ---- cache to CSV ----
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")

        return df, csv_path

    except Exception:
        logging.exception("fetch_drug_targets_chembl failed for '%s'", drug_name)
        return pd.DataFrame(), None
