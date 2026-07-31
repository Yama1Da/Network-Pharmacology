import os
import logging

import pandas as pd
import requests
import streamlit as st

OPENTARGETS_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def _safe_round(value, ndigits: int = 4) -> float:
    """None 安全取分：Open Targets 部分靶点的指标字段可能为 null，统一转 0.0 再四舍五入。"""
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return 0.0


@st.cache_data
def get_disease_targets(disease_name: str, top_n: int = 200) -> pd.DataFrame:
    """
    通过 Open Targets Platform GraphQL API 获取疾病的关联靶点及各项关联指标。

    请求逻辑分两步：
      1. search 查询：按疾病英文名模糊搜索，取排名第一的 disease 实体的
         疾病标识符（EFO / MONDO 等）。
      2. disease 查询：用该标识符查询 associatedTargets，获取关联分数、
         新颖性、各数据类型分数以及靶点的基因组/类别信息，提取前 top_n 条。

    Parameters
    ----------
    disease_name : str
        疾病英文名，如 "lung cancer"。
    top_n : int
        返回的靶点数量上限，默认 200。

    Returns
    -------
    pd.DataFrame
        每行一个靶点，包含基因符号、ENSG ID、基因全名、染色体位置、
        靶点类别、关联分数、新颖性及各数据类型分数等指标；
        API 请求失败或未找到该疾病时返回空 DataFrame。
    """
    try:
        # ========== 第一步：search 查询，根据疾病名获取疾病 ID ==========
        # 这个 search 端点会在整个平台（疾病、药物、靶点、研究等实体）中做模糊搜索。
        # 返回结果按相关性排序，hits[0] 是匹配度最高的记录。
        # 我们只挑 entity == "disease" 的记录，取第一个的 id 作为疾病标识符
        # （可能是 EFO_xxx、MONDO_xxx 或 Orphanet 等，disease 查询都接受）。
        search_query = """
        query SearchDisease($queryString: String!) {
            search(queryString: $queryString) {
                hits {
                    id
                    entity
                    name
                }
            }
        }
        """
        resp = requests.post(
            OPENTARGETS_API_URL,
            json={"query": search_query, "variables": {"queryString": disease_name.strip()}},
            timeout=30,
        )
        resp.raise_for_status()  # 非 2xx 状态码会抛出异常
        payload = resp.json()

        # GraphQL 层错误（如查询语法问题）会出现在 payload["errors"] 中，而不是抛异常
        if "errors" in payload or "data" not in payload:
            logging.error("Open Targets search query failed: %s", payload.get("errors"))
            return pd.DataFrame()

        disease_id = None
        for hit in payload["data"]["search"]["hits"]:
            if hit.get("entity") == "disease":
                disease_id = hit["id"]
                break
        if disease_id is None:
            logging.warning("No disease found for '%s'", disease_name)
            return pd.DataFrame()

        # ========== 第二步：disease 查询，获取关联靶点及指标 ==========
        # disease(efoId: ...) 返回疾病本体信息；
        # associatedTargets 是"疾病-靶点"关联列表，官方已按关联分数降序排列。
        # 每个关联（AssociatedTarget）携带的指标：
        #   - score          : 总关联分数（0~1，越大关联越强）
        #   - novelty        : 新颖性（0~1，越大表示该关联越少被人研究）
        #   - datatypeScores : 各数据类型贡献分数（遗传关联、体细胞突变、
        #                      文献、动物模型、通路、临床证据等）
        #   - target         : 靶点详情（ENSG ID、基因符号、全名、生物类型、
        #                      基因组位置、靶点类别）
        # page 参数控制分页：index 是页码（从 0 开始），size 是每页条数，
        # 这里 size 直接设为 top_n（默认 200），一次取回关联分数最高的前 200 个靶点。
        disease_query = """
        query DiseaseTargets($diseaseId: String!, $size: Int!) {
            disease(efoId: $diseaseId) {
                id
                name
                associatedTargets(page: {index: 0, size: $size}) {
                    rows {
                        score
                        novelty
                        datatypeScores {
                            id
                            score
                        }
                        target {
                            id
                            approvedSymbol
                            approvedName
                            biotype
                            genomicLocation {
                                chromosome
                                start
                                end
                            }
                            targetClass {
                                label
                            }
                        }
                    }
                }
            }
        }
        """
        resp = requests.post(
            OPENTARGETS_API_URL,
            json={"query": disease_query, "variables": {"diseaseId": disease_id, "size": top_n}},
            timeout=30,
        )
        resp.raise_for_status()
        payload = resp.json()

        if "errors" in payload or "data" not in payload:
            logging.error("Open Targets disease query failed: %s", payload.get("errors"))
            return pd.DataFrame()

        rows = payload["data"]["disease"]["associatedTargets"]["rows"]
        # 按关联分数降序排序（API 通常已排序，这里再排一次更保险），只保留前 top_n 条
        rows_sorted = sorted(rows, key=lambda r: r.get("score", 0.0), reverse=True)[:top_n]

        # 把每条关联记录整理成一行 DataFrame
        records: list[dict] = []
        for r in rows_sorted:
            tgt = r.get("target", {})
            gloc = tgt.get("genomicLocation", {}) or {}
            # targetClass 是列表（可能有多个类别标签），取第一个
            tc = tgt.get("targetClass") or []
            # datatypeScores 形如 [{id: "genetic_association", score: 0.93}, ...]，
            # 转成 {类型名: 分数} 字典，便于按数据类型名取分
            dt_scores = {s.get("id"): s.get("score") for s in r.get("datatypeScores", []) if s.get("id")}
            records.append({
                "approvedSymbol": tgt.get("approvedSymbol"),      # 官方基因符号
                "ensembl_id": tgt.get("id"),                       # ENSG 编号
                "approvedName": tgt.get("approvedName"),           # 基因全名
                "biotype": tgt.get("biotype"),                     # 生物类型
                "chromosome": gloc.get("chromosome"),              # 所在染色体
                "genomic_start": gloc.get("start"),                # 基因组起始位置
                "genomic_end": gloc.get("end"),                    # 基因组结束位置
                "targetClass": tc[0].get("label") if tc else None, # 靶点类别（取第一个）
                "association_score": _safe_round(r.get("score", 0.0)),           # 总关联分数
                "novelty": _safe_round(r.get("novelty", 0.0)),                   # 新颖性
                "genetic_association": _safe_round(dt_scores.get("genetic_association")),  # 遗传关联分数
                "somatic_mutation": _safe_round(dt_scores.get("somatic_mutation")),        # 体细胞突变分数
                "literature": _safe_round(dt_scores.get("literature")),                    # 文献分数
                "animal_model": _safe_round(dt_scores.get("animal_model")),                # 动物模型分数
                "affected_pathway": _safe_round(dt_scores.get("affected_pathway")),        # 通路分数
                "clinical": _safe_round(dt_scores.get("clinical")),                        # 临床证据分数
            })

        df = pd.DataFrame(records)
        logging.info("Disease '%s' (%s): %d targets fetched", disease_name, disease_id, len(df))
        return df

    except Exception:
        # 网络错误、超时、JSON 解析失败等统一在此兜底：
        # 打印完整错误日志（含堆栈），返回空 DataFrame，不让页面崩溃
        logging.exception("get_disease_targets failed for '%s'", disease_name)
        return pd.DataFrame()


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

    Uses ALL activity records of the drug to identify target proteins
    (no count limit — return as many targets as ChEMBL has), extracts
    gene symbols from target component synonyms, and sorts the result
    by pChEMBL value in descending order.

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
    cache_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache"
    )
    os.makedirs(cache_dir, exist_ok=True)
    csv_path = os.path.join(cache_dir, f"ChEMBL_{drug_name}.csv")

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

        # ---- Step 2: get ALL activities (no pChEMBL filter, no count limit) ----
        acts = list(new_client.activity.filter(
            molecule_chembl_id=chembl_id,
        ))

        # ---- Step 3: get target details (cache by target_id) ----
        from collections import OrderedDict
        target_cache: dict[str, dict | None] = {}
        for act in acts:
            tid = act.get("target_chembl_id")
            if tid and tid not in target_cache:
                try:
                    target_cache[tid] = new_client.target.get(tid)
                except Exception:
                    target_cache[tid] = None

        def _extract_gene_name(t: dict) -> str | None:
            """Extract gene symbol from a target dict."""
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
                continue

            gene_name = _extract_gene_name(t)
            if not gene_name:
                continue

            pref_name = t.get("pref_name", "") or ""
            organism = t.get("organism", "") or ""
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
