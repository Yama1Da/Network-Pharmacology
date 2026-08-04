import os
import logging
import hashlib
from datetime import datetime

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

    缓存
    ----
    结果持久化为 CSV（Static/data/cache/disease/OpenTargets_{disease_name}.csv），
    命中缓存时直接读取本地文件，不访问远端 API。
    """
    cache_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache", "disease"
    )
    os.makedirs(cache_dir, exist_ok=True)
    csv_path = os.path.join(cache_dir, f"OpenTargets_{disease_name}.csv")

    # ---- cache hit: load local ----
    if os.path.exists(csv_path):
        logging.info("Loading cached disease targets for '%s'", disease_name)
        return pd.read_csv(csv_path)

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
        # ---- cache to CSV ----
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")
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
            os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache", "drug"
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

    当前为【全靶点模式】：获取药物在 ChEMBL 中的全部活性记录
    （不限数量、不过滤 pChEMBL 值），再提取基因符号，并按 pChEMBL 降序排列。
    靶点更全，但会包含 ALB、CYP450 等弱关联/非特异性结合靶点。

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
    cache_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache", "drug"
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

        # ---- Step 2: 获取全部活性数据（全靶点模式） ----
        # 【全靶点模式】不筛选 pChEMBL 值，取药物全部活性记录，靶点更全。
        # 后续如需切换为【高置信度模式】，在此 filter 中加回 pchembl_value__isnull=False，
        # 即可只保留有 pChEMBL 定量值的靶点（更少但更可靠）。
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


def list_cache_files(category: str) -> pd.DataFrame:
    """列出 Static/data/cache/<category> 目录下缓存文件的文件名与创建时间。

    Parameters
    ----------
    category : str
        子目录名，如 "drug" 或 "disease"。

    Returns
    -------
    pd.DataFrame
        含 "File name"、"Created time" 两列；目录不存在或为空时返回空表。
    """
    cache_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache", category
    )
    if not os.path.isdir(cache_dir):
        return pd.DataFrame(columns=["File name", "Created time"])
    rows = []
    for name in sorted(os.listdir(cache_dir)):
        path = os.path.join(cache_dir, name)
        if os.path.isfile(path):
            created = datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append({"File name": name, "Created time": created})
    return pd.DataFrame(rows, columns=["File name", "Created time"])


# ==================== PPI 网络模块 ====================


def get_common_targets(drug_name: str, disease_name: str, top_n: int = 1000) -> list[str]:
    """
    计算药物靶点与疾病靶点的交集基因列表（PPI 网络输入）。

    复用 fetch_drug_targets_chembl 与 get_disease_targets（自带 CSV 缓存），
    取两表 approvedSymbol 的交集，按字母序排序并截断到 top_n（默认对齐 STRING 查询上限 1000）。

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
    drug_df, _ = fetch_drug_targets_chembl(drug_name)
    disease_df = get_disease_targets(disease_name)
    drug_set = set(drug_df["approvedSymbol"].dropna()) if not drug_df.empty else set()
    disease_set = set(disease_df["approvedSymbol"].dropna()) if not disease_df.empty else set()
    return sorted(drug_set & disease_set)[:top_n]


def _ppi_cache_key(genes: list[str], species: str, score: float, network_type: str) -> str:
    """根据基因列表与查询参数生成 PPI 缓存文件名（md5）：相同输入命中同一缓存文件。"""
    payload = (
        "\r".join(sorted(set(genes)))
        + f"|{species}|{int(score * 1000)}|{network_type}"
    )
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def _string_post(endpoint: str, params: dict) -> requests.Response | None:
    """对 STRING API 发送 POST 请求，带超时与重试；最终失败返回 None。

    STRING 对"提交的标识符全部无法解析"的请求返回 404，这里视为空结果而非错误。
    """
    import time
    import config

    url = f"{config.STRING_BASE_URL}/{endpoint}"
    for attempt in range(config.API_MAX_RETRIES):
        try:
            resp = requests.post(url, data=params, timeout=config.API_TIMEOUT)
            if resp.status_code == 404:
                logging.warning("STRING %s: 404 (no identifiers resolved)", endpoint)
                return None
            resp.raise_for_status()
            return resp
        except Exception:
            logging.exception(
                "STRING %s request failed (attempt %d/%d)",
                endpoint, attempt + 1, config.API_MAX_RETRIES,
            )
            if attempt < config.API_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)
    return None


def fetch_string_ppi(
    genes: list[str],
    species: str = "9606",
    score: float = 0.7,
    network_type: str = "full",
) -> pd.DataFrame:
    """
    通过 STRING API 获取基因列表的蛋白互作边表。

    请求逻辑分两步（STRING 官方最佳实践）：
      1. get_string_ids：把基因符号映射为 STRING 蛋白 ID，同时过滤掉
         STRING 不认识的基因（映射成功的才算数）；
      2. network：用映射后的 STRING ID 获取相互作用边，按 required_score
         （= score×1000）过滤低置信度连接。

    结果缓存在 Static/data/cache/ppi 下（按 基因列表+物种+阈值+网络类型 生成
    md5 文件名），命中缓存时直接读取本地文件，不访问远端 API。

    Parameters
    ----------
    genes : list[str]
        基因符号列表（如 ["TP53", "BRCA1"]），会自动去重排序并截断到 1000。
    species : str
        NCBI taxon id，如 "9606"（人类）。
    score : float
        置信度阈值（0~1），映射为 STRING 的 required_score（0~1000）。
    network_type : str
        "full"（功能关联，默认）或 "physical"（仅物理结合）。

    Returns
    -------
    pd.DataFrame
        每行一条相互作用边，含 protein_A / protein_B / score 及各证据分数列；
        请求失败、全部基因不被识别或没有符合阈值的边时返回空 DataFrame。
    """
    import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes:
        return pd.DataFrame()
    if len(genes) > config.STRING_QUERY_LIMIT:
        genes = genes[:config.STRING_QUERY_LIMIT]

    os.makedirs(config.PPI_CACHE_DIR, exist_ok=True)
    cache_file = os.path.join(
        config.PPI_CACHE_DIR, f"StringPPI_{_ppi_cache_key(genes, species, score, network_type)}.csv"
    )

    # ---- cache hit: load local ----
    if os.path.exists(cache_file):
        logging.info("Loading cached STRING PPI for %d genes", len(genes))
        return pd.read_csv(cache_file)

    # ---- Step 1: 标识符映射（基因符号 → STRING ID） ----
    id_resp = _string_post("get_string_ids", {
        "identifiers": "\r".join(genes),
        "species": species,
        "caller_identity": config.STRING_CALLER_IDENTITY,
    })
    if id_resp is None:
        return pd.DataFrame()
    try:
        id_rows = id_resp.json()
    except Exception:
        logging.exception("STRING get_string_ids response parse failed")
        return pd.DataFrame()

    string_ids = [r.get("stringId") for r in id_rows if r.get("stringId")]
    if not string_ids:
        logging.warning("STRING: none of %d genes were recognized", len(genes))
        return pd.DataFrame()

    # ---- Step 2: 获取相互作用边 ----
    net_resp = _string_post("network", {
        "identifiers": "\r".join(string_ids),
        "species": species,
        "required_score": int(score * 1000),
        "network_type": network_type,
        "caller_identity": config.STRING_CALLER_IDENTITY,
    })
    if net_resp is None:
        return pd.DataFrame()
    try:
        rows = net_resp.json()
    except Exception:
        logging.exception("STRING network response parse failed")
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    if df.empty:
        logging.info("STRING: no edges above score %s for %d genes", score, len(genes))
        return df

    # 只保留输入基因内部的边（大小写不敏感），整理列名
    gene_set = {g.upper() for g in genes}
    df = df[
        df["preferredName_A"].str.upper().isin(gene_set)
        & df["preferredName_B"].str.upper().isin(gene_set)
    ].rename(columns={"preferredName_A": "protein_A", "preferredName_B": "protein_B"})

    keep = ["protein_A", "protein_B", "score", "nscore", "fscore", "pscore",
            "ascore", "escore", "dscore", "tscore"]
    df = df[[c for c in keep if c in df.columns]]

    # ---- cache to CSV ----
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    logging.info("STRING PPI for %d genes: %d edges cached", len(genes), len(df))
    return df


def build_ppi_graph(edges_df: pd.DataFrame, all_genes: list[str]) -> "nx.Graph":
    """
    用 STRING 边表构建 networkx 无向图。

    - 全部输入基因都作为节点（无边的基因为孤立节点，degree=0）；
    - 去除自环与重复边（同一对基因保留置信度最高的一条）；
    - 边附带 score（STRING 原始分 0~1000）与 weight（归一化 0~1）。
    """
    import networkx as nx

    G = nx.Graph()
    G.add_nodes_from(all_genes)
    if edges_df.empty:
        return G

    edges = edges_df[edges_df["protein_A"] != edges_df["protein_B"]].copy()
    # 统一 A/B 方向后按键去重，保留 score 最高的边
    edges["_key"] = edges.apply(
        lambda r: tuple(sorted((str(r["protein_A"]), str(r["protein_B"])))), axis=1
    )
    edges = edges.sort_values("score", ascending=False).drop_duplicates(subset="_key")
    for _, r in edges.iterrows():
        s = float(r["score"])
        G.add_edge(r["protein_A"], r["protein_B"], score=s, weight=s / 1000.0)
    return G


def compute_ppi_metrics(graph: "nx.Graph") -> pd.DataFrame:
    """
    计算网络节点的拓扑指标：degree / betweenness / closeness。

    孤立节点（degree=0）的中心性为 0；结果按 degree 降序排列。
    """
    import networkx as nx

    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph)
    closeness = nx.closeness_centrality(graph)

    df = pd.DataFrame({
        "gene": list(graph.nodes()),
        "degree": [degree[n] for n in graph.nodes()],
        "betweenness": [betweenness.get(n, 0.0) for n in graph.nodes()],
        "closeness": [closeness.get(n, 0.0) for n in graph.nodes()],
    })
    return df.sort_values("degree", ascending=False).reset_index(drop=True)


def mark_hub_genes(metrics_df: pd.DataFrame, min_degree: int = 3) -> pd.DataFrame:
    """
    标记 Hub 基因：degree 大于全部节点 degree 的中位数，且 degree ≥ min_degree（兜底阈值）。

    返回原表并新增 is_hub 布尔列。
    """
    df = metrics_df.copy()
    threshold = df["degree"].median()
    df["is_hub"] = (df["degree"] > threshold) & (df["degree"] >= min_degree)
    return df


def render_pyvis_html(graph: "nx.Graph", metrics_df: pd.DataFrame) -> str:
    """
    用 PyVis 生成交互式 PPI 网络 HTML。

    - 节点大小 ∝ degree，Hub 基因金色高亮、普通节点淡蓝；
    - 边宽 ∝ 置信度分数；
    - 物理引擎迭代次数取 config.PHYSICS_ITERATIONS，防止节点过度拥挤。

    Returns
    -------
    str
        HTML 内容字符串（可直接用 st.iframe 内嵌或保存为 .html 文件）。
    """
    import config
    from pyvis.network import Network

    hub_set = set(metrics_df.loc[metrics_df["is_hub"], "gene"])
    degree_map = dict(zip(metrics_df["gene"], metrics_df["degree"]))
    max_deg = max(metrics_df["degree"].max(), 1)

    nt = Network(height="600px", width="100%", directed=False, notebook=False)
    for node in graph.nodes():
        deg = degree_map.get(node, 0)
        size = 10 + int(40 * deg / max_deg)
        color = "#F5C518" if node in hub_set else "#64B5F6"
        nt.add_node(
            node, label=node, size=size, color=color,
            title=f"{node}<br>degree: {deg}",
        )
    for u, v, data in graph.edges(data=True):
        s = float(data.get("score", 0))
        nt.add_edge(u, v, value=s / 1000, title=f"score: {s}")

    nt.set_options(f"""
    {{
      "physics": {{
        "solver": "forceAtlas2Based",
        "stabilization": {{"iterations": {config.PHYSICS_ITERATIONS}}},
        "forceAtlas2Based": {{
          "gravitationalConstant": -50,
          "springLength": 100,
          "springConstant": 0.08
        }}
      }},
      "interaction": {{"hover": true}}
    }}
    """)
    return nt.generate_html()


def render_static_plot(graph: "nx.Graph", metrics_df: pd.DataFrame) -> "matplotlib.figure.Figure":
    """
    用 networkx + matplotlib 绘制静态 PPI 网络图（论文用）。

    节点颜色映射 degree（附 colorbar）、大小 ∝ degree；只标注 Hub 基因标签避免拥挤。
    """
    import networkx as nx
    import matplotlib.pyplot as plt
    from matplotlib import cm

    degrees = dict(graph.degree())
    node_degrees = [degrees[n] for n in graph.nodes()]
    max_deg = max(node_degrees, default=1)
    n_nodes = max(len(graph.nodes()), 1)

    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(graph, k=2.2 / (n_nodes ** 0.5), iterations=200, seed=42)

    nx.draw_networkx_edges(graph, pos, ax=ax, alpha=0.25, edge_color="#9e9e9e", width=0.8)
    nodes = nx.draw_networkx_nodes(
        graph, pos, ax=ax,
        node_size=[80 + 300 * d / max_deg for d in node_degrees],
        node_color=node_degrees,
        cmap=cm.viridis, vmin=0, vmax=max_deg,
        edgecolors="white", linewidths=0.6,
    )
    hub_set = set(metrics_df.loc[metrics_df["is_hub"], "gene"])
    nx.draw_networkx_labels(
        graph, pos, ax=ax,
        labels={n: n for n in graph.nodes() if n in hub_set}, font_size=8,
    )
    cb = fig.colorbar(nodes, ax=ax, shrink=0.8)
    cb.set_label("Degree")
    ax.set_title("PPI Network (node color = degree)", fontsize=14)
    ax.axis("off")
    fig.tight_layout()
    return fig


def save_ppi_results(
    graph: "nx.Graph",
    static_fig: "matplotlib.figure.Figure",
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    html_content: str,
    prefix: str = "PPI",
) -> dict:
    """
    将一次 PPI 分析的全部结果保存到 DOWNLOADS_DIR/<prefix>/ 目录（替代式）。

    【替代式保存】：保存前会删除同前缀的旧目录（含历史时间戳目录），
    同一 prefix（同一药物+疾病）始终只保留最新一次分析，不会堆积。

    目录内容：
      - {prefix}_nodes.csv      节点拓扑指标表（含 is_hub 列）
      - {prefix}_edges.csv      边表（protein_A / protein_B / score ...）
      - {prefix}_network.html   交互式网络
      - {prefix}_graph.graphml  可编辑网络文件（供 Cytoscape / Gephi 精修）
      - {prefix}_network.svg    静态图矢量版（论文出版用）

    Parameters
    ----------
    graph : networkx.Graph
        分析构建的网络图，用于 GraphML 导出。
    static_fig : matplotlib.figure.Figure
        静态网络图，用于 SVG 导出。
    nodes_df : pd.DataFrame
        节点拓扑指标表。
    edges_df : pd.DataFrame
        边表。
    html_content : str
        交互式网络的 HTML 内容。
    prefix : str
        文件名前缀，建议用 药物_疾病 标识本次分析。

    Returns
    -------
    dict
        {"nodes", "edges", "html", "graphml", "svg", "directory"} 各文件路径与目录路径。
    """
    import glob
    import shutil
    import networkx as nx
    import config

    os.makedirs(config.DOWNLOADS_DIR, exist_ok=True)

    # 替代式保存：删除同前缀的旧目录（含历史时间戳目录），只保留本次最新分析
    for old in glob.glob(os.path.join(config.DOWNLOADS_DIR, f"{prefix}_*")):
        if os.path.isdir(old):
            shutil.rmtree(old)
    subdir = os.path.join(config.DOWNLOADS_DIR, prefix)
    if os.path.isdir(subdir):
        shutil.rmtree(subdir)
    os.makedirs(subdir, exist_ok=True)

    paths = {
        "nodes": os.path.join(subdir, f"{prefix}_nodes.csv"),
        "edges": os.path.join(subdir, f"{prefix}_edges.csv"),
        "html": os.path.join(subdir, f"{prefix}_network.html"),
        "graphml": os.path.join(subdir, f"{prefix}_graph.graphml"),
        "svg": os.path.join(subdir, f"{prefix}_network.svg"),
    }
    nodes_df.to_csv(paths["nodes"], index=False, encoding="utf-8-sig")
    edges_df.to_csv(paths["edges"], index=False, encoding="utf-8-sig")
    with open(paths["html"], "w", encoding="utf-8") as f:
        f.write(html_content)
    nx.write_graphml(graph, paths["graphml"])
    static_fig.savefig(paths["svg"], format="svg", bbox_inches="tight")

    paths["directory"] = subdir
    return paths


def make_ppi_zip(file_paths: dict) -> bytes:
    """
    将一次分析的结果文件（nodes/edges/html/graphml/svg）打包为 zip 字节。

    供 st.download_button 直接下载；文件在 zip 内使用各自的文件名，
    配合独立子目录，用户解压后即可获得该次分析的全部成品。
    """
    import io
    import zipfile

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name in ("nodes", "edges", "html", "graphml", "svg"):
            path = file_paths.get(name)
            if path and os.path.exists(path):
                zf.write(path, arcname=os.path.basename(path))
    return buf.getvalue()
