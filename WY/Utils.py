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
    # 通过Open Targets API搜索药物名称，返回其ChEMBL ID（如果找到），否则返回None
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
        import time
        import config
        from collections import OrderedDict
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


#TTL时间修改
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


def list_download_dirs() -> pd.DataFrame:
    """列出 DOWNLOADS_DIR 下的子目录（一次分析一组）及其创建时间。

    Returns
    -------
    pd.DataFrame
        含 "Directory"、"Created time" 两列；目录不存在或无子目录时返回空表。
    """
    import config

    if not os.path.isdir(config.DOWNLOADS_DIR):
        return pd.DataFrame(columns=["Directory", "Created time"])
    rows = []
    for name in sorted(os.listdir(config.DOWNLOADS_DIR)):
        path = os.path.join(config.DOWNLOADS_DIR, name)
        if os.path.isdir(path):
            created = datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append({"Directory": name, "Created time": created})
    return pd.DataFrame(rows, columns=["Directory", "Created time"])


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
    # 从药物和疾病靶点表中提取基因符号（approvedSymbol），并转换为集合
    drug_set = set(drug_df["approvedSymbol"].dropna()) if not drug_df.empty else set()
    disease_set = set(disease_df["approvedSymbol"].dropna()) if not disease_df.empty else set()
    # 取药物和疾病靶点的交集，按字母序排序并截断到 top_n
    return sorted(drug_set & disease_set)[:top_n]


def _ppi_cache_key(genes: list[str], species: str, score: float, network_type: str) -> str:
    # 根据基因列表与查询参数生成 PPI 缓存文件名（md5）：相同输入命中同一缓存文件。
    #该函数作为后续存储为CSV文件的辅助函数，确保每个查询结果对应一个唯一的文件名
    payload = (
        "\r".join(sorted(set(genes)))
        + f"|{species}|{int(score * 1000)}|{network_type}"
    )
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def _read_cache_csv(cache_file: str) -> pd.DataFrame:
    # 读取本地缓存 CSV（兼容旧版带 # 元数据头的md5文件；新文件为纯数据 CSV）
    return pd.read_csv(cache_file, encoding="utf-8-sig", comment="#")


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
    label: str | None = None,
) -> pd.DataFrame:
    """
    通过 STRING API 获取基因列表的蛋白互作边表。

    请求逻辑分两步（STRING 官方最佳实践）：
      1. get_string_ids：把基因符号映射为 STRING 蛋白 ID，同时过滤掉
         STRING 不认识的基因（映射成功的才算数）；
      2. network：用映射后的 STRING ID 获取相互作用边，按 required_score
         （= score×1000）过滤低置信度连接。

    结果缓存在 Static/data/cache/ppi 下，文件名为 StringPPI_{label}_{md5}.csv
    （label 为来源标识，如 药物名_疾病名；无 label 时退化为 StringPPI_{md5}.csv），
    缓存首行为 # 元数据注释头（记录 label/物种/阈值等），命中缓存直接读本地，不访问远端 API。

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
    label : str | None
        来源标识（如 "aspirin_cancer" 或 "manual"），用于缓存文件名前缀，便于用户识别。

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
    md5 = _ppi_cache_key(genes, species, score, network_type)
    name = f"StringPPI_{label}_{md5}.csv" if label else f"StringPPI_{md5}.csv"
    cache_file = os.path.join(config.PPI_CACHE_DIR, name)

    # ---- cache hit: load local ----
    if os.path.exists(cache_file):
        logging.info("Loading cached STRING PPI for %d genes", len(genes))
        return _read_cache_csv(cache_file)

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
    prefix: str,
    enrichment_df: pd.DataFrame | None = None,
    enrichment_bar_fig: "matplotlib.figure.Figure | None" = None,
    enrichment_dot_fig: "matplotlib.figure.Figure | None" = None,
) -> dict:
    """
    将一次 PPI 分析的全部结果保存到 DOWNLOADS_DIR/<prefix>_data/ 目录（替代式）。

    目录名 = <prefix>_data（prefix 如 药物名_疾病名，例如 aspirin_cancer_data），
    文件名 = <prefix>_xxx（如 aspirin_cancer_nodes.csv）。

    【替代式保存】：保存前会删除同 prefix 的旧目录（含历史时间戳目录），
    同一 prefix（同一药物+疾病）始终只保留最新一次分析，不会堆积。

    目录内容：
      - {prefix}_nodes.csv      节点拓扑指标表（含 is_hub 列）
      - {prefix}_edges.csv      边表（protein_A / protein_B / score ...）
      - {prefix}_network.html   交互式网络
      - {prefix}_graph.graphml  可编辑网络文件（供 Cytoscape / Gephi 精修）
      - {prefix}_network.svg    静态图矢量版（论文出版用）
      - {prefix}_enrichment.csv           富集结果表（有富集数据时）
      - {prefix}_enrichment_bar.png/.svg  富集柱状图（有富集数据时）
      - {prefix}_enrichment_dot.png/.svg  富集气泡图（有富集数据时）

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
        文件名/目录名前缀，如 "aspirin_cancer"（目录会追加 "_data" 后缀）。
    enrichment_df : pd.DataFrame | None
        富集结果表（统一列）；None 表示本次分析无富集数据。
    enrichment_bar_fig / enrichment_dot_fig : Figure | None
        富集柱状图 / 气泡图；None 表示无富集数据。

    Returns
    -------
    dict
        各文件路径与目录路径（含可选富集文件路径）。
    """
    import glob
    import shutil
    import networkx as nx
    import config

    os.makedirs(config.DOWNLOADS_DIR, exist_ok=True)

    # 替代式保存：删除同前缀的旧目录（含历史时间戳目录），只保留本次最新分析
    for old in glob.glob(os.path.join(config.DOWNLOADS_DIR, f"{prefix}_data*")):
        if os.path.isdir(old):
            shutil.rmtree(old)
    subdir = os.path.join(config.DOWNLOADS_DIR, f"{prefix}_data")
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

    # 富集文件（可选）
    if enrichment_df is not None and not enrichment_df.empty:
        paths["enrichment_csv"] = os.path.join(subdir, f"{prefix}_enrichment.csv")
        enrichment_df.to_csv(paths["enrichment_csv"], index=False, encoding="utf-8-sig")
    if enrichment_bar_fig is not None:
        paths["enrichment_bar_png"] = os.path.join(subdir, f"{prefix}_enrichment_bar.png")
        paths["enrichment_bar_svg"] = os.path.join(subdir, f"{prefix}_enrichment_bar.svg")
        enrichment_bar_fig.savefig(paths["enrichment_bar_png"], format="png", bbox_inches="tight")
        enrichment_bar_fig.savefig(paths["enrichment_bar_svg"], format="svg", bbox_inches="tight")
    if enrichment_dot_fig is not None:
        paths["enrichment_dot_png"] = os.path.join(subdir, f"{prefix}_enrichment_dot.png")
        paths["enrichment_dot_svg"] = os.path.join(subdir, f"{prefix}_enrichment_dot.svg")
        enrichment_dot_fig.savefig(paths["enrichment_dot_png"], format="png", bbox_inches="tight")
        enrichment_dot_fig.savefig(paths["enrichment_dot_svg"], format="svg", bbox_inches="tight")

    paths["directory"] = subdir
    return paths


def make_ppi_zip(file_paths: dict) -> bytes:
    """
    将一次分析的全部结果文件打包为 zip 字节（nodes/edges/html/graphml/svg
    及富集 CSV/图等所有已落盘的文件，排除 directory 键）。

    供 st.download_button 直接下载；文件在 zip 内使用各自的文件名，
    配合独立子目录，用户解压后即可获得该次分析的全部成品。
    """
    import io
    import zipfile

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, path in file_paths.items():
            if name == "directory":
                continue
            if path and os.path.exists(path):
                zf.write(path, arcname=os.path.basename(path))
    return buf.getvalue()


def make_dir_zip(directory: str) -> bytes:
    """
    将某个目录下的全部文件打包为 zip 字节（供管理页按目录整体下载）。

    zip 内使用各文件的文件名（不带目录前缀），返回空 zip 当目录不存在或为空。
    """
    import io
    import zipfile

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        if os.path.isdir(directory):
            for name in sorted(os.listdir(directory)):
                path = os.path.join(directory, name)
                if os.path.isfile(path):
                    zf.write(path, arcname=name)
    return buf.getvalue()


# ==================== 富集分析模块 ====================

# gseapy (Enrichr) 基因集库名与 g:Profiler sources 的映射关系：
#   g:Profiler source    →  Enrichr 库名
#   GO:BP                →  GO_Biological_Process_2021
#   GO:MF                →  GO_Molecular_Function_2021
#   GO:CC                →  GO_Cellular_Component_2021
#   KEGG                 →  KEGG_2021_Human
ENRICHR_LIB_MAP = {
    "GO:BP": "GO_Biological_Process_2021",
    "GO:MF": "GO_Molecular_Function_2021",
    "GO:CC": "GO_Cellular_Component_2021",
    "KEGG": "KEGG_2021_Human",
}

# g:Profiler 物种名 → Enrichr 物种名。
# 注意：Enrichr 只支持 human / mouse 等，大鼠（rnorvegicus）无对应库，
# 因此大鼠场景下降级 Enrichr 会被跳过，仅依赖 g:Profiler。
_ENRICHR_ORGANISM_MAP = {
    "hsapiens": "human",
    "mmusculus": "mouse",
}


def _enrichment_cache_key(genes: list[str], organism: str, sources, pvalue: float) -> str:
    """根据基因列表与查询参数生成富集缓存文件名（md5）：相同输入命中同一缓存文件。"""
    payload = (
        "\r".join(sorted(set(genes)))
        + f"|{organism}|{','.join(sources)}|{pvalue}"
    )
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def fetch_enrichment_gprofiler(
    genes: list[str],
    organism: str = "hsapiens",
    sources: tuple = ("GO:BP", "GO:MF", "GO:CC", "KEGG"),
    pvalue: float = 0.05,
    label: str | None = None,
) -> pd.DataFrame:
    """
    通过 g:Profiler 官方客户端（gprofiler-official）获取功能富集结果。

    结果缓存在 Static/data/cache/enrichment 下，文件名为 gProfiler_{label}_{md5}.csv
    （label 为来源标识，如 药物名_疾病名），缓存首行为 # 元数据注释头；
    命中缓存时直接读取本地文件，不访问远端 API。

    Parameters
    ----------
    genes : list[str]
        输入基因符号列表。
    organism : str
        g:Profiler 物种名，如 "hsapiens"。
    sources : tuple
        g:Profiler 数据源，如 ("GO:BP", "GO:MF", "GO:CC", "KEGG")。
    pvalue : float
        显著性阈值（g:Profiler 的 user_threshold，默认即 g_SCS 校正后阈值）。
    label : str | None
        来源标识（如 "aspirin_cancer" 或 "manual"），用于缓存文件名前缀，便于用户识别。

    Returns
    -------
    pd.DataFrame
        统一列富集结果（source / term_id / term_name / p_value /
        adjusted_p_value / intersection_size / query_size / genes）；
        请求失败或无显著结果时返回空 DataFrame。
    """
    import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes:
        return pd.DataFrame()
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    md5 = _enrichment_cache_key(genes, organism, sources, pvalue)
    name = f"gProfiler_{label}_{md5}.csv" if label else f"gProfiler_{md5}.csv"
    cache_file = os.path.join(config.ENRICHMENT_CACHE_DIR, name)
    if os.path.exists(cache_file):
        logging.info("Loading cached g:Profiler enrichment for %d genes", len(genes))
        return _read_cache_csv(cache_file)

    import time
    import socket

    # 请求保护：gprofiler-official 内部 requests.post 未传 timeout、也无重试，
    # 这里在调用层统一补上【请求间隔 + 失败重试 + socket 层超时兜底】。
    # socket.setdefaulttimeout 兜底库内部的 requests 调用（未显式传 timeout 时生效），
    # try/finally 确保调用结束后恢复，不影响其他请求。
    from gprofiler import GProfiler

    res = None
    for attempt in range(config.API_MAX_RETRIES):
        time.sleep(config.API_REQUEST_DELAY)  # 每次请求前等待，避免快速连发被限流
        socket.setdefaulttimeout(config.API_TIMEOUT)
        try:
            gp = GProfiler(return_dataframe=True)
            res = gp.profile(
                query=genes,
                organism=organism,
                sources=list(sources),
                user_threshold=pvalue,
            )
            break  # 成功
        except Exception:
            logging.exception(
                "g:Profiler enrichment failed (attempt %d/%d)",
                attempt + 1, config.API_MAX_RETRIES,
            )
            if attempt < config.API_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)  # 失败退避后重试
        finally:
            socket.setdefaulttimeout(None)

    if res is None or len(res) == 0:
        logging.info("g:Profiler: no significant terms (p < %s) for %d genes", pvalue, len(genes))
        return pd.DataFrame()

    df = _normalize_enrichment_gprofiler(res)
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    logging.info("g:Profiler enrichment for %d genes: %d terms cached", len(genes), len(df))
    return df


def _normalize_enrichment_gprofiler(df: pd.DataFrame) -> pd.DataFrame:
    """把 g:Profiler 返回的 DataFrame 整理为统一列。

    g:Profiler 默认已做 g_SCS 多重检验校正，p_value 即校正后 p 值，
    因此 adjusted_p_value 与 p_value 相同。
    """
    df = df.reset_index()  # return_dataframe=True 时 term_id 常在 index 中
    out = pd.DataFrame({
        "source": df["source"].astype(str),
        "term_id": df["term_id"].astype(str),
        "term_name": df["term_name"].astype(str),
        "p_value": pd.to_numeric(df["p_value"], errors="coerce"),
        "adjusted_p_value": pd.to_numeric(df["p_value"], errors="coerce"),
        "intersection_size": pd.to_numeric(df["intersection_size"], errors="coerce"),
        "query_size": pd.to_numeric(df["query_size"], errors="coerce"),
        "genes": df["intersection"].fillna(""),
    })
    out = out.dropna(subset=["term_id", "p_value"]).sort_values("p_value")
    return out.reset_index(drop=True)


def fetch_enrichment_enrichr(
    genes: list[str],
    libs: list[str],
    organism: str = "human",
    pvalue: float = 0.05,
    label: str | None = None,
) -> pd.DataFrame:
    """
    通过 gseapy（Enrichr 官方 wrapper）获取功能富集结果（作为 g:Profiler 的降级方案）。

    结果同样缓存在 Static/data/cache/enrichment 下，文件名为 Enrichr_{label}_{md5}.csv，
    缓存首行为 # 元数据注释头。

    Parameters
    ----------
    genes : list[str]
        输入基因符号列表。
    libs : list[str]
        Enrichr 基因集库名（见 ENRICHR_LIB_MAP 的映射）。
    organism : str
        Enrichr 物种名（"human" / "mouse"）。
    pvalue : float
        显著性阈值（cutoff）。
    label : str | None
        来源标识（如 "aspirin_cancer" 或 "manual"），用于缓存文件名前缀，便于用户识别。

    Returns
    -------
    pd.DataFrame
        统一列富集结果；请求失败或无显著结果时返回空 DataFrame。
    """
    import config

    genes = sorted({g.strip() for g in genes if g and g.strip()})
    if not genes or not libs:
        return pd.DataFrame()
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    md5 = _enrichment_cache_key(genes, organism, libs, pvalue)
    name = f"Enrichr_{label}_{md5}.csv" if label else f"Enrichr_{md5}.csv"
    cache_file = os.path.join(config.ENRICHMENT_CACHE_DIR, name)
    if os.path.exists(cache_file):
        logging.info("Loading cached Enrichr enrichment for %d genes", len(genes))
        return _read_cache_csv(cache_file)

    import time
    import socket

    # 请求保护：gseapy 内部 session 虽自带 5xx 重试，但部分请求 timeout=None
    # 无超时、且无请求间隔；这里在调用层统一补上【请求间隔 + 失败重试 + socket 层超时兜底】。
    # socket.setdefaulttimeout 兜底库内部的 requests 调用（未显式传 timeout 时生效），
    # try/finally 确保调用结束后恢复，不影响其他请求。
    import gseapy as gp

    res = None
    for attempt in range(config.API_MAX_RETRIES):
        time.sleep(config.API_REQUEST_DELAY)  # 每次请求前等待，避免快速连发被限流
        socket.setdefaulttimeout(config.API_TIMEOUT)
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=list(libs),
                organism=organism,
                outdir=None,  # 不写盘，结果由我们自行缓存
                cutoff=pvalue,
            )
            res = enr.results
            break  # 成功
        except Exception:
            logging.exception(
                "Enrichr enrichment failed (attempt %d/%d)",
                attempt + 1, config.API_MAX_RETRIES,
            )
            if attempt < config.API_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)  # 失败退避后重试
        finally:
            socket.setdefaulttimeout(None)

    if res is None or len(res) == 0:
        return pd.DataFrame()

    df = _normalize_enrichment_enrichr(res)
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    return df


def _normalize_enrichment_enrichr(df: pd.DataFrame) -> pd.DataFrame:
    """把 Enrichr 返回的 DataFrame 整理为统一列（并把库名反查回 g:Profiler 的 source 命名）。"""
    import re

    def _parse_overlap(s) -> tuple[int, int]:
        try:
            a, b = str(s).split("/")
            return int(a), int(b)
        except (ValueError, TypeError):
            return 0, 0

    def _extract_term_id(term) -> str:
        # Enrichr 的 Term 形如 "peptidyl-tyrosine modification (GO:0018212)"，提取末尾括号内的 ID
        m = re.search(r"\(([^()]+)\)\s*$", str(term))
        return m.group(1) if m else str(term)

    def _strip_suffix(term) -> str:
        return re.sub(r"\s*\([^()]*\)\s*$", "", str(term)).strip() or str(term)

    parsed = df["Overlap"].map(_parse_overlap)
    source_rev = {v: k for k, v in ENRICHR_LIB_MAP.items()}  # 库名 → g:Profiler source 命名
    out = pd.DataFrame({
        "source": df["Gene_set"].map(source_rev).fillna(df["Gene_set"]),
        "term_id": df["Term"].map(_extract_term_id),
        "term_name": df["Term"].map(_strip_suffix),
        "p_value": pd.to_numeric(df["P-value"], errors="coerce"),
        "adjusted_p_value": pd.to_numeric(df["Adjusted P-value"], errors="coerce"),
        "intersection_size": parsed.map(lambda t: t[0]),
        "query_size": parsed.map(lambda t: t[1]),
        "genes": df["Genes"].fillna(""),
    })
    out = out.dropna(subset=["term_id", "p_value"]).sort_values("p_value")
    return out.reset_index(drop=True)


def run_enrichment(
    genes: list[str],
    organism: str = "hsapiens",
    pvalue: float = 0.05,
    sources: tuple | None = None,
    label: str | None = None,
) -> tuple[pd.DataFrame, str | None]:
    """
    富集分析主入口：g:Profiler 优先，失败降级 Enrichr。

    Parameters
    ----------
    genes : list[str]
        输入基因符号列表（通常为药物∩疾病交集靶点）。
    organism : str
        g:Profiler 物种名（hsapiens / mmusculus / rnorvegicus）。
    pvalue : float
        显著性阈值。
    sources : tuple | None
        g:Profiler 数据源，默认取 config.ENRICHMENT_SOURCES（GO三分类 + KEGG）。
    label : str | None
        来源标识（如 "aspirin_cancer" 或 "manual"），用于缓存文件名前缀，便于用户识别。

    Returns
    -------
    tuple[pd.DataFrame, str | None]
        (统一列富集结果, 实际使用的数据源 "gProfiler" / "Enrichr")；
        两者都失败时返回 (空表, None)。
    """
    import config

    sources = sources or config.ENRICHMENT_SOURCES

    # 主：g:Profiler
    df = fetch_enrichment_gprofiler(genes, organism=organism, sources=sources, pvalue=pvalue, label=label)
    if not df.empty:
        return df, "gProfiler"

    # 备：Enrichr（不支持大鼠，跳过）
    enrichr_organism = _ENRICHR_ORGANISM_MAP.get(organism)
    if enrichr_organism:
        libs = [ENRICHR_LIB_MAP[s] for s in sources if s in ENRICHR_LIB_MAP]
        df = fetch_enrichment_enrichr(genes, libs=libs, organism=enrichr_organism, pvalue=pvalue, label=label)
        if not df.empty:
            return df, "Enrichr"
    return pd.DataFrame(), None


def render_enrichment_plots(df: pd.DataFrame) -> tuple:
    """
    绘制富集结果图：柱状图（每个 source 取前 N 条，-log10 校正 p）+ 气泡图（全局前 N 条）。

    Returns
    -------
    tuple
        (柱状图 Figure, 气泡图 Figure)；df 无有效行时返回 (None, None)。
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import config

    df = df[df["adjusted_p_value"].notna() & (df["adjusted_p_value"] > 0)].copy()
    if df.empty:
        return None, None
    df["neg_log10_p"] = -np.log10(df["adjusted_p_value"])

    source_colors = {
        "GO:BP": "#4C72B0", "GO:MF": "#DD8452",
        "GO:CC": "#55A868", "KEGG": "#C44E52",
    }

    # ---- 柱状图：每个 source 取校正 p 最小的前 N 条 ----
    top = (
        df.sort_values("adjusted_p_value")
        .groupby("source", as_index=False)
        .head(config.ENRICHMENT_TOP_PER_SOURCE)
    )
    fig_bar, ax_bar = plt.subplots(figsize=(9, max(5, 0.28 * len(top))))
    colors = [source_colors.get(s, "#7f7f7f") for s in top["source"]]
    ax_bar.barh(top["term_name"], top["neg_log10_p"], color=colors, height=0.6)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("-log10(adjusted p-value)")
    ax_bar.set_title("Functional Enrichment (top terms per source)")
    ax_bar.tick_params(axis="y", labelsize=7)
    fig_bar.tight_layout()

    # ---- 气泡图：全局校正 p 最小的前 N 条 ----
    bubble = df.sort_values("adjusted_p_value").head(config.PATHWAY_DISPLAY_NUM)
    sizes = 40 + 120 * (bubble["intersection_size"] / max(bubble["intersection_size"].max(), 1))
    fig_dot, ax_dot = plt.subplots(figsize=(10, max(6, 0.3 * len(bubble))))
    sc = ax_dot.scatter(
        bubble["neg_log10_p"], bubble["term_name"],
        s=sizes, c=bubble["adjusted_p_value"], cmap="viridis_r",
        edgecolors="grey", alpha=0.85,
    )
    cb = fig_dot.colorbar(sc, ax=ax_dot)
    cb.set_label("adjusted p-value")
    ax_dot.set_xlabel("-log10(adjusted p-value)")
    ax_dot.set_title("Enrichment bubble plot")
    ax_dot.tick_params(axis="y", labelsize=7)
    fig_dot.tight_layout()

    return fig_bar, fig_dot
