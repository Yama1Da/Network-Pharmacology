import os
import logging
import logging.handlers
import hashlib
from datetime import datetime

import pandas as pd
import requests
import streamlit as st

OPENTARGETS_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"


def _opentargets_post(query: str, variables: dict, timeout: int = 30) -> dict | None:
    """向 Open Targets GraphQL 端点发送 POST，带失败重试；最终失败返回 None。

    网络抖动 / 临时 5xx 时按退避重试（config.LIGHT_API_MAX_RETRIES 次），
    避免一次失败就让整个查询返回空表。GraphQL 层错误（payload["errors"]）
    不属于网络异常，由调用方自行判断。
    """
    import time
    import config

    for attempt in range(config.LIGHT_API_MAX_RETRIES):
        try:
            resp = requests.post(
                OPENTARGETS_API_URL,
                json={"query": query, "variables": variables},
                timeout=timeout,
            )
            resp.raise_for_status()  # 非 2xx 状态码会抛出异常
            return resp.json()
        except Exception:
            logging.exception(
                "Open Targets POST failed (attempt %d/%d)",
                attempt + 1, config.LIGHT_API_MAX_RETRIES,
            )
            if attempt < config.LIGHT_API_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)
    return None


# ==================== 日志系统 ====================


class _DailyRotatingFileHandler(logging.handlers.TimedRotatingFileHandler):
    """按天轮转的文件 handler（每天零点切文件）。

    继承 TimedRotatingFileHandler 并修复其清理逻辑：原版 getFilesToDelete()
    会把"当前正在写入的文件"（如 app.log）也计入待删除列表，在 Windows 上
    删除被占用的文件会抛 PermissionError，导致两年后的过期清理失效。
    这里在删除前排除当前文件，使 backupCount 清理只删真正的旧轮转文件。
    """

    def getFilesToDelete(self) -> list[str]:
        files = super().getFilesToDelete()
        current = os.path.basename(self.baseFilename)
        return [f for f in files if os.path.basename(f) != current]


def setup_logging() -> None:
    """初始化日志系统：控制台 + 每日文件，保留两年，幂等。

    - 日志目录 config.LOG_DIR（WY/Static/log）
    - 【按天轮转】TimedRotatingFileHandler 每天零点自动把日志写入新的
      轮转文件 app.log.YYYY-MM-DD（服务器不重启也能正确切分到当天文件）
    - 【保留两年】backupCount = LOG_RETENTION_DAYS（730），超过自动删除
      最旧的轮转文件，取代手写的 _cleanup_old_logs
    - 【轻量版】root logger 级别设为 INFO：屏蔽 matplotlib/PIL 等第三方库的
      DEBUG 噪音（字体查找、图片流），仅录项目自身 INFO 及以上日志；
      若后续需要精细控制（项目 DEBUG + 第三方 WARNING），再升级为
      "项目专用 logger（logging.getLogger('WY')）+ root 降级" 的方案
    - 幂等：Streamlit 每次交互都会重跑脚本，用 hasHandlers 守卫避免重复挂 handler
    - 目录创建失败时降级为仅控制台输出，不抛异常
    """
    import config

    logger = logging.getLogger()
    if logger.hasHandlers():
        return  # 已初始化过

    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(config.LOG_FORMAT, config.LOG_DATE_FORMAT)

    # 控制台 handler（实时排障）
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)
    logger.addHandler(console)

    # 文件 handler（按天轮转 + 保留两年自动清理）
    try:
        os.makedirs(config.LOG_DIR, exist_ok=True)
        log_file = os.path.join(config.LOG_DIR, "app.log")
        file_handler = _DailyRotatingFileHandler(
            log_file,
            when="midnight",                              # 每天零点切文件
            backupCount=config.LOG_RETENTION_DAYS,        # 保留两年，超期自动删除
            encoding="utf-8",
        )
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except OSError:
        # 目录/文件创建失败：降级为仅控制台，不影响应用
        logger.warning("日志文件创建失败，降级为仅控制台输出: %s", config.LOG_DIR)


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
        payload = _opentargets_post(search_query, {"queryString": disease_name.strip()})
        if payload is None:
            return pd.DataFrame()

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
        payload = _opentargets_post(disease_query, {"diseaseId": disease_id, "size": top_n})
        if payload is None:
            return pd.DataFrame()

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


def _acquire_disgenet_rate_limit_slot() -> None:
    """DisGeNET 主动限流：发请求前节流，保证低于官方配额（每分钟/每天）。

    用时间戳列表（滑动窗口）持久化到 config.DISGENET_RATE_LIMIT_STATE_FILE，
    Streamlit 每次 rerun 重跑脚本也能共享计数。窗口已满时：
      - 每分钟窗口：等待至最早一次请求滑出窗口（最多等 DISGENET_RATE_LIMIT_MAX_WAIT 秒）；
      - 每天窗口：直接抛异常，由调用方降级返回空表。
    """
    import json
    import time
    import config

    now = time.time()
    minute_window = 60.0
    day_window = 24 * 3600.0

    state_file = config.DISGENET_RATE_LIMIT_STATE_FILE
    timestamps: list[float] = []
    if os.path.exists(state_file):
        try:
            with open(state_file, "r", encoding="utf-8") as f:
                timestamps = [float(t) for t in json.load(f)]
        except Exception:
            logging.warning("Failed to load DisGeNET rate-limit state; resetting")

    # 只保留一天内的请求记录
    timestamps = [t for t in timestamps if now - t <= day_window]

    # 每日硬上限：等待可能数小时，直接拒绝
    if len(timestamps) >= config.DISGENET_RATE_LIMIT_PER_DAY:
        raise RuntimeError(
            "DisGeNET daily rate limit (%d) reached" % config.DISGENET_RATE_LIMIT_PER_DAY
        )

    # 每分钟窗口：满了则等到最早一次请求滑出窗口
    recent = [t for t in timestamps if t > now - minute_window]
    if len(recent) >= config.DISGENET_RATE_LIMIT_PER_MINUTE:
        wait = recent[0] + minute_window - now + 0.05
        if wait > config.DISGENET_RATE_LIMIT_MAX_WAIT:
            raise RuntimeError(
                "DisGeNET per-minute rate limit (%d) reached"
                % config.DISGENET_RATE_LIMIT_PER_MINUTE
            )
        logging.warning("DisGeNET per-minute limit reached, waiting %.1fs", wait)
        time.sleep(wait)
        now = time.time()
        timestamps = [t for t in timestamps if now - t <= day_window]

    timestamps.append(now)

    # 原子写入（先写临时文件再替换），避免并发写坏状态文件
    os.makedirs(os.path.dirname(state_file), exist_ok=True)
    tmp_file = state_file + ".tmp"
    with open(tmp_file, "w", encoding="utf-8") as f:
        json.dump(timestamps, f)
    os.replace(tmp_file, state_file)


def fetch_disease_targets_disgenet(disease_name: str) -> pd.DataFrame:
    """从 DisGeNET 新版 REST API 获取疾病关联基因（gene-disease association, GDA）。

    调用 gda summary 端点，按 disease 参数查询。注意当前 API.env 里的
    DisGeNET key 为 TRIAL 账号，存在以下限制：
      - 每页最多 30 条、不支持分页；
      - 最多查询 10 个标识符；
      - 请求频率限制严格（频繁调用会返回 429）。
    因此本函数只取第一页结果，并在遇到 429 时做有限次等待重试。

    Parameters
    ----------
    disease_name : str
        疾病英文名，如 "Alzheimer's Disease"。

    Returns
    -------
    pd.DataFrame
        列: approvedSymbol、entrez_id、ensembl_id、disease_name、
        disease_umls_cui、disgenet_score、normalized_score、num_pmids 等；
        无数据或失败时返回空表。

    缓存
    ----
    结果持久化为 CSV（Static/data/cache/disease/DisGeNET_{disease_name}.csv）。
    """
    import time
    import config

    cache_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache", "disease"
    )
    os.makedirs(cache_dir, exist_ok=True)
    csv_path = os.path.join(cache_dir, f"DisGeNET_{disease_name}.csv")

    # ---- cache hit: load local ----
    if os.path.exists(csv_path):
        logging.info("Loading cached DisGeNET targets for '%s'", disease_name)
        return pd.read_csv(csv_path)

    if not config.DisGeNET_API_KEY:
        logging.error("DisGeNET_API_KEY is empty; check API.env")
        return pd.DataFrame()

    headers = {"Authorization": config.DisGeNET_API_KEY, "accept": "application/json"}
    # DisGeNET gda summary 端点按疾病查询使用 `disease` 参数（值通常为疾病名或 UMLS CUI）
    params = {"disease": disease_name.strip()}

    try:
        # 每次发请求前先过主动限流（每分钟/每天），再请求；429 时按响应头有限重试。
        resp = None
        for _ in range(config.ENRICHMENT_MAX_RETRIES + 1):
            _acquire_disgenet_rate_limit_slot()
            resp = requests.get(
                config.DisGeNET_URL, params=params, headers=headers, timeout=config.API_TIMEOUT
            )
            if resp.status_code != 429:
                break

            retry_after = int(
                resp.headers.get("x-rate-limit-retry-after-seconds", config.API_RETRY_DELAY)
            )
            # 若需等待过久（超过最大等待阈值）则直接放弃，避免阻塞界面
            if retry_after > config.API_RETRY_DELAY * config.ENRICHMENT_MAX_RETRIES * 10:
                logging.error(
                    "DisGeNET rate limit requires %ss wait; aborting to avoid blocking", retry_after
                )
                return pd.DataFrame()
            logging.warning("DisGeNET rate-limited (429), retrying after %ss", retry_after)
            time.sleep(retry_after)

        resp.raise_for_status()
        payload = resp.json()

        if payload.get("status") != "OK":
            logging.error("DisGeNET query failed: %s", payload.get("payload"))
            return pd.DataFrame()

        rows = payload.get("payload", [])
        if not isinstance(rows, list) or not rows:
            logging.warning("DisGeNET returned no associations for '%s'", disease_name)
            return pd.DataFrame()

        records: list[dict] = []
        for r in rows:
            ensembl_ids = r.get("geneEnsemblIDs") or []
            records.append({
                "approvedSymbol": r.get("symbolOfGene"),
                "entrez_id": r.get("geneNcbiID"),
                "ensembl_id": ensembl_ids[0] if ensembl_ids else None,
                "disease_name": r.get("diseaseName"),
                "disease_umls_cui": r.get("diseaseUMLSCUI"),
                "disgenet_score": r.get("score"),
                "normalized_score": r.get("normalized_score"),
                "num_pmids": r.get("numPMIDs"),
                "year_initial": r.get("yearInitial"),
                "year_final": r.get("yearFinal"),
            })

        df = (
            pd.DataFrame(records)
            .drop_duplicates(subset=["approvedSymbol"])
            .reset_index(drop=True)
        )
        df.to_csv(csv_path, index=False, encoding="utf-8-sig")
        logging.info("DisGeNET '%s': %d targets fetched", disease_name, len(df))
        return df

    except Exception:
        logging.exception("fetch_disease_targets_disgenet failed for '%s'", disease_name)
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

    cache_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "Static", "data", "cache", "drug"
    )
    os.makedirs(cache_dir, exist_ok=True)
    csv_path = os.path.join(cache_dir, f"DGIdb_{drug_name}.csv")

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
        # 跳过更新过程中产生的备份/临时文件（.bak / .tmp），避免残留虚影显示在页面上
        if name.endswith((".bak", ".tmp")):
            continue
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


def list_log_files() -> pd.DataFrame:
    """列出 LOG_DIR 下的日志文件及其创建/修改时间。

    兼容两种命名：
      - 旧格式：app_YYYY-MM-DD.log
      - TimedRotating 轮转格式：app.log（当前活跃）、app.log.YYYY-MM-DD（历史）

    返回含 "File name"、"Created time"、"Modified time" 三列；
    其中 app.log 为当前正在写入的文件（即"今天"的日志）。

    Returns
    -------
    pd.DataFrame
        含 "File name"、"Created time"、"Modified time" 列；
        目录不存在或无日志文件时返回空表。
    """
    import config

    if not os.path.isdir(config.LOG_DIR):
        return pd.DataFrame(columns=["File name", "Created time", "Modified time"])
    rows = []
    for name in sorted(os.listdir(config.LOG_DIR)):
        path = os.path.join(config.LOG_DIR, name)
        is_log = name.endswith(".log") or name.startswith("app.log.")
        if os.path.isfile(path) and is_log:
            created = datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")
            modified = datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append({"File name": name, "Created time": created, "Modified time": modified})
    return pd.DataFrame(rows, columns=["File name", "Created time", "Modified time"])


# ==================== PPI 网络模块 ====================


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
    merged_path = os.path.join(config.MERGED_CACHE_DIR, f"merged_drug_{drug_name}.csv")
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
    merged_path = os.path.join(config.MERGED_CACHE_DIR, f"merged_disease_{disease_name}.csv")
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


def _is_cache_fresh(cache_file: str, ttl_seconds: int | None = None) -> bool:
    """缓存文件是否仍有效：存在且修改时间未超过 TTL。

    ttl_seconds 缺省（None）时取 config.CACHE_TTL_SECONDS。
    """
    import time
    import config

    if not os.path.exists(cache_file):
        return False
    if ttl_seconds is None:
        ttl_seconds = config.CACHE_TTL_SECONDS
    try:
        return (time.time() - os.path.getmtime(cache_file)) < ttl_seconds
    except OSError:
        return False


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

    # ---- cache hit: load local (若超过 TTL 视为过期，重新请求) ----
    if _is_cache_fresh(cache_file):
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


# ---------- g:Profiler 熔断（circuit breaker） ----------
# g:Profiler 服务器不可达时，单次连接可能空等 20s+ 才超时（socket.setdefaulttimeout
# 对库内部 requests 调用不生效），重试次数越多拖得越久。因此失败时写一个带时间戳的
# 标记文件，TTL 内直接跳过 g:Profiler 走 Enrichr；TTL 过期后自动重试，服务恢复即复用。


def _gprofiler_circuit_file() -> str:
    """g:Profiler 熔断标记文件路径（内容为上次失败时间戳）。"""
    import config

    return os.path.join(config.ENRICHMENT_CACHE_DIR, ".gprofiler_unreachable")


def _gprofiler_circuit_open() -> bool:
    """熔断是否开启：失败标记存在且未超过 TTL。"""
    import config
    import time

    path = _gprofiler_circuit_file()
    if not os.path.exists(path):
        return False
    try:
        with open(path, "r", encoding="utf-8") as f:
            failed_at = float(f.read().strip())
    except Exception:
        return False
    return (time.time() - failed_at) < config.GPROFILER_CIRCUIT_TTL_SECONDS


def _gprofiler_circuit_open_trip() -> None:
    """记录 g:Profiler 本次失败时间（开启熔断）。"""
    import time

    path = _gprofiler_circuit_file()
    try:
        with open(path, "w", encoding="utf-8") as f:
            f.write(str(time.time()))
    except Exception:
        logging.exception("Failed to write g:Profiler circuit breaker file")


def _gprofiler_circuit_close() -> None:
    """g:Profiler 请求成功时清除熔断标记。"""
    path = _gprofiler_circuit_file()
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception:
        logging.exception("Failed to remove g:Profiler circuit breaker file")


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
    if _is_cache_fresh(cache_file):
        logging.info("Loading cached g:Profiler enrichment for %d genes", len(genes))
        return _read_cache_csv(cache_file)

    # 熔断：gProfiler 近期失败过则直接跳过，不再空等连接超时（走 Enrichr 兜底）
    if _gprofiler_circuit_open():
        logging.info("g:Profiler circuit breaker open — skipping, Enrichr fallback")
        return pd.DataFrame()

    import time
    import socket

    # 请求保护：gprofiler-official 内部 requests.post 未传 timeout、也无重试，
    # 这里在调用层统一补上【请求间隔 + 失败重试 + socket 层超时兜底】。
    # socket.setdefaulttimeout 兜底库内部的 requests 调用（未显式传 timeout 时生效），
    # try/finally 确保调用结束后恢复，不影响其他请求。
    from gprofiler import GProfiler

    res = None
    for attempt in range(config.ENRICHMENT_MAX_RETRIES):
        time.sleep(config.API_REQUEST_DELAY)  # 每次请求前等待，避免快速连发被限流
        socket.setdefaulttimeout(config.API_TIMEOUT)
        try:
            gp = GProfiler(return_dataframe=True)
            res = gp.profile(
                query=genes,
                organism=organism,
                sources=list(sources),
                user_threshold=pvalue,
                no_evidences=False,  # 需要返回交集基因（intersections），供 genes 列使用
            )
            break  # 成功
        except Exception:
            logging.exception(
                "g:Profiler enrichment failed (attempt %d/%d)",
                attempt + 1, config.ENRICHMENT_MAX_RETRIES,
            )
            _gprofiler_circuit_open_trip()  # 记录失败时间，开启熔断
            if attempt < config.ENRICHMENT_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)  # 失败退避后重试
        finally:
            socket.setdefaulttimeout(None)

    if res is None or len(res) == 0:
        logging.info("g:Profiler: no significant terms (p < %s) for %d genes", pvalue, len(genes))
        return pd.DataFrame()

    _gprofiler_circuit_close()  # 请求成功，关闭熔断
    df = _normalize_enrichment_gprofiler(res)
    df.to_csv(cache_file, index=False, encoding="utf-8-sig")
    logging.info("g:Profiler enrichment for %d genes: %d terms cached", len(genes), len(df))
    return df


def _normalize_enrichment_gprofiler(df: pd.DataFrame) -> pd.DataFrame:
    """把 g:Profiler 返回的 DataFrame 整理为统一列。

    g:Profiler 默认已做 g_SCS 多重检验校正，p_value 即校正后 p 值，
    因此 adjusted_p_value 与 p_value 相同。

    兼容 gprofiler-official 0.3.x 与 1.x 的列名差异：
    0.3.x 返回 term_id / term_name / intersection；
    1.x 返回 native / name / intersections（需 no_evidences=False）。
    """
    df = df.reset_index()  # return_dataframe=True 时 term_id 常在 index 中

    def _pick(*names) -> pd.Series:
        """按优先级取第一个存在的列，都不存在时返回全 NA 序列。"""
        for n in names:
            if n in df.columns:
                return df[n]
        return pd.Series(pd.NA, index=df.index, dtype=object)

    out = pd.DataFrame({
        "source": df["source"].astype(str),
        "term_id": _pick("term_id", "native"),
        "term_name": _pick("term_name", "name"),
        "p_value": pd.to_numeric(df["p_value"], errors="coerce"),
        "adjusted_p_value": pd.to_numeric(df["p_value"], errors="coerce"),
        "intersection_size": pd.to_numeric(df["intersection_size"], errors="coerce"),
        "query_size": pd.to_numeric(df["query_size"], errors="coerce"),
        "genes": _pick("intersection", "intersections").fillna(""),
    })
    out = out.dropna(subset=["term_id", "p_value"]).sort_values("p_value")
    out["term_id"] = out["term_id"].astype(str)
    out["term_name"] = out["term_name"].astype(str)
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
    if _is_cache_fresh(cache_file):
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
    for attempt in range(config.ENRICHMENT_MAX_RETRIES):
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
                attempt + 1, config.ENRICHMENT_MAX_RETRIES,
            )
            if attempt < config.ENRICHMENT_MAX_RETRIES - 1:
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


# ==================== 分子对接模块（AutoDock Vina v1.2.7） ====================
# 配体 pdbqt 用 meeko（确定性）；受体 pdbqt 用 openbabel -xr（去水/加氢 pH7.4/AD typing）；
# 配体 3D SDF 统一 pubchempy.get_sdf(record_type="3d")。
# 缓存策略：受体/配体 pdbqt 与对接盒是确定性中间物 → 安全缓存（TTL 见
# config.DOCKING_CACHE_TTL_DAYS）；vina 输出（best pose）是随机搜索产物 → 只归档不缓存。


def _docking_fingerprint(path: str) -> str:
    """文件内容 MD5 指纹（前 12 位），用于 pdbqt 缓存命中判断。

    结构变化（不同 PDB 版本 / 不同 SDF 内容）指纹即变化，保证缓存不串货。
    """
    try:
        with open(path, "rb") as f:
            return hashlib.md5(f.read()).hexdigest()[:12]
    except OSError:
        return ""


def _unique_path(path: str) -> str:
    """防覆盖：目标路径已存在时追加时间戳（必要时再加序号），绝不覆写固定名。

    对接成果/元数据文件名都应经此函数生成，避免 open(..., "w") 直接覆盖旧结果。
    """
    if not os.path.exists(path):
        return path
    stem, ext = os.path.splitext(path)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = f"{stem}_{ts}{ext}"
    n = 1
    while os.path.exists(candidate):
        candidate = f"{stem}_{ts}_{n}{ext}"
        n += 1
    return candidate


def _strip_receptor_pdb(pdb_path: str, out_path: str) -> None:
    """文本层清洗受体 PDB：仅保留蛋白 ATOM/TER/END。

    移除 HETATM（水、配体、离子——apo 受体对接不需要）并去重 altloc
    （只保留主构象，即 altLoc 为空或 'A'）。清洗失败抛 RuntimeError。
    """
    seen = set()
    lines = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM":
                if line[16] not in (" ", "A"):
                    continue  # 跳过 altloc 冗余构象
                key = (line[21], line[22:26].strip(), line[17:20], line[12:16].strip())
                if key in seen:
                    continue
                seen.add(key)
                lines.append(line.rstrip("\n"))
            elif rec in ("TER", "END"):
                lines.append(line.rstrip("\n"))
    if not lines:
        raise RuntimeError(f"no ATOM records in receptor PDB: {pdb_path}")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _pdb_coords(pdb_path: str) -> tuple[list, list]:
    """解析 PDB 原子坐标（单位 Å）。

    Returns
    -------
    tuple[list, list]
        (蛋白 ATOM 坐标列表, 非水 HETATM 坐标列表)；每项为 (x, y, z) float 元组。
    """
    protein, het = [], []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if line[17:20].strip() in ("HOH", "H2O", "WAT", "DOD"):
                continue  # 水分子不参与结合位点定位
            try:
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            (het if rec == "HETATM" else protein).append(coord)
    return protein, het


# -------------------- 配体链 --------------------


def fetch_ligand_sdf_3d(cid: int, out_dir: str | None = None) -> str:
    """通过 pubchempy 获取含 3D 坐标的配体 SDF 并落盘。

    Parameters
    ----------
    cid : int
        PubChem Compound ID（配体来源方案 A 的现成 3D 结构）。
    out_dir : str | None
        输出目录，默认 config.DOCKING_RAW_CACHE_DIR（原始文件，随 TTL 清理）。

    Returns
    -------
    str
        SDF 文件路径；PubChem 无此 CID / 无 3D 结构时返回空串。
    """
    import config
    import pubchempy

    out_dir = out_dir or config.DOCKING_RAW_CACHE_DIR
    os.makedirs(out_dir, exist_ok=True)
    # 确定性中间物：同 CID 同名幂等覆盖（不使用 _unique_path，避免缓存堆积）。
    # 已存在且非空时直接复用，不重复下载；过期由 clean_docking_cache 按 TTL 清理。
    path = os.path.join(out_dir, f"ligand_{cid}.sdf")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        logging.info("Reusing cached 3D SDF for CID %d: %s", cid, path)
        return path
    sdf_text = pubchempy.get_sdf(cid, namespace="cid", record_type="3d")
    if not sdf_text:
        logging.warning("pubchempy: no 3D SDF for CID %s", cid)
        return ""
    with open(path, "w", encoding="utf-8") as f:
        f.write(sdf_text)
    logging.info("Fetched 3D SDF for CID %d -> %s", cid, path)
    return path


def prepare_ligand_pdbqt(sdf_path: str, out_path: str | None = None) -> str:
    """用 meeko 将配体 3D SDF 转为对接用 pdbqt（确定性中间物，可缓存）。

    Parameters
    ----------
    sdf_path : str
        含 3D 坐标的 SDF（meeko 要求分子带 conformer）。
    out_path : str | None
        输出路径；默认写入 config.DOCKING_LIGAND_CACHE_DIR，文件名带内容指纹。

    Returns
    -------
    str
        生成的 pdbqt 文件路径。

    Raises
    ------
    RuntimeError
        rdkit 读 SDF 失败 / meeko 原子类型分配或写出失败时抛出。
    """
    import config
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem

    mol = None
    try:
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False)  # 默认 sanitize，保证键序/形式电荷正确
    except Exception:
        logging.warning("ligand SDF sanitize failed, fallback to raw parse: %s", sdf_path)
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"rdkit failed to read ligand SDF: {sdf_path}")
    if mol.GetNumConformers() == 0:
        raise RuntimeError(f"ligand SDF has no 3D conformer: {sdf_path}（需 record_type=3d）")
    mol = Chem.AddHs(mol)  # 补显式氢，保证原子类型分配完整

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
    if not is_ok:
        raise RuntimeError(f"meeko failed to write ligand pdbqt: {error_msg}")

    if not out_path:
        fingerprint = _docking_fingerprint(sdf_path)
        out_path = os.path.join(config.DOCKING_LIGAND_CACHE_DIR, f"ligand_{fingerprint}.pdbqt")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(pdbqt_string)
    logging.info("Prepared ligand pdbqt: %s", out_path)
    return out_path


# -------------------- 受体链 --------------------


def fetch_receptor_pdb_candidates(gene_symbol: str, max_hits: int = 50) -> pd.DataFrame:
    """按基因符号查询 RCSB PDB 候选蛋白结构。

    实现说明：RCSB Search API 已不再对 gene/organism 类 attribute 开放 text 检索
    （"search is not enabled"），故改用 full-text 精确短语匹配基因符号，再经
    RCSB Data API 并发补充每个条目的分辨率/方法/标题。

    Returns
    -------
    pd.DataFrame
        列：pdb_id / resolution / method / title；查询失败或无结果返回空表。
    """
    import config

    # 1) RCSB Search API：full-text 精确短语匹配基因符号
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": f'"{gene_symbol.strip().upper()}"'},
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    try:
        resp = requests.post(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            json=query, timeout=config.API_TIMEOUT,
            headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"},
        )
        resp.raise_for_status()
        hits = resp.json().get("result_set", [])[:max_hits]
    except Exception as exc:
        # 搜索失败必须向上抛（带原因），供页面提示；200 空 result_set 才是"无结构"
        logging.exception("RCSB search failed for gene %s", gene_symbol)
        raise RuntimeError(f"RCSB 搜索失败（{gene_symbol}）：{exc}") from exc

    # 2) RCSB Data API：并发补充元数据（分辨率/方法/标题）
    def _fetch_meta(entry_id: str) -> tuple[str, float | None, str, str]:
        try:
            d = requests.get(
                f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}",
                timeout=15, headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"},
            ).json()
            res = d.get("rcsb_entry_info", {}).get("resolution_combined")
            res = res[0] if isinstance(res, list) and res else res
            method = ", ".join(m.get("method", "") for m in (d.get("exptl") or [])) or ""
            title = (d.get("struct") or {}).get("title", "") or ""
            return entry_id, res, method, title
        except Exception:
            return entry_id, None, "", ""

    ids = [hit.get("identifier", "") for hit in hits if hit.get("identifier")]
    meta = {}
    if ids:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=8) as pool:
            for eid, res, method, title in pool.map(_fetch_meta, ids):
                meta[eid] = (res, method, title)

    rows = []
    for eid in ids:
        res, method, title = meta.get(eid, (None, "", ""))
        rows.append({"pdb_id": eid, "resolution": res, "method": method, "title": title})
    return pd.DataFrame(rows)


def prepare_receptor_pdbqt(pdb_path: str, out_path: str | None = None) -> str:
    """将受体 PDB 转为对接用 pdbqt（确定性中间物，可缓存）。

    处理链：文本层清洗（去水/去 HETATM 配体/altloc 去重）→
    openbabel 受体模式（加氢 pH7.4、AutoDock 原子类型分配）。

    Parameters
    ----------
    pdb_path : str
        受体原始 PDB 文件路径（如 RCSB 下载的 {pdb_id}.pdb）。
    out_path : str | None
        输出路径；默认写入 config.DOCKING_RECEPTOR_CACHE_DIR，文件名带内容指纹。

    Returns
    -------
    str
        生成的 pdbqt 文件路径。

    Raises
    ------
    RuntimeError
        openbabel 读 PDB 或写 pdbqt 失败时抛出。
    """
    import config
    from openbabel import openbabel as ob

    os.makedirs(config.DOCKING_RAW_CACHE_DIR, exist_ok=True)
    cleaned = os.path.join(
        config.DOCKING_RAW_CACHE_DIR,
        f"receptor_clean_{os.path.splitext(os.path.basename(pdb_path))[0]}.pdb",
    )
    try:
        _strip_receptor_pdb(pdb_path, cleaned)
        source = cleaned
    except Exception:
        logging.exception("receptor PDB cleaning failed, fallback to raw PDB")
        source = pdb_path

    if not out_path:
        fingerprint = _docking_fingerprint(source)
        out_path = os.path.join(config.DOCKING_RECEPTOR_CACHE_DIR, f"receptor_{fingerprint}.pdbqt")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdbqt")
    conv.SetOptions("r", ob.OBConversion.OUTOPTIONS)  # 受体模式：加氢 / 去水 / AD typing
    mol = ob.OBMol()
    if not conv.ReadFile(mol, source):
        raise RuntimeError(f"openbabel failed to read receptor PDB: {source}")
    if not conv.WriteFile(mol, out_path):
        raise RuntimeError(f"openbabel failed to write receptor pdbqt: {out_path}")
    logging.info("Prepared receptor pdbqt: %s", out_path)
    return out_path


# -------------------- 对接盒 --------------------


def infer_box_from_cofactor(pdb_path: str, padding: float = 8.0) -> dict:
    """依据 PDB 中的配体/辅因子（HETATM，非水）坐标推断对接盒。

    PDB 含 HETATM 时以其质心/包围盒定位结合位点；无 HETATM（apo 结构）
    时退化为蛋白整体包围盒。返回值为科研决策点参考，仍需用户人工确认。

    Parameters
    ----------
    pdb_path : str
        PDB 文件路径。
    padding : float
        结合位点包围盒外扩距离（Å），默认 8.0（覆盖口袋残基）。

    Returns
    -------
    dict
        {"center": (x, y, z), "size": (sx, sy, sz)}，单位 Å，保留 4 位小数。
    """
    protein, het = _pdb_coords(pdb_path)
    if not protein and not het:
        raise RuntimeError(f"no coordinates found in PDB: {pdb_path}")
    source = het if het else protein  # 优先配体/辅因子定位口袋
    xs, ys, zs = zip(*source)
    center = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    size = (
        max(xs) - min(xs) + 2 * padding,
        max(ys) - min(ys) + 2 * padding,
        max(zs) - min(zs) + 2 * padding,
    )
    return {
        "center": tuple(round(v, 4) for v in center),
        "size": tuple(round(v, 4) for v in size),
    }


# -------------------- 执行链 --------------------


def run_vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    box: dict,
    out_path: str | None = None,
) -> "subprocess.CompletedProcess":
    """subprocess 调用 AutoDock Vina 执行对接（随机搜索）。

    Parameters
    ----------
    receptor_pdbqt / ligand_pdbqt : str
        预处理后的 pdbqt 路径。
    box : dict
        {"center": (x,y,z), "size": (sx,sy,sz)}，由 infer_box_from_cofactor 得到。
    out_path : str | None
        vina 输出 pdbqt 路径；默认写入 config.MD_RESULT_DIR 且自动防覆盖。

    Returns
    -------
    subprocess.CompletedProcess
        包含 stdout/stderr（vina 打分 log 在 stdout，供前端展示）。

    Raises
    ------
    RuntimeError
        vina 超时（config.DOCKING_TIMEOUT）或非零退出时抛出。
    """
    import subprocess
    import config

    if not os.path.isfile(config.VINA_EXE):
        raise RuntimeError(
            f"未找到 AutoDock Vina 可执行文件：{config.VINA_EXE}。"
            "请设置环境变量 VINA_DIR 指向 vina 所在目录，或将 vina 加入 PATH。")

    out_path = out_path or _unique_path(os.path.join(config.MD_RESULT_DIR, "vina_out.pdbqt"))
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    cmd = [
        config.VINA_EXE,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(config.EXHAUSTIVENESS),
        "--num_modes", str(config.NUM_MODES),
        "--cpu", str(config.DOCKING_CPU),
        "--seed", str(config.DOCKING_SEED),
        "--out", out_path,
    ]
    logging.info("Running vina: %s", " ".join(cmd))
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=config.DOCKING_TIMEOUT,
            encoding="utf-8", errors="replace",
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"vina timed out after {config.DOCKING_TIMEOUT}s")
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "")[-500:]
        raise RuntimeError(f"vina failed (rc={proc.returncode}): {tail}")
    return proc


def _split_models_python(out_pdbqt: str, out_dir: str) -> list[str]:
    """纯 Python 按 MODEL/ENDMDL 拆分 vina 多模型输出（vina_split 不可用时的降级）。

    文件名格式：{输入名}_{序号}.pdbqt（vina_split 的替代，不依赖外部 exe）。
    """
    stem = os.path.splitext(os.path.basename(out_pdbqt))[0]
    blocks = []
    current = None
    with open(out_pdbqt, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                current = [line]
            elif line.startswith("ENDMDL") and current is not None:
                current.append(line)
                blocks.append(current)
                current = None
            elif current is not None:
                current.append(line)
    files = []
    for i, block in enumerate(blocks, 1):
        path = os.path.join(out_dir, f"{stem}_{i}.pdbqt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("".join(block))
        files.append(path)
    return files


def split_vina(out_pdbqt: str, out_dir: str | None = None) -> list[str]:
    """将 vina 多模型输出拆分为独立 pose 文件（每个结合模式一个 pdbqt）。

    优先调用 vina_split exe；若 exe 缺失/损坏（如 0 字节）或执行失败，
    自动降级为纯 Python 按 MODEL 拆分，保证流程不被外部工具阻塞。

    Parameters
    ----------
    out_pdbqt : str
        run_vina 的输出文件。
    out_dir : str | None
        拆分产物输出目录；默认与 out_pdbqt 同目录。

    Returns
    -------
    list[str]
        拆分出的 pose 文件路径（按文件名排序）。

    Raises
    ------
    RuntimeError
        vina_split 与纯 Python 拆分均失败时抛出。
    """
    import subprocess
    import glob
    import config

    out_dir = out_dir or os.path.dirname(out_pdbqt)
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(out_pdbqt))[0]

    # 优先：vina_split exe（须是有效 PE 文件，非 0 字节）
    split_exe = config.VINA_SPLIT_EXE
    if os.path.isfile(split_exe) and os.path.getsize(split_exe) > 0:
        cmd = [split_exe, "--input", out_pdbqt]
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, cwd=out_dir,
                timeout=config.DOCKING_TIMEOUT, encoding="utf-8", errors="replace",
            )
            if proc.returncode == 0:
                files = sorted(glob.glob(os.path.join(out_dir, f"{stem}*")))
                if files:
                    return files
        except OSError:
            logging.warning("vina_split unavailable (%s), fallback to python split", split_exe)

    # 降级：纯 Python 拆分
    files = _split_models_python(out_pdbqt, out_dir)
    if not files:
        raise RuntimeError(f"failed to split vina output: {out_pdbqt}")
    return files


def parse_vina_scores(out_pdbqt: str) -> pd.DataFrame:
    """解析 vina 输出 pdbqt 中各结合模式的打分与 RMSD。

    Vina 在每个 MODEL 块写入 "REMARK VINA RESULT <affinity> <rmsd_lb> <rmsd_ub>"，
    本函数按此提取并按亲和力升序排序（mode 1 即最佳 pose）。

    Returns
    -------
    pd.DataFrame
        列：mode / affinity / rmsd_lb / rmsd_ub；无可解析模式时返回空表。
    """
    modes = []
    current = None
    with open(out_pdbqt, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                current = {}
            elif line.startswith("REMARK VINA RESULT") and current is not None:
                parts = line.split()
                if len(parts) >= 5:
                    current["affinity"] = float(parts[3])
                    current["rmsd_lb"] = float(parts[4])
                    current["rmsd_ub"] = float(parts[5]) if len(parts) > 5 else None
            elif line.startswith("ENDMDL") and current:
                modes.append(current)
                current = None
    df = pd.DataFrame(modes)
    if df.empty:
        return pd.DataFrame(columns=["mode", "affinity", "rmsd_lb", "rmsd_ub"])
    df.insert(0, "mode", range(1, len(df) + 1))
    return df.sort_values("affinity").reset_index(drop=True)


# -------------------- 元数据 / 归档 / 清理 --------------------


def write_docking_metadata(meta_path: str, fields: dict) -> str:
    """写对接元数据文本（# 注释行，对齐富集/PPI 缓存的 # 头风格）。

    Parameters
    ----------
    meta_path : str
        目标路径（自动防覆盖，同名时追加时间戳）。
    fields : dict
        键值对，按序写为 "# key: value"；时间戳自动追加到末尾。

    Returns
    -------
    str
        实际写入的文件路径。
    """
    lines = [f"# {key}: {value}" for key, value in fields.items()]
    lines.append(f"# timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    path = _unique_path(meta_path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    logging.info("Wrote docking metadata: %s", path)
    return path


def star_docking_result(result_dir: str) -> str:
    """将成果目录完整复制到收藏区（MD_STAR_DIR/star），长期保留（仅可手动删除）。

    收藏副本与原成果互不影响：原成果仍受 prune_docking_results 版本控制，
    收藏副本独立存在于 star/ 下，不受自动删除影响。

    Parameters
    ----------
    result_dir : str
        本次对接成果目录（MD_RESULT_DIR 下，时间戳命名）。

    Returns
    -------
    str
        收藏区中的目录路径（自动防覆盖）。
    """
    import shutil
    import config

    os.makedirs(config.MD_STAR_DIR, exist_ok=True)
    dest = _unique_path(os.path.join(config.MD_STAR_DIR, os.path.basename(result_dir)))
    shutil.copytree(result_dir, dest)
    logging.info("Starred docking result: %s -> %s", result_dir, dest)
    return dest


def prune_docking_results(result_dir: str, max_versions: int | None = None) -> list[str]:
    """成果版本控制：同条件（drug_disease_gene）目录超过上限时，直接删除最旧版本。

    目录名格式约定：{prefix}_{YYYYMMDD_HHMMSS}，prefix 如 aspirin_cancer_BRCA1。

    Parameters
    ----------
    result_dir : str
        本次对接成果目录（MD_RESULT_DIR 下）。
    max_versions : int | None
        最多保留版本数，默认 config.DOCKING_MAX_VERSIONS（3）。

    Returns
    -------
    list[str]
        本次实际删除的旧版目录路径列表；无需删除返回空列表。
    """
    import re
    import shutil
    import glob
    import config

    max_versions = max_versions or config.DOCKING_MAX_VERSIONS
    parent = os.path.dirname(result_dir)
    basename = os.path.basename(result_dir)
    prefix = re.sub(r"_\d{8}_\d{6}$", "", basename)  # 剥离末尾时间戳段
    if prefix == basename:
        return []  # 目录名不符合命名约定，不做版本控制
    siblings = sorted(
        d for d in glob.glob(os.path.join(parent, prefix + "_*"))
        if os.path.isdir(d)
    )
    removed = []
    while len(siblings) > max_versions:
        oldest = siblings.pop(0)
        shutil.rmtree(oldest, ignore_errors=True)
        removed.append(oldest)
        logging.info("Pruned old docking result: %s", oldest)
    return removed


def clean_docking_cache(ttl_days: int | None = None) -> int:
    """清理 docking 缓存目录下超过 TTL 的中间物文件。

    TTL 分级：
      - raw 子目录（原始 SDF/PDB/vina 输出，可随时重新下载生成）→
        config.DOCKING_RAW_TTL_DAYS（24 小时）
      - ligand / receptor / box（加工后的 pdbqt 等确定性产物）→
        config.DOCKING_CACHE_TTL_DAYS（7 天）

    Parameters
    ----------
    ttl_days : int | None
        显式指定时对全部目录统一使用该 TTL（便于手动全清，如传 0）；
        默认 None 按上面的分级策略清理。

    Returns
    -------
    int
        删除的文件数。
    """
    import time
    import config

    if ttl_days is not None:
        cache_ttl = raw_ttl = ttl_days
    else:
        cache_ttl = config.DOCKING_CACHE_TTL_DAYS
        raw_ttl = config.DOCKING_RAW_TTL_DAYS
    removed = 0
    if not os.path.isdir(config.DOCKING_CACHE_DIR):
        return 0
    raw_dir = config.DOCKING_RAW_CACHE_DIR
    for root, _, files in os.walk(config.DOCKING_CACHE_DIR):
        for name in files:
            path = os.path.join(root, name)
            ttl = raw_ttl if root == raw_dir else cache_ttl
            cutoff = time.time() - ttl * 86400
            try:
                if os.path.getmtime(path) < cutoff:
                    os.remove(path)
                    removed += 1
            except OSError:
                logging.warning("clean_docking_cache: cannot remove %s", path)
    if removed:
        logging.info("Docking cache cleaned: %d files removed", removed)
    return removed
