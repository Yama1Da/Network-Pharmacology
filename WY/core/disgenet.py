"""DisGeNET 数据源：疾病关联基因（GDA），含主动限流（每分钟/每天滑动窗口）。"""

import logging
import os

import pandas as pd
import requests

from .paths import BASE_DIR
from .text_utils import safe_name


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

    cache_dir = os.path.join(BASE_DIR, "Static", "data", "cache", "disease")
    os.makedirs(cache_dir, exist_ok=True)
    csv_path = os.path.join(cache_dir, f"DisGeNET_{safe_name(disease_name)}.csv")

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
