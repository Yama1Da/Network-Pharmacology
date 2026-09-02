"""DisGeNET 数据源：查询疾病的关联基因（GDA），并按官方配额做主动限流。"""

import logging
import os

import pandas as pd
import requests

from .caching import open_csv_cache, save_csv_cache
from .text_utils import safe_name


def _acquire_disgenet_rate_limit_slot() -> None:
    """DisGeNET 主动限流：发请求前检查官方配额（每分钟/每天），避免超限被拒。

    每分钟配额用满时短暂等待；每日配额用满时直接报错（由调用方返回空表）。
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

    # 每日配额用满：等待可能长达数小时，直接拒绝
    if len(timestamps) >= config.DISGENET_RATE_LIMIT_PER_DAY:
        raise RuntimeError(
            "DisGeNET daily rate limit (%d) reached" % config.DISGENET_RATE_LIMIT_PER_DAY
        )

    # 每分钟配额用满：等到最早一次请求移出窗口后再继续
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

    # 先写临时文件再替换，避免同时写入弄坏记录文件
    os.makedirs(os.path.dirname(state_file), exist_ok=True)
    tmp_file = state_file + ".tmp"
    with open(tmp_file, "w", encoding="utf-8") as f:
        json.dump(timestamps, f)
    os.replace(tmp_file, state_file)


def fetch_disease_targets_disgenet(disease_name: str) -> pd.DataFrame:
    """从 DisGeNET 获取疾病的关联基因（gene-disease association, GDA）。

    当前 API key 为试用（TRIAL）账号：每页最多 30 条、不支持分页，
    因此只取第一页结果；遇到限流（429）时做有限次等待重试。
    失败或无数据时返回空表；结果存为本地 CSV，已有缓存则直接复用。
    """
    import time
    import config

    csv_path = os.path.join(config.DISEASE_CACHE_DIR, f"DisGeNET_{safe_name(disease_name)}.csv")

    # 已有缓存则直接读取本地文件
    cached = open_csv_cache(csv_path, label=f"DisGeNET targets for '{disease_name}'")
    if cached is not None:
        return cached

    if not config.DisGeNET_API_KEY:
        logging.error("DisGeNET_API_KEY is empty; check API.env")
        return pd.DataFrame()

    headers = {"Authorization": config.DisGeNET_API_KEY, "accept": "application/json"}
    # gda summary 端点按疾病名（或 UMLS CUI）查询
    params = {"disease": disease_name.strip()}

    try:
        # 发请求前先过限流检查；遇到 429 限流时按服务器指示等待后重试
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
            # 若需等待过久则直接放弃，避免界面一直卡住
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
        save_csv_cache(df, csv_path)
        logging.info("DisGeNET '%s': %d targets fetched", disease_name, len(df))
        return df

    except Exception:
        logging.exception("fetch_disease_targets_disgenet failed for '%s'", disease_name)
        return pd.DataFrame()
