"""DisGeNET 数据源：查询疾病的关联基因（GDA），并按官方配额做主动限流。"""

import logging
import os
import threading

import pandas as pd

from .caching import load_json_file, open_named_csv_cache, save_csv_cache, save_json_file
from .http_client import _get_session
from .text_utils import remember_name_mapping

# 限流记录文件的进程内互斥锁:并发获取靶点时防止同时读写把时间戳记乱
_rate_limit_lock = threading.Lock()


def _acquire_disgenet_rate_limit_slot() -> None:
    """DisGeNET 主动限流：发请求前检查官方配额（每分钟/每天/每月），避免超限被拒。
    每分钟配额用满时短暂等待；每天/每月配额用满时直接报错（由调用方返回空表）。本地限额值已按学术计划官方配额(50/分、25,000/天、500,000/月)留余量,见 config。
    """
    import time
    import config

    with _rate_limit_lock:  # 读-改-写整段互斥:配额统计和写入必须视为一个事务
        now = time.time()
        minute_window = 60.0
        day_window = 24 * 3600.0
        month_window = 30 * day_window

        state_file = config.DISGENET_RATE_LIMIT_STATE_FILE
        timestamps = load_json_file(state_file, [])
        try:
            timestamps = [float(t) for t in timestamps]
        except (TypeError, ValueError):
            logging.warning("Failed to load DisGeNET rate-limit state; resetting")

        # 只保留一个月内的请求记录(同时供月/日/分钟三级配额统计)
        timestamps = [t for t in timestamps if now - t <= month_window]

        # 每月/每天配额用满：等待可能长达数小时~数周，直接拒绝
        if len(timestamps) >= config.DISGENET_RATE_LIMIT_PER_MONTH:
            raise RuntimeError(
                "DisGeNET monthly rate limit (%d) reached"
                % config.DISGENET_RATE_LIMIT_PER_MONTH
            )
        if len([t for t in timestamps if now - t <= day_window]) >= config.DISGENET_RATE_LIMIT_PER_DAY:
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
            timestamps = [t for t in timestamps if now - t <= month_window]

        timestamps.append(now)

        # 原子写入(临时文件+替换),并发写不会弄坏记录文件
        save_json_file(state_file, timestamps)


def _resolve_disgenet_query(disease_name: str) -> str:
    """把疾病输入转成 DisGeNET 认的查询串。
    DisGeNET 的 disease 参数只认带词表前缀的编号(UMLS_/MONDO_ 等),自由文本疾病名恒返回 0 行;裸 CUI(如 C0002395)需补 UMLS_ 前缀。自由文本先借 Open Targets 的疾病搜索解析成编号(解析结果有 HTTP 缓存)。
    """
    import re

    from .opentargets import find_disease_id

    query = disease_name.strip().replace("\u2019", "'")
    if re.fullmatch(r"C\d{6,8}", query):
        return f"UMLS_{query}"  # 裸 UMLS CUI
    if re.match(r"^[A-Z][A-Za-z0-9]*_[A-Za-z0-9_]+$", query):
        return query  # 已带词表前缀
    resolved = find_disease_id(query)
    if resolved:
        logging.info("DisGeNET: '%s' resolved to '%s' via Open Targets search", query, resolved)
        return resolved
    logging.warning("DisGeNET: could not resolve '%s' to a vocabulary id", query)
    return query


def fetch_disease_targets_disgenet(disease_name: str, need_rows: int | None = None) -> pd.DataFrame:
    """从 DisGeNET 获取疾病的关联基因（gene-disease association, GDA）;失败或无数据时返回空表,结果存为本地 CSV,已有缓存则直接复用。
    学术(Academic)账号只能访问 CURATED 来源(CLINVAR、UNIPROT、ORPHANET 等),每页固定 100 条;按 config.DISGENET_PAGES
    抓取前若干页(page_number 分页,响应里剩余条数不足一页时提前停止);遇到限流（429）时做有限次等待重试。
    need_rows 为本次分析至少需要的条数(参与调用方语义,原始缓存行数不足时自动重新抓取)。
    """
    import time
    import config

    csv_path, cached = open_named_csv_cache(
        config.DISEASE_CACHE_DIR, "DisGeNET_", disease_name,
        label=f"DisGeNET targets for '{disease_name}'",
        remember=remember_name_mapping,  # 记录安全文件名→原始名,供管理页还原
    )
    # 隐藏哨兵文件记录上次抓取进度:complete=已确认抓全(疾病关联总数不足页数上限),pages=已抓页数
    meta_path = os.path.join(os.path.dirname(csv_path), "." + os.path.basename(csv_path) + ".meta")
    meta = load_json_file(meta_path, {}) or {}
    if cached is not None:
        if need_rows is None:
            return cached
        if meta.get("complete") or int(meta.get("pages") or 0) >= int(config.DISGENET_PAGES):
            return cached
        if len(cached) >= need_rows:
            return cached

    if not config.DisGeNET_API_KEY:
        logging.error("DisGeNET_API_KEY is empty; check API.env")
        return pd.DataFrame()

    headers = {"Authorization": config.DisGeNET_API_KEY, "accept": "application/json"}
    # 弯撇号(’)统一成直撇号,疾病名解析成 DisGeNET 认的带前缀编号
    base_params = {"disease": _resolve_disgenet_query(disease_name)}
    # 超时收紧为(连接 10s, 读 30s):官方接口偶尔长时间不返回数据,
    # 用短超时把最坏等待控制住,避免拖住整个分析页面
    request_timeout = (10, 30)

    def _request_page(page: int) -> list | None:
        """抓取单页并返回 payload 行列表;429 按官方提示退避重试;任何失败(网络/HTTP/状态码)返回 None。"""
        params = {**base_params, "page_number": page}
        resp = None
        try:
            for _ in range(config.ENRICHMENT_MAX_RETRIES + 1):
                _acquire_disgenet_rate_limit_slot()
                resp = _get_session().get(
                    config.DisGeNET_URL, params=params, headers=headers, timeout=request_timeout
                )
                if resp.status_code != 429:
                    break
                retry_after = int(
                    resp.headers.get("x-rate-limit-retry-after-seconds", config.API_RETRY_DELAY)
                )
                # 若需等待过久则直接放弃，避免界面一直卡住
                if retry_after > config.API_RETRY_DELAY * config.ENRICHMENT_MAX_RETRIES * 10:
                    logging.error(
                        "DisGeNET rate limit requires %ss wait; aborting to avoid blocking",
                        retry_after
                    )
                    return None
                logging.warning("DisGeNET rate-limited (429), retrying after %ss", retry_after)
                time.sleep(retry_after)
            resp.raise_for_status()
            payload = resp.json()
        except Exception:
            logging.exception("DisGeNET request failed (%s page %d)", disease_name, page)
            return None
        if payload.get("status") != "OK":
            logging.error("DisGeNET query failed: %s", payload.get("payload"))
            return None
        return payload.get("payload") or []

    try:
        all_rows: list = []
        complete = False   # 是否确认翻到了结果末尾(疾病关联数不足页数上限)
        pages_done = 0
        for page in range(max(int(config.DISGENET_PAGES), 1)):
            rows = _request_page(page)
            if rows is None:
                # 中途某页失败:保留已抓到的页提前收工,不清空重来(首页失败才按整体失败处理)
                if all_rows:
                    logging.warning("DisGeNET page %d failed for '%s'; keeping %d rows "
                                    "from earlier pages", page, disease_name, len(all_rows))
                break
            all_rows.extend(rows)
            pages_done = page + 1
            # 空页或不满一页说明已是最后一页,继续翻页只会浪费配额
            if len(rows) < config.DISGENET_PAGE_SIZE:
                complete = True
                break
        if not all_rows:
            logging.warning("DisGeNET returned no associations for '%s'", disease_name)
            return pd.DataFrame()

        records: list[dict] = []
        for r in all_rows:
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
            # 跨页可能重复:先按分数降序排稳再去重,保证每个基因保留其最高分记录
            .assign(_s=lambda d: pd.to_numeric(d["disgenet_score"], errors="coerce"))
            .sort_values("_s", ascending=False, kind="stable")
            .drop(columns="_s")
            .drop_duplicates(subset=["approvedSymbol"])
            .reset_index(drop=True)
        )
        save_csv_cache(df, csv_path)
        save_json_file(meta_path, {"complete": complete, "pages": pages_done})
        logging.info("DisGeNET '%s': %d targets fetched (%d pages, complete=%s)",
                     disease_name, len(df), pages_done, complete)
        return df

    except Exception:
        logging.exception("fetch_disease_targets_disgenet failed for '%s'", disease_name)
        return pd.DataFrame()
