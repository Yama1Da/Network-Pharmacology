"""PK-DB 镜像层:把 pk-db.com 的开放药代动力学数据批量搬到本地(二期方案 §3.1 第 1 层)。

分层纪律(架构红线):本模块只搬运——TSV 原文逐字保存,不换算单位、不过滤剂量、
不做任何清洗;清洗/换算在归一化层(core/reachability.py,PK_NORM_VERSION 控制)。
分层后改清洗规则只需重跑归一化,不重新联网。

数据路径(2026-11-17 探针定案,"数字化源表直下"路线):
  1. GET {PKDB_API_BASE}/pkdata/studies/?page=N   研究清单(803 个:sid/substances/pmid)
  2. GET {PKDB_API_BASE}/pkdata/studies/<sid>/    研究摘要
  3. GET {PKDB_API_BASE}/_studies/<sid>/          全量嵌套(files = 数字化源表 URL 清单)
  4. GET files 清单里的绝对 URL                    下载 .tsv(跳过 .png/.pdf/.xlsx 并记账;
                                                   URL 含 media 域,无需拼接 base——PKDB_MEDIA_BASE 已删)
探针记录:outputs 结构化端点是黑盒(列表空参 count=0、详情 404、多种过滤参数
  名均失败),源表直下为实测可行路径——不要再去尝试 outputs 列表端点。
  源表即药代参数本体:`Tab*.tsv` 含 measurement_type=cmax/auc_end/clearance 行,
  且带 intervention 列回链给药方案表(剂型/途径/剂量),供归一化层做"治疗剂量
  过滤"与暴露-活性配对。

可靠性三件套(2026-11-17 定案,当日实跑 3 个失败全定性后加):
  ① 并发搬运——下载是纯 I/O 等待,线程池重叠回包等待提速 2.5–3 倍;http_client
    的 0.2s 礼貌间隔有全局锁,请求发起节奏与串行完全一致(加 worker 不加压);
  ② Range 断点续传兜底——上游对个别较大源表做"确定性"chunked 截断(实测
    .Lenuzza2016_Fig1.tsv 恒断于 65269/310696 字节,每次重试死在同一位置),
    全量 GET 永远拿不全;实测服务器支持 Range(206),从断点续传拼接即拿全;
  ③ 研究级 404 归"缺失"——索引有名字、文件集从未发布的悬空引用(如
    PKDB00965),重试无意义:不算镜像失败、不进补试,单独立账。

每研究落盘 {CACHE_DIR}/{PKDB_CACHE_SUBDIR}/<sid>.json:
  研究摘要(sid/name/date/licence/access/substances/pmid/title) + files 清单 +
  tsv:[{file,url,rows}] (rows 为字符串二维表,逐字保存) + skipped + fetched_at。
licence/access 逐字保留(许可审计);数值一律存原始字符串,换算在归一化层。
另存 _meta.json(构建统计,供管理页状态)与 substance_map.json(物质同义词/
  ChEBI 注释,归一化层对齐 ChEMBL ID 的桥)。
API 路由:/api/pk-mirror{,/build} 在 WY/api/pk_mirror.py(jobs 提交-轮询)。
"""

import csv
import io
import json
import logging
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

from WY import config
from WY.core.http_client import request_with_retry, throttle

logger = logging.getLogger(__name__)

_TSV_EXT = ".tsv"
_SKIP_EXT = (".png", ".pdf", ".jpg", ".xlsx")  # xlsx 需 openpyxl,暂不镜像,记入 skipped


class MirrorIndexError(RuntimeError):
    """研究清单拉取失败(无法继续镜像,整批终止)。"""


class StudyMissingError(RuntimeError):
    """研究文件集上游缺失(summary 或 _studies 端点 404,悬空引用):重试无意义,
    不算镜像失败、不进补试,单独立账进 _meta.json 的 missing 清单。"""


def _raw_dir() -> str:
    """镜像落盘目录:调用时读 config(沙箱/测试改 CACHE_DIR 即生效)。"""
    return os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR)


def _meta_path() -> str:
    return os.path.join(_raw_dir(), "_meta.json")


def _get_json(url: str, label: str) -> tuple[dict | None, int | None]:
    """拉 JSON。返回 (payload, status):网络重试耗尽 → (None, None);
    404(上游悬空引用,短路不重试) → (None, 404);成功 → (payload, 200)。"""
    resp = request_with_retry(
        url, method="get", max_retries=config.LIGHT_API_MAX_RETRIES,
        timeout=config.API_TIMEOUT, label=label, throttle_key="pkdb",
        use_cache=True, retry_delay=config.PKDB_RETRY_DELAY,
        no_retry_status=(404,),
        cache_ttl_hours=24,   # PK-DB 缓存 1 天(2026-11-17 定稿):补试近零成本,>1 天的"更新"能拿到策展修订
    )
    if resp is None:
        return None, None
    if resp.status_code == 404:
        return None, 404
    return resp.json(), int(resp.status_code)


def _stream_with_resume(url: str, label: str, max_parts: int = 30) -> str | None:
    """流式下载 + Range 断点续传兜底(2026-11-17 定案):上游对个别较大源表做
    确定性 chunked 截断(实测 .Lenuzza2016_Fig1.tsv 恒断于 65269 字节),全量
    GET 永远不完整、逐次重试只会死在同一位置;实测服务器支持 Range(206),
    从断点续传拼接即可拿全。返回 UTF-8 全文;超轮限/非 200|206 → None。
    只在上游缓存路径重试耗尽后才走这里(正常文件不受影响),故不走 sqlite
    HTTP 缓存——残缺响应不能进缓存,整文件本地持久化在 <sid>.json 里。

    两个实现细节(2026-11-17 实跑教训,首版续传空手而归的根因):
    - 强制 Accept-Encoding: identity——requests 默认发 gzip,服务端若压缩,
      Range 偏移是对压缩流的,而 iter_content 给的是解码后字节,拼接必错位;
    - 自适应缩块——urllib3 的 iter_content 不足一块的断前尾部随 IncompleteRead
      一起丢弃(截断点 65269 < 64KB 首块 ⇒ 首版每轮零进展、20 轮空转),
      零进展轮把块 16384→…→1 收缩,块小于"到截断点的余量"后断前字节即可全部入袋。"""
    import requests

    parts: list[bytes] = []
    got = 0
    chunk = 16384
    for round_no in range(1, max_parts + 1):
        throttle("pkdb")   # 续传请求同样守 0.2s 礼貌间隔
        headers = {"Accept-Encoding": "identity"}
        if got:
            headers["Range"] = f"bytes={got}-"
        try:
            with requests.get(url, headers=headers, timeout=config.API_TIMEOUT,
                              stream=True) as resp:
                if resp.status_code == 200 and got:
                    # 服务器无视 Range 回了整段:丢弃旧断片从头收
                    parts, got = [], 0
                elif resp.status_code not in (200, 206):
                    logger.warning("%s: 续传响应异常 status=%s", label, resp.status_code)
                    return None
                got_before = got
                for piece in resp.iter_content(chunk_size=chunk):
                    parts.append(piece)
                    got += len(piece)
                # 读到尾即拿全(200 一次拿全 / 206 收到段尾)
                if round_no > 1:
                    logger.info("%s: Range 续传完成,共 %d 字节(%d 轮)", label, got, round_no)
                return b"".join(parts).decode("utf-8", errors="replace")
        except Exception as exc:
            if got == got_before and chunk > 1:
                # 本轮零进展:截断点落在本轮首块内部,不足一块的断前尾部被整块丢弃
                chunk = max(1, chunk // 8)
                logger.warning("%s: 零进展轮,块缩至 %d 字节(%s 于 %d 字节)",
                               label, chunk, type(exc).__name__, got)
            else:
                logger.warning("%s: 传输中断于 %d 字节(%s),第 %d 轮续传",
                               label, got, type(exc).__name__, round_no + 1)
    logger.error("%s: Range 续传 %d 轮仍未拿全", label, max_parts)
    return None


def _get_text(url: str, label: str) -> tuple[str | None, int | None]:
    """下载文本。返回 (text, status):404(上游引用悬空,重试无意义) → (None, 404);
    成功 → (text, 200);缓存路径重试耗尽 → Range 续传兜底,仍败 → (None, None)。"""
    resp = request_with_retry(
        url, method="get", max_retries=config.API_MAX_RETRIES,
        timeout=config.API_TIMEOUT, label=label, throttle_key="pkdb",
        use_cache=True, no_retry_status=(404,),
        retry_delay=config.PKDB_RETRY_DELAY, cache_ttl_hours=24,
    )
    if resp is not None:
        if resp.status_code == 404:
            return None, 404
        # 显式按 UTF-8 解码:requests 对无 charset 的 text/* 默认按 latin-1,
        # 作者名/化学品名的非 ASCII 字符(如 König)会变乱码——镜像要逐字保存原文。
        return resp.content.decode("utf-8", errors="replace"), int(resp.status_code)
    # 重试耗尽多为确定性截断(IncompleteRead),重试永不过关:转 Range 续传
    text = _stream_with_resume(url, label)
    if text is None:
        return None, None
    return text, 200


def _fetch_index(progress=None) -> list[dict]:
    """全量研究清单(pkdata/studies 分页;翻页跟 next_page_url,不自己拼页码)。"""
    studies: list[dict] = []
    url: str | None = f"{config.PKDB_API_BASE}/pkdata/studies/"
    pages = 0
    while url:
        payload, _ = _get_json(url, "PK-DB study index")
        if payload is None:
            raise MirrorIndexError("PK-DB 研究清单拉取失败(重试后仍不可达)")
        page = payload.get("data", {})
        studies.extend(page.get("data", []))
        pages += 1
        if progress:
            progress(text=f"研究清单:已取 {len(studies)} 个(第 {pages} 页)")
        url = payload.get("next_page_url")
    return studies


def _collect_media(full: dict) -> tuple[list[str], list[dict]]:
    """从 _studies 嵌套的 files 清单分流:tsv 待下载,其余(png/pdf/xlsx)记入 skipped。"""
    download: list[str] = []
    skipped: list[dict] = []
    seen: set[str] = set()
    for url in full.get("files", []) or []:
        if not url or url in seen:
            continue
        seen.add(url)
        name = url.rsplit("/", 1)[-1]
        lower = url.lower()
        if lower.endswith(_TSV_EXT):
            download.append(url)
        elif lower.endswith(_SKIP_EXT):
            skipped.append({"file": name, "reason": "非文本源表(png/pdf/xlsx),不镜像"})
        else:
            skipped.append({"file": name, "reason": "未知类型,不镜像"})
    return download, skipped


def _parse_tsv(text: str) -> list[list[str]]:
    """TSV → 字符串二维表(逐字保存;不转数值、不填 NA——清洗在归一化层)。"""
    return [row for row in csv.reader(io.StringIO(text), delimiter="\t")]


def _fetch_substance_map(progress=None) -> list[dict]:
    """物质映射(sid/同义词/ChEBI 注释):归一化层对齐 ChEMBL ID 的桥。"""
    nodes: list[dict] = []
    url: str | None = f"{config.PKDB_API_BASE}/info_nodes/?ntype=substance"
    while url:
        payload, _ = _get_json(url, "PK-DB substance map")
        if payload is None:
            raise RuntimeError("PK-DB 物质映射拉取失败")
        page = payload.get("data", {})
        nodes.extend(page.get("data", []))
        if progress:
            progress(text=f"物质映射:已取 {len(nodes)} 个")
        url = payload.get("next_page_url")
    return nodes


def _fetch_study_payload(sid: str) -> dict:
    """单研究全量搬运:摘要 + 嵌套文件清单 + 逐个下载 .tsv 解析为二维表。
    摘要或文件清单 404(上游悬空引用) → StudyMissingError,由调用方归缺失立账。"""
    summary, s_status = _get_json(f"{config.PKDB_API_BASE}/pkdata/studies/{sid}/",
                                  f"PK-DB study summary {sid}")
    full, f_status = _get_json(f"{config.PKDB_API_BASE}/_studies/{sid}/",
                               f"PK-DB study files {sid}")
    if summary is None or full is None:
        if 404 in (s_status, f_status):
            raise StudyMissingError(f"{sid}: 研究摘要或文件集上游缺失(404)")
        raise RuntimeError(f"{sid}: 研究摘要/文件清单拉取失败(网络)")
    download, skipped = _collect_media(full)
    tsv_entries: list[dict] = []
    for url in download:
        name = url.rsplit("/", 1)[-1]
        text, status = _get_text(url, f"PK-DB tsv {name}")
        if text is None:
            if status == 404:
                # 上游引用了已删除/未发布的文件:重试无意义,记缺失不算研究失败
                skipped.append({"file": name, "reason": "404 上游缺失"})
                continue
            raise RuntimeError(f"{sid}: 源表下载失败(网络) {name}")
        tsv_entries.append({"file": name, "url": url, "rows": _parse_tsv(text)})
    return {
        "sid": sid,
        "name": summary.get("name", sid),
        "date": summary.get("date", ""),
        "licence": summary.get("licence", ""),
        "access": summary.get("access", ""),
        "substances": summary.get("substances", []),
        "reference_pmid": summary.get("reference_pmid", ""),
        "reference_title": summary.get("reference_title", ""),
        "files": full.get("files", []),
        "tsv": tsv_entries,
        "skipped": skipped,
        "fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def _atomic_write_json(path: str, obj) -> None:
    """原子替换写盘(tmp + os.replace,照 core/caching.py 模式),防并发读写撞坏。"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False)
    os.replace(tmp, path)


def build_mirror(force: bool = False, progress=None,
                 max_studies: int | None = None) -> dict:
    """全量镜像:研究清单 → 并发搬运(线程池,礼貌间隔不变) → 落盘。

    并发(2026-11-17 定案):下载是纯 I/O 等待,worker 只重叠回包等待;0.2s
    礼貌间隔有全局锁,请求发起节奏与串行完全一致,对 PK-DB 不加压。主循环与
    失败补试共用 _run_batch,记账同锁。
    断点续传:已存在的 <sid>.json 在 force=False 时跳过(顺序先筛掉,不进池);
    单研究失败记入 failures 不断批(失败清单随 _meta.json 持久化,重跑可补齐);
    研究级 404 归 missing 立账,不算失败、不进补试。
    返回统计 dict(同时持久化到 _meta.json,供管理页状态查询)。
    """
    started = time.time()
    os.makedirs(_raw_dir(), exist_ok=True)
    studies = _fetch_index()
    if max_studies:
        studies = studies[:max_studies]
    sids = [str(item.get("sid", "") or "") for item in studies]
    sids = [s for s in sids if s]
    total = len(sids)

    lock = threading.Lock()
    # 计数续传累计(2026-11-17 修复):tsv/skipped 原本只计本次新建,断点续传
    # 一跑(全部命中已存文件)就把上一轮累计覆盖清零,状态页出现"645 研究 ·
    # 0 源表"。force 全量重建从零起算;续传时接上上一轮 meta 的累计值。
    # 上一轮 meta 缺失(首次构建被中断)时从 0 起算——完成一轮后恢复精确。
    prev_tsv = prev_skipped = 0
    if not force and os.path.exists(_meta_path()):
        try:
            with open(_meta_path(), encoding="utf-8") as f:
                _prev = json.load(f)
            prev_tsv = int(_prev.get("tsv_files") or 0)
            prev_skipped = int(_prev.get("skipped_files") or 0)
        except (OSError, ValueError):
            logger.warning("PK-DB 镜像上一轮 _meta.json 不可读,源表计数从 0 重新累计")
    stats = {"built": 0, "skipped_existing": 0, "tsv": prev_tsv,
             "skipped_files": prev_skipped,
             "done": 0, "failures": [], "missing": []}

    def _attempt(sid: str) -> None:
        """单研究搬运(下载+落盘);缺失抛 StudyMissingError、失败抛异常,
        由 _run_batch 统一记账。"""
        payload = _fetch_study_payload(sid)
        _atomic_write_json(os.path.join(_raw_dir(), f"{sid}.json"), payload)
        with lock:
            stats["built"] += 1
            stats["tsv"] += len(payload["tsv"])
            stats["skipped_files"] += len(payload["skipped"])

    def _run_batch(batch: list[str], *, frac: float | None,
                   text_prefix: str) -> None:
        """一批研究并发搬运:完成顺序不定,进度按完成数推进(不能按下标)。
        frac=None 时进度条随完成数爬升(主循环);给定 frac 时定格(补试)。"""
        if not batch:
            return
        workers = max(1, int(getattr(config, "PKDB_DOWNLOAD_WORKERS", 4)))
        with ThreadPoolExecutor(max_workers=workers,
                                thread_name_prefix="pkdb-mirror") as pool:
            futures = {pool.submit(_attempt, sid): sid for sid in batch}
            done_local = 0   # 仅 as_completed 循环线程(调用方主线程)读写,无需锁
            for fut in as_completed(futures):
                sid = futures[fut]
                try:
                    fut.result()
                except StudyMissingError as exc:
                    with lock:
                        stats["missing"].append({"sid": sid, "reason": str(exc)})
                    logger.info("PK-DB 镜像研究 %s 上游缺失(404),归 missing 不算失败", sid)
                except Exception as exc:
                    with lock:
                        stats["failures"].append({"sid": sid,
                                                  "error": f"{type(exc).__name__}: {exc}"})
                    logger.exception("PK-DB 镜像研究 %s 失败", sid)
                finally:
                    with lock:
                        stats["done"] += 1
                        d = stats["done"]
                    done_local += 1
                if progress:
                    # 主批按全局完成数/总数(含续传跳过);补试批按本批进度——
                    # 此时全局 done 已 ≥ total,再当分母会显示"804/803"(2026-11-17 审查)
                    shown = (f"{d}/{total}" if frac is None
                             else f"{done_local}/{len(batch)}")
                    p_frac = (min(d / max(total, 1) * 0.95, 0.95)
                              if frac is None else frac)
                    progress(fraction=p_frac, text=f"{text_prefix} {shown}: {sid}")

    todo: list[str] = []
    for sid in sids:
        if not force and os.path.exists(os.path.join(_raw_dir(), f"{sid}.json")):
            with lock:
                stats["skipped_existing"] += 1
                stats["done"] += 1
                d = stats["done"]
            if progress:
                progress(fraction=min(d / max(total, 1) * 0.95, 0.95),
                         text=f"续传跳过 {d}/{total}: {sid}")
            continue
        todo.append(sid)
    _run_batch(todo, frac=None, text_prefix="镜像研究")

    # 失败补试:网络抖动/上游瞬时故障多数一轮内恢复;连续两轮补试仍败才最终
    # 记失败。缺失(missing)不进补试——404 是永久状态,重试无意义。失败研究
    # 无落盘文件,下次构建断点续传仍会自动补(不会被跳过)。
    retry_passes_used = 0
    for pass_no in (1, 2):
        with lock:
            pending = [f["sid"] for f in stats["failures"]]
            stats["failures"] = []
        if not pending:
            break
        retry_passes_used += 1
        if progress:
            progress(fraction=0.96,
                     text=f"失败补试第 {pass_no} 轮({len(pending)} 个研究)")
        _run_batch(pending, frac=0.96, text_prefix=f"补试第 {pass_no} 轮")
    with lock:
        failures = list(stats["failures"])
        missing = list(stats["missing"])
    substance_nodes = 0
    substance_map_ok = True
    try:
        if progress:
            progress(fraction=0.97, text="物质映射")
        substance_map = _fetch_substance_map(progress=progress)
        substance_nodes = len(substance_map)
        _atomic_write_json(os.path.join(_raw_dir(), "substance_map.json"),
                           {"fetched_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                            "nodes": substance_map})
    except Exception as exc:
        substance_map_ok = False
        failures.append({"sid": "_substance_map", "error": f"{type(exc).__name__}: {exc}"})
        logger.exception("PK-DB 物质映射拉取失败")
    # 在物质映射记账之后取值(2026-11-17 审查修正:此前提前取值,映射失败会
    # 漏计——failures 清单含 _substance_map 但 studies_failed 不增,前端据此
    # 不出失败警告与重试入口,静默半失败)
    failed = len(failures)
    meta = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "studies_indexed": total,
        "studies_built": stats["built"],
        "studies_skipped_existing": stats["skipped_existing"],
        "studies_failed": failed,
        "studies_missing": len(missing),
        "retry_passes": retry_passes_used,
        "failures": failures,
        "missing": missing,
        "tsv_files": stats["tsv"],
        "skipped_files": stats["skipped_files"],
        "substance_nodes": substance_nodes,
        "substance_map_ok": substance_map_ok,
        "duration_sec": round(time.time() - started, 1),
    }
    _atomic_write_json(_meta_path(), meta)
    if progress:
        done_text = f"完成:新建 {stats['built']},续传跳过 {stats['skipped_existing']}"
        if failed:
            # 收尾警告(2026-11-17):补试双轮仍失败的计数要明确提示用户(404 缺失不计)
            done_text += f",仍有 {failed} 个研究失败(404 缺失不计,可点重试)"
        if missing:
            done_text += f";{len(missing)} 个研究上游文件集缺失(404)"
        progress(fraction=1.0, text=done_text)
    return meta


def mirror_status() -> dict:
    """数据管理页"PK-DB 镜像"卡片的状态查询:built/研究数/建库时间/是否过 TTL。
    只扫描目录与读 _meta.json,不联网。"""
    raw_dir = _raw_dir()
    ttl = int(config.PK_MIRROR_TTL_DAYS)
    if not os.path.isdir(raw_dir):
        return {"built": False, "studies": 0, "age": None, "stale": True,
                "ttl_days": ttl, "built_at": ""}
    sids = [f for f in os.listdir(raw_dir) if f.startswith("PKDB") and f.endswith(".json")]
    meta: dict = {}
    if os.path.exists(_meta_path()):
        try:
            with open(_meta_path(), encoding="utf-8") as f:
                meta = json.load(f)
        except Exception:
            logger.exception("PK-DB 镜像 _meta.json 读取失败")
    built_at = str(meta.get("built_at", "") or "")
    if not built_at and sids:
        newest = max(os.path.getmtime(os.path.join(raw_dir, f)) for f in sids)
        built_at = datetime.fromtimestamp(newest, timezone.utc).isoformat(timespec="seconds")
    age = None
    stale = True
    if built_at:
        try:
            started_at = datetime.fromisoformat(built_at)
            if started_at.tzinfo is None:
                started_at = started_at.replace(tzinfo=timezone.utc)
            age = round((datetime.now(timezone.utc) - started_at).total_seconds() / 86400, 1)
            stale = age > ttl
        except ValueError:
            logger.exception("PK-DB 镜像建库时间解析失败: %s", built_at)
    return {"built": len(sids) > 0, "studies": len(sids),
            "tsv_files": meta.get("tsv_files"), "failures": meta.get("studies_failed", 0),
            "missing": meta.get("studies_missing", 0),
            "built_at": built_at[:10], "age": age, "stale": stale, "ttl_days": ttl}
