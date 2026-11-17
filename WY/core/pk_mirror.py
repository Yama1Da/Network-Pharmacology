"""PK-DB 镜像层:把 pk-db.com 的开放药代动力学数据批量搬到本地(二期方案 §3.1 第 1 层)。

分层纪律(架构红线):本模块只搬运——TSV 原文逐字保存,不换算单位、不过滤剂量、
不做任何清洗;清洗/换算在归一化层(core/reachability.py,PK_NORM_VERSION 控制)。
分层后改清洗规则只需重跑归一化,不重新联网。

数据路径(2026-11-17 探针定案,"数字化源表直下"路线):
  1. GET {PKDB_API_BASE}/pkdata/studies/?page=N   研究清单(803 个:sid/substances/pmid)
  2. GET {PKDB_API_BASE}/pkdata/studies/<sid>/    研究摘要
  3. GET {PKDB_API_BASE}/_studies/<sid>/          全量嵌套(files = 数字化源表 URL 清单)
  4. GET {PKDB_MEDIA_BASE}/<file>                 下载 .tsv(跳过 .png/.pdf/.xlsx 并记账)
探针记录:outputs 结构化端点是黑盒(列表空参 count=0、详情 404、多种过滤参数
  名均失败),源表直下为实测可行路径——不要再去尝试 outputs 列表端点。
  源表即药代参数本体:`Tab*.tsv` 含 measurement_type=cmax/auc_end/clearance 行,
  且带 intervention 列回链给药方案表(剂型/途径/剂量),供归一化层做"治疗剂量
  过滤"与暴露-活性配对。

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
import time
from datetime import datetime, timezone

from WY import config
from WY.core.http_client import request_with_retry

logger = logging.getLogger(__name__)

_TSV_EXT = ".tsv"
_SKIP_EXT = (".png", ".pdf", ".jpg", ".xlsx")  # xlsx 需 openpyxl,暂不镜像,记入 skipped


class MirrorIndexError(RuntimeError):
    """研究清单拉取失败(无法继续镜像,整批终止)。"""


def _raw_dir() -> str:
    """镜像落盘目录:调用时读 config(沙箱/测试改 CACHE_DIR 即生效)。"""
    return os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR)


def _meta_path() -> str:
    return os.path.join(_raw_dir(), "_meta.json")


def _get_json(url: str, label: str):
    resp = request_with_retry(
        url, method="get", max_retries=config.LIGHT_API_MAX_RETRIES,
        timeout=config.API_TIMEOUT, label=label, throttle_key="pkdb",
        use_cache=True,
    )
    if resp is None:
        return None
    return resp.json()


def _get_text(url: str, label: str) -> str | None:
    resp = request_with_retry(
        url, method="get", max_retries=config.LIGHT_API_MAX_RETRIES,
        timeout=config.API_TIMEOUT, label=label, throttle_key="pkdb",
        use_cache=True,
    )
    if resp is None:
        return None
    # 显式按 UTF-8 解码:requests 对无 charset 的 text/* 默认按 latin-1,
    # 作者名/化学品名的非 ASCII 字符(如 König)会变乱码——镜像要逐字保存原文。
    return resp.content.decode("utf-8", errors="replace")


def _fetch_index(progress=None) -> list[dict]:
    """全量研究清单(pkdata/studies 分页;翻页跟 next_page_url,不自己拼页码)。"""
    studies: list[dict] = []
    url: str | None = f"{config.PKDB_API_BASE}/pkdata/studies/"
    pages = 0
    while url:
        payload = _get_json(url, "PK-DB study index")
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
        payload = _get_json(url, "PK-DB substance map")
        if payload is None:
            raise RuntimeError("PK-DB 物质映射拉取失败")
        page = payload.get("data", {})
        nodes.extend(page.get("data", []))
        if progress:
            progress(text=f"物质映射:已取 {len(nodes)} 个")
        url = payload.get("next_page_url")
    return nodes


def _fetch_study_payload(sid: str) -> dict:
    """单研究全量搬运:摘要 + 嵌套文件清单 + 逐个下载 .tsv 解析为二维表。"""
    summary = _get_json(f"{config.PKDB_API_BASE}/pkdata/studies/{sid}/",
                        f"PK-DB study summary {sid}")
    full = _get_json(f"{config.PKDB_API_BASE}/_studies/{sid}/",
                     f"PK-DB study files {sid}")
    if summary is None or full is None:
        raise RuntimeError(f"{sid}: 研究摘要/文件清单拉取失败")
    download, skipped = _collect_media(full)
    tsv_entries: list[dict] = []
    for url in download:
        name = url.rsplit("/", 1)[-1]
        text = _get_text(url, f"PK-DB tsv {name}")
        if text is None:
            raise RuntimeError(f"{sid}: 源表下载失败 {name}")
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
    """全量镜像:研究清单 → 逐研究(摘要+嵌套+源表) → 落盘。

    断点续传:已存在的 <sid>.json 在 force=False 时跳过;单研究失败记入
    failures 不断批(失败清单随 _meta.json 持久化,重跑 force=True 可补齐)。
    返回统计 dict(同时持久化到 _meta.json,供管理页状态查询)。
    """
    started = time.time()
    os.makedirs(_raw_dir(), exist_ok=True)
    studies = _fetch_index()
    if max_studies:
        studies = studies[:max_studies]
    total = len(studies)
    built = skipped_existing = failed = 0
    tsv_files = skipped_files = 0
    failures: list[dict] = []
    for i, item in enumerate(studies, start=1):
        sid = str(item.get("sid", "") or "")
        if progress:
            progress(fraction=(i - 1) / max(total, 1) * 0.95,
                     text=f"镜像研究 {i}/{total}: {sid}")
        if not sid:
            continue
        path = os.path.join(_raw_dir(), f"{sid}.json")
        if not force and os.path.exists(path):
            skipped_existing += 1
            continue
        try:
            payload = _fetch_study_payload(sid)
            _atomic_write_json(path, payload)
            built += 1
            tsv_files += len(payload["tsv"])
            skipped_files += len(payload["skipped"])
        except Exception as exc:
            failed += 1
            failures.append({"sid": sid, "error": f"{type(exc).__name__}: {exc}"})
            logger.exception("PK-DB 镜像研究 %s 失败", sid)
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
    meta = {
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "studies_indexed": total,
        "studies_built": built,
        "studies_skipped_existing": skipped_existing,
        "studies_failed": failed,
        "failures": failures,
        "tsv_files": tsv_files,
        "skipped_files": skipped_files,
        "substance_nodes": substance_nodes,
        "substance_map_ok": substance_map_ok,
        "duration_sec": round(time.time() - started, 1),
    }
    _atomic_write_json(_meta_path(), meta)
    if progress:
        progress(fraction=1.0, text=f"完成:新建 {built},续传跳过 {skipped_existing},失败 {failed}")
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
            "built_at": built_at[:10], "age": age, "stale": stale, "ttl_days": ttl}
