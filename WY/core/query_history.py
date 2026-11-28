"""查询历史台账(2026-11-28 负责人提议):全部双源整合成功过的药/病查询名,
单 JSON 单一真相源,供分析页/协同页选择器下拉复用。

生命周期与 merged 缓存**解耦**:设置页改筛选参数会 purge 整个 merged 目录
(设计内,见 api/settings._PURGE_MERGED_KEYS),台账保留——历史 ≠ 缓存,
选择器下拉因此不缩水。transfer scope "query_history" 可导出/导入;
导入按名字合并(union),绝不覆盖本机历史。

写入点唯一:core.targets._merge_source_pair 双源整合成功落盘后 record()——
台账里的每个名字都是"双源整合命中过"的查询名,语义与 merged 一致
(负责人定稿:不用"跑过但缓存未命中"的名字混入)。
文件: cache/query_history.json(刻意不在 merged/ 子目录,防 purge 误删)。
首次读取自动种子:现有 merged 缓存名 + 分析运行目录 analysis_meta 跑过名
(一次性;旧版运行无 meta 不计,诚实边界)。
"""

import json
import os
import threading
from datetime import datetime

_LOCK = threading.Lock()


def _path() -> str:
    from WY import config

    return os.path.join(config.CACHE_DIR, "query_history.json")


def _empty() -> dict:
    """空台账:drugs/diseases 两侧,名字 → {first_at, last_at}。"""
    return {"drugs": {}, "diseases": {}}


def _scan_run_meta_names() -> tuple[set[str], set[str]]:
    """分析运行目录 analysis_meta.json 的跑过名(种子源②);手动模式空名不计。"""
    import glob

    from WY import config

    drugs: set[str] = set()
    diseases: set[str] = set()
    root = os.path.join(config.DOWNLOADS_DIR, "ppi-results")
    if not os.path.isdir(root):
        return drugs, diseases
    for meta_path in glob.glob(os.path.join(root, "*", "analysis_meta.json")):
        try:
            with open(meta_path, encoding="utf-8") as f:
                meta = json.load(f)
        except (OSError, ValueError):
            continue
        drug = str(meta.get("drug") or "").strip()
        disease = str(meta.get("disease") or "").strip()
        if drug:
            drugs.add(drug)
        if disease:
            diseases.add(disease)
    return drugs, diseases


def _ensure_seeded_unlocked() -> dict:
    """首次读取种子(调用方须已持 _LOCK):台账文件不存在时,用现有 merged
    缓存名 ∪ 运行 meta 跑过名一次性建账;建不出任何名字则不落文件(下次再试)。"""
    from .targets import list_known_targets

    if os.path.isfile(_path()):
        return _load_unlocked()
    merged_drugs = set(list_known_targets("drug"))
    merged_diseases = set(list_known_targets("disease"))
    meta_drugs, meta_diseases = _scan_run_meta_names()
    now = datetime.now().isoformat(timespec="seconds")
    data = _empty()
    for name in merged_drugs | meta_drugs:
        data["drugs"][name] = {"first_at": now, "last_at": now}
    for name in merged_diseases | meta_diseases:
        data["diseases"][name] = {"first_at": now, "last_at": now}
    if data["drugs"] or data["diseases"]:
        _save_unlocked(data)
    return data


def _load_unlocked() -> dict:
    """读台账;缺失/损坏按空台账处理(以本次写入为准,不抛异常打断分析)。"""
    try:
        with open(_path(), encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, ValueError):
        return _empty()
    out = _empty()
    for side in ("drugs", "diseases"):
        v = data.get(side)
        if isinstance(v, dict):
            out[side] = {str(k): dict(rec) for k, rec in v.items()
                         if isinstance(rec, dict)}
    return out


def _save_unlocked(data: dict) -> None:
    data["updated_at"] = datetime.now().isoformat(timespec="seconds")
    from .caching import save_json_file

    # indent=2 + sort_keys:台账是人要看的登记簿(区别于 ui_settings 等机器状态文件),
    # 缩进分行 + 键序稳定,人工核对/逐字节 diff 都友好
    save_json_file(_path(), data, indent=2, sort_keys=True)


def load() -> dict:
    """读台账快照(种子语义同 names):drugs/diseases 名字 → {first_at, last_at}。"""
    with _LOCK:
        return _ensure_seeded_unlocked()


def record(kind: str, name: str) -> None:
    """双源整合成功落盘后记一笔(kind: "drug"|"disease");同名字刷新 last_at。
    并发安全(targets 取数线程并发写同一文件);损坏文件按空台账起底。"""
    name = str(name or "").strip()
    if not name or kind not in ("drug", "disease"):
        return
    now = datetime.now().isoformat(timespec="seconds")
    with _LOCK:
        data = _ensure_seeded_unlocked()
        side = data[kind + "s"]
        rec = side.get(name)
        if isinstance(rec, dict):
            rec["last_at"] = now
        else:
            side[name] = {"first_at": now, "last_at": now}
        _save_unlocked(data)


def names(kind: str) -> list[str]:
    """台账查询名列表(kind: "drug"|"disease"),字典序;首次调用自动种子。"""
    if kind not in ("drug", "disease"):
        return []
    with _LOCK:
        data = _ensure_seeded_unlocked()
        return sorted(data.get(kind + "s", {}))


def merge_blob(blob: dict) -> int:
    """导入合并(transfer scope query_history):包内台账按名字并入本地,
    已有则刷新 last_at,缺 first_at 补 now;返回处理的名字总数。
    永不覆盖/删除本地已有条目——台账只增不减。"""
    now = datetime.now().isoformat(timespec="seconds")
    n = 0
    with _LOCK:
        data = _ensure_seeded_unlocked()
        for side, incoming in (("drugs", blob.get("drugs")),
                               ("diseases", blob.get("diseases"))):
            if not isinstance(incoming, dict):
                continue
            for name, rec in incoming.items():
                name = str(name).strip()
                if not name:
                    continue
                existing = data[side].get(name)
                if isinstance(existing, dict):
                    existing["last_at"] = now
                else:
                    first = str(rec.get("first_at") or now) \
                        if isinstance(rec, dict) else now
                    data[side][name] = {"first_at": first, "last_at": now}
                n += 1
        _save_unlocked(data)
    return n


def export_blob() -> dict | None:
    """transfer 导出用:台账内容(含 updated_at);两侧全空返回 None(无物可导)。"""
    with _LOCK:
        data = _ensure_seeded_unlocked()
        if data.get("drugs") or data.get("diseases"):
            return data
        return None
