"""分子对接多步会话(无状态 API 的服务端状态存储)。

TS 前端按步骤调用 /api/docking/* 端点;会话字典镜像原 Streamlit dock 会话
(键名保持一致便于对照迁移),进程内 dict + TTL 清理,单机单用户。
数据源/配体/受体/盒子/结果的相互失效规则与原页面 _reset_docking_downstream 一致。"""

import time
import uuid

_SESSIONS: dict[str, dict] = {}
_TTL = 4 * 3600  # 会话保留 4 小时,过期在新建时清扫


def create() -> str:
    _sweep()
    sid = uuid.uuid4().hex[:12]
    # _touched 必须在创建时落地:sweep 按 d.get("_touched", 0) 判过期,
    # 缺了它新会话会被下一次 create 立即清掉(2026-10 API 冒烟发现)
    _SESSIONS[sid] = {"_touched": time.time()}
    return sid


def get(sid: str) -> dict | None:
    dock = _SESSIONS.get(sid)
    if dock is not None:
        dock["_touched"] = time.time()   # 摸一次续命:按最后访问计时
    return dock


def drop(sid: str) -> None:
    _SESSIONS.pop(sid, None)


def reset_downstream(dock: dict) -> None:
    """数据源/基因/结构变更后,清空下游步骤的旧结果(原页面同名函数的会话版):
    残留的旧配体/受体/盒会静默覆盖新推断值,必须整段失效。"""
    for k in ("lig_pdbqt", "cid", "lig_auto_done", "gene", "rec_pdbqt", "pdb_id",
              "pdb_path", "box", "result", "result_dir", "starred", "receptor_source",
              "gene_structures", "candidate_pre", "batch_results",
              "gene_sel", "pdb_sel", "pose_sel", "conclusion",
              "drug_name", "disease_name", "lig_query"):
        dock.pop(k, None)


def _sweep() -> None:
    now = time.time()
    for sid in [s for s, d in _SESSIONS.items()
                if d.get("_touched", 0) < now - _TTL]:
        del _SESSIONS[sid]
