"""进程内长任务注册表:TS API 的分钟级操作(来源刷新、互作组构建)的提交-轮询模型。

单机单用户定位,不引入 Celery/Redis:任务在 daemon 工作线程执行,进度经回调写入
注册表,HTTP 层只暴露 submit/get 两个动作。任务完成后结果保留 _TTL 秒供轮询取走,
过期自动回收。FastAPI 的同步端点跑在线程池里,与任务线程共享 core 的线程安全约定
(各数据源自带锁保护),与 Streamlit 页面的并发语义一致。
"""

import logging
import threading
import time
import traceback
import uuid
from typing import Callable

_LOCK = threading.Lock()
_JOBS: dict[str, dict] = {}
_TTL = 3600  # 已完成任务保留 1 小时,过期回收


def submit(fn: Callable[[Callable], object], description: str = "") -> str:
    """提交长任务并立即返回 job_id。fn 接收一个进度回调 progress(fraction=None,
    text="", **extra),任务返回值即 job 的 result。"""
    job_id = uuid.uuid4().hex[:12]
    with _LOCK:
        _JOBS[job_id] = {
            "status": "running", "description": description,
            "progress": {"fraction": None, "text": ""},
            "result": None, "error": None, "created": time.time(),
        }
    threading.Thread(target=_run, args=(job_id, fn), daemon=True,
                     name=f"npap-job-{job_id}").start()
    return job_id


def _run(job_id: str, fn: Callable[[Callable], object]) -> None:
    def progress(fraction: float | None = None, text: str = "", **extra) -> None:
        with _LOCK:
            job = _JOBS.get(job_id)
            if job is not None:
                p = dict(job["progress"])
                if fraction is not None:
                    p["fraction"] = max(0.0, min(1.0, float(fraction)))
                if text:
                    p["text"] = text
                p.update(extra)
                job["progress"] = p

    try:
        result = fn(progress)
        with _LOCK:
            if job_id in _JOBS:
                _JOBS[job_id].update(status="done", result=result,
                                     progress={"fraction": 1.0,
                                               "text": _JOBS[job_id]["progress"]["text"]})
    except Exception as exc:
        logging.exception("job %s (%s) failed", job_id,
                          _JOBS.get(job_id, {}).get("description", ""))
        with _LOCK:
            if job_id in _JOBS:
                _JOBS[job_id].update(status="error", error=f"{type(exc).__name__}: {exc}",
                                     trace=traceback.format_exc(limit=8))
    _gc()


def get(job_id: str) -> dict | None:
    """任务快照:status(running/done/error/unknown) + progress + result/error。"""
    with _LOCK:
        job = _JOBS.get(job_id)
        if job is None:
            return None
        return {**job, "job_id": job_id}


def _gc() -> None:
    now = time.time()
    with _LOCK:
        for jid in [j for j, job in _JOBS.items()
                    if job["status"] in ("done", "error") and now - job["created"] > _TTL]:
            del _JOBS[jid]
