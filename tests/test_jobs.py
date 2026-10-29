"""core.jobs 的测试:长任务注册表的提交-轮询模型(完成/出错/进度钳制/TTL 回收)。

阶段三第 3 批新写(0% 覆盖小模块)。任务表经 monkeypatch 隔离,任务线程为
daemon,轮询带超时防挂死。
"""
import time

import pytest

from WY.core import jobs


def _wait_status(job_id: str, timeout: float = 5.0) -> dict:
    """轮询直到任务离开 running 状态(或超时失败)。"""
    deadline = time.time() + timeout
    while time.time() < deadline:
        snap = jobs.get(job_id)
        if snap and snap["status"] != "running":
            return snap
        time.sleep(0.02)
    raise AssertionError(f"job {job_id} did not finish within {timeout}s")


def test_submit_done_with_result_and_progress():
    """正常任务:返回值进入 result,状态 done,进度最终钉在 1.0。"""
    def work(progress):
        progress(fraction=0.5, text="half")
        return 42
    jid = jobs.submit(work, "测试任务")
    snap = _wait_status(jid)
    assert snap["status"] == "done"
    assert snap["result"] == 42
    assert snap["progress"]["fraction"] == 1.0
    assert snap["description"] == "测试任务"


def test_progress_fraction_clamped():
    """进度回调:分数钳制在 [0,1](1.7 → 1.0)。用事件在任务中途定格观察。"""
    import threading

    started = threading.Event()
    proceed = threading.Event()

    def work(progress):
        progress(fraction=1.7, text="over")
        started.set()
        proceed.wait(timeout=5)
        progress(fraction=-0.5)

    jid = jobs.submit(work, "clamp")
    assert started.wait(timeout=5)
    mid = jobs.get(jid)["progress"]
    assert mid["fraction"] == 1.0 and mid["text"] == "over"  # 超界钳到 1.0
    proceed.set()
    snap = _wait_status(jid)
    assert snap["status"] == "done"


def test_submit_error_captured():
    """任务抛异常:状态 error,错误消息含异常类型,不让线程把进程带崩。"""
    def bad(progress):
        raise ValueError("boom")

    jid = jobs.submit(bad, "会炸的任务")
    snap = _wait_status(jid)
    assert snap["status"] == "error"
    assert "ValueError" in snap["error"] and "boom" in snap["error"]
    assert "trace" in snap


def test_get_unknown_job():
    assert jobs.get("never") is None


def test_gc_reaps_only_expired_finished_jobs(monkeypatch):
    """TTL 回收:只清理超过保留期的已完成/出错任务;运行中的与新鲜的保留。"""
    fresh: dict = {}
    monkeypatch.setattr(jobs, "_JOBS", fresh)
    now = time.time()
    fresh["old_done"] = {"status": "done", "created": now - 7200}
    fresh["fresh_done"] = {"status": "done", "created": now - 60}
    fresh["old_running"] = {"status": "running", "created": now - 7200}
    fresh["old_error"] = {"status": "error", "created": now - 7200}
    jobs._gc()
    assert set(fresh) == {"fresh_done", "old_running"}
