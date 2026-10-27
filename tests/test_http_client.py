"""core.http_client 的测试:同站最小请求间隔(throttle)的限速/不限速/并发排队。

迁移自 selfcheck_docking.py(18 项中的 1 项)。限流配置与状态表经 monkeypatch
隔离并自动还原,不影响其他测试。
"""
import threading
import time

from WY import config
import WY.core.http_client as hc
from WY.core.http_client import throttle


def test_throttle_spacing(monkeypatch):
    """对同一网站的请求有最小时间间隔;未配置的请求不限速;并发调用也正确排队。"""
    orig = config.HOST_MIN_REQUEST_INTERVAL
    monkeypatch.setattr(config, "HOST_MIN_REQUEST_INTERVAL",
                        {**(orig or {}), "selftest": 0.3})
    # 限流状态表整体换新:测试内的记时残留随 monkeypatch 自动还原,不泄漏给其他测试
    monkeypatch.setattr(hc, "_next_allowed_at", {})

    # 串行 3 次 → 至少跨越 2 个完整间隔
    t0 = time.monotonic()
    for _ in range(3):
        throttle("selftest")
    el = time.monotonic() - t0
    assert el >= 0.6, f"3 calls in {el:.3f}s, expected >= 0.6 (2 intervals)"

    # 未知键 = 不限速
    t0 = time.monotonic()
    for _ in range(3):
        throttle("unknown_host_x")
    assert time.monotonic() - t0 < 0.05

    # 4 个线程同时请求 → 总耗时至少 3 个间隔(不会一起放行)
    barrier = threading.Barrier(4)
    t0 = time.monotonic()
    threads = [threading.Thread(target=lambda: (barrier.wait(), throttle("selftest")))
               for _ in range(4)]
    [t.start() for t in threads]
    [t.join() for t in threads]
    el = time.monotonic() - t0
    assert el >= 0.9, f"4 concurrent calls in {el:.3f}s, expected >= 0.9 (3 intervals)"
