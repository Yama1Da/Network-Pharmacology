"""core.http_client 的测试:同站最小请求间隔(throttle)与统一重试入口
(request_with_retry)的成败/重试/豁免语义。

第 1 批阶段三扩充(重试语义原自检未覆盖);请求一律打桩,不联网。
"""
import threading
import time

import requests

from WY import config
import WY.core.http_client as hc
from WY.core.http_client import throttle


class _FakeResp:
    """最小假响应:只需要 status_code 与 raise_for_status。"""

    def __init__(self, status: int):
        self.status_code = status

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"http {self.status_code}")


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


def test_get_session_is_singleton():
    """共享会话:两次获取是同一个对象(懒加载双重检查锁),且为带落盘缓存的会话。"""
    s1 = hc.get_session()
    s2 = hc.get_session()
    assert s1 is s2
    assert hasattr(s1, "cache")  # requests_cache 可用:sqlite 落盘缓存会话


def test_retry_success_first_try(monkeypatch):
    """一次成功:直接返回响应,只发一次请求。"""
    calls = {"n": 0}

    def fake_post(url, **kw):
        calls["n"] += 1
        return _FakeResp(200)

    monkeypatch.setattr(requests, "post", fake_post)
    resp = hc.request_with_retry("http://x", max_retries=3, timeout=5,
                                 label="test", retry_delay=0)
    assert resp is not None and resp.status_code == 200
    assert calls["n"] == 1


def test_retry_recovers_after_failures(monkeypatch):
    """前两次失败第三次成功:自动重试并返回成功响应。"""
    calls = {"n": 0}

    def flaky_post(url, **kw):
        calls["n"] += 1
        if calls["n"] < 3:
            raise ConnectionError("simulated outage")
        return _FakeResp(200)

    monkeypatch.setattr(requests, "post", flaky_post)
    resp = hc.request_with_retry("http://x", max_retries=3, timeout=5,
                                 label="test", retry_delay=0)
    assert resp is not None and calls["n"] == 3


def test_retry_exhausted(monkeypatch):
    """全部失败:默认返回 None 不抛异常;raise_on_exhausted=True 抛最后一次异常。"""
    calls = {"n": 0}

    def boom(url, **kw):
        calls["n"] += 1
        raise ConnectionError("simulated outage")

    monkeypatch.setattr(requests, "post", boom)
    resp = hc.request_with_retry("http://x", max_retries=2, timeout=5,
                                 label="test", retry_delay=0)
    assert resp is None and calls["n"] == 2
    try:
        hc.request_with_retry("http://x", max_retries=2, timeout=5,
                              label="test", retry_delay=0, raise_on_exhausted=True)
        raise AssertionError("raise_on_exhausted must re-raise")
    except ConnectionError:
        pass


def test_no_retry_status_passthrough(monkeypatch):
    """no_retry_status 豁免(如 STRING 404=基因全未识别属正常业务):原样返回,不重试不抛错。"""
    calls = {"n": 0}

    def not_found(url, **kw):
        calls["n"] += 1
        return _FakeResp(404)

    monkeypatch.setattr(requests, "get", not_found)
    resp = hc.request_with_retry("http://x", method="get", max_retries=3, timeout=5,
                                 label="test", retry_delay=0, no_retry_status=(404,))
    assert resp is not None and resp.status_code == 404
    assert calls["n"] == 1, "豁免状态码不得触发重试"


def test_error_status_retries(monkeypatch):
    """服务器错误(500):raise_for_status 视为失败,重试至耗尽返回 None。"""
    calls = {"n": 0}

    def server_error(url, **kw):
        calls["n"] += 1
        return _FakeResp(500)

    monkeypatch.setattr(requests, "get", server_error)
    resp = hc.request_with_retry("http://x", method="get", max_retries=2, timeout=5,
                                 label="test", retry_delay=0)
    assert resp is None and calls["n"] == 2
