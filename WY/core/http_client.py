"""远端 API 的 HTTP 公共层：端点常量 + 带重试的请求封装 + 主机级主动限速。"""

import logging
import threading

import requests

OPENTARGETS_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# 主机级限速状态（进程内）：键 = config.HOST_MIN_REQUEST_INTERVAL 的主机名，
# 值 = 该主机下一次允许发起请求的 monotonic 时间。
_throttle_lock = threading.Lock()
_next_allowed_at: dict[str, float] = {}


def throttle(host_key: str) -> None:
    """主机级最小请求间隔（主动限速）。

    并发场景下在锁内预约时间槽、锁外睡眠：任意主机的请求始终以不小于
    config.HOST_MIN_REQUEST_INTERVAL[host_key] 的间隔串行发出，线程不会
    同时睡醒后连发。未知主机键 = 不限速，直接返回。
    """
    import time

    import config

    interval = float(config.HOST_MIN_REQUEST_INTERVAL.get(host_key, 0.0))
    if interval <= 0:
        return
    with _throttle_lock:
        now = time.monotonic()
        prev = _next_allowed_at.get(host_key, 0.0)
        wait = max(0.0, prev - now)
        _next_allowed_at[host_key] = max(now, prev) + interval
    if wait > 0:
        time.sleep(wait)


def _opentargets_post(query: str, variables: dict, timeout: int = 30) -> dict | None:
    """向 Open Targets GraphQL 端点发送 POST，带失败重试；最终失败返回 None。

    网络抖动 / 临时 5xx 时按退避重试（config.LIGHT_API_MAX_RETRIES 次），
    避免一次失败就让整个查询返回空表。GraphQL 层错误（payload["errors"]）
    不属于网络异常，由调用方自行判断。
    """
    import time
    import config

    for attempt in range(config.LIGHT_API_MAX_RETRIES):
        try:
            resp = requests.post(
                OPENTARGETS_API_URL,
                json={"query": query, "variables": variables},
                timeout=timeout,
            )
            resp.raise_for_status()  # 非 2xx 状态码会抛出异常
            return resp.json()
        except Exception:
            logging.exception(
                "Open Targets POST failed (attempt %d/%d)",
                attempt + 1, config.LIGHT_API_MAX_RETRIES,
            )
            if attempt < config.LIGHT_API_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)
    return None
