"""访问远程数据库的公共设置:统一的请求重试和访问频率控制。"""

import logging
import threading

import requests

OPENTARGETS_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

# 记录每个数据源下一次允许访问的时间,
# 用于控制访问频率,避免请求过快被对方拒绝。
_throttle_lock = threading.Lock()
_next_allowed_at: dict[str, float] = {}


def throttle(host_key: str) -> None:
    """控制对同一数据源的请求间隔,避免访问过快被对方限流。

    没有设置间隔要求的数据源不限速,直接通过。
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
    """向 Open Targets 数据库发送查询;网络失败时自动重试几次,仍失败则返回 None。

    查询本身报错(返回结果带 errors)不算网络失败,由调用方自行处理。
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
            resp.raise_for_status()  # 返回异常状态码(如服务器报错)时视为请求失败
            return resp.json()
        except Exception:
            logging.exception(
                "Open Targets POST failed (attempt %d/%d)",
                attempt + 1, config.LIGHT_API_MAX_RETRIES,
            )
            if attempt < config.LIGHT_API_MAX_RETRIES - 1:
                time.sleep(config.API_RETRY_DELAY)
    return None
