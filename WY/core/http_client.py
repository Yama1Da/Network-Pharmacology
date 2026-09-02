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


def request_with_retry(
    url: str,
    *,
    method: str = "post",
    max_retries: int,
    timeout: int,
    label: str,
    json: dict | None = None,
    data: dict | None = None,
    headers: dict | None = None,
    throttle_key: str | None = None,
    retry_delay: float | None = None,
    no_retry_status: tuple = (),
    raise_on_exhausted: bool = False,
) -> requests.Response | None:
    """带重试的 HTTP 请求,项目内所有数据源请求的统一入口。

    每次失败都记录日志并等待后重试;全部失败时默认返回 None,
    raise_on_exhausted=True 则抛出最后一次异常(调用方需要错误详情时使用)。
    no_retry_status 中的状态码不视为失败,直接原样返回(如 STRING 的 404 表示
    基因全部未识别,属正常业务情况)。throttle_key 传入时每次尝试前先过限速。
    """
    import time

    import config

    if retry_delay is None:
        retry_delay = config.API_RETRY_DELAY
    last_exc: Exception | None = None
    for attempt in range(max_retries):
        try:
            if throttle_key:
                throttle(throttle_key)
            # 按方法分别走 requests.post/get(而不是 requests.request),
            # 保持对 requests.post/get 打补丁的测试与工具兼容。
            if method == "get":
                resp = requests.get(url, timeout=timeout, headers=headers)
            else:
                resp = requests.post(url, json=json, data=data, timeout=timeout, headers=headers)
            if resp.status_code in no_retry_status:
                return resp
            resp.raise_for_status()  # 返回异常状态码(如服务器报错)时视为请求失败
            return resp
        except Exception as exc:
            last_exc = exc
            logging.exception("%s failed (attempt %d/%d)", label, attempt + 1, max_retries)
            if attempt < max_retries - 1:
                time.sleep(retry_delay)
    if raise_on_exhausted and last_exc is not None:
        raise last_exc
    return None


def _opentargets_post(query: str, variables: dict, timeout: int = 30) -> dict | None:
    """向 Open Targets 数据库发送查询;网络失败时自动重试几次,仍失败则返回 None。

    查询本身报错(返回结果带 errors)不算网络失败,由调用方自行处理。
    """
    import config

    resp = request_with_retry(
        OPENTARGETS_API_URL, json={"query": query, "variables": variables},
        max_retries=config.LIGHT_API_MAX_RETRIES, timeout=timeout,
        label="Open Targets POST",
    )
    return resp.json() if resp is not None else None
