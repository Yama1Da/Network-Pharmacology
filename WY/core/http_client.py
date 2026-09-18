"""访问远程数据库的公共设置:统一的请求重试和访问频率控制。"""

import logging
import threading

import requests

# 记录每个数据源下一次允许访问的时间,
# 用于控制访问频率,避免请求过快被对方拒绝。
_throttle_lock = threading.Lock()
_next_allowed_at: dict[str, float] = {}

# 共享 HTTP 会话(懒加载):带 sqlite 落盘缓存,供 DGIdb/Open Targets/DisGeNET 复用
_session: requests.Session | None = None
_session_lock = threading.Lock()


def get_session() -> requests.Session:
    """返回共享 HTTP 会话;requests_cache 可用时为带 TTL 的落盘缓存会话。
    只缓存 GET/POST 的 200 响应(失败/限流不缓存,重试逻辑不受影响);POST 按完整请求体作缓存键,同一查询命中、不同查询互不干扰。
    """
    global _session
    if _session is None:
        with _session_lock:  # 双重检查加锁:并发首请求只创建一个会话,避免先创建的 sqlite 连接泄漏
            if _session is None:
                try:
                    import os

                    import config
                    from requests_cache import CachedSession

                    os.makedirs(os.path.dirname(config.HTTP_CACHE_FILE), exist_ok=True)
                    _session = CachedSession(
                        config.HTTP_CACHE_FILE, backend="sqlite",
                        expire_after=config.HTTP_CACHE_TTL_HOURS * 3600,
                        allowable_methods=("GET", "POST"), allowable_codes=(200,),
                    )
                except Exception:
                    logging.exception("requests_cache unavailable, falling back to plain session")
                    _session = requests.Session()
    return _session


def clear_http_cache() -> None:
    """清空 HTTP 响应缓存(管理页强制更新前调用,确保抓到新数据)。"""
    try:
        get_session().cache.clear()
    except Exception:
        logging.exception("Failed to clear HTTP cache")


def throttle(host_key: str) -> None:
    """控制对同一数据源的请求间隔,避免访问过快被对方限流;没有设置间隔要求的数据源不限速,直接通过。"""
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
    params: dict | None = None,
    headers: dict | None = None,
    throttle_key: str | None = None,
    retry_delay: float | None = None,
    no_retry_status: tuple = (),
    raise_on_exhausted: bool = False,
    use_cache: bool = False,
    cache_ttl_hours: float | None = None,
) -> requests.Response | None:
    """带重试的 HTTP 请求,项目内所有数据源请求的统一入口;每次失败记录日志并等待后重试,全部失败默认返回 None
    (raise_on_exhausted=True 则抛出最后一次异常)。no_retry_status 中的状态码不视为失败直接原样返回
    (如 STRING 的 404 表示基因全部未识别,属正常业务情况);use_cache=True 走共享缓存会话,cache_ttl_hours 覆盖该请求有效期。"""
    import time

    import config

    if retry_delay is None:
        retry_delay = config.API_RETRY_DELAY
    last_exc: Exception | None = None
    for attempt in range(max_retries):
        try:
            if throttle_key:
                throttle(throttle_key)
            # use_cache=True 走共享缓存会话,否则用 requests 本体;
            # 按方法分别走 .get/.post(而不是 .request),
            # 保持对 requests.post/get 打补丁的测试与工具兼容。
            http = get_session() if use_cache else requests
            # 按请求缓存有效期仅对真正的缓存会话生效(降级为普通会话时忽略)
            extra: dict = {}
            if use_cache and cache_ttl_hours is not None and hasattr(http, "cache"):
                extra["expire_after"] = float(cache_ttl_hours) * 3600
            if method == "get":
                resp = http.get(url, params=params, timeout=timeout, headers=headers, **extra)
            else:
                resp = http.post(url, json=json, data=data, params=params,
                                 timeout=timeout, headers=headers, **extra)
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
