"""远端 API 的 HTTP 公共层：端点常量 + 带重试的请求封装。"""

import logging

import requests

OPENTARGETS_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"


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
