"""分析结果的本地缓存读写:刚算过的结果短期内存起来,重复分析时直接复用。"""

import os

import pandas as pd


def _read_cache_csv(cache_file: str) -> pd.DataFrame:
    # 读取本地缓存的结果表(跳过以 # 开头的旧格式说明行)
    return pd.read_csv(cache_file, encoding="utf-8-sig", comment="#")


def _is_cache_fresh(cache_file: str, ttl_seconds: int | None = None) -> bool:
    """缓存文件是否仍有效:文件存在、且保存时间距今没超过有效期。

    超过有效期视为过期,需要重新请求。有效期默认取全局配置。
    """
    import time
    import config

    if not os.path.exists(cache_file):
        return False
    if ttl_seconds is None:
        ttl_seconds = config.CACHE_TTL_SECONDS
    try:
        return (time.time() - os.path.getmtime(cache_file)) < ttl_seconds
    except OSError:
        return False
