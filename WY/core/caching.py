"""通用 CSV 缓存读写与 TTL 新鲜度判断（PPI / 富集等带 # 元数据头的缓存共用）。"""

import os

import pandas as pd


def _read_cache_csv(cache_file: str) -> pd.DataFrame:
    # 读取本地缓存 CSV（兼容旧版带 # 元数据头的md5文件；新文件为纯数据 CSV）
    return pd.read_csv(cache_file, encoding="utf-8-sig", comment="#")


def _is_cache_fresh(cache_file: str, ttl_seconds: int | None = None) -> bool:
    """缓存文件是否仍有效：存在且修改时间未超过 TTL。

    ttl_seconds 缺省（None）时取 config.CACHE_TTL_SECONDS。
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
