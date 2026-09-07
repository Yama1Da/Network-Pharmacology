"""分析结果的本地缓存读写:刚算过的结果短期内存起来,重复分析时直接复用。"""

import logging
import os
from collections.abc import Callable

import pandas as pd

from .text_utils import safe_name


def _read_cache_csv(cache_file: str) -> pd.DataFrame:
    # 读取本地缓存的结果表(跳过以 # 开头的旧格式说明行)
    return pd.read_csv(cache_file, encoding="utf-8-sig", comment="#")


def _is_cache_fresh(cache_file: str, ttl_seconds: int | None = None) -> bool:
    """缓存文件是否仍有效:文件存在、且保存时间距今没超过有效期（过期需重新请求,有效期默认取全局配置）。"""
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


def open_csv_cache(cache_file: str, *, fresh: bool = False, label: str = "") -> pd.DataFrame | None:
    """命中缓存则读取返回;未命中返回 None,由调用方自行请求数据。
    fresh=True 时按 TTL 判定有效性(超期视为未命中),False 时只要文件存在即命中;label 用于命中日志。
    """
    hit = _is_cache_fresh(cache_file) if fresh else os.path.exists(cache_file)
    if not hit:
        return None
    logging.info("Loading cached %s", label or os.path.basename(cache_file))
    return _read_cache_csv(cache_file)


def open_named_csv_cache(
    cache_dir: str, prefix: str, name: str, *,
    label: str = "", remember: Callable[[str], None] | None = None,
) -> tuple[str, pd.DataFrame | None]:
    """统一各数据源"按规范名拼缓存路径+读缓存"的样板:路径固定为
    <cache_dir>/<prefix><safe_name(name)>.csv,与既有各处的拼接方式逐字一致
    (绝不改变文件名,否则会作废用户现有缓存);命中返回读到的 DataFrame,未命中返回 None。
    remember 为可选的名称映射回调(如 remember_name_mapping),在读缓存前以原始名调用一次。"""
    csv_path = os.path.join(cache_dir, f"{prefix}{safe_name(name)}.csv")
    if remember is not None:
        remember(name)
    return csv_path, open_csv_cache(csv_path, label=label)


def save_csv_cache(df: pd.DataFrame, cache_file: str, *, skip_empty: bool = True,
                   header_lines: list[str] | None = None) -> bool:
    """把结果表落盘缓存;默认空表不写(避免把失败/空结果当成有效缓存),返回是否实际写出,目录不存在时自动创建。
    header_lines 为可选缓存键行(逐行以 "# " 前缀写在文件最前部),供同名覆盖式缓存在命中时校验内容对应当前查询;
    读取端会跳过 # 行,互不干扰。"""
    if skip_empty and (df is None or df.empty):
        return False
    os.makedirs(os.path.dirname(cache_file) or ".", exist_ok=True)
    with open(cache_file, "w", encoding="utf-8-sig", newline="") as f:
        for line in header_lines or []:
            f.write(f"# {line}\n")
        df.to_csv(f, index=False)
    return True


def read_cache_key(cache_file: str) -> str | None:
    """读取缓存文件首行 "# <键>" 形式的缓存键(save_csv_cache 的 header_lines 写入);没有该行(旧格式文件)返回 None。"""
    try:
        with open(cache_file, "r", encoding="utf-8-sig") as f:
            for line in f:
                line = line.rstrip("\r\n")
                if not line.strip():
                    continue
                if line.startswith("#"):
                    return line[1:].strip()
                return None
    except OSError:
        return None
    return None
