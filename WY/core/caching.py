"""分析结果的本地缓存读写:刚算过的结果短期内存起来,重复分析时直接复用。"""

import logging
import os
import threading
from collections.abc import Callable

import pandas as pd

# 状态/缓存文件按路径的写锁:同一文件的并发写串行化,不同文件互不阻塞。
# JSON/文本/CSV 写入共用:没有锁的话,两个写者会互写同一个 .tmp 文件,
# Windows 上 os.replace 撞打开句柄直接 PermissionError,POSIX 上会交错写坏内容。
_path_locks: dict[str, threading.Lock] = {}
_locks_guard = threading.Lock()


def _path_lock(path: str) -> threading.Lock:
    with _locks_guard:
        return _path_locks.setdefault(os.path.abspath(path), threading.Lock())


def _atomic_write(tmp_file: str, write_fn: Callable, encoding: str = "utf-8") -> None:
    """先写 .tmp 再原子替换的公共骨架:写过程中途失败(异常/中断)时清理 .tmp,
    避免残留的半截文件后来被"文件存在"判定误当成有效缓存。"""
    try:
        with open(tmp_file, "w", encoding=encoding, newline="") as f:
            write_fn(f)
        os.replace(tmp_file, tmp_file[: -len(".tmp")])
    except BaseException:
        try:
            os.remove(tmp_file)
        except OSError:
            pass
        raise


def load_json_file(path: str, default=None):
    """读取 JSON 状态文件;不存在或损坏时返回 default(损坏按不存在处理,由调用方决定是否重建)。"""
    import json

    try:
        with _path_lock(path), open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        return default
    except Exception:
        logging.exception("Failed to load JSON state file: %s", path)
        return default


def save_json_file(path: str, value, *, indent: int | None = None,
                   sort_keys: bool = False) -> None:
    """把 JSON 值写入状态文件:先写 .tmp 再原子替换,并发写不会弄坏文件;目录不存在时自动创建。"""
    import json

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with _path_lock(path):
        _atomic_write(path + ".tmp",
                      lambda f: json.dump(value, f, ensure_ascii=False,
                                          indent=indent, sort_keys=sort_keys))


def save_text_atomic(path: str, text: str) -> None:
    """文本文件(如下载的 PDB/SDF)原子写:先写 .tmp 再原子替换,中断不会留下半截文件
    被后续"文件存在即复用"的判定永久信任;目录不存在时自动创建。按路径加锁,
    并发写同一文件时不会互写 .tmp(见 _path_locks 注释)。"""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with _path_lock(path):
        _atomic_write(path + ".tmp", lambda f: f.write(text))


def _read_cache_csv(cache_file: str) -> pd.DataFrame | None:
    # 读取本地缓存的结果表(跳过以 # 开头的旧格式说明行)。
    # 空结果缓存(查无数据的药物/疾病也落盘)只有 BOM 没有表头,按空表读回而不是报错;
    # 内容残缺/损坏的缓存删除后按未命中(None)返回,由调用方重新抓取覆盖。
    try:
        return pd.read_csv(cache_file, encoding="utf-8-sig", comment="#")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    except (pd.errors.ParserError, UnicodeDecodeError):
        logging.exception("Cache file unreadable, removing: %s", cache_file)
        try:
            os.remove(cache_file)
        except OSError:
            pass
        return None


def is_cache_fresh(cache_file: str, ttl_seconds: int | None = None) -> bool:
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
    hit = is_cache_fresh(cache_file) if fresh else os.path.exists(cache_file)
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
    from .text_utils import safe_name  # 函数内导入:text_utils 反向依赖本模块的 JSON 助手,避免循环导入

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
    # 先写 .tmp 再原子替换:并发读者不会读到半截 CSV(半截文件曾直接抛 ParserError);
    # 按路径加锁,并发写同一文件时不会互写 .tmp
    def _write(f) -> None:
        for line in header_lines or []:
            f.write(f"# {line}\n")
        df.to_csv(f, index=False)

    with _path_lock(cache_file):
        _atomic_write(cache_file + ".tmp", _write, encoding="utf-8-sig")
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


def clean_expired_cache_files(directory: str, *, exts: tuple[str, ...] = (".csv", ".json"),
                              ttl_seconds: int | None = None) -> int:
    """删除目录里超过缓存有效期(CACHE_TTL_SECONDS)的缓存文件,返回删除数;目录不存在为 0。
    TTL 在各命中判定里只挡复用、不删文件,长期使用目录只增不减;这些结果都已永久保存在
    分析结果目录(<prefix>_data)中,过期文件直接删除即可,下次分析会重算并覆盖。"""
    import time

    import config

    if not os.path.isdir(directory):
        return 0
    if ttl_seconds is None:
        ttl_seconds = config.CACHE_TTL_SECONDS
    cutoff = time.time() - ttl_seconds
    removed = 0
    for name in os.listdir(directory):
        if not name.endswith(tuple(exts)):
            continue
        path = os.path.join(directory, name)
        try:
            if os.path.isfile(path) and os.path.getmtime(path) < cutoff:
                os.remove(path)
                removed += 1
        except OSError:
            logging.exception("Failed to remove expired cache %s", path)
    if removed:
        logging.info("Expired caches removed in %s: %d file(s)", directory, removed)
    return removed
