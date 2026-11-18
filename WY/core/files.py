"""结果文件与目录管理：列出缓存/下载/日志文件，并支持 zip 打包下载。"""

import logging
import os
from datetime import datetime

import pandas as pd


def _mtime_str(path: str) -> str:
    return datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")


def _list_to_table(directory: str, *, name_col: str, dirs: bool,
                   skip_suffixes: tuple = (), name_pred=None,
                   exts: tuple = ()) -> pd.DataFrame:
    """列出 directory 下满足条件的条目名与最后修改时间;目录不存在时返回空表,隐藏文件(以 . 开头,如熔断状态标记)一律不列出。
    用 mtime 而非 getctime:后者在 Windows 上是创建时间、在 Linux 上是 inode 变更时间,跨平台语义不一致。"""
    cols = [name_col, "Modified time"]
    if not os.path.isdir(directory):
        return pd.DataFrame(columns=cols)
    rows = []
    for name in sorted(os.listdir(directory)):
        if name.startswith("."):
            continue
        if name.endswith(skip_suffixes):
            continue
        if exts and not name.endswith(exts):
            continue
        if name_pred is not None and not name_pred(name):
            continue
        path = os.path.join(directory, name)
        if os.path.isdir(path) == dirs:
            rows.append({name_col: name, "Modified time": _mtime_str(path)})
    return pd.DataFrame(rows, columns=cols)


def list_cache_files(category: str, *, exts: tuple = ()) -> pd.DataFrame:
    """列出指定类别缓存目录（如 "drug" 或 "disease"）下的文件名与最后修改时间；目录为空时返回空表。"""
    from WY import config

    return _list_to_table(os.path.join(config.CACHE_DIR, category), name_col="File name",
                          dirs=False, skip_suffixes=(".bak", ".tmp"), exts=exts)


def list_download_dirs() -> pd.DataFrame:
    """分析产物目录的运行清单(2026-11-18 聚合后口径):PPI 运行目录聚合在
    ppi-results/ 一层下,此处走进去按 run 列出(名称带 "ppi-results/" 前缀);
    synergy-results 由协同分析结果页专管,不出现在本清单;GEOresult 由 ML 管理
    页专管,顶层项保留;MDresult/star 本页从来不管,排除。"""
    from WY import config

    downloads = config.DOWNLOADS_DIR
    excluded_top = {"MDresult", "star", "synergy-results"}
    category = {"ppi-results"}
    rows: list[dict] = []
    if not os.path.isdir(downloads):
        return pd.DataFrame(columns=["Directory", "Modified time"])
    for name in sorted(os.listdir(downloads)):
        path = os.path.join(downloads, name)
        if not os.path.isdir(path) or name in excluded_top or name.startswith("."):
            continue
        if name in category:
            for child in sorted(os.listdir(path)):
                cpath = os.path.join(path, child)
                if os.path.isdir(cpath) and not child.startswith("."):
                    rows.append({"Directory": f"{name}/{child}",
                                 "Modified time": _mtime_str(cpath)})
        else:
            rows.append({"Directory": name, "Modified time": _mtime_str(path)})
    return pd.DataFrame(rows, columns=["Directory", "Modified time"])


def migrate_results_layout() -> None:
    """downloads 顶层归类迁移(2026-11-18 聚合,幂等,由 API 启动兜底调用一次):
    PPI 运行目录(*_data)移入 ppi-results/,协同运行目录(synergy_*)移入
    synergy-results/;GEOresult/MDresult/star/readme.md.md 与两个类别目录本身不动;
    未知目录保守不搬(不误移用户手工放置的内容);同名目标已存在时跳过不覆盖。"""
    from WY import config

    downloads = config.DOWNLOADS_DIR
    if not os.path.isdir(downloads):
        return
    keep = {"GEOresult", "MDresult", "star", "ppi-results", "synergy-results"}
    for name in sorted(os.listdir(downloads)):
        path = os.path.join(downloads, name)
        if not os.path.isdir(path) or name in keep or name.startswith("."):
            continue
        if name.endswith("_data"):
            dest = os.path.join(downloads, "ppi-results", name)
        elif name.startswith("synergy_"):
            dest = os.path.join(downloads, "synergy-results", name)
        else:
            continue
        try:
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            if not os.path.exists(dest):
                os.replace(path, dest)
                logging.info("results layout migration: %s -> %s", name, dest)
        except OSError:
            logging.exception("results layout migration failed for %s", name)


def list_log_files() -> pd.DataFrame:
    """列出日志文件的文件名与最后修改时间；无日志文件时返回空表。
    app.log 为当前正在写入的日志（即"今天"的），app.log.日期 为历史日志。
    """
    from WY import config

    return _list_to_table(
        config.LOG_DIR, name_col="File name", dirs=False,
        name_pred=lambda n: n.endswith(".log") or n.startswith("app.log."),
    )


def _zip_bytes(writer) -> bytes:
    """在内存 zip 中执行 writer(zf) 并返回字节流,供页面直接下载。"""
    import io
    import zipfile

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        writer(zf)
    return buf.getvalue()


def make_ppi_zip(file_paths: dict) -> bytes:
    """把一次分析产生的全部结果文件（PPI 网络、富集结果等已落盘的文件）打包为 zip,供页面直接下载。"""
    def write(zf):
        for name, path in file_paths.items():
            if name != "directory" and path and os.path.exists(path):
                zf.write(path, arcname=os.path.basename(path))
    return _zip_bytes(write)


def make_bytes_zip(files: dict) -> bytes:
    """把内存中的若干文件（文件名 → 字节）打包为 zip,供不落盘的数据(如主页检索结果)直接下载。"""
    def write(zf):
        for name, data in files.items():
            zf.writestr(name, data)
    return _zip_bytes(write)


def make_dir_zip(directory: str) -> bytes:
    """把某个目录下的全部文件打包为 zip,供管理页按目录整体下载;目录不存在或为空时返回空 zip。"""
    def write(zf):
        if os.path.isdir(directory):
            for name in sorted(os.listdir(directory)):
                path = os.path.join(directory, name)
                if os.path.isfile(path):
                    zf.write(path, arcname=name)
    return _zip_bytes(write)


def make_dirs_zip(dir_paths: dict) -> bytes:
    """把多个目录打包进同一个 zip(每个目录作为 zip 内的一层文件夹,保留子目录结构),供管理页批量下载;空字典返回空 zip。"""
    def write(zf):
        for name, directory in dir_paths.items():
            if not os.path.isdir(directory):
                continue
            for root, _, files in os.walk(directory):
                for fn in files:
                    fp = os.path.join(root, fn)
                    # 相对各自父目录取弧名,得到 "目录名/文件名(相对路径)"
                    zf.write(fp, arcname=os.path.relpath(fp, os.path.dirname(directory)))
    return _zip_bytes(write)
