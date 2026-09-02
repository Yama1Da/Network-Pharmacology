"""结果文件与目录管理：列出缓存/下载/日志文件，并支持 zip 打包下载。"""

import os
from datetime import datetime

import pandas as pd


def list_cache_files(category: str) -> pd.DataFrame:
    """列出指定类别缓存目录（如 "drug" 或 "disease"）下的文件名与创建时间；目录为空时返回空表。"""
    import config

    cache_dir = os.path.join(config.CACHE_DIR, category)
    if not os.path.isdir(cache_dir):
        return pd.DataFrame(columns=["File name", "Created time"])
    rows = []
    for name in sorted(os.listdir(cache_dir)):
        # 跳过更新过程产生的备份/临时文件，不在页面上显示
        if name.endswith((".bak", ".tmp")):
            continue
        path = os.path.join(cache_dir, name)
        if os.path.isfile(path):
            created = datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append({"File name": name, "Created time": created})
    return pd.DataFrame(rows, columns=["File name", "Created time"])


def list_download_dirs() -> pd.DataFrame:
    """列出下载目录下的各次分析子目录及其创建时间；无子目录时返回空表。"""
    import config

    if not os.path.isdir(config.DOWNLOADS_DIR):
        return pd.DataFrame(columns=["Directory", "Created time"])
    rows = []
    for name in sorted(os.listdir(config.DOWNLOADS_DIR)):
        path = os.path.join(config.DOWNLOADS_DIR, name)
        if os.path.isdir(path):
            created = datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append({"Directory": name, "Created time": created})
    return pd.DataFrame(rows, columns=["Directory", "Created time"])


def list_log_files() -> pd.DataFrame:
    """列出日志文件的文件名、创建与修改时间；无日志文件时返回空表。

    app.log 为当前正在写入的日志（即"今天"的），app.log.日期 为历史日志。
    """
    import config

    if not os.path.isdir(config.LOG_DIR):
        return pd.DataFrame(columns=["File name", "Created time", "Modified time"])
    rows = []
    for name in sorted(os.listdir(config.LOG_DIR)):
        path = os.path.join(config.LOG_DIR, name)
        is_log = name.endswith(".log") or name.startswith("app.log.")
        if os.path.isfile(path) and is_log:
            created = datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")
            modified = datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append({"File name": name, "Created time": created, "Modified time": modified})
    return pd.DataFrame(rows, columns=["File name", "Created time", "Modified time"])


def make_ppi_zip(file_paths: dict) -> bytes:
    """
    把一次分析产生的全部结果文件（PPI 网络、富集结果等已落盘的文件）
    打包为 zip，供页面直接下载，解压后即得该次分析的全部成品。
    """
    import io
    import zipfile

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, path in file_paths.items():
            if name == "directory":
                continue
            if path and os.path.exists(path):
                zf.write(path, arcname=os.path.basename(path))
    return buf.getvalue()


def make_dir_zip(directory: str) -> bytes:
    """
    把某个目录下的全部文件打包为 zip，供管理页按目录整体下载；
    目录不存在或为空时返回空 zip。
    """
    import io
    import zipfile

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        if os.path.isdir(directory):
            for name in sorted(os.listdir(directory)):
                path = os.path.join(directory, name)
                if os.path.isfile(path):
                    zf.write(path, arcname=name)
    return buf.getvalue()
