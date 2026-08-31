"""文件 / 目录管理：缓存列表、下载目录列表、日志文件列表、zip 打包。"""

import os
from datetime import datetime

import pandas as pd


def list_cache_files(category: str) -> pd.DataFrame:
    """列出 Static/data/cache/<category> 目录下缓存文件的文件名与创建时间。

    Parameters
    ----------
    category : str
        子目录名，如 "drug" 或 "disease"。

    Returns
    -------
    pd.DataFrame
        含 "File name"、"Created time" 两列；目录不存在或为空时返回空表。
    """
    import config

    cache_dir = os.path.join(config.CACHE_DIR, category)
    if not os.path.isdir(cache_dir):
        return pd.DataFrame(columns=["File name", "Created time"])
    rows = []
    for name in sorted(os.listdir(cache_dir)):
        # 跳过更新过程中产生的备份/临时文件（.bak / .tmp），避免残留虚影显示在页面上
        if name.endswith((".bak", ".tmp")):
            continue
        path = os.path.join(cache_dir, name)
        if os.path.isfile(path):
            created = datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append({"File name": name, "Created time": created})
    return pd.DataFrame(rows, columns=["File name", "Created time"])


def list_download_dirs() -> pd.DataFrame:
    """列出 DOWNLOADS_DIR 下的子目录（一次分析一组）及其创建时间。

    Returns
    -------
    pd.DataFrame
        含 "Directory"、"Created time" 两列；目录不存在或无子目录时返回空表。
    """
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
    """列出 LOG_DIR 下的日志文件及其创建/修改时间。

    兼容两种命名：
      - 旧格式：app_YYYY-MM-DD.log
      - TimedRotating 轮转格式：app.log（当前活跃）、app.log.YYYY-MM-DD（历史）

    返回含 "File name"、"Created time"、"Modified time" 三列；
    其中 app.log 为当前正在写入的文件（即"今天"的日志）。

    Returns
    -------
    pd.DataFrame
        含 "File name"、"Created time"、"Modified time" 列；
        目录不存在或无日志文件时返回空表。
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
    将一次分析的全部结果文件打包为 zip 字节（nodes/edges/html/graphml/svg
    及富集 CSV/图等所有已落盘的文件，排除 directory 键）。

    供 st.download_button 直接下载；文件在 zip 内使用各自的文件名，
    配合独立子目录，用户解压后即可获得该次分析的全部成品。
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
    将某个目录下的全部文件打包为 zip 字节（供管理页按目录整体下载）。

    zip 内使用各文件的文件名（不带目录前缀），返回空 zip 当目录不存在或为空。
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
