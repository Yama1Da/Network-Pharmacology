"""结果文件与目录管理：列出缓存/下载/日志文件，并支持 zip 打包下载。"""

import os
from datetime import datetime

import pandas as pd


def _ctime_str(path: str) -> str:
    return datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")


def _list_to_table(directory: str, *, name_col: str, dirs: bool,
                   skip_suffixes: tuple = (), name_pred=None,
                   with_mtime: bool = False, exts: tuple = ()) -> pd.DataFrame:
    """列出 directory 下满足条件的条目名与创建时间;目录不存在时返回空表。

    dirs=True 选子目录,False 选文件;skip_suffixes 跳过指定后缀(如更新过程的备份/临时文件);
    name_pred 为额外的文件名过滤;exts 非空时只保留这些扩展名的文件;
    隐藏文件(以 . 开头,如熔断状态标记)一律不列出;with_mtime=True 时附带修改时间列。
    """
    cols = [name_col, "Created time"] + (["Modified time"] if with_mtime else [])
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
            row = {name_col: name, "Created time": _ctime_str(path)}
            if with_mtime:
                row["Modified time"] = datetime.fromtimestamp(
                    os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")
            rows.append(row)
    return pd.DataFrame(rows, columns=cols)


def list_cache_files(category: str, *, exts: tuple = ()) -> pd.DataFrame:
    """列出指定类别缓存目录（如 "drug" 或 "disease"）下的文件名与创建时间；目录为空时返回空表。

    exts 非空时只列出这些扩展名的文件(如富集目录传 (".csv",) 排除非结果的杂项文件)。
    """
    import config

    return _list_to_table(os.path.join(config.CACHE_DIR, category), name_col="File name",
                          dirs=False, skip_suffixes=(".bak", ".tmp"), exts=exts)


def list_download_dirs() -> pd.DataFrame:
    """列出下载目录下的各次分析子目录及其创建时间；无子目录时返回空表。"""
    import config

    return _list_to_table(config.DOWNLOADS_DIR, name_col="Directory", dirs=True)


def list_log_files() -> pd.DataFrame:
    """列出日志文件的文件名、创建与修改时间；无日志文件时返回空表。

    app.log 为当前正在写入的日志（即"今天"的），app.log.日期 为历史日志。
    """
    import config

    return _list_to_table(
        config.LOG_DIR, name_col="File name", dirs=False, with_mtime=True,
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
    """
    把一次分析产生的全部结果文件（PPI 网络、富集结果等已落盘的文件）
    打包为 zip，供页面直接下载，解压后即得该次分析的全部成品。
    """
    def write(zf):
        for name, path in file_paths.items():
            if name != "directory" and path and os.path.exists(path):
                zf.write(path, arcname=os.path.basename(path))
    return _zip_bytes(write)


def make_dir_zip(directory: str) -> bytes:
    """
    把某个目录下的全部文件打包为 zip，供管理页按目录整体下载；
    目录不存在或为空时返回空 zip。
    """
    def write(zf):
        if os.path.isdir(directory):
            for name in sorted(os.listdir(directory)):
                path = os.path.join(directory, name)
                if os.path.isfile(path):
                    zf.write(path, arcname=name)
    return _zip_bytes(write)
