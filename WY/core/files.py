"""结果文件与目录管理：列出缓存/下载/日志文件，并支持 zip 打包下载。"""

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
    """列出下载目录下的各次分析子目录及其最后修改时间；无子目录时返回空表。"""
    from WY import config

    return _list_to_table(config.DOWNLOADS_DIR, name_col="Directory", dirs=True)


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


def make_fig_zip(png_path: str, svg_path: str) -> bytes:
    """把一张图的 PNG + SVG 两个文件打包为小 zip,供下载面板"一图一钮"同时带走位图与矢量版
    （Streamlit 一个按钮只能携带一个文件）。"""
    def write(zf):
        for path in (png_path, svg_path):
            if path and os.path.exists(path):
                zf.write(path, arcname=os.path.basename(path))
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
