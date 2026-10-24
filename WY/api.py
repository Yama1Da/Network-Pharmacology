"""FastAPI 应用：TS+Vue 前端的 API 层（阶段 1：Logs 页）。

挂载 /api/* 端点与 WY/web/dist 静态托管（SPA 回退），由 `npap serve` 启动，
与 Streamlit 页面共享同一套 core 与缓存。端点只是 core 函数的薄壳：
业务语义（日志级别持久化、文件列表）的唯一定义仍在 core，本文件不做业务决策。
"""

import os

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from WY import config
from WY.core.files import list_log_files
from WY.core.logging_setup import set_logging_level, setup_logging
from WY.core.management import save_ui_settings

setup_logging()

app = FastAPI(title="NPAP API", docs_url=None, redoc_url=None)

_ALLOWED_LEVELS = ("INFO", "DEBUG")
_LOG_TAIL_MAX = 2 * 1024 * 1024  # 单次最多回传 2MB 尾部,防止历史大日志拖垮响应


def _log_file_path(name: str) -> str:
    """校验文件名并解析为日志目录内路径:必须是纯文件名且在 list_log_files 列表内,
    堵住路径穿越(../)与任意文件读取。"""
    if not name or name != os.path.basename(name) or "\\" in name or "/" in name:
        raise HTTPException(status_code=400, detail="invalid log file name")
    listed = set(list_log_files()["File name"].tolist())
    if name not in listed:
        raise HTTPException(status_code=404, detail="no such log file")
    return os.path.join(config.LOG_DIR, name)


@app.get("/api/logs/files")
def logs_files() -> dict:
    """日志文件清单:文件名/修改时间/大小/是否当前(app.log)。"""
    df = list_log_files()
    files = []
    for _, row in df.iterrows():
        path = os.path.join(config.LOG_DIR, row["File name"])
        files.append({
            "name": row["File name"],
            "modified": str(row["Modified time"]),
            "size": os.path.getsize(path) if os.path.exists(path) else 0,
            "is_current": row["File name"] == "app.log",
        })
    return {"files": files}


@app.get("/api/logs/content")
def logs_content(name: str, max_bytes: int = _LOG_TAIL_MAX) -> dict:
    """单个日志文件内容;超过 max_bytes 只回传尾部(掐掉截断的半行,从完整行开始)。"""
    path = _log_file_path(name)
    max_bytes = max(1_000, min(max_bytes, _LOG_TAIL_MAX))
    size = os.path.getsize(path) if os.path.exists(path) else 0
    truncated = size > max_bytes
    with open(path, "rb") as f:
        if truncated:
            f.seek(-max_bytes, os.SEEK_END)
        raw = f.read()
    text = raw.decode("utf-8", errors="replace")
    if truncated:
        text = text.split("\n", 1)[-1]
    return {"name": name, "content": text, "truncated": truncated, "size": size}


@app.get("/api/logs/level")
def logs_level_get() -> dict:
    return {"level": config.LOG_LEVEL}


class LevelBody(BaseModel):
    level: str


@app.post("/api/logs/level")
def logs_level_set(body: LevelBody) -> dict:
    """切换日志级别:与 Streamlit Logs 页同一条链路——改运行时 → 落盘 ui_settings.json
    → 根 logger 即时生效,重启后仍保留。"""
    level = body.level.strip().upper()
    if level not in _ALLOWED_LEVELS:
        raise HTTPException(status_code=400,
                            detail=f"level must be one of {_ALLOWED_LEVELS}")
    config.LOG_LEVEL = level
    save_ui_settings({"LOG_LEVEL": level})
    set_logging_level(level)
    return {"level": level}


# ======== 前端静态托管(SPA):WY/web/dist 由 npm run build 产出 ========
_DIST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web", "dist")

if os.path.isdir(os.path.join(_DIST, "assets")):
    app.mount("/assets", StaticFiles(directory=os.path.join(_DIST, "assets")), name="assets")

    @app.get("/{full_path:path}", include_in_schema=False)
    async def spa(full_path: str) -> FileResponse:
        if ".." in full_path:
            raise HTTPException(status_code=404)
        candidate = os.path.join(_DIST, full_path)
        if full_path and os.path.isfile(candidate):
            return FileResponse(candidate)
        return FileResponse(os.path.join(_DIST, "index.html"))
else:
    @app.get("/", include_in_schema=False)
    async def no_frontend() -> dict:
        return {"hint": "前端未构建:cd WY/web && npm install && npm run build,再重启 npap serve"}
