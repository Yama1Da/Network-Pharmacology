"""Logs(阶段 1 端点,原样保留)。"""

from fastapi import APIRouter
import os
from fastapi import HTTPException
from pydantic import BaseModel
from WY import config
from WY.core.files import list_log_files
from WY.core.logging_setup import set_logging_level
from WY.core.management import save_ui_settings
from .common import (_ALLOWED_LEVELS, _LOG_TAIL_MAX)


router = APIRouter()


def _log_file_path(name: str) -> str:
    """校验文件名并解析为日志目录内路径:必须是纯文件名且在 list_log_files 列表内,
    堵住路径穿越(../)与任意文件读取。"""
    if not name or name != os.path.basename(name) or "\\" in name or "/" in name:
        raise HTTPException(status_code=400, detail="invalid log file name")
    if name not in _listed_log_names():
        raise HTTPException(status_code=404, detail="no such log file")
    return os.path.join(config.LOG_DIR, name)


def _listed_log_names() -> set:
    return set(list_log_files()["File name"].tolist())


@router.get("/api/logs/files")
def logs_files() -> dict:
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


@router.get("/api/logs/content")
def logs_content(name: str, max_bytes: int = _LOG_TAIL_MAX) -> dict:
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


@router.get("/api/logs/level")
def logs_level_get() -> dict:
    return {"level": config.LOG_LEVEL}


class LevelBody(BaseModel):
    level: str


@router.post("/api/logs/level")
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
