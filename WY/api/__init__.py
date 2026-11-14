"""FastAPI 应用：TS+Vue 前端的 API 层。

阶段 1(Logs) + 阶段 2(信息管理类六页)。端点只是 core 函数的薄壳:
业务语义(来源刷新编排、日志级别持久化、集合清单)的唯一定义在 core,
本文件做的是 HTTP 形态(路径/校验/打包/长任务提交)。
路由按页面拆分在各子模块(WY/api/<页>.py),本 __init__ 只做装配与 SPA 托管。
由 `npap serve` 启动,与 Streamlit 页面共享同一套 core 与缓存。
"""
import logging
import os
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from WY.core import transfer
from WY.core.logging_setup import setup_logging


setup_logging()

# 对接缓存清理与富集缓存改名:原在 Streamlit 页进入时做一次,API 进程启动时做一次
try:
    from WY.core.docking import clean_docking_cache
    from WY.core.enrichment import migrate_legacy_cache_names
    clean_docking_cache()
    migrate_legacy_cache_names()
except Exception:
    logging.exception("docking cache cleanup at startup failed")

# 数据迁移暂存清扫:进程崩溃/重启残留的导入暂存与过期导出包,启动时兜底删除
try:
    transfer.cleanup_transfer_staging()
except Exception:
    logging.exception("transfer staging cleanup at startup failed")

app = FastAPI(title="NPAP API", docs_url=None, redoc_url=None)


# 各页路由:include 顺序与原单文件定义顺序一致,SPA 通配路由最后注册
# 注:上方 `from WY.core import jobs, transfer` 已把名字 transfer 占为 WY.core.transfer,
# `from . import transfer` 会被同名属性短路,故该子模块须显式路径导入并别名。
from . import collections
import WY.api.transfer as transfer_routes
from . import sources
from . import interactome
from . import logs
from . import ppi
from . import stats
from . import settings
from . import retrieve
from . import docking
from . import ml

app.include_router(collections.router)
app.include_router(transfer_routes.router)
app.include_router(sources.router)
app.include_router(interactome.router)
app.include_router(logs.router)
app.include_router(ppi.router)
app.include_router(stats.router)
app.include_router(settings.router)
app.include_router(retrieve.router)
app.include_router(docking.router)
app.include_router(ml.router)



# 拆包后 __file__ 在 WY/api/,dist 在上一级 WY/web/dist
_DIST = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                     "web", "dist")

if os.path.isdir(os.path.join(_DIST, "assets")):
    app.mount("/assets", StaticFiles(directory=os.path.join(_DIST, "assets")), name="assets")

    @app.get("/{full_path:path}", include_in_schema=False)
    async def spa(full_path: str) -> FileResponse:
        # realpath 包含校验(与 /api/artifacts 同款,2026-11 P1-A):仅挡 ".."
        # 挡不住 Windows 盘符绝对路径——os.path.join(_DIST, "C:/...") 会整段
        # 丢弃 dist 前缀,任意文件可读
        base = os.path.realpath(_DIST)
        full = os.path.realpath(os.path.join(base, full_path))
        if not (full == base or full.startswith(base + os.sep)):
            raise HTTPException(status_code=404)
        if full_path and os.path.isfile(full):
            return FileResponse(full)
        # index.html 不缓存:构建产物文件名带哈希,浏览器缓存旧 index 会引用已删除的旧资产
        return FileResponse(os.path.join(_DIST, "index.html"),
                            headers={"Cache-Control": "no-cache"})
else:
    @app.get("/", include_in_schema=False)
    async def no_frontend() -> dict:
        return {"hint": "前端未构建:cd WY/web && npm install && npm run build,再重启 npap serve"}
