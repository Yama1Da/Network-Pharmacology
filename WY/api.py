"""FastAPI 应用：TS+Vue 前端的 API 层。

阶段 1(Logs) + 阶段 2(信息管理类六页)。端点只是 core 函数的薄壳:
业务语义(来源刷新编排、日志级别持久化、集合清单)的唯一定义在 core,
本文件做的是 HTTP 形态(路径/校验/打包/长任务提交)。
由 `npap serve` 启动,与 Streamlit 页面共享同一套 core 与缓存。
"""

import os
import shutil

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from WY import config
from WY.core import jobs
from WY.core.enrichment import parse_enrichment_filename
from WY.core.files import (list_cache_files, list_download_dirs, list_log_files,
                           make_dirs_zip, make_ppi_zip)
from WY.core.interactome import (build_interactome, clean_old_version_caches,
                                 interactome_age_days, is_built, is_stale,
                                 require_meta)
from WY.core.logging_setup import set_logging_level, setup_logging
from WY.core.management import (enrichment_display_name, save_ui_settings,
                                update_all_cache_files)
from WY.core.sources import disease_spec, drug_spec
from WY.core.sources import file_path as spec_file_path

setup_logging()

app = FastAPI(title="NPAP API", docs_url=None, redoc_url=None)

_ALLOWED_LEVELS = ("INFO", "DEBUG")
_LOG_TAIL_MAX = 2 * 1024 * 1024  # 单次最多回传 2MB 尾部,防止历史大日志拖垮响应

_ZIP_NAMES = {
    "drug_sources": "Drug_sources.zip",
    "disease_sources": "Disease_sources.zip",
    "ppi": "PPI_files.zip",
    "enrichment": "Enrichment_files.zip",
    "stars": "Starred_results.zip",
    "downloads": "Data_files.zip",
}


# ======== 通用文件集合:六类管理页的清单/打包/删除 ========
# mode "files" 条目为文件(os.remove 删除,make_ppi_zip 打包);
# mode "dirs" 条目为目录(rmtree 删除,make_dirs_zip 打包)。

def _list_mtime(path: str) -> str:
    import datetime

    try:
        return datetime.datetime.fromtimestamp(os.path.getmtime(path)).strftime(
            "%Y-%m-%d %H:%M:%S")
    except OSError:
        return ""


def _items_for(kind: str) -> tuple[str, list[dict], object]:
    """返回 (mode, 条目清单, 路径解析器)。条目含 name/modified/size 与各 kind 的附加字段。"""
    if kind in ("drug_sources", "disease_sources"):
        spec = drug_spec() if kind == "drug_sources" else disease_spec()
        ns = "drug" if kind == "drug_sources" else "disease"
        source_df = list_cache_files(ns)
        merged_df = list_cache_files("merged")
        merged_df = merged_df[merged_df["File name"].str.startswith(f"merged_{ns}_")]
        names = [*source_df["File name"].tolist(), *merged_df["File name"].tolist()]

        def resolver(name: str) -> str:
            return spec_file_path(name, spec)

        items = [{"name": n, "modified": _list_mtime(resolver(n)),
                  "size": os.path.getsize(resolver(n)) if os.path.exists(resolver(n)) else 0,
                  "is_merged": n.startswith(spec["merged_prefix"]),
                  "source": next((p.rstrip("_") for p in spec["source_prefixes"]
                                  if n.startswith(p)), "merged")}
                 for n in names]
        return "files", items, resolver

    if kind == "ppi":
        df = list_cache_files("ppi")

        def resolver(name: str) -> str:
            return os.path.join(config.PPI_CACHE_DIR, name)

        items = [{"name": n, "modified": m, "size": os.path.getsize(resolver(n))
                  if os.path.exists(resolver(n)) else 0}
                 for n, m in zip(df["File name"].tolist(), df["Modified time"].tolist())]
        return "files", items, resolver

    if kind == "enrichment":
        df = list_cache_files("enrichment")

        def resolver(name: str) -> str:
            return os.path.join(config.ENRICHMENT_CACHE_DIR, name)

        items = []
        for n, m in zip(df["File name"].tolist(), df["Modified time"].tolist()):
            info = parse_enrichment_filename(n) or {}
            items.append({
                "name": n, "modified": m,
                "size": os.path.getsize(resolver(n)) if os.path.exists(resolver(n)) else 0,
                "display": enrichment_display_name(n),
                "source": info.get("source", ""),
                "side": info.get("side") or "",
                "gsea": info.get("kind") == "gsea",
            })
        return "files", items, resolver

    if kind == "stars":
        star_dir = config.MD_STAR_DIR
        names = sorted(n for n in os.listdir(star_dir)
                       if os.path.isdir(os.path.join(star_dir, n))) \
            if os.path.isdir(star_dir) else []
        return "dirs", [{"name": n, "modified": _list_mtime(os.path.join(star_dir, n))}
                        for n in names], lambda n: os.path.join(star_dir, n)

    if kind == "downloads":
        df = list_download_dirs()
        excluded = {os.path.basename(config.MD_RESULT_DIR),
                    os.path.basename(config.MD_STAR_DIR)}
        df = df[~df["Directory"].isin(excluded)]
        return "dirs", [{"name": d, "modified": m}
                        for d, m in zip(df["Directory"].tolist(),
                                        df["Modified time"].tolist())], \
            lambda n: os.path.join(config.DOWNLOADS_DIR, n)

    raise HTTPException(status_code=404, detail=f"unknown collection '{kind}'")


@app.get("/api/collections/{kind}")
def collection_items(kind: str) -> dict:
    mode, items, _ = _items_for(kind)
    return {"kind": kind, "mode": mode, "items": items}


class NamesBody(BaseModel):
    names: list[str]


def _validated_paths(kind: str, names: list[str]) -> list[tuple[str, str]]:
    """把选中名解析为 (name, 路径),并确认都在集合清单内(防越权删除/打包)。"""
    mode, items, resolver = _items_for(kind)
    known = {it["name"] for it in items}
    out = []
    for n in names:
        if n != os.path.basename(n) or n not in known:
            raise HTTPException(status_code=400, detail=f"invalid name: {n!r}")
        out.append((n, resolver(n)))
    return out


@app.post("/api/collections/{kind}/zip")
def collection_zip(kind: str, body: NamesBody) -> Response:
    mode, _, resolver = _items_for(kind)
    pairs = _validated_paths(kind, body.names)
    if not pairs:
        raise HTTPException(status_code=400, detail="empty selection")
    if mode == "files":
        data = make_ppi_zip({n: p for n, p in pairs})
    else:
        data = make_dirs_zip({n: p for n, p in pairs})
    return Response(
        content=data,
        media_type="application/zip",
        headers={"Content-Disposition": f'attachment; filename="{_ZIP_NAMES[kind]}"'})


@app.post("/api/collections/{kind}/delete")
def collection_delete(kind: str, body: NamesBody) -> dict:
    mode, _, resolver = _items_for(kind)
    pairs = _validated_paths(kind, body.names)
    deleted, failed = 0, []
    for n, path in pairs:
        try:
            if mode == "dirs":
                shutil.rmtree(path, ignore_errors=True)
                if os.path.isdir(path):  # rmtree 静默失败(Windows 目录被占用)按失败处理
                    raise OSError(path)
            else:
                os.remove(path)
        except OSError:
            failed.append(n)
        else:
            deleted += 1
    return {"deleted": deleted, "failed": failed}


# ======== 来源刷新(药物/疾病管理页):单文件与批量,经 job 提交-轮询 ========

class RefreshBody(BaseModel):
    filenames: list[str] | None = None   # None=全部更新;列表=更新选中
    hours: int | None = None             # 缓存 TTL;缺省按页面惯例(drug 12/1,disease 默认)


def _refresh_hours(kind: str, filenames: list[str] | None) -> int | None:
    if kind == "drug":
        return 12 if filenames is None else 1
    return None


@app.post("/api/sources/{kind}/refresh")
def sources_refresh(kind: str, body: RefreshBody) -> dict:
    if kind not in ("drug", "disease"):
        raise HTTPException(status_code=404, detail=f"unknown source kind '{kind}'")
    spec = drug_spec() if kind == "drug" else disease_spec()
    hours = body.hours if body.hours is not None else _refresh_hours(kind, body.filenames)
    filenames = body.filenames

    def task(progress) -> dict:
        if filenames is None:
            results = update_all_cache_files(spec, cache_hours=hours, on_progress=progress)
        else:
            results = update_all_cache_files(spec, cache_hours=hours, filenames=filenames,
                                             on_progress=progress)
        return {"rows": results.to_dict(orient="records")}

    job_id = jobs.submit(task, f"{kind} sources refresh")
    return {"job_id": job_id}


# ======== 背景互作组:状态与一键重建(Data_Management 的另一半) ========

@app.get("/api/interactome")
def interactome_status() -> dict:
    if not is_built("Human"):
        return {"built": False, "version": config.STRING_VERSION}
    meta = require_meta("Human")
    age = interactome_age_days("Human")
    return {"built": True,
            "nodes": meta["n_nodes"], "edges": meta["n_edges"],
            "version": meta["string_version"], "confidence": meta["confidence"],
            "built_at": str(meta.get("built_at", ""))[:10],
            "age": age if age is not None else None,
            "stale": is_stale("Human")}


@app.post("/api/interactome/build")
def interactome_build() -> dict:
    def task(progress) -> dict:
        meta_new = build_interactome(progress_cb=progress)
        cleaned = clean_old_version_caches("Human")
        return {"nodes": meta_new["n_nodes"], "edges": meta_new["n_edges"],
                "cleaned": cleaned}

    return {"job_id": jobs.submit(task, "interactome build")}


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str) -> dict:
    job = jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="no such job")
    return job


# ======== Logs(阶段 1 端点,原样保留) ========

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


@app.get("/api/logs/files")
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


@app.get("/api/logs/content")
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
