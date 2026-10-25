"""FastAPI 应用：TS+Vue 前端的 API 层。

阶段 1(Logs) + 阶段 2(信息管理类六页)。端点只是 core 函数的薄壳:
业务语义(来源刷新编排、日志级别持久化、集合清单)的唯一定义在 core,
本文件做的是 HTTP 形态(路径/校验/打包/长任务提交)。
由 `npap serve` 启动,与 Streamlit 页面共享同一套 core 与缓存。
"""

import os
import re
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
from WY.core.management import (enrichment_display_name, save_api_env,
                                save_ui_settings, update_all_cache_files)
from WY.core.sources import disease_spec, drug_spec
from WY.core.text_utils import safe_name
from WY.core.analysis import run_ppi_analysis
from WY.core.sources import file_path as spec_file_path

setup_logging()

# 对接缓存清理与富集缓存改名:原在 Streamlit 页进入时做一次,API 进程启动时做一次
try:
    from WY.core.docking import clean_docking_cache
    from WY.core.enrichment import migrate_legacy_cache_names
    clean_docking_cache()
    migrate_legacy_cache_names()
except Exception:
    logging.exception("docking cache cleanup at startup failed")

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


# ======== PPI 分析(阶段 3a:Data_Network 页) ========
# 产物直出:分析的全部图表/表格由 run_ppi_analysis 经 save_ppi_results 落盘到
# downloads/<prefix>_data/,本层只做序列化与文件伺服,不重绘任何图。
# 交互区重跑(模块/GSEA)是无状态的:从落盘的 nodes/edges CSV 重建图,与网页展示一致。


class _JobProgress:
    """把 core.analysis.AnalysisProgress 协议适配到 core.jobs 的进度回调:
    阶段名/状态进 extra 字段,前端据此渲染阶段图标行。"""

    def __init__(self, cb):
        self._cb = cb

    def set(self, key: str, state: str) -> None:
        self._cb(text=f"{key}:{state}", step=key, state=state)

    def advance(self, fraction: float) -> None:
        self._cb(fraction=fraction)


def _df_rows(df, limit: int | None = None) -> list[dict]:
    if df is None or len(df) == 0:
        return []
    if limit is not None:
        df = df.head(limit)
    return df.to_dict(orient="records")


def _serialize_ppi(ppi: dict) -> dict:
    """把 PpiAnalysisResult 压成 JSON 安全的结果:图表/表格引用落盘产物的
    downloads 相对路径,大表截断。"""
    saved = ppi["saved"]

    def rel(path) -> str | None:
        if not path:
            return None
        return os.path.relpath(path, config.DOWNLOADS_DIR).replace(os.sep, "/")

    artifacts = {k: rel(v) for k, v in saved.items()
                 if isinstance(v, str) and k != "directory"}

    views = ["interactive", "static"]
    if ppi.get("string_fig") is not None:
        views.append("string")
    if ppi.get("cyto_fig") is not None:
        views.append("cytoscape")

    expansion_df = ppi.get("expansion_df")
    expansion = None
    if expansion_df is not None:
        expansion = {
            "seeds": int((expansion_df["source"] == "seed").sum()),
            "expanded": int((expansion_df["source"] == "expanded").sum()),
        }

    gsea = {}
    for side, info in (ppi.get("gsea") or {}).items():
        gsea[side] = {"reason": info.get("reason", ""), "stats": info.get("stats", {}),
                      "rows": _df_rows(info.get("df"), 100)}

    def _cap_rwr(res: dict | None) -> dict | None:
        if not res:
            return res
        capped = {k: v for k, v in res.items() if k not in ("seeds", "top_nonseed")}
        capped["seeds"] = (res.get("seeds") or [])[:50]
        capped["top_nonseed"] = (res.get("top_nonseed") or [])[:100]
        return capped

    rwr = {side: _cap_rwr(res) for side, res in (ppi.get("rwr") or {}).items()}

    modules = ppi.get("modules")
    modules_out = None
    if modules is not None:
        modules_out = {"summary": _df_rows(modules.get("summary_df")),
                       "failures": modules.get("failures", []),
                       "no_sig": modules.get("no_sig", [])}

    enr_df = ppi.get("enrichment_df")
    return {
        "prefix": ppi["prefix"],
        "directory_rel": rel(saved["directory"]),
        "params": list(ppi["params"]),
        "hub_rule": ppi.get("hub_rule", ""),
        "hub_method_eff": ppi.get("hub_method_eff", ""),
        "genes_count": len(ppi["genes"]),
        "nodes": len(ppi["metrics_df"]),
        "edges": len(ppi["edges_df"]),
        "hubs": len(ppi["hub_df"]),
        "expansion": expansion,
        "hub_table": _df_rows(ppi["hub_df"], 50),
        "views": views,
        "artifacts": artifacts,
        "enrichment": {
            "source": ppi.get("enrichment_source", ""),
            "reason": ppi.get("enrichment_reason", ""),
            "rows": _df_rows(enr_df, 200),
        },
        "gsea": gsea,
        "proximity": ppi.get("proximity"),
        "rwr": rwr,
        "modules": modules_out,
    }


class AnalyzeBody(BaseModel):
    mode: str                      # "names" | "genes"
    drug: str = ""
    disease: str = ""
    genes_text: str = ""
    score: float = 0.45
    network_type: str = "full"
    rwr_expand: int = 0
    hub_top_n: int = 10
    hub_method: str = "degree"


@app.post("/api/ppi/analyze")
def ppi_analyze(body: AnalyzeBody) -> dict:
    if body.mode not in ("names", "genes"):
        raise HTTPException(status_code=400, detail="mode must be names|genes")
    if body.mode == "names" and (not body.drug.strip() or not body.disease.strip()):
        raise HTTPException(status_code=400, detail="drug and disease are required in names mode")

    def task(progress) -> dict:
        result = run_ppi_analysis(
            body.mode, body.drug.strip(), body.disease.strip(), body.genes_text.strip(),
            float(body.score),
            network_type=body.network_type, rwr_expand=int(body.rwr_expand),
            hub_top_n=int(body.hub_top_n), hub_method=body.hub_method,
            progress=_JobProgress(progress))
        return _serialize_ppi(result)

    return {"job_id": jobs.submit(task, f"ppi analyze {body.drug}/{body.disease}")}


@app.get("/api/artifacts")
def artifact(rel: str, download: int = 0) -> FileResponse:
    """伺服 downloads/ 下分析产物(PNG/SVG/HTML/CSV/graphml/json)。
    rel 为 downloads 相对路径,realpath 校验防穿越。"""
    base = os.path.realpath(config.DOWNLOADS_DIR)
    full = os.path.realpath(os.path.join(base, rel))
    if not (full == base or full.startswith(base + os.sep)) or not os.path.isfile(full):
        raise HTTPException(status_code=404, detail="no such artifact")
    media = {
        ".png": "image/png", ".svg": "image/svg+xml", ".html": "text/html",
        ".csv": "text/csv", ".graphml": "application/xml", ".json": "application/json",
    }.get(os.path.splitext(full)[1].lower(), "application/octet-stream")
    fname = os.path.basename(full)
    disposition = "attachment" if download else "inline"
    return FileResponse(full, media_type=media,
                        headers={"Content-Disposition": f'{disposition}; filename="{fname}"'})


class RerunModulesBody(BaseModel):
    directory_rel: str   # downloads 相对路径(分析目录)
    prefix: str


@app.post("/api/ppi/rerun-modules")
def ppi_rerun_modules(body: RerunModulesBody) -> dict:
    """交互区手动重跑模块检测:从落盘 nodes/edges CSV 无状态重建图,与展示一致。"""
    base = os.path.realpath(config.DOWNLOADS_DIR)
    directory = os.path.realpath(os.path.join(base, body.directory_rel))
    if not directory.startswith(base + os.sep) or not os.path.isdir(directory):
        raise HTTPException(status_code=404, detail="no such analysis directory")

    def task(progress) -> dict:
        import pandas as pd

        from WY.core.network_algorithms import compute_module_analysis
        from WY.core.ppi import build_ppi_graph, save_module_computation

        progress(fraction=0.05, text="loading graph from artifacts")
        nodes = pd.read_csv(os.path.join(directory, f"{body.prefix}_nodes.csv"),
                            encoding="utf-8-sig", comment="#")
        edges = pd.read_csv(os.path.join(directory, f"{body.prefix}_edges.csv"),
                            encoding="utf-8-sig", comment="#")
        graph = build_ppi_graph(edges, nodes["gene"].astype(str).tolist())
        progress(fraction=0.2, text="louvain + module enrichment")
        organism = config.SPECIES_MAP["Human"]["gprofiler"]
        computed = compute_module_analysis(graph, organism, label_prefix=body.prefix)
        save_module_computation(directory, body.prefix, computed)
        return {"summary": _df_rows(computed["summary_df"]),
                "failures": computed["failures"], "no_sig": computed["no_sig"]}

    return {"job_id": jobs.submit(task, "ppi rerun modules")}


class RerunGseaBody(BaseModel):
    side: str            # "drug" | "disease"
    drug: str
    disease: str
    prefix: str


@app.post("/api/ppi/rerun-gsea")
def ppi_rerun_gsea(body: RerunGseaBody) -> dict:
    if body.side not in ("drug", "disease"):
        raise HTTPException(status_code=400, detail="side must be drug|disease")

    def task(progress) -> dict:
        from WY.core.analysis import gsea_prerank_side

        progress(fraction=0.1, text=f"gsea prerank ({body.side})")
        gdf, rstats, reason = gsea_prerank_side(body.side, body.drug, body.disease, body.prefix)
        return {"reason": reason, "stats": rstats, "rows": _df_rows(gdf, 100)}

    return {"job_id": jobs.submit(task, "ppi rerun gsea")}


# ======== 设置(App.py 齿轮职责的 TS 化):密钥/工具路径/STRING 版本/口袋阈值 ========

class ApiEnvBody(BaseModel):
    changes: dict


@app.get("/api/settings")
def settings_get() -> dict:
    """设置快照:密钥只回布尔(不回值),路径明文,STRING 版本与已知清单,口袋阈值。"""
    return {
        "keys": {"disgenet": bool(config.DisGeNET_API_KEY),
                 "ncbi": bool(config.NCBI_API_KEY)},
        "vina_dir": config.VINA_DIR,
        "obabel_dir": config.OBABEL_DIR,
        "cli_results_dir": config.CLI_RESULTS_DIR,
        "string_version": config.STRING_VERSION,
        "known_versions": list(config.STRING_KNOWN_VERSIONS),
        "pocket_cutoff": config.POCKET_CUTOFF,
    }


@app.post("/api/settings/api-env")
def settings_api_env(body: ApiEnvBody) -> dict:
    """密钥/工具路径写入:走 core.management.save_api_env(白名单+消毒+set_key+运行时回填)。"""
    return {"applied": save_api_env(body.changes)}


class StringVersionBody(BaseModel):
    version: str


@app.post("/api/settings/string-version")
def settings_string_version(body: StringVersionBody) -> dict:
    """STRING 版本切换:格式校验 → 下载站探活 → 激活(core.activate_string_version,
    版本常量/派生 URL/互作组运行时缓存一并处理)。"""
    from WY.core.interactome import string_version_exists
    from WY.core.management import activate_string_version

    target = config.normalized_string_version(body.version)
    if not target:
        raise HTTPException(status_code=400, detail=f"bad version format: {body.version!r}")
    if target == config.STRING_VERSION:
        return {"version": target, "changed": False}
    if not string_version_exists(target):
        raise HTTPException(status_code=400, detail=f"version {target} not found on download site")
    activate_string_version(target)
    return {"version": target, "changed": True}


class UiSettingsBody(BaseModel):
    changes: dict


_ALLOWED_UI_KEYS = ("POCKET_CUTOFF",)  # 段 2 设置页先放开这一项,其余齿轮参数后续逐个开放


@app.post("/api/settings/ui")
def settings_ui(body: UiSettingsBody) -> dict:
    changes = {k: v for k, v in body.changes.items() if k in _ALLOWED_UI_KEYS}
    if not changes:
        raise HTTPException(status_code=400,
                            detail=f"no allowed keys (allowed: {_ALLOWED_UI_KEYS})")
    for k, v in changes.items():
        setattr(config, k, v)
    save_ui_settings(changes)
    return {"pocket_cutoff": config.POCKET_CUTOFF}


# ======== 主页四库检索(App.py 检索职责的 TS 化) ========

class RetrieveBody(BaseModel):
    drug: str
    disease: str


_VENN_CACHE: dict[str, dict] = {}   # job_id → {"png": bytes}(任务线程写入,端点读取)
_VENN_ORDER: list[str] = []


def _render_venn_png(drug_set: set, disease_set: set, drug_name: str, disease_name: str) -> bytes:
    """venn2 服务端渲染(与 App.py 主页同款;线程安全:用 Figure API 而非 pyplot 全局状态)。
    图内文字保持英文:matplotlib 默认字体不支持中文,中文会显示成方框。"""
    import io

    from matplotlib.figure import Figure
    from matplotlib_venn import venn2

    fig = Figure(figsize=(6, 4))
    ax = fig.add_subplot()
    v = venn2(
        subsets=(len(drug_set - disease_set), len(disease_set - drug_set),
                 len(drug_set & disease_set)),
        set_labels=(f"{drug_name} targets", f"{disease_name} targets"),
        ax=ax,
    )
    if v and v.get_label_by_id("11"):
        v.get_label_by_id("11").set_fontsize(12)  # 交集区域字号
    ax.set_title(f"{drug_name} vs {disease_name} — Venn Diagram", fontsize=14)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    return buf.getvalue()


@app.post("/api/targets/retrieve")
def targets_retrieve(body: RetrieveBody) -> dict:
    """四库并发检索(ChEMBL/DGIdb/OpenTargets/DisGeNET),写入来源与 merged 缓存;
    分析页的 targets 阶段将直接命中本次缓存。完成后附带主页 venn 图(内存缓存,按 job_id 伺服)。"""
    if not body.drug.strip() or not body.disease.strip():
        raise HTTPException(status_code=400, detail="drug and disease are required")

    venn_holder: dict = {}

    def task(progress) -> dict:
        from WY.core.targets import get_targets_concurrent

        def cb(done: list, total: int) -> None:
            progress(fraction=(len(done) / total) if total else None,
                     text=f"四库检索中 · 已完成 {len(done)}/{total}", done=done)

        drug_df, disease_df, failed = get_targets_concurrent(
            body.drug.strip(), body.disease.strip(), progress_cb=cb)
        drug_set = (set(drug_df["approvedSymbol"].dropna().astype(str))
                    if drug_df is not None and not drug_df.empty else set())
        disease_set = (set(disease_df["approvedSymbol"].dropna().astype(str))
                       if disease_df is not None and not disease_df.empty else set())
        if drug_set or disease_set:
            venn_holder["png"] = _render_venn_png(drug_set, disease_set,
                                                  body.drug.strip(), body.disease.strip())
        common = sorted(drug_set & disease_set)
        return {"drug_count": len(drug_set), "disease_count": len(disease_set),
                "common_count": len(common),
                "common_genes": common[:200],
                "failed": list(failed)}

    job_id = jobs.submit(task, f"targets retrieve {body.drug}/{body.disease}")
    _VENN_CACHE[job_id] = venn_holder
    _VENN_ORDER.append(job_id)
    if len(_VENN_ORDER) > 10:  # 只保留最近 10 次检索的 venn,防内存增长
        old = _VENN_ORDER.pop(0)
        _VENN_CACHE.pop(old, None)
    return {"job_id": job_id}


@app.get("/api/targets/venn")
def targets_venn(job_id: str) -> Response:
    holder = _VENN_CACHE.get(job_id)
    png = (holder or {}).get("png")
    if not png:
        raise HTTPException(status_code=404, detail="no venn for this job")
    return Response(content=png, media_type="image/png")


# ======== 分子对接(阶段 3b:Molecular_Docking 页) ========
# 多步向导的无状态 API:会话状态存 core.docking_session(镜像原 dock 会话键),
# 每个端点对应一步;网络重的步骤(结构批量检查/受体准备/vina)经 job 提交-轮询。
# 3D 预览走产物直出:py3Dmol 在服务端渲染为 HTML,前端 iframe 内嵌。

from functools import lru_cache

from WY.core import docking_session


def _embed_centered(html: str) -> str:
    """py3Dmol 生成的 HTML 注入 iframe 适配样式:去 body 边距、禁滚动条、
    视图在框内水平垂直居中(用户拖转/滚轮缩放不受影响)。"""
    style = ("<style>html,body{margin:0;padding:0;overflow:hidden;height:100%}"
             "body{display:flex;align-items:center;justify-content:center;"
             "background:#ffffff}</style>")
    if "<head>" in html:
        return html.replace("<head>", "<head>" + style, 1)
    return style + html


@lru_cache(maxsize=256)
def _pdb_candidates_cached(gene: str, allow_af: bool):
    """单基因受体候选查询的进程内缓存(原页面 @st.cache_data(ttl=1h) 的等价物)。"""
    from WY.core.docking import fetch_receptor_pdb_candidates

    return fetch_receptor_pdb_candidates(gene, max_hits=20, allow_alphafold=allow_af)


def _session(sid: str) -> dict:
    dock = docking_session.get(sid)
    if dock is None:
        raise HTTPException(status_code=404, detail="docking session expired, create a new one")
    return dock


@app.post("/api/docking/session")
def docking_create_session() -> dict:
    return {"sid": docking_session.create()}


@app.get("/api/docking/sources")
def docking_sources() -> dict:
    """Step 1 的两个候选来源清单:PPI 分析目录(展示名 = 药物 → 疾病)与 ORA 富集文件。"""
    from WY.core.docking import split_result_prefix
    from WY.core.enrichment import parse_enrichment_filename
    from WY.core.files import list_cache_files, list_download_dirs
    from WY.core.management import enrichment_display_name

    dirs = list_download_dirs()
    ppi_dirs = dirs[dirs["Directory"].str.endswith("_data", na=False)] if not dirs.empty else dirs
    ppi_out = []
    for d in ppi_dirs["Directory"].tolist():
        prefix = d[:-5] if d.endswith("_data") else d
        drug, disease = split_result_prefix(prefix)
        ppi_out.append({"name": d, "display": f"{drug} → {disease}" or d})

    enr = list_cache_files("enrichment", exts=(".csv",))
    if not enr.empty:
        is_ora = enr["File name"].map(
            lambda n: (parse_enrichment_filename(n) or {}).get("kind") == "ora")
        enr = enr[is_ora].sort_values("Modified time", ascending=False)
    enr_out = [{"name": row["File name"],
                "display": f"{enrichment_display_name(row['File name'])} · {row['Modified time']}"}
               for _, row in enr.iterrows()]
    return {"ppi_dirs": ppi_out, "enr_files": enr_out}


class DockSourceBody(BaseModel):
    mode: str            # "ppi" | "enr" | "manual"
    name: str = ""       # ppi 目录名 / enr 文件名
    manual_genes: str = ""


@app.post("/api/docking/{sid}/source")
def docking_source(sid: str, body: DockSourceBody) -> dict:
    dock = _session(sid)
    import pandas as pd

    from WY.core.docking import split_result_prefix
    from WY.core.enrichment import enrichment_filename_label
    from WY.core.files import list_cache_files
    from WY.core.ppi import read_ppi_dataset
    from WY.core.text_utils import genes_from_dataframe

    hub_genes: list[str] = []
    if body.mode == "ppi":
        dock["ppi_dir"] = body.name
        dock["hub_genes"], dock["enrichment_genes"] = read_ppi_dataset(body.name)
        prefix = body.name[:-5] if body.name.endswith("_data") else body.name
        dock["label"] = prefix
        dock["drug"], dock["disease"] = split_result_prefix(prefix)
    elif body.mode == "enr":
        dock["enr_file"] = body.name
        df = pd.read_csv(os.path.join(config.ENRICHMENT_CACHE_DIR, body.name),
                         encoding="utf-8-sig", comment="#")
        dock["enrichment_genes"] = sorted(genes_from_dataframe(df))
        dock["hub_genes"] = []
        label = enrichment_filename_label(body.name)
        dock["label"] = label
        dock["drug"], dock["disease"] = split_result_prefix(label)
    elif body.mode == "manual":
        dock["genes"] = sorted({g.strip().upper() for g in re.split(r"[,;\s]+", body.manual_genes)
                                if g.strip()})
        dock["label"] = ""
        dock["drug"] = ""
        dock["disease"] = ""
        dock["hub_genes"], dock["enrichment_genes"] = [], []
    else:
        raise HTTPException(status_code=400, detail="mode must be ppi|enr|manual")

    docking_session.reset_downstream(dock)
    if body.mode != "manual":
        hub_set = set(dock.get("hub_genes") or [])
        enr_set = set(dock.get("enrichment_genes") or [])
        dock["candidate_pre"] = sorted(hub_set & enr_set) if hub_set else sorted(enr_set)
    return {
        "label": dock.get("label", ""), "drug": dock.get("drug", ""),
        "disease": dock.get("disease", ""),
        "hub_count": len(dock.get("hub_genes") or []),
        "enr_count": len(dock.get("enrichment_genes") or []),
        "candidates": dock.get("candidate_pre") or dock.get("genes") or [],
    }


class DockStructBody(BaseModel):
    use_struct: bool = True
    allow_af: bool = False


@app.post("/api/docking/{sid}/check-structures")
def docking_check_structures(sid: str, body: DockStructBody) -> dict:
    """Step 1 结构可用性批量筛选(RCSB 4 路并发 + GraphQL 批量归因):经 job 提交。"""
    dock = _session(sid)
    candidate = dock.get("candidate_pre") or []
    if not candidate:
        raise HTTPException(status_code=400, detail="no candidates: select a source first")
    if body.use_struct:
        from WY.core.docking import fetch_receptor_candidates_batch

        def task(progress) -> dict:
            progress(fraction=0.15, text="RCSB 批量检索中…")
            results = fetch_receptor_candidates_batch(list(candidate), max_hits=20,
                                                      allow_alphafold=body.allow_af)
            gene_structures = {g: len(df) for g, df in results.items() if not df.empty}
            dock["batch_results"] = results
            dock["gene_structures"] = gene_structures
            dock["genes"] = sorted(gene_structures)
            hub_set = set(dock.get("hub_genes") or [])
            no_struct = [g for g in candidate if g not in gene_structures]
            return {"genes": dock["genes"], "no_struct": no_struct,
                    "hub_no_struct": sorted(hub_set & set(no_struct)),
                    "structures": {g: gene_structures[g] for g in dock["genes"]}}
        return {"job_id": jobs.submit(task, "docking structure check")}
    dock["genes"] = list(candidate)
    return {"job_id": None, "direct": {"genes": dock["genes"], "no_struct": [],
                                       "hub_no_struct": [], "structures": {}}}


class DockLigandBody(BaseModel):
    source: str          # "pubchem" | "smiles"
    query: str = ""
    smiles: str = ""


@app.post("/api/docking/{sid}/ligand")
def docking_ligand(sid: str, body: DockLigandBody) -> dict:
    """Step 2 配体获取:PubChem 名/CID → 3D SDF → pdbqt,或 SMILES → 3D → pdbqt。
    同步端点(线程池执行);失败以 error 字段返回(与原页面分类一致),不抛 500。"""
    dock = _session(sid)
    from datetime import datetime

    from WY.core.docking import (PubChemUnavailable, fetch_ligand_sdf_3d,
                                 prepare_ligand_pdbqt, resolve_cid, smiles_to_sdf_3d)

    if body.source == "pubchem":
        query = body.query.strip()
        if not query:
            return {"ok": False, "error": "未提供药物名或 CID"}
        try:
            cid = resolve_cid(query)
        except PubChemUnavailable as exc:
            return {"ok": False, "error": f"PubChem 暂不可用({exc}),可稍后重试"}
        if cid is None:
            return {"ok": False, "error": f"PubChem 查无化合物:{query}"}
        try:
            sdf_path = fetch_ligand_sdf_3d(cid)
        except Exception as exc:
            return {"ok": False, "error": f"配体获取失败:{exc}"}
        if not sdf_path:
            return {"ok": False, "error": f"CID {cid} 无 3D 结构,请换一个条目或用 SMILES"}
        dock["cid"] = cid
        try:
            dock["lig_pdbqt"] = prepare_ligand_pdbqt(sdf_path)
        except Exception as exc:
            return {"ok": False, "error": f"配体准备失败:{exc}"}
    elif body.source == "smiles":
        if not body.smiles.strip():
            return {"ok": False, "error": "未提供 SMILES"}
        try:
            os.makedirs(config.DOCKING_RAW_CACHE_DIR, exist_ok=True)
            sdf_path = smiles_to_sdf_3d(
                body.smiles,
                os.path.join(config.DOCKING_RAW_CACHE_DIR,
                             f"ligand_smiles_{datetime.now():%Y%m%d_%H%M%S}.sdf"))
            dock["cid"] = None
            dock["lig_pdbqt"] = prepare_ligand_pdbqt(sdf_path)
        except Exception as exc:
            return {"ok": False, "error": f"配体生成失败:{exc}"}
    else:
        raise HTTPException(status_code=400, detail="source must be pubchem|smiles")
    dock["lig_auto_done"] = True
    return {"ok": True, "cid": dock.get("cid"), "lig_pdbqt": dock.get("lig_pdbqt")}


def _candidate_rows(cand) -> list[dict]:
    """候选结构表 → JSON 行(含核酸防呆标注)。"""
    flag_notes = {"dna_only": "纯核酸结构,不含蛋白,不能对接",
                  "has_nucleic": "含核酸组分,注意确认靶点"}
    rows = []
    for _, r in cand.iterrows():
        flag = str(r.get("flag") or "")
        rows.append({"pdb_id": r["pdb_id"], "resolution": str(r.get("resolution", "")),
                     "method": str(r.get("method", "")), "title": str(r.get("title", "")),
                     "source": str(r.get("source", "")), "flag": flag,
                     "note": flag_notes.get(flag, "")})
    return rows


class DockReceptorQueryBody(BaseModel):
    gene: str
    allow_af: bool = False


@app.post("/api/docking/{sid}/receptor/candidates")
def docking_receptor_candidates(sid: str, body: DockReceptorQueryBody) -> dict:
    """Step 3:该基因的受体候选列表;优先复用 Step 1 批量结果,查不到再单基因查(带缓存)。"""
    dock = _session(sid)
    if body.gene not in (dock.get("genes") or []):
        raise HTTPException(status_code=400, detail=f"gene {body.gene} not in confirmed genes")
    try:
        cand = (dock.get("batch_results") or {}).get(body.gene)
        if cand is None:
            cand = _pdb_candidates_cached(body.gene, bool(body.allow_af and config.ALPHAFOLD_ENABLED))
    except Exception as exc:
        return {"error": f"RCSB 检索失败:{exc}", "rows": []}
    # 换基因:失效下游(与原页面一致)
    if dock.get("gene") != body.gene:
        dock["gene"] = body.gene
        for _k in ("rec_pdbqt", "pdb_path", "pdb_id", "box", "result", "result_dir",
                   "starred", "receptor_source", "pdb_sel", "pose_sel", "conclusion"):
            dock.pop(_k, None)
    return {"rows": _candidate_rows(cand) if cand is not None else [],
            "af_enabled": bool(config.ALPHAFOLD_ENABLED)}


class DockReceptorPrepareBody(BaseModel):
    gene: str
    pdb_id: str


@app.post("/api/docking/{sid}/receptor/prepare")
def docking_receptor_prepare(sid: str, body: DockReceptorPrepareBody) -> dict:
    """下载受体 PDB(或 AlphaFold 兜底)+ 转 pdbqt:网络+本地转换,经 job 提交。"""
    dock = _session(sid)
    if body.gene not in (dock.get("genes") or []):
        raise HTTPException(status_code=400, detail=f"gene {body.gene} not in confirmed genes")
    if body.pdb_id.startswith("AFDB:") and not config.ALPHAFOLD_ENABLED:
        raise HTTPException(status_code=400, detail="AlphaFold fallback is disabled")

    def task(progress) -> dict:
        from WY.core.docking import (download_pdb, fetch_alphafold_receptor,
                                     prepare_receptor_pdbqt)

        progress(fraction=0.2, text="下载受体结构…")
        if body.pdb_id.startswith("AFDB:"):
            pdb_path, af_acc = fetch_alphafold_receptor(
                body.gene, taxon=int(config.SPECIES_MAP["Human"]["ncbi_taxon"]))
            dock["pdb_id"] = f"AFDB:{af_acc}"
            dock["receptor_source"] = "alphafold"
        else:
            pdb_path = download_pdb(body.pdb_id)
            dock["pdb_id"] = body.pdb_id
            dock["receptor_source"] = "experiment"
        progress(fraction=0.7, text="受体转 pdbqt…")
        dock["pdb_path"] = pdb_path
        dock["rec_pdbqt"] = prepare_receptor_pdbqt(pdb_path)
        dock["pdb_sel"] = body.pdb_id
        dock.pop("result", None)
        return {"pdb_id": dock["pdb_id"], "receptor_source": dock["receptor_source"],
                "pdb_path": pdb_path}

    return {"job_id": jobs.submit(task, f"docking receptor {body.gene}/{body.pdb_id}")}


class DockBoxBody(BaseModel):
    infer: bool = False
    center: list[float] | None = None
    size: list[float] | None = None


@app.post("/api/docking/{sid}/box")
def docking_box(sid: str, body: DockBoxBody) -> dict:
    dock = _session(sid)
    from WY.core.docking import clamp_box, infer_box_from_cofactor, parse_hetatm

    if not dock.get("pdb_path"):
        raise HTTPException(status_code=400, detail="receptor not prepared")
    if body.infer:
        box = infer_box_from_cofactor(dock["pdb_path"])
        clamped = False
        if (dock.get("receptor_source") == "alphafold"
                and max(box["size"]) > config.DOCKING_BOX_MAX_SIZE_AF):
            box = clamp_box(box, config.DOCKING_BOX_MAX_SIZE_AF)
            clamped = True
        dock["box"] = box
    else:
        if not body.center or not body.size or len(body.center) != 3 or len(body.size) != 3:
            raise HTTPException(status_code=400, detail="center/size must be [x,y,z]")
        dock["box"] = {"center": tuple(body.center), "size": tuple(body.size)}
    het = parse_hetatm(dock["pdb_path"])
    ligs = sorted({f"{r} (chain {c}{n})" for r, c, n, _ in het}) if het else []
    from WY.core.docking import residues_in_box

    residues = residues_in_box(dock["pdb_path"], dock["box"]) if dock.get("box") else []
    return {"box": {"center": list(dock["box"]["center"]), "size": list(dock["box"]["size"])},
            "clamped": clamped if body.infer else False,
            "inferred_ligands": ligs[:5] if body.infer else None,
            "residues": residues}


@app.get("/api/docking/{sid}/box-3d")
def docking_box_3d(sid: str) -> Response:
    """Step 4 盒子 3D 预览:服务端 py3Dmol 渲染为自包含 HTML,前端 iframe 内嵌。"""
    dock = _session(sid)
    from WY.core.docking import cofactor_hetatm_pdb

    if not (dock.get("pdb_path") and dock.get("box")):
        raise HTTPException(status_code=400, detail="receptor/box not ready")
    try:
        import py3Dmol
    except ImportError:
        raise HTTPException(status_code=503, detail="py3Dmol not installed")
    pdb_path = dock["pdb_path"]
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    atom_lines = [ln for ln in text.splitlines() if ln.startswith(("ATOM", "TER", "END"))]
    het_pdb = cofactor_hetatm_pdb(pdb_path, receptor_text=text)
    c, s = dock["box"]["center"], dock["box"]["size"]
    view = py3Dmol.view(width=980, height=560)
    view.addModel("\n".join(atom_lines) + "\n", "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
    if het_pdb:
        view.addModel(het_pdb + "\n", "pdb")
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "greyCarbon"}})
    view.addBox({"center": {"x": c[0], "y": c[1], "z": c[2]},
                 "dimensions": {"w": s[0], "h": s[1], "d": s[2]},
                 "color": "orange", "opacity": 0.6})
    view.zoomTo()
    try:
        html = view._make_html()
    except AttributeError:
        html = view.write()
    return Response(content=_embed_centered(html), media_type="text/html")


class DockRunBody(BaseModel):
    pass


@app.post("/api/docking/{sid}/run")
def docking_run(sid: str, body: DockRunBody) -> dict:
    dock = _session(sid)
    from datetime import datetime

    if not (dock.get("rec_pdbqt") and dock.get("lig_pdbqt") and dock.get("box")):
        raise HTTPException(status_code=400, detail="ligand/receptor/box not ready")

    def task(progress) -> dict:
        from WY.core.docking import parse_vina_scores, run_vina, split_vina

        progress(fraction=None, text=f"vina 运行中(exhaustiveness={config.EXHAUSTIVENESS},"
                                     f" 最长 {config.DOCKING_TIMEOUT}s)…")
        out_pdbqt = os.path.join(config.DOCKING_RAW_CACHE_DIR,
                                 f"vina_out_{datetime.now():%Y%m%d_%H%M%S}.pdbqt")
        proc = run_vina(dock["rec_pdbqt"], dock["lig_pdbqt"], dock["box"], out_path=out_pdbqt)
        progress(fraction=0.9, text="解析构象与打分…")
        poses = split_vina(out_pdbqt, out_dir=os.path.dirname(out_pdbqt))
        scores = parse_vina_scores(out_pdbqt)
        # 新一轮对接:清掉上一轮残留(与原页面一致)
        for _k in ("pose_sel", "conclusion", "result_dir", "starred"):
            dock.pop(_k, None)
        dock["result"] = {"log_tail": (proc.stdout or "")[-1500:], "out": out_pdbqt,
                          "poses": poses, "scores_rows": scores.to_dict(orient="records")}
        return {"log_tail": dock["result"]["log_tail"],
                "scores": dock["result"]["scores_rows"], "poses": poses}

    return {"job_id": jobs.submit(task, "docking vina run")}


class DockSaveBody(BaseModel):
    pose_mode: int
    conclusion: str = ""
    drug: str = ""
    disease: str = ""


@app.post("/api/docking/{sid}/save")
def docking_save(sid: str, body: DockSaveBody) -> dict:
    dock = _session(sid)
    import pandas as pd

    from WY.core.docking import docking_result_prefix, save_docking_result

    if not dock.get("result"):
        raise HTTPException(status_code=400, detail="no vina result to save")
    drug = body.drug.strip() or dock.get("drug", "")
    disease = body.disease.strip() or dock.get("disease", "")
    gene = dock.get("gene", "gene")
    prefix = docking_result_prefix(drug, disease, gene)
    scores = pd.DataFrame(dock["result"]["scores_rows"])
    result_dir = save_docking_result(
        prefix, gene, dock.get("pdb_id", ""), dock.get("cid"), dock["box"],
        dock["result"]["out"], scores, body.conclusion,
        receptor_source=dock.get("receptor_source", "experiment"),
        pose_mode=int(body.pose_mode))
    dock["result_dir"] = result_dir
    dock["drug_name"], dock["disease_name"] = drug, disease
    dock.pop("starred", None)
    return {"result_dir": result_dir}


class DockStarBody(BaseModel):
    on: bool


@app.post("/api/docking/{sid}/star")
def docking_star(sid: str, body: DockStarBody) -> dict:
    dock = _session(sid)
    import glob
    import shutil as _shutil

    if not dock.get("result_dir"):
        raise HTTPException(status_code=400, detail="save the result first")
    if body.on:
        from WY.core.docking import star_docking_result

        star_docking_result(dock["result_dir"])
        dock["starred"] = True
        return {"starred": True}
    star_base = os.path.basename(dock["result_dir"])
    for p in glob.glob(os.path.join(config.MD_STAR_DIR, star_base) + "*"):
        if os.path.isdir(p):
            _shutil.rmtree(p, ignore_errors=True)
    dock.pop("starred", None)
    return {"starred": False}


@app.get("/api/docking/{sid}/view3d")
def docking_view3d(sid: str, pose: int, cutoff: float = 4.0, focus: int = 0) -> Response:
    """Step 6 3D 预览:受体卡通 + 构象绿棍 + 口袋残基橙棍标注 + 相互作用虚线 +
    共晶配体参照;服务端 py3Dmol 渲染,图例由前端展示。"""
    dock = _session(sid)
    if not dock.get("result") or not dock.get("pdb_path"):
        raise HTTPException(status_code=400, detail="no docking result")
    try:
        import py3Dmol
    except ImportError:
        raise HTTPException(status_code=503, detail="py3Dmol not installed")
    from WY.core.docking import (cofactor_hetatm_pdb, extract_model,
                                 interactions_near_ligand, pose_pdbqt_to_pdb,
                                 residues_near_ligand)

    out_stem = os.path.splitext(os.path.basename(dock["result"]["out"]))[0]
    pose_path = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{out_stem}_pose_{int(pose)}.pdbqt")
    extract_model(dock["result"]["out"], int(pose), pose_path)
    pose_pdb = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{out_stem}_pose_view.pdb")
    try:
        pose_pdbqt_to_pdb(pose_path, pose_pdb)
    except Exception:
        raise HTTPException(status_code=500, detail="pose pdbqt→pdb conversion failed")
    with open(dock["pdb_path"], encoding="utf-8", errors="replace") as f:
        receptor_text = f.read()
    with open(pose_pdb, encoding="utf-8") as f:
        ligand_pdb = f.read()
    het_pdb = cofactor_hetatm_pdb(dock["pdb_path"], receptor_text=receptor_text)
    pocket = residues_near_ligand(dock["pdb_path"], pose_path, cutoff=float(cutoff))
    interactions = interactions_near_ligand(dock["pdb_path"], pose_path)

    view = py3Dmol.view(width=1100, height=560)
    view.addModel(receptor_text, "pdb")
    view.addModel(ligand_pdb, "pdb")
    cartoon: dict = {"color": "spectrum"}
    if focus:
        cartoon["opacity"] = 0.3
    view.setStyle({"model": 0}, {"cartoon": cartoon})
    view.setStyle({"model": 1}, {"stick": {"colorscheme": "greenCarbon"}})
    if het_pdb:
        view.addModel(het_pdb + "\n", "pdb")
        view.setStyle({"model": 2}, {"stick": {"colorscheme": "greyCarbon",
                                               "opacity": 0.5, "radius": 0.15}})
    by_chain: dict[str, list] = {}
    for r in pocket or []:
        by_chain.setdefault(r["chain"], []).append(r["resi"])
    for chain, resis in by_chain.items():
        sel: dict = {"model": 0, "resi": resis}
        if chain:
            sel["chain"] = chain
        view.setStyle(sel, {"cartoon": cartoon,
                            "stick": {"colorscheme": "orangeCarbon", "radius": 0.2}})
    for r in pocket or []:
        x, y, z = r["pos"]
        view.addLabel(f"{r['resn']}-{r['resi']}", {
            "position": {"x": x, "y": y, "z": z},
            "fontSize": 11, "fontColor": "black",
            "showBackground": True, "backgroundColor": "white",
            "backgroundOpacity": 0.7,
        })
    for it in interactions or []:
        view.addCylinder({"dashed": True, "radius": 0.06, "dashLength": 0.3, "gapLength": 0.22,
                          "color": {"hbond": "#0f62fe", "salt": "#c026d3",
                                    "hydro": "#808080"}.get(it["kind"], "#777777"),
                          "start": {"x": it["pos1"][0], "y": it["pos1"][1], "z": it["pos1"][2]},
                          "end": {"x": it["pos2"][0], "y": it["pos2"][1], "z": it["pos2"][2]}})
        if it["kind"] != "hydro":
            mid = [(a + b) / 2 for a, b in zip(it["pos1"], it["pos2"])]
            view.addLabel(f'{it["dist"]:.1f} Å', {
                "position": {"x": mid[0], "y": mid[1], "z": mid[2]},
                "fontSize": 10, "fontColor": "#333333",
                "showBackground": True, "backgroundColor": "white",
                "backgroundOpacity": 0.6,
            })
    if focus:
        view.zoomTo({"model": 1})
        view.zoom(0.55)
    else:
        view.zoomTo()
    try:
        html = view._make_html()
    except AttributeError:
        html = view.write()
    return Response(content=_embed_centered(html), media_type="text/html")


@app.get("/api/docking/{sid}/interactions")
def docking_interactions(sid: str, pose: int) -> dict:
    dock = _session(sid)
    from WY.core.docking import extract_model, interactions_near_ligand

    if not dock.get("result") or not dock.get("pdb_path"):
        raise HTTPException(status_code=400, detail="no docking result")
    out_stem = os.path.splitext(os.path.basename(dock["result"]["out"]))[0]
    pose_path = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{out_stem}_pose_{int(pose)}.pdbqt")
    extract_model(dock["result"]["out"], int(pose), pose_path)
    return {"rows": interactions_near_ligand(dock["pdb_path"], pose_path)}


class DockSave3dBody(BaseModel):
    pose: int
    caption: str = ""


@app.post("/api/docking/{sid}/save3d")
def docking_save3d(sid: str, body: DockSave3dBody) -> dict:
    """保存交互式 3D 视图 HTML 到结果目录(caption 由前端传入,离线可用)。"""
    dock = _session(sid)
    from WY.core.docking import save_interactive_3d_html

    if not (dock.get("result") and dock.get("pdb_path") and dock.get("result_dir")):
        raise HTTPException(status_code=400, detail="save the result first")
    out_stem = os.path.splitext(os.path.basename(dock["result"]["out"]))[0]
    pose_path = os.path.join(config.DOCKING_RAW_CACHE_DIR,
                             f"{out_stem}_pose_{int(body.pose)}.pdbqt")
    if not os.path.exists(pose_path):
        from WY.core.docking import extract_model

        extract_model(dock["result"]["out"], int(body.pose), pose_path)
    file_name = f"3d_view_{safe_name(dock.get('gene', 'target'))}_pose{int(body.pose)}.html"
    out_path, offline_ok = save_interactive_3d_html(
        dock["pdb_path"], pose_path,
        os.path.join(dock["result_dir"], file_name), body.caption)
    rel = os.path.relpath(out_path, config.DOWNLOADS_DIR).replace(os.sep, "/")
    return {"path": out_path, "artifact_rel": rel, "offline_ok": offline_ok}


@app.get("/api/docking/{sid}/zip")
def docking_zip(sid: str) -> Response:
    dock = _session(sid)
    from WY.core.files import make_dir_zip

    if not dock.get("result_dir"):
        raise HTTPException(status_code=400, detail="save the result first")
    data = make_dir_zip(dock["result_dir"])
    return Response(content=data, media_type="application/zip",
                    headers={"Content-Disposition":
                             f'attachment; filename="{os.path.basename(dock["result_dir"])}.zip"'})


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
        # index.html 不缓存:构建产物文件名带哈希,浏览器缓存旧 index 会引用已删除的旧资产
        return FileResponse(os.path.join(_DIST, "index.html"),
                            headers={"Cache-Control": "no-cache"})
else:
    @app.get("/", include_in_schema=False)
    async def no_frontend() -> dict:
        return {"hint": "前端未构建:cd WY/web && npm install && npm run build,再重启 npap serve"}
