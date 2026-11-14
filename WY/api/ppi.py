"""PPI 分析(阶段 3a:Data_Network 页)。"""

from fastapi import APIRouter
import os
from fastapi import HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel
from WY import config
from WY.core import jobs
from WY.core.analysis import run_ppi_analysis
from .common import (_df_rows)


router = APIRouter()


# 产物直出:分析的全部图表/表格由 run_ppi_analysis 经 save_ppi_results 落盘到
# downloads/<prefix>_data/,本层只做序列化与文件伺服,不重绘任何图。
# 交互区重跑(模块/GSEA)是无状态的:从落盘的 nodes/edges CSV 重建图,与网页展示一致。


class _JobProgress:
    """把 core.analysis.AnalysisProgress 协议适配到 core.jobs 的进度回调:
    阶段名/状态进 extra 字段,前端据此渲染阶段图标行。
    富集/GSEA/邻近度+RWR/模块检测四个任务在线程池并发执行,这里维护完整的
    step→state 映射并整体上报(steps),前端可同时点亮多个运行中的阶段。"""

    def __init__(self, cb):
        self._cb = cb
        self._states: dict[str, str] = {}

    def set(self, key: str, state: str) -> None:
        self._states[key] = state
        self._cb(text=f"{key}:{state}", steps=dict(self._states))

    def advance(self, fraction: float) -> None:
        self._cb(fraction=fraction)




def _modules_payload(modules: dict | None) -> dict | None:
    """模块检测结果 → JSON 安全载荷:汇总表 + 模块富集长表(截断前 100 行)+ 单基因隐藏数。
    long_df 在 core 侧已按 module_id, p_value 排序;完整结果由 modules_csv /
    module_enrichment_csv 产物文件承载,前端下载区可直接取用。"""
    if modules is None:
        return None
    long_df = modules.get("long_df")
    has_enr = long_df is not None and len(long_df) > 0
    return {
        "summary": _df_rows(modules.get("summary_df")),
        "failures": modules.get("failures", []),
        "no_sig": modules.get("no_sig", []),
        "hidden_single": sum(1 for m in modules.get("modules", []) if m.get("size", 0) < 2),
        "enr_rows": _df_rows(long_df, 100) if has_enr else [],
        "enr_total": int(len(long_df)) if has_enr else 0,
        "enr_modules": int(long_df["module_id"].nunique()) if has_enr else 0,
    }


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
    modules_out = _modules_payload(modules)

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


@router.post("/api/ppi/analyze")
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


@router.get("/api/artifacts")
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


@router.post("/api/ppi/rerun-modules")
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
        return _modules_payload(computed)

    return {"job_id": jobs.submit(task, "ppi rerun modules")}


class RerunGseaBody(BaseModel):
    side: str            # "drug" | "disease"
    drug: str
    disease: str
    prefix: str


@router.post("/api/ppi/rerun-gsea")
def ppi_rerun_gsea(body: RerunGseaBody) -> dict:
    if body.side not in ("drug", "disease"):
        raise HTTPException(status_code=400, detail="side must be drug|disease")

    def task(progress) -> dict:
        from WY.core.analysis import gsea_prerank_side

        progress(fraction=0.1, text=f"gsea prerank ({body.side})")
        gdf, rstats, reason = gsea_prerank_side(body.side, body.drug, body.disease, body.prefix)
        return {"reason": reason, "stats": rstats, "rows": _df_rows(gdf, 100)}

    return {"job_id": jobs.submit(task, "ppi rerun gsea")}
