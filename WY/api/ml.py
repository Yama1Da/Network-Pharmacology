"""GEO 表达谱 + ML 共识筛选(阶段 3c:Expression_ML_Screening + 结果管理)。"""

from fastapi import APIRouter
import logging
import os
import time
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel
from WY import config
from WY.core import jobs, transfer
from WY.core.ml_pipeline import STAGES, run_screen_chain
from WY.core.ppi import (base_disease_of, list_analyzed_diseases,
                        list_analyzed_drugs,
                        list_ranked_csv_candidates)
from .common import (_df_rows)


router = APIRouter()


# 全链编排唯一定义在 core.ml_pipeline.run_screen_chain(CLI/页面/API 共用);
# 本层只做 HTTP 形态:检索/预检(同步)、总览与运行(job)、结果管理与下载。


@router.get("/api/ml/diseases")
def ml_diseases() -> dict:
    from WY.core.ppi import list_analyzed_diseases, list_analyzed_drugs

    return {"diseases": list_analyzed_diseases(), "drugs": list_analyzed_drugs()}


class MlSearchBody(BaseModel):
    disease: str
    limit: int = 15


@router.post("/api/ml/search")
def ml_search(body: MlSearchBody) -> dict:
    """GEO 数据集检索(同步:NCBI eutils,数秒)。"""
    from WY.core.geo import search_geo_datasets

    if not body.disease.strip():
        raise HTTPException(status_code=400, detail="disease is required")
    try:
        df = search_geo_datasets(body.disease.strip(), limit=max(5, min(50, int(body.limit))))
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"search failed: {exc}")
    return {"results": df.to_dict(orient="records") if not df.empty else [],
            "disease": body.disease.strip()}


def _ml_progress(cb):
    """ml_pipeline progress(stage, message) → job 进度(阶段号 + 附加字段)。"""
    stages = STAGES

    def _p(stage, message=""):
        if stage in stages:
            idx = stages.index(stage)
            cb(fraction=idx / len(stages),
               text=f"阶段 {idx + 1}/{len(stages)}:{STAGE_ZH.get(stage, stage)}",
               step=stage, state="running")
        else:
            cb(text=str(message)[:120])
    return _p


STAGE_ZH = {"fetch": "获取数据集", "group": "样本分组", "preprocess": "预处理",
            "deg": "差异表达", "ml": "ML 共识", "shortlist": "shortlist", "save": "落盘"}


class MlProbeBody(BaseModel):
    gse: str
    platform: str = ""   # 检索结果带出的平台号(芯片型才探注释)
    type_: str = ""      # 检索结果带出的类型("芯片"/"测序")


def _empty_matrix_assess(gse: str) -> dict:
    """空矩阵类数据集的可用性判定(series matrix 不存在或为元数据空表):
    counts 补充文件决定生死 —— 有则运行时自动回退 DESeq2,可用;无则确定性不可用。"""
    from WY.core.geo import has_counts_supplement

    counts = has_counts_supplement(gse)
    return {"matrix_empty": True, "counts_available": counts,
            "symbol_level": None, "annot": None, "gpl": "",
            "group_fields": [], "n_samples": None, "usable": counts}


def _probe_one(gse: str, platform: str = "", type_: str = "") -> dict:
    """单 GSE 预检判定(ml_probe 与 probe-batch 共用):
    体积/速度 + 可用性(已缓存→权威;未缓存→元数据条件版)+ counts 补充 + 注释探活。

    批检附加:counts_confirmed —— counts 可用性是否经本地缓存实锤
    (True=缓存盘命中确凿;False=已确认不可用;None=文件名级判定,未下载无法确认)。
    """
    from WY.core.expression import (assess_expression_dataset, normalize_gpl,
                                    probe_platform_annot)
    from WY.core.geo import _geo_cache_path, fetch_geo_dataset, \
        probe_series_matrix

    out: dict = {"gse": gse, "r": None, "err": "", "assess": None,
                 "annot_missing": False, "gpl": "", "counts_confirmed": None}
    try:
        out["r"] = probe_series_matrix(gse)
    except Exception as exc:
        out["err"] = str(exc)
        # 404 = series matrix 不存在:元数据型提交(RNA-seq 常态)或提交者未提供。
        # 这是确定性事实而非探测失败 —— 与"小体积空表"同路径体检 counts,
        # 给出绿/红结论,不留灰态(2026-11-07 GSE250283 实测:404+注释缺失
        # 却亮灰"待运行时确认",用户被迫跑完总览才知道不可用)
        if "404" in str(exc) or "Not Found" in str(exc):
            out["assess"] = _empty_matrix_assess(gse)
    if type_.strip() == "芯片" and platform.strip():
        try:
            out["annot"] = probe_platform_annot(platform.strip())
            out["gpl"] = normalize_gpl(platform.strip())
            out["annot_missing"] = out["annot"] == "missing"
        except Exception:
            pass
    if os.path.isfile(_geo_cache_path(gse)):
        try:
            out["assess"] = assess_expression_dataset(fetch_geo_dataset(gse))
            if out["assess"].get("matrix_empty"):
                # 空矩阵分支的 counts_available 由 has_counts_supplement 文件名级判定;
                # 缓存盘上若已有 counts 矩阵文件则为实锤可用(运行时不会重下)
                out["counts_confirmed"] = bool(out["assess"].get("counts_available"))
        except Exception:
            pass
    elif out["r"] and out["r"].get("size_bytes", 10 ** 9) < 20 * 1024:
        # 小体积 series matrix(≈3KB)几乎必然是元数据型空表:预检 counts 可用性
        out["assess"] = _empty_matrix_assess(gse)
        # counts_confirmed 保持 None:补充文件未落盘,文件名级判定待运行时实锤
    return out


def _probe_status(p: dict) -> str:
    """把 _probe_one 的原料归纳为四态徽章:err/unusable/ok/unknown。

    判定优先级:assess 的可用/不可用结论 > 注释缺失(探针级矩阵的确定性
    不可用)> 探测本身失败(网络层,真未知)> 其余。"""
    a = p.get("assess") or {}
    if a:
        if a.get("usable") or a.get("counts_available"):
            return "ok"
        if a.get("matrix_empty") and not a.get("counts_available"):
            return "unusable"
    if p.get("annot_missing"):
        return "unusable"
    if p.get("err"):
        return "err"
    return "unknown"


@router.post("/api/ml/probe")
def ml_probe(body: MlProbeBody) -> dict:
    """GSE 预检(同步,小 Range 请求,每 GSE 前端调用一次):
    体积/速度/ETA + 可用性判定(已缓存→权威;未缓存→元数据条件版)+ counts 补充 + 注释探活。"""
    gse = body.gse.strip().upper()
    if not gse:
        raise HTTPException(status_code=400, detail="gse is required")
    return _probe_one(gse, body.platform, body.type_)


class MlProbeBatchBody(BaseModel):
    class Row(BaseModel):
        gse: str
        platform: str = ""
        type_: str = ""

    rows: list[Row]


@router.post("/api/ml/probe-batch")
def ml_probe_batch(body: MlProbeBatchBody) -> dict:
    """检索表批量体检(同步,行间 0.3s 礼貌间隔):
    每行跑与单检同款的判定,返回轻量徽章原料 {gse, err, assess, counts_confirmed, status}。
    上限 50 行 = 检索条上限(limit stepper max),50 行最坏 1~3 分钟;
    不体积/测速 —— 徽章只需要判定,不需要 ETA。"""
    if not body.rows:
        raise HTTPException(status_code=400, detail="rows is required")
    if len(body.rows) > 50:
        raise HTTPException(status_code=400, detail="rows limited to 50")
    results = []
    for i, row in enumerate(body.rows):
        gse = row.gse.strip().upper()
        if not gse:
            continue
        p = _probe_one(gse, row.platform, row.type_)
        results.append({"gse": gse, "err": p["err"], "assess": p["assess"],
                        "counts_confirmed": p.get("counts_confirmed"),
                        "status": _probe_status(p)})
        if i < len(body.rows) - 1:
            time.sleep(0.3)   # 行间礼貌间隔:探测走裸 requests,不经统一节流
    return {"results": results}


class MlOverviewBody(BaseModel):
    gse: str


@router.post("/api/ml/overview")
def ml_overview(body: MlOverviewBody) -> dict:
    """分组总览(经 job:可能触发数据集下载,数秒到数分钟):
    分组字段分布 + counts 列清单(空矩阵类)+ 可用性权威判定。"""
    import pandas as pd

    from WY.core.expression import (assess_expression_dataset, counts_picker_fields,
                                    group_picker_fields)
    from WY.core.geo import fetch_gse_count_matrix, fetch_geo_dataset, group_fields_overview

    gse = body.gse.strip().upper()
    if not gse:
        raise HTTPException(status_code=400, detail="gse is required")

    def task(progress) -> dict:
        progress(fraction=0.1, text="下载/读取数据集…")
        ds = fetch_geo_dataset(gse)
        progress(fraction=0.75, text="分组字段总览…")
        overview = group_fields_overview(ds)
        a = assess_expression_dataset(ds)
        picker = group_picker_fields(ds)
        counts_picker = None
        counts_match = None
        counts_rows: list[dict] = []
        if ds.issues and ds.table.empty:
            try:
                cnt = fetch_gse_count_matrix(gse)
                count_columns = [str(c) for c in cnt["counts"].columns]
                counts_picker = counts_picker_fields(ds, count_columns)
                grp_cols = [c for c in ds.samples.columns
                            if c.startswith("characteristics")
                            and 2 <= ds.samples[c].astype(str).nunique() <= 8
                            and not any(k in c.lower() for k in
                                        ("age", "sex", "gender", "date", "batch"))]
                title = (ds.samples["title"].astype(str)
                         if "title" in ds.samples.columns else None)
                for c in count_columns:
                    row = {"sample": c}
                    if title is not None:
                        m = ds.samples[title == c]
                        if not m.empty:
                            for gc in grp_cols:
                                row[gc] = str(m.iloc[0][gc])[:64]
                    counts_rows.append(row)
            except Exception as exc:
                logging.exception("counts overview failed for %s", gse)
                counts_rows = []
                counts_picker = None
                return {"overview": overview.to_dict(orient="records"),
                        "counts_rows": counts_rows, "assess": a,
                        "group_picker": {"fields": picker}, "counts_picker": None,
                        "counts_match": None,
                        "issues": list(ds.issues), "counts_error": str(exc)}
            # 列名↔样本标题匹配率:提交者常在 counts 文件里用另一套缩写
            # (实测 GSE353436:列名 SH1/NCM1 vs 标题 "AGS gastric cancer cells, NC, …"),
            # 全部对不上时分组信息不可自动提取,前端据此明示而非静默无选择器
            counts_match = {"matched": sum(1 for row in counts_rows if len(row) > 1),
                            "total": len(count_columns)}
        return {"overview": overview.to_dict(orient="records"),
                "counts_rows": counts_rows, "assess": a, "group_picker": {"fields": picker},
                "counts_picker": {"fields": counts_picker}, "counts_match": counts_match,
                "issues": list(ds.issues)}

    return {"job_id": jobs.submit(task, f"ml overview {gse}")}


class MlRankedBody(BaseModel):
    disease: str = ""
    drug: str = ""


@router.post("/api/ml/ranked")
def ml_ranked(body: MlRankedBody) -> dict:
    """平台排序 CSV 候选(按疾病上下文与药物过滤,规则与原页面一致)。"""
    from WY.core.ppi import base_disease_of, list_ranked_csv_candidates

    ranked = list_ranked_csv_candidates()
    ctx = base_disease_of(body.disease).lower() if body.disease.strip() else ""
    f_drug = body.drug.strip().lower()
    shown = [rf for rf in ranked
             if (not ctx or ctx == rf["base_disease"].lower()
                 or ctx in rf["base_disease"].lower())
             and (not f_drug or f_drug in rf["drug"].lower())]
    return {"files": [{"label": rf["label"], "base_disease": rf["base_disease"],
                       "drug": rf["drug"]} for rf in shown]}


class MlRunBody(BaseModel):
    mode: str                     # "gse" | "local"
    gse: str = ""
    case_pattern: str = ""
    control_pattern: str = ""
    expr: str = ""
    labels: str = ""
    case_label: str = ""
    ranked_labels: list[str] = []
    top_k: int | None = None
    prefilter_top: int | None = None
    bootstrap_rounds: int | None = None
    map_ids: bool = False         # --expr 模式 Ensembl/探针级矩阵经 mygene 映射(需联网)
    perm_rounds: int = 0          # 标签置换检验轮次(0=关;离线对照,长任务,UI 暂不暴露)
    cv_repeats: int | None = None  # 分层 CV 重复次数(None → config.ML_CV_REPEATS)


def _rel_downloads(path: str) -> str | None:
    try:
        real = os.path.realpath(path)
        base = os.path.realpath(config.DOWNLOADS_DIR)
        if real.startswith(base + os.sep):
            return os.path.relpath(real, base).replace(os.sep, "/")
    except OSError:
        pass
    return None


@router.post("/api/ml/run")
def ml_run(body: MlRunBody) -> dict:
    """ML 筛选全链(经 job:DEG + ML 共识,本机实测数分钟)。"""
    if body.mode not in ("gse", "local"):
        raise HTTPException(status_code=400, detail="mode must be gse|local")
    if body.mode == "gse" and not (body.gse.strip() and body.case_pattern.strip()
                                   and body.control_pattern.strip()):
        raise HTTPException(status_code=400, detail="gse/case/control required in gse mode")
    if body.mode == "local" and not (body.expr.strip() and body.labels.strip()):
        raise HTTPException(status_code=400, detail="expr/labels required in local mode")
    # ranked 标签 → 路径(白名单:只接受 list_ranked_csv_candidates 给出的标签)
    candidates = {rf["label"]: rf["path"] for rf in list_ranked_csv_candidates()}
    ranked_csvs = {}
    for lbl in body.ranked_labels:
        path = candidates.get(lbl)
        if not path:
            raise HTTPException(status_code=400, detail=f"unknown ranked label: {lbl}")
        ranked_csvs[os.path.splitext(os.path.basename(path))[0]] = path

    def task(progress) -> dict:
        common = dict(prefilter_top=int(body.prefilter_top or 500),
                      bootstrap_rounds=int(body.bootstrap_rounds or 30),
                      ranked_csvs=ranked_csvs, top_k=int(body.top_k or 50),
                      out_base=config.GEO_RESULT_DIR, out_subdir=None,
                      progress=_ml_progress(progress),
                      map_ids=bool(body.map_ids),
                      perm_rounds=max(0, int(body.perm_rounds or 0)),
                      cv_repeats=body.cv_repeats)
        if body.mode == "local":
            res = run_screen_chain(expr=body.expr.strip(), labels=body.labels.strip(),
                                   case_label=body.case_label.strip() or None, **common)
        else:
            res = run_screen_chain(gse=body.gse.strip().upper(),
                                   case_pattern=body.case_pattern.strip(),
                                   control_pattern=body.control_pattern.strip(), **common)
        out_dir = res["out_dir"]
        artifacts = {}
        if os.path.isdir(out_dir):
            for f in sorted(os.listdir(out_dir)):
                fp = os.path.join(out_dir, f)
                if os.path.isfile(fp):
                    artifacts[f] = _rel_downloads(fp)
        roc_png_rel = _rel_downloads(res.get("roc_png") or "")
        curve_rel = _rel_downloads(os.path.join(out_dir, "consensus_curve.png"))
        return {
            "gse": res["gse"] or res["tag"], "tag": res["tag"],
            "out_dir": out_dir, "out_dir_rel": _rel_downloads(out_dir),
            "artifacts": artifacts,
            "roc_png_rel": roc_png_rel,
            "consensus_curve_rel": curve_rel,
            "elapsed": res["summary"]["elapsed_total_s"],
            "meta": res["meta"],
            "geo_fallback": bool(res["summary"].get("geo_cache_fallback")),
            "geo_cache_date": res["summary"].get("geo_cache_date", ""),
            "funnel": _df_rows(res["funnel"]),
            "consensus": _df_rows(res["consensus"].rename_axis("gene").reset_index(), 300),
            "consensus_total": int(len(res["consensus"])),
            "shortlist": _df_rows(res["shortlist"]),
            "shortlist_count": int(res["shortlist"]["in_shortlist"].sum())
            if "in_shortlist" in res["shortlist"] else None,
            "deg": _df_rows(res["deg"], 200),
        }

    return {"job_id": jobs.submit(task, f"ml run {body.gse or 'local'}")}


@router.get("/api/ml/history")
def ml_history() -> dict:
    """最近 3 次 run 的摘要(原页面历史区)。"""
    from WY.core.ml_results import list_ml_runs

    runs = []
    for r in list_ml_runs()[:3]:
        s = r["summary"] or {}
        m = s.get("ml_meta", {})
        runs.append({"name": r["name"], "auc": m.get("cv_auc"),
                     "n_samples": s.get("n_samples"),
                     "shortlist_count": s.get("shortlist_count"),
                     "ml_version": s.get("ml_version"),
                     "stale": [{"type": x["type"],
                                "old": x.get("old"), "new": x.get("new"),
                                "name": x.get("name")} for x in r["stale"]]})
    return {"runs": runs}


def _ml_run_rows() -> list[dict]:
    from WY.core.ml_results import list_ml_runs

    rows = []
    for r in list_ml_runs():
        s = r["summary"] or {}
        m = s.get("ml_meta", {})
        rows.append({"name": r["name"], "path": r["path"],
                     "dir_rel": _rel_downloads(r["path"]),
                     "size_bytes": r["size_bytes"], "mtime": r["mtime"],
                     "gse": r["gse"], "source": r["source"],
                     "n_samples": r["n_samples"], "auc": r["auc"],
                     "shortlist_count": r["shortlist_count"],
                     "ml_version": r["ml_version"],
                     "stale": [{"type": x["type"], "old": x.get("old"),
                                "new": x.get("new"), "name": x.get("name")}
                               for x in r["stale"]]})
    return rows


@router.get("/api/ml/runs")
def ml_runs() -> dict:
    """全部 ML run 清单(结果管理):新→旧,含陈旧原因与体积。"""
    return {"runs": _ml_run_rows()}


class MlRunsBody(BaseModel):
    paths: list[str]


@router.post("/api/ml/runs/delete")
def ml_runs_delete(body: MlRunsBody) -> dict:
    from WY.core.ml_results import delete_ml_runs

    return {"deleted": delete_ml_runs(body.paths)}


@router.post("/api/ml/runs/zip")
def ml_runs_zip(body: MlRunsBody) -> Response:
    from WY.core.ml_results import zip_ml_runs

    try:
        data = zip_ml_runs(body.paths)
    except ValueError as exc:
        # 登记根外/非 run 目录一律视为请求错误,400 带被拒路径明细(P1-B;
        # 与 delete 的"跳过+warning"语义不同:zip 的调用方是前端勾选列表)
        raise HTTPException(status_code=400, detail=str(exc))
    return Response(content=data, media_type="application/zip",
                    headers={"Content-Disposition": "attachment; filename=\"ML_runs.zip\""})


class MlCleanupBody(BaseModel):
    keep: int = 6


@router.post("/api/ml/runs/cleanup")
def ml_runs_cleanup(body: MlCleanupBody) -> dict:
    from WY.core.ml_results import cleanup_keep_per_gse

    deleted = cleanup_keep_per_gse(keep=max(1, int(body.keep)))
    return {"deleted": deleted}


class MlForceBody(BaseModel):
    gse: str


@router.post("/api/ml/force-refresh")
def ml_force_refresh(body: MlForceBody) -> dict:
    """GEO 缓存强制重下(降级场景②,经 job)。"""
    from WY.core.geo import fetch_geo_dataset

    gse = body.gse.strip().upper()
    if not gse:
        raise HTTPException(status_code=400, detail="gse is required")

    def task(progress) -> dict:
        progress(fraction=None, text=f"强制重下 {gse}…")
        fetch_geo_dataset(gse, force_refresh=True)
        return {"ok": True}

    return {"job_id": jobs.submit(task, f"geo force refresh {gse}")}
