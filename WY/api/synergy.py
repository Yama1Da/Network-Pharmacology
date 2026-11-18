"""分离度协同预测 API:数据流 = 药/疾病名 → 四库合并靶点(缓存) → top-N 证据切片 →
core.analyze_pair(jobs 提交-轮询,单对分钟级,§1.3 定案 3)。

错误约定(前端按稳定令牌映射 web.synergy.* 文案,核心与 API 均不产出人言文案,
避开后端单语言;清单 §1.3 致队友第 5 条同款纪律):
  - "SynergyPrecheckError: too_few_targets_a/b"   入 LCC 靶点 <5(E3 拒算)
  - "InteractomeNotBuiltError"                    背景互作组未构建
  - "SynergyNoTargetsError: synergy:no_targets_a/b/disease"  四库无数据/名称不匹配
进度:job.progress.stage ∈ targets/separation/proximity(fraction 由 core 零分布
循环按轮推进),前端按 stage 映射文案。

页级参数(top_n/n_perm)随请求传入,Settings 的 SYNERGY_* 仅作默认值(页内覆盖,
2026-11-18 定案);切片在本模块任务体内完成,契约 meta.run 记账(v1.1)。

归档:每次成功分析把契约 JSON 原样写入
<DOWNLOADS_DIR>/synergy_<A>_x_<B>_<时间戳>/result.json——自动出现在"分析产物
目录"页(可打包/删除),前端据 job 结果里的 directory_rel 提供直接下载;结果自带
meta.run 与互作组版本,同参数重跑可复现(核心确定性已单测验证)。
"""

import os
import threading
import time

import pandas as pd
from fastapi import APIRouter, HTTPException, Response
from pydantic import BaseModel

from WY import config
from WY.core import jobs
from WY.core.caching import save_json_file
from WY.core.synergy import (
    analyze_pair,
    delete_runs,
    list_results,
    read_run,
    zip_runs,
)
from WY.core.targets import (
    get_disease_targets_merged,
    get_drug_targets_merged,
    list_known_targets,
)
from WY.core.text_utils import safe_name

router = APIRouter()

# 证据切片排序规则(冒烟启发版):两源分数各自归一到 [0,1] 后取最大,降序取前 N。
# 正式口径归 §1.6 DrugComb 基准扫描定案;规则 id 进契约 meta.run.sort 记账。
SORT_RULE = "dual_source_max_norm"

# 最近一次分析任务的登记(页面重进恢复进度展示用,2026-11-18):与 pk-mirror 的
# 单飞不同,协同允许多对并行(互不落盘冲突),这里只记最新供前端找回轮询
_analyze_job_lock = threading.Lock()
_analyze_job_id: str | None = None


class SynergyNoTargetsError(RuntimeError):
    """四库均无该药/疾病的靶点(名称不匹配或库无收录),msg 内带稳定令牌。"""


class AnalyzeBody(BaseModel):
    drug_a: str
    drug_b: str
    disease: str = ""            # 空 = 无疾病锚(proximity 三字段 null,几何照算)
    top_n: int | None = 100      # None=全量;切片发生在任务体内(定案 #7,核心不截断)
    n_perm: int | None = None    # None=配置默认;一杆驱动 s_AB 零分布与三次 proximity


def _top_n_slice(df: pd.DataFrame, n: int | None) -> list[str]:
    """dual_source_max_norm 切片(score_pchembl/score_dgidb 各自归一取最大,降序前 N);
    n=None 原样全量;分数全缺的行垫 0 殿后。"""
    if df is None or df.empty:
        return []
    if n is None:
        return sorted({str(g) for g in df["approvedSymbol"].dropna()})

    def _norm(s: pd.Series) -> pd.Series:
        m = s.max()
        return (s / m).fillna(0.0) if pd.notna(m) and m > 0 else s.fillna(0.0) * 0.0

    a = _norm(pd.to_numeric(df["score_pchembl"], errors="coerce"))
    b = _norm(pd.to_numeric(df["score_dgidb"], errors="coerce"))
    best = pd.concat([a, b], axis=1).max(axis=1)
    order = best.sort_values(ascending=False).head(int(n))
    return sorted({str(g) for g in df.loc[order.index, "approvedSymbol"].dropna()})


def _gene_list(df: pd.DataFrame) -> list[str]:
    return sorted({str(g) for g in df["approvedSymbol"].dropna()}) if df is not None and not df.empty else []


@router.get("/api/synergy/results")
def synergy_results_list() -> dict:
    """协同分析运行清单(结果管理页,新→旧,行级摘要)。"""
    return {"runs": list_results()}


class SynergyRunsBody(BaseModel):
    runs: list[str]


@router.post("/api/synergy/results/delete")
def synergy_results_delete(body: SynergyRunsBody) -> dict:
    """删除勾选的运行目录(逐个校验,被拒项进 skipped)。"""
    return delete_runs(body.runs)


@router.post("/api/synergy/results/zip")
def synergy_results_zip(body: SynergyRunsBody) -> Response:
    """勾选运行目录打包 ZIP(批量下载);越权名 400 带明细。"""
    try:
        data = zip_runs(body.runs)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    filename = (f"synergy_results_{len(body.runs)}runs.zip"
                if len(body.runs) > 1 else f"{body.runs[0]}.zip")
    return Response(data, media_type="application/zip",
                    headers={"Content-Disposition": f'attachment; filename="{filename}"'})


@router.get("/api/synergy/results/detail")
def synergy_results_detail(run: str) -> dict:
    """单次运行的契约 JSON(结果详情渲染/复现回填)。"""
    try:
        return read_run(run)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except OSError as exc:
        raise HTTPException(status_code=404, detail=str(exc))


@router.get("/api/synergy/library")
def synergy_library() -> dict:
    """选择器数据源:本地 merged 缓存已收录的药物/疾病查询名(只读本地,不联网)。
    前端 datalist 下拉 + 保留自由手输。"""
    return {"drugs": list_known_targets("drug"), "diseases": list_known_targets("disease")}


@router.post("/api/synergy/analyze")
def synergy_analyze(body: AnalyzeBody) -> dict:
    drug_a, drug_b, disease = body.drug_a.strip(), body.drug_b.strip(), body.disease.strip()
    if not drug_a or not drug_b:
        raise HTTPException(status_code=400, detail="drug_a and drug_b are required")
    if body.top_n is not None and body.top_n < 1:
        raise HTTPException(status_code=400, detail="top_n must be a positive integer or null")
    if body.n_perm is not None and body.n_perm < 10:
        raise HTTPException(status_code=400, detail="n_perm must be >= 10 or null")

    def task(progress) -> dict:
        progress(fraction=0.02, stage="targets")
        d_genes = _gene_list(get_disease_targets_merged(disease)) if disease else []
        if disease and not d_genes:
            raise SynergyNoTargetsError(f"synergy:no_targets_disease: {disease}")
        ga = _top_n_slice(get_drug_targets_merged(drug_a), body.top_n)
        progress(fraction=0.06, stage="targets")
        gb = _top_n_slice(get_drug_targets_merged(drug_b), body.top_n)
        if not ga:
            raise SynergyNoTargetsError(f"synergy:no_targets_a: {drug_a}")
        if not gb:
            raise SynergyNoTargetsError(f"synergy:no_targets_b: {drug_b}")

        def core_progress(fraction=None, **kw) -> None:
            stage = "separation" if (fraction or 0.0) < 0.6 else "proximity"
            progress(fraction=fraction, stage=stage)

        contract = analyze_pair(
            ga, gb, d_genes, drug_a=drug_a, drug_b=drug_b,
            disease_name=disease, disease_source="platform_targets",
            n_perm=body.n_perm, proximity_n_perm=body.n_perm,
            top_n=body.top_n, sort_rule=SORT_RULE, progress=core_progress,
        )

        # 归档进 synergy-results/(协同分析结果页专管,downloads 清单已排除该类别)
        dir_name = (f"synergy_{safe_name(drug_a)}_x_{safe_name(drug_b)}"
                    f"_{time.strftime('%Y%m%d_%H%M%S')}")
        out_dir = os.path.join(config.DOWNLOADS_DIR, "synergy-results", dir_name)
        os.makedirs(out_dir, exist_ok=True)
        save_json_file(os.path.join(out_dir, "result.json"), contract, indent=1)
        progress(fraction=1.0, stage="done")
        return {**contract, "directory_rel": dir_name}

    job_id = jobs.submit(task, f"synergy {drug_a} x {drug_b}")
    with _analyze_job_lock:
        global _analyze_job_id
        _analyze_job_id = job_id
    return {"job_id": job_id}


@router.get("/api/synergy/job")
def synergy_job_status() -> dict:
    """最近一次协同分析任务的状态(页面重进时恢复进度展示;跨页离开任务照跑)。"""
    with _analyze_job_lock:
        job_id = _analyze_job_id
    if job_id is None:
        return {"running": False}
    job = jobs.get(job_id)
    if job is None:
        return {"running": False, "job_id": job_id, "status": "unknown"}
    return {"running": job["status"] == "running", "job_id": job_id,
            "status": job["status"], "progress": job["progress"]}
