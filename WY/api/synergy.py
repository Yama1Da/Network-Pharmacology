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

import threading
import time

from fastapi import APIRouter, HTTPException, Response
from pydantic import BaseModel

from WY.core import jobs
from WY.core.synergy import (
    SORT_RULE,
    analyze_pair,
    archive_run,
    delete_runs,
    list_results,
    mount_reachability,
    read_run,
    slice_targets_top_n,
    zip_runs,
)
from WY.core.targets import (
    get_disease_targets_merged,
    get_drug_targets_merged,
    list_known_targets,
)

router = APIRouter()

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


def _gene_list(df) -> list[str]:
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
    """选择器数据源(只读本地,不联网),两路取并集:
    ① merged 缓存命中名(list_known_targets);
    ② 查询历史台账(query_history.names——双源整合成功过的全部查询名,
    与 merged 缓存生命周期解耦,设置页 purge 后下拉不缩水)。
    前端 datalist 下拉 + 保留自由手输。"""
    from WY.core.query_history import names as qh_names

    drugs = sorted(set(list_known_targets("drug")) | set(qh_names("drug")))
    diseases = sorted(set(list_known_targets("disease")) | set(qh_names("disease")))
    return {"drugs": drugs, "diseases": diseases}


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
        ma = get_drug_targets_merged(drug_a)  # v1.2:留完整明细 df 供 E6 挂载
        ga = slice_targets_top_n(ma, body.top_n)  # 切片只喂核心(定案 #7)
        progress(fraction=0.06, stage="targets")
        mb = get_drug_targets_merged(drug_b)
        gb = slice_targets_top_n(mb, body.top_n)
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
        contract = mount_reachability(contract, ma, mb)  # v1.2 E6 挂载(明细 df 复用)

        # 归档进 synergy-results/(协同分析结果页专管,downloads 清单已排除该类别);
        # archive 内先外移 reachability 明细 CSV,之后契约才是纯 JSON(轮询回放安全)
        dir_name = archive_run(contract, drug_a, drug_b)
        progress(fraction=1.0, stage="done")
        return {**contract, "directory_rel": dir_name}

    job_id = jobs.submit(task, f"synergy {drug_a} x {drug_b}")
    with _analyze_job_lock:
        global _analyze_job_id
        _analyze_job_id = job_id
    return {"job_id": job_id}


@router.get("/api/synergy/job")
def synergy_job_status() -> dict:
    """最近一次协同分析任务的状态(跨页离开任务照跑,页面重进恢复):
    running 照常吐 job_id 供续接;已完成但在回填窗口内(jobs.RESULT_REPLAY_WINDOW_S)
    也吐——前端 resume 第一拍命中 done 即回填结果面板,同当场跑完的体验;
    超窗不吐,空表单,结果走协同分析结果页。"""
    with _analyze_job_lock:
        job_id = _analyze_job_id
    if job_id is None:
        return {"running": False}
    job = jobs.get(job_id)
    if job is None:
        return {"running": False, "job_id": job_id, "status": "unknown"}
    if job["status"] != "running" \
            and time.time() - job["created"] > jobs.RESULT_REPLAY_WINDOW_S:
        return {"running": False}
    return {"running": job["status"] == "running", "job_id": job_id,
            "status": job["status"], "progress": job["progress"]}
