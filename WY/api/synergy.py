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
"""

import threading

import pandas as pd
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from WY.core import jobs
from WY.core.synergy import analyze_pair
from WY.core.targets import (
    get_disease_targets_merged,
    get_drug_targets_merged,
    list_known_targets,
)

router = APIRouter()

# 证据切片排序规则(冒烟启发版):两源分数各自归一到 [0,1] 后取最大,降序取前 N。
# 正式口径归 §1.6 DrugComb 基准扫描定案;规则 id 进契约 meta.run.sort 记账。
SORT_RULE = "dual_source_max_norm"


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

        return analyze_pair(
            ga, gb, d_genes, drug_a=drug_a, drug_b=drug_b,
            disease_name=disease, disease_source="platform_targets",
            n_perm=body.n_perm, proximity_n_perm=body.n_perm,
            top_n=body.top_n, sort_rule=SORT_RULE, progress=core_progress,
        )

    job_id = jobs.submit(task, f"synergy {drug_a} x {drug_b}")
    return {"job_id": job_id}
