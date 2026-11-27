"""疾病模块多方法一致性(二期方案 §4.4 工程落点:API 层,数据网络页面板后端)。

core/disease_modules.py 的三方法 + Borda 共识在网页侧的入口:POST 发起 jobs
异步任务,前端经通用 /api/jobs/{job_id} 轮询(同 rerun-modules 模式)。产物落
当前分析运行目录({prefix}_consensus_targets.csv / _consensus_methods.csv /
_consensus_significance.csv),运行目录全量文件机制让它们自动进下载产物列表。

分层与性能红线:
- 与主管线 run_ppi_analysis 完全解耦——DIAMOnD 真实图单次扩张分钟级,置换检验
  = perms × 三方法各一次完整扩张,绝不允许进每次分析的主链路;
- 网页默认 n_perms=0(跳过置换检验),交互场景上限 200;RWR_NULL_PERMS=200 全量
  档位留给 CLI 与基准层(模块 docstring 成本注记同源);
- 种子来源(方案 §4.1:种子集 = 疾病靶点集 T_D):names 模式取
  get_disease_targets_merged 疾病侧四源整合全集(分析刚跑过,缓存必命中);
  genes 模式取运行目录 {prefix}_query_genes.csv(gene 列,与本次分析完全同一份,
  表单文本可能已被用户改过,以落盘产物为准),该文件在未启用 RWR 扩展的分析里
  不落盘(save_query_genes 契约),此时兜底解析请求携带的 genes_text。
数据侧失败一律返回 {"ok": False, "reason": 稳定枚举}(前端映射 i18n),不抛
HTTP 500;HTTPException 只用于请求形态错误(目录不存在/参数越界)。
"""

import os

import pandas as pd
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from WY import config
from WY.core import jobs
from WY.core.disease_modules import CONSENSUS_SCHEMA
from .common import _df_rows, _resolve_under

router = APIRouter()

# 网页交互场景的参数边界(N 放宽到小实验规模;permutation 上限钉死:DIAMOnD
# 零分布 = perms × 逐步扩张,真实图分钟级起,200 全量属 CLI/基准层档位)
_N_MIN, _N_MAX = 10, 300
_PERMS_MAX = 200
# 载荷行数上限(前端展示前 100;完整表始终由 CSV 产物承载)
_ROWS_CAP = 200
# 交集外基因清单上限(计数恒为真值,清单截断防载荷膨胀)
_OUT_OF_LCC_CAP = 30


class ConsensusBody(BaseModel):
    directory_rel: str            # downloads 相对路径(分析运行目录,产物落这里)
    prefix: str
    mode: str                     # "names" | "genes"(种子来源见模块 docstring)
    drug: str = ""
    disease: str = ""
    genes_text: str = ""          # genes 模式 query_genes.csv 缺失时的兜底解析源
    n: int | None = None          # 模块扩张规模 N,None = config.MODULE_EXPAND_N
    k: int | None = None          # 一致性标记 k,None = config.CONSENSUS_MIN_METHODS
    n_perms: int | None = 0       # 置换检验次数,网页默认 0 = 跳过


def _consensus_rows(df: pd.DataFrame, limit: int | None = None) -> list[dict]:
    """consensus 表 → JSON 安全行:秩/共识分的 NaN 显式转 None(Starlette
    JSONResponse allow_nan=False,NaN 直传会 500;同 retrieve.py r_value→null
    约定),布尔逐列转 Python 原生类型(防 numpy 标量漏进 JSON)。"""
    if limit is not None:
        df = df.head(limit)
    out: list[dict] = []
    for rec in df.to_dict(orient="records"):
        def _rank(key: str) -> int | None:
            v = rec.get(key)
            return None if pd.isna(v) else int(v)

        score = rec.get("consensus_score")
        out.append({
            "gene": str(rec["gene"]),
            "rwr_rank": _rank("rwr_rank"),
            "diamond_rank": _rank("diamond_rank"),
            "neighbor_rank": _rank("neighbor_rank"),
            "n_methods": int(rec["n_methods"]),
            "consensus_score": None if pd.isna(score) else float(score),
            "consistent": bool(rec["consistent"]),
            "is_seed": bool(rec["is_seed"]),
        })
    return out


def _resolve_seeds(body: ConsensusBody) -> list[str]:
    """种子集解析(见模块 docstring 的来源契约);解析不到基因抛 ValueError,
    稳定枚举文案由调用方映射进 reason。"""
    if body.mode == "genes":
        directory = _resolve_under(config.DOWNLOADS_DIR, body.directory_rel)
        qpath = (os.path.join(directory, f"{body.prefix}_query_genes.csv")
                 if directory else "")
        genes: list[str] = []
        if qpath and os.path.isfile(qpath):
            df = pd.read_csv(qpath, encoding="utf-8-sig", comment="#")
            col = "gene" if "gene" in df.columns else df.columns[0]
            genes = [str(g).strip().upper() for g in df[col].dropna()]
            genes = [g for g in genes if g]
        if not genes and body.genes_text.strip():
            # 兜底:未启用 RWR 扩展的分析没有 query_genes.csv(见 docstring),
            # 解析请求携带的基因文本(与手动模式表单同一套分词规则)
            tokens = body.genes_text.replace(";", ",").replace("\n", ",").split(",")
            genes = [t.strip().upper() for t in tokens if t.strip()]
        if genes:
            return sorted(set(genes))
        raise ValueError("no_seeds")
    disease = body.disease.strip()
    if not disease:
        raise ValueError("disease_required")
    from WY.core.targets import get_disease_targets_merged

    df = get_disease_targets_merged(disease)
    if df is None or df.empty or "approvedSymbol" not in df.columns:
        raise ValueError("no_disease_targets")
    genes = sorted({str(g).strip().upper() for g in df["approvedSymbol"].dropna()
                    if str(g).strip()})
    if not genes:
        raise ValueError("no_disease_targets")
    return genes


@router.post("/api/modules/consensus")
def modules_consensus(body: ConsensusBody) -> dict:
    """发起三方法一致性任务;前端经 /api/jobs/{job_id} 轮询(结果含 rows 前
    _ROWS_CAP 行 + coverage + 可选显著性),产物 CSV 落运行目录进下载列表。"""
    if body.mode not in ("names", "genes"):
        raise HTTPException(status_code=400, detail="mode must be names|genes")
    if not body.prefix.strip():
        raise HTTPException(status_code=400, detail="prefix required")
    directory = _resolve_under(config.DOWNLOADS_DIR, body.directory_rel)
    if directory is None or not os.path.isdir(directory):
        raise HTTPException(status_code=404, detail="no such analysis directory")
    if body.n is not None and not (_N_MIN <= body.n <= _N_MAX):
        raise HTTPException(status_code=400, detail=f"n must be {_N_MIN}..{_N_MAX}")
    if body.k is not None and not (1 <= body.k <= 3):
        raise HTTPException(status_code=400, detail="k must be 1..3")
    perms = int(body.n_perms or 0)
    if not (0 <= perms <= _PERMS_MAX):
        raise HTTPException(status_code=400, detail=f"n_perms must be 0..{_PERMS_MAX}")

    def task(progress) -> dict:
        from WY.core.disease_modules import (consensus_targets, export_consensus_csv,
                                             module_significance)
        from WY.core.interactome import is_built
        from WY.core.network_algorithms import _lcc_context

        if not is_built():
            return {"ok": False, "reason": "interactome_not_built"}
        progress(fraction=0.05, text="resolving seeds")
        try:
            seeds = _resolve_seeds(body)
        except ValueError as exc:
            return {"ok": False, "reason": str(exc)}
        progress(fraction=0.15,
                 text=f"three-method expansion (rwr/diamond/neighbor), seeds={len(seeds)}")
        ctx = _lcc_context("Human")
        res = consensus_targets(ctx, seeds, n=body.n, k=body.k)
        if not res.get("ok"):
            return {"ok": False, "reason": res.get("error") or "rwr_failed"}

        written: list[str] = []
        cons_path = export_consensus_csv(
            res["consensus"],
            os.path.join(directory, f"{body.prefix}_consensus_targets.csv"))
        if cons_path:
            written.append(os.path.basename(cons_path))
        # 三方法成员长表(method,rank,gene;秩从 1 起,陷阱#5 同 CLI 口径)
        methods_df = pd.DataFrame(
            [{"method": m, "rank": r, "gene": g}
             for m in ("rwr", "diamond", "neighbor")
             for r, g in enumerate(res["modules"][m], 1)])
        methods_path = os.path.join(directory, f"{body.prefix}_consensus_methods.csv")
        if not methods_df.empty:
            methods_df.to_csv(methods_path, index=False, encoding="utf-8-sig")
            written.append(os.path.basename(methods_path))

        sig_payload: dict | None = None
        if perms > 0:
            progress(fraction=0.6,
                     text=f"significance: {perms} degree-matched perms x 3 methods")
            sig = module_significance(ctx, seeds, n=body.n, n_perms=perms)
            sig_df = pd.DataFrame([{"method": m, **v}
                                   for m, v in sig["methods"].items()])
            sig_path = os.path.join(directory,
                                    f"{body.prefix}_consensus_significance.csv")
            if not sig_df.empty:
                sig_df.to_csv(sig_path, index=False, encoding="utf-8-sig")
                written.append(os.path.basename(sig_path))
            sig_payload = {m: {"p": v["p"], "n_perms": v["n_perms"],
                               "observed_cohesion": v["observed_cohesion"],
                               "null_mean": v["null_mean"], "null_sd": v["null_sd"]}
                           for m, v in sig["methods"].items()}

        progress(fraction=0.95, text="saving outputs")
        cov = res.get("coverage") or {}
        out_lcc = sorted(set(cov.get("seeds_out_of_lcc") or []))
        n_eff = (int(sig["n"]) if sig_payload else
                 (body.n if body.n is not None else int(config.MODULE_EXPAND_N)))
        cons = res["consensus"]
        consistent_list = sorted(
            cons.loc[cons["consistent"].astype(bool), "gene"].astype(str))
        return {
            "ok": True, "schema": res.get("schema", CONSENSUS_SCHEMA),
            "params": {"n": n_eff, "k": int(res.get("k") or 2),
                       "n_perms": perms if perms > 0 else None},
            "coverage": {"seeds_total": int(len(seeds)),
                         "seeds_in_lcc": int(cov.get("seeds_in_lcc") or 0),
                         "seeds_out_of_lcc_count": len(out_lcc),
                         "seeds_out_of_lcc": out_lcc[:_OUT_OF_LCC_CAP]},
            "module_sizes": res.get("module_sizes"),
            "consistent_genes": len(consistent_list),
            "consistent_list": consistent_list,
            "rows_total": int(len(cons)),
            "rows": _consensus_rows(cons, _ROWS_CAP),
            "significance": sig_payload,
            "files_written": written,
        }

    return {"job_id": jobs.submit(task, "modules consensus")}
