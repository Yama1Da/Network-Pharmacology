#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""消融 benchmark(ML 优化期 B1):九臂同台,回答"六算法融合是否真的有用"。

同一折内协议(core.ml_screening.crossval_roc,经 select_fn/select_seed 钩子注入),
只换"折内选基因"那一步,臂间共享同一折划分(配对比较)。九臂:

    fusion(默认生产路径,兼作自检:其 AUC 必须与正常 run 的 cv_auc 一致)
    lasso / enet / ridge / rf_boruta / svm_rfe / xgboost   单算法臂(只跑该选择器)
    deg_only        折内重跑 Welch t 按 FDR 序取 K(裸统计对照)
    random          折内池均匀抽 K;固定划分、只换抽签整 CV 重跑 N 次(零分布)

三条公平性护栏(实现已内置,解读时请核对报告):
    1. 全臂共享同一折划分(同 seed)——配对比较的前提;
    2. 全臂走同一折内口径(cv_matrix 全基因矩阵折内重筛)——否则比的是泄漏红利;
    3. K 全臂统一(生产同款 _top_k 规则;数据准备与 core.ml_pipeline 同构)。

用法(在 WY 目录下,二选一数据源):
    # 本地表达矩阵(GSE253502 口径,推荐主战场)
    python benchmark/ml_ablation_benchmark.py \
        --expr results/gse253502_big/expr.csv \
        --labels results/gse253502_big/labels.csv --case-label T2DM --tag gse253502
    # GEO series matrix(缓存命中时离线)
    python benchmark/ml_ablation_benchmark.py \
        --gse GSE281291 --case-pat T2DM --ctrl-pat healthy
    # 合成数据冒烟(埋信号=前 10 基因,验证管线接线,不出结论)
    python benchmark/ml_ablation_benchmark.py --synthetic

产出 <结果根>/ml_ablation/<tag>_<时间戳>/:
    ablation.csv     臂 × 指标(auc/cv_mean±sd/ΔAUC vs fusion/折审计)
    random_null.csv  random 臂 N 次 AUC 零分布(经验 p 用)
    summary.json + report.md(三结局判定指引内置)

注意:3v3 这类 LOOCV 小样本的 AUC 只能 ±1/K 网格跳变,臂间差异仅方向参考;
结论以 ≥每组 5-6 样本、5 折×5 重复的数据集为准(如 GSE253502)。
"""

import argparse
import json
import os
import sys
import time
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# random 臂零分布抽签次数:经验 p 粒度 ~1/(N+1)≈0.005,与方案 §B1 定稿一致
_N_RANDOM_DRAWS = 200
_RANDOM_SEED_OFFSET = 55021


def _config():
    from WY import config
    return config


# ---------------- 数据准备(与 core/ml_pipeline 同构,防漂移逐行对照) ----------------


def _build_matrices(args) -> tuple[pd.DataFrame, pd.Series, str]:
    """返回 (候选池 X: 样本×基因, 标签 y, 全基因宇宙 cv_universe: 样本×基因)。

    与 ml_pipeline.run_screen_chain 的准备段逐行同构(DEG→预筛选基因、重复列
    均值聚合、仅取有标签样本);差异只有一处:这里不进 _standardize_features
    ——fold_mode 下 crossval_roc 只用 X 的索引排划分,数值全部取自 cv_universe
    并由折内 scaler 标准化,X 的取值本就不参与。
    """
    from core.expression import (differential_expression, map_probe_to_symbol,
                                 preprocess)
    from core.geo import assign_groups, fetch_geo_dataset
    from core.ml_pipeline import load_local_expr
    from core.ml_screening import deg_prefilter

    if args.expr:
        matrix, labels, tag, _ = load_local_expr(
            args.expr, args.labels, args.case_label, map_ids=args.map_ids)
        log_branch = "local(启发式 log 判定在 load_local_expr 内)"
    else:
        gse = args.gse.strip().upper()
        tag = gse
        ds = fetch_geo_dataset(gse, force_refresh=False)
        if ds.issues:
            sys.exit(f"{gse}: {';'.join(ds.issues)}——B1 暂不支持 counts 回退,"
                     "请换 series matrix 数据集或用 --expr 本地矩阵")
        labels = assign_groups(ds, args.case_pat, args.ctrl_pat)
        matrix, log_branch = preprocess(ds.table)
        matrix = map_probe_to_symbol(matrix, ds.gpl)
        print(f"数据集 {gse}: 样本 {len(ds.samples)} | 平台 {ds.gpl} | {log_branch}")
    n_case = int((labels == 1).sum())
    n_ctrl = int((labels == 0).sum())
    print(f"分组 case={n_case} control={n_ctrl}")

    deg = differential_expression(matrix, labels)
    deg_top, _ = deg_prefilter(deg)
    genes = deg_top["gene"].tolist()

    labeled = list(labels[labels == 1].index) + list(labels[labels == 0].index)
    # 与 ml_pipeline 同款:重复基因列(符号级多探针/映射多对一)均值聚合,
    # 否则 X 重复列会在 ML 阶段以难懂的库错误崩溃(2026-11 P0-C)
    X = matrix[labeled].T[genes]
    if X.columns.duplicated().any():
        X = X.T.groupby(level=0).mean().T
    cv_universe = matrix[labeled].T
    if cv_universe.columns.duplicated().any():
        cv_universe = cv_universe.T.groupby(level=0).mean().T
    return X, labels.loc[X.index], cv_universe, tag


def _build_synthetic() -> tuple[pd.DataFrame, pd.Series, pd.DataFrame, str]:
    """合成冒烟:10v10、100 基因、前 10 个植入 +2.0 信号(与 tests 同口径)。"""
    rng = np.random.default_rng(42)
    genes = [f"G{i:03d}" for i in range(100)]
    X = rng.normal(0, 1, (20, 100))
    X[10:, :10] += 2.0
    X = pd.DataFrame(X, columns=genes)
    y = pd.Series([0] * 10 + [1] * 10)
    return X, y, X.copy(), "synthetic"


# ---------------- 九臂 select_fn(签名 (Xtr, ytr, sel_seed, k) → 基因列表) ----------------


def _single(fn):
    """单算法臂:只跑该选择器的重要性,按秩取前 k(不跑全六再取一,省 6 倍成本)。
    rank 后必须 sort_values 再 head——rank 返回的 Series 仍按原列序排列,池的
    原序恰是 FDR 序,漏 sort 会把所有单算法臂静默降级成 deg_only 的克隆
    (2026-11-14 实测踩坑:七臂 AUC 分毫不差即是此症的指纹)。"""
    def _f(Xtr, ytr, sel_seed, k):
        imp, _ = fn(Xtr, ytr, sel_seed)
        return list(imp.rank(ascending=False, method="average")
                    .sort_values().head(k).index)
    return _f


def _arm_deg_only(Xtr, ytr, sel_seed, k):
    """裸 DEG 对照:折内重跑 Welch t 按 (FDR, p) 序取 K——与折内圈池同方法,
    对照"ML 层相对裸统计加了什么"。"""
    from core.expression import differential_expression
    deg = differential_expression(Xtr.T, pd.Series(np.asarray(ytr), index=Xtr.index))
    return list(deg.sort_values(["FDR", "p"]).head(k)["gene"])


def _arm_random(Xtr, ytr, sel_seed, k):
    """随机臂:折内池均匀抽 K(不放回)。随机流由 sel_seed 驱动——固定划分、
    只换抽签整 CV 重跑 N 次,即得同划分下的选基因零分布。"""
    rng = np.random.default_rng(sel_seed)
    cols = list(Xtr.columns)
    idx = rng.choice(len(cols), size=min(k, len(cols)), replace=False)
    return [str(cols[i]) for i in idx]


def _build_arms() -> dict:
    """臂名 → select_fn(None = 生产默认融合路径)。顺序即 ablation.csv 行序。"""
    import core.ml_screening as ml

    from functools import partial

    logistic = ml._sel_logistic
    return {
        "fusion": None,
        "lasso": _single(partial(logistic, solver="saga", l1_ratio=1.0)),
        "enet": _single(partial(logistic, solver="saga", l1_ratio=0.5)),
        "ridge": _single(partial(logistic, solver="lbfgs")),
        "rf_boruta": _single(partial(ml._sel_rf_boruta,
                                     iterations=ml._BORUTA_ITER_BOOT)),
        "svm_rfe": _single(ml._sel_svm_rfe),
        "xgboost": _single(ml._sel_xgboost),
        "deg_only": _arm_deg_only,
        "random": _arm_random,
    }


# ---------------- 运行器 ----------------


def _run_arm(ml, X, y, cv_universe, arm_name, fn, cfg, n_random, random_seeds):
    """单臂执行:常规臂一次 crossval_roc;random 臂 N 次(同划分换抽签)取分布。"""
    t0 = time.time()
    if arm_name != "random":
        auc, y_true, scores, mode, aucs, folds = ml.crossval_roc(
            X, y, cfg.ML_SEED, cv_matrix=cv_universe,
            repeats=int(cfg.ML_CV_REPEATS), select_fn=fn)
        cv_mean = float(np.nanmean(aucs)) if aucs and not all(np.isnan(aucs)) else None
        cv_sd = (float(np.nanstd(aucs, ddof=1))
                 if len(aucs) > 1 and not all(np.isnan(aucs)) else None)
        return {"arm": arm_name, "auc": None if np.isnan(auc) else round(auc, 4),
                "cv_mode": mode, "cv_mean": None if cv_mean is None else round(cv_mean, 4),
                "cv_sd": None if cv_sd is None else round(cv_sd, 4),
                "folds_done": folds["folds_done"], "folds_total": folds["folds_total"],
                "n_draws": 1, "elapsed_s": round(time.time() - t0, 1),
                "random_aucs": None}
    aucs = []
    mode = None
    folds_total = folds_done = 0
    for i, s in enumerate(random_seeds):
        auc, _, _, mode, _, folds = ml.crossval_roc(
            X, y, cfg.ML_SEED, cv_matrix=cv_universe,
            repeats=int(cfg.ML_CV_REPEATS), select_fn=fn, select_seed=int(s))
        aucs.append(None if np.isnan(auc) else round(float(auc), 4))
        folds_total = max(folds_total, folds["folds_total"])
        folds_done += folds["folds_done"]
        if (i + 1) % 50 == 0:
            print(f"  random {i + 1}/{n_random} 次完成(累计 {time.time() - t0:.0f}s)")
    return {"arm": "random", "auc": None, "cv_mode": mode,
            "cv_mean": round(float(np.nanmean(aucs)), 4),
            "cv_sd": round(float(np.nanstd(aucs, ddof=1)), 4),
            "folds_done": int(folds_done), "folds_total": int(folds_total),
            "n_draws": int(n_random), "elapsed_s": round(time.time() - t0, 1),
            "random_aucs": aucs}


def _write_report(out_dir: str, s: dict) -> None:
    lines = [
        "# 消融 benchmark 报告(B1)",
        f"- 数据集: `{s['tag']}` | n={s['n_samples']} | K=生产同款 _top_k | "
        f"random 零分布 {s['n_random']} 次(同划分换抽签)",
        "- 公平性护栏:全臂同折划分、同折内口径、同 K;fusion 臂=生产默认路径",
        "",
        "| 臂 | AUC(池化) | cv_mean±sd | ΔAUC vs fusion | 经验 p vs random | 折审计 |",
        "|---|---|---|---|---|---|",
    ]
    fusion = next((r for r in s["arms"] if r["arm"] == "fusion"), None)
    rand = next((r for r in s["arms"] if r["arm"] == "random"), None)
    null_aucs = [a for a in (rand or {}).get("random_aucs") or [] if a is not None]
    for r in s["arms"]:
        dauc = (None if fusion is None or r["auc"] is None or fusion["auc"] is None
                else round(r["auc"] - fusion["auc"], 4))
        p_emp = None
        if null_aucs and r["auc"] is not None:
            p_emp = round((1 + sum(a >= r["auc"] for a in null_aucs)) / (1 + len(null_aucs)), 4)
        audit = f"{r['folds_done']}/{r['folds_total']}" + (
            f"×{r['n_draws']}次" if r["arm"] == "random" else "")
        lines.append(f"| {r['arm']} | {r['auc']} | {r['cv_mean']}±{r['cv_sd']} | "
                     f"{dauc} | {p_emp} | {audit} |")
    lines += ["", "## 结局判定指引(跑前已签,勿事后改口径)",
              "- 融合显著最高 → 「融合具有判别增益」;",
              "- 融合 ≈ 最优单算法但跨数据集更稳 → 主张改写为**稳定性价值**"
              "(用 cross_dataset_robustness 的 Jaccard 佐证);",
              "- 融合 ≈ 裸 DEG → 诚实结论:ML 层未加判别信息(降级/明示/移除三选一)——"
              "**这是科学结果,不是失败**;",
              "- 3v3 LOOCV 数据集的 AUC 只能 ±1/K 网格跳变,仅作方向参考,结论以"
              " ≥每组 5 样本、重复分层 CV 的数据集为准。"]
    with open(os.path.join(out_dir, "report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="ML 九臂消融 benchmark(B1)")
    ap.add_argument("--expr", help="本地表达矩阵 CSV(样本×基因或基因×样本自动识别)")
    ap.add_argument("--labels", help="两列 CSV:样本ID, 分组")
    ap.add_argument("--case-label", default=None, help="病例组取值(缺省自动推断)")
    ap.add_argument("--map-ids", action="store_true",
                    help="行标识非符号级时经 mygene 映射(需联网)")
    ap.add_argument("--gse", help="GEO series matrix 数据集(缓存命中可离线)")
    ap.add_argument("--case-pat", default=None, help="病例组正则(GEO 模式)")
    ap.add_argument("--ctrl-pat", default=None, help="对照组正则(GEO 模式)")
    ap.add_argument("--synthetic", action="store_true", help="合成数据冒烟(不出结论)")
    ap.add_argument("--tag", default=None, help="产物目录名(缺省取数据集标识)")
    ap.add_argument("--n-random", type=int, default=None,
                    help="random 臂零分布次数(默认 200)")
    ap.add_argument("--out", default=None, help="输出根目录(默认 <结果根>/ml_ablation)")
    args = ap.parse_args()
    modes = sum(bool(a) for a in (args.expr, args.gse, args.synthetic))
    if modes != 1:
        sys.exit("数据源三选一:--expr(+--labels)/ --gse(+--case-pat --ctrl-pat)/ --synthetic")
    if args.expr and not args.labels:
        sys.exit("--expr 需同时提供 --labels(两列 CSV:样本ID, 分组)")
    if args.gse and not (args.case_pat and args.ctrl_pat):
        sys.exit("--gse 需同时提供 --case-pat 与 --ctrl-pat")

    cfg = _config()
    n_random = int(args.n_random if args.n_random is not None else _N_RANDOM_DRAWS)
    rng = np.random.default_rng(int(cfg.ML_SEED) + _RANDOM_SEED_OFFSET)
    random_seeds = [int(v) for v in rng.integers(1, 10 ** 9, size=n_random)]

    if args.synthetic:
        X, y, cv_universe, tag = _build_synthetic()
    else:
        X, y, cv_universe, tag = _build_matrices(args)
    tag = args.tag or tag
    k_eff = _ml()._top_k(min(int(cfg.ML_DEG_PREFILTER_TOP_N), cv_universe.shape[1]))
    print(f"== B1 消融 ==\n数据集 {tag} | n={len(X)} | 宇宙基因 {cv_universe.shape[1]} "
          f"| K={k_eff}(生产同款 _top_k(预筛{int(cfg.ML_DEG_PREFILTER_TOP_N)})) "
          f"| random {n_random} 次")

    arms = _build_arms()
    ml = _ml()
    rows = []
    for arm_name, fn in arms.items():
        print(f"[{arm_name}] 运行中…")
        r = _run_arm(ml, X, y, cv_universe, arm_name, fn, cfg, n_random, random_seeds)
        rows.append(r)
        print(f"[{arm_name}] auc={r['auc']} cv_mean={r['cv_mean']}±{r['cv_sd']} "
              f"折 {r['folds_done']}/{r['folds_total']} 耗时 {r['elapsed_s']}s")

    fusion_auc = next((r["auc"] for r in rows if r["arm"] == "fusion"), None)
    null_aucs = next((r["random_aucs"] for r in rows if r["arm"] == "random"), []) or []
    for r in rows:
        r["delta_vs_fusion"] = (None if r["auc"] is None or fusion_auc is None
                                else round(r["auc"] - fusion_auc, 4))
        r["p_vs_random"] = (round((1 + sum(a >= r["auc"] for a in null_aucs if a is not None))
                                  / (1 + len([a for a in null_aucs if a is not None])), 4)
                            if (null_aucs and r["auc"] is not None) else None)

    out_base = args.out or os.path.join(_default_root(), "ml_ablation")
    out_dir = os.path.join(out_base, f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}")
    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame([{k: v for k, v in r.items() if k != "random_aucs"} for r in rows]
                 ).to_csv(os.path.join(out_dir, "ablation.csv"), index=False, encoding="utf-8-sig")
    rand_row = next((r for r in rows if r["arm"] == "random"), None)
    pd.DataFrame({"draw": range(1, len(null_aucs) + 1),
                  "auc": rand_row["random_aucs"]}).to_csv(
        os.path.join(out_dir, "random_null.csv"), index=False, encoding="utf-8-sig")
    summary = {"module": "ml-ablation-benchmark", "b1_version": 1,
               "tag": tag, "n_samples": int(len(X)),
               "universe_n": int(cv_universe.shape[1]),
               "k_rule": "生产同款 _top_k", "n_random": n_random,
               "seed": int(cfg.ML_SEED),
               "arms": [{k: v for k, v in r.items() if k != "random_aucs"} for r in rows],
               "random_aucs": null_aucs, "args": vars(args)}
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    _write_report(out_dir, summary)
    print(f"落盘 → {out_dir}(ablation.csv / random_null.csv / report.md / summary.json)")


def _ml():
    import core.ml_screening as ml
    return ml


def _default_root() -> str:
    from core.ml_pipeline import default_results_root
    return default_results_root()


if __name__ == "__main__":
    main()
