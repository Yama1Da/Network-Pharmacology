#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""跨数据集稳健性验证(线性 benchmark 板块 3)。

对同一疾病的两个独立 GEO 数据集跑同一 ML 共识流程,比较共识基因 top-K 的
Jaccard 重叠,回答"换个数据集结论还成立吗"。

用法(在 WY 目录下):
    python scripts/cross_dataset_robustness.py \\
        --gse-a GSE1297 --case-a "group: Severe|group: Incipient|group: Moderate" \\
        --control-a "group: Control" \\
        --gse-b GSE63067 --case-b "steatosis|non-alcoholic" \\
        --control-b "healthy" \\
        --top-k 50

产出 results/cross_robustness/ 下:
    consensus_A.csv / consensus_B.csv  各数据集的共识打分表
    overlap.csv                        共同基因清单与各项指标
    summary.json                       参数快照与 Jaccard 等关键数值
"""

import argparse
import json
import os
import sys
import time
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

from core.geo import assign_groups, fetch_geo_dataset, group_fields_overview
from core.expression import differential_expression, map_probe_to_symbol, preprocess
from core.ml_screening import deg_prefilter, run_ml_screening


def run_one_dataset(gse, case_pattern, control_pattern, top_genes, bootstrap, seed):
    """对一个数据集跑 DEG → ML 共识,返回共识打分 DataFrame 和 DEG 表。"""
    print(f"\n{'='*60}")
    print(f"数据集 {gse}")
    print(f"{'='*60}")
    ds = fetch_geo_dataset(gse)
    if ds.issues:
        raise RuntimeError(f"{gse}: {';'.join(ds.issues)}")
    print(f"  样本 {len(ds.samples)} | 矩阵 {ds.table.shape} | 平台 {ds.gpl}")

    labels = assign_groups(ds, case_pattern, control_pattern)
    n_case, n_ctrl = int((labels == 1).sum()), int((labels == 0).sum())
    print(f"  分组: case={n_case} control={n_ctrl}")

    matrix, log_branch = preprocess(ds.table)
    matrix = map_probe_to_symbol(matrix, ds.gpl)
    print(f"  映射后基因数: {len(matrix)}")

    deg = differential_expression(matrix, labels)
    deg_top, _ = deg_prefilter(deg, top_n=top_genes)
    genes = deg_top["gene"].tolist()
    labeled = list(labels[labels == 1].index) + list(labels[labels == 0].index)
    X = matrix[labeled].T[genes]
    print(f"  ML 输入: {len(labeled)} 样本 x {len(genes)} 基因")

    res = ml_screening(X, labels.loc[X.index], bootstrap, seed)
    m = res["meta"]
    print(f"  {m['cv_mode']} AUC={m['cv_auc']} | 耗时 {m['elapsed_s']}s")
    return res["consensus"], deg, {"gse": gse, "n_samples": m["n_samples"],
        "n_case": m["n_case"], "n_control": m["n_control"],
        "cv_mode": m["cv_mode"], "cv_auc": m["cv_auc"],
        "n_genes": m["n_genes"]}


def ml_screening(X, y, bootstrap, seed):
    """调用 core.ml_screening 并确保结果包含共识打分。"""
    from core.ml_screening import run_ml_screening
    return run_ml_screening(X, y, bootstrap_rounds=bootstrap, seed=seed)


def main():
    ap = argparse.ArgumentParser(description="跨数据集稳健性验证(板块 3)")
    ap.add_argument("--gse-a", required=True, help="数据集 A 的 GSE 编号")
    ap.add_argument("--case-a", required=True, help="数据集 A 的病例组正则")
    ap.add_argument("--control-a", required=True, help="数据集 A 的对照组正则")
    ap.add_argument("--gse-b", required=True, help="数据集 B 的 GSE 编号")
    ap.add_argument("--case-b", required=True, help="数据集 B 的病例组正则")
    ap.add_argument("--control-b", required=True, help="数据集 B 的对照组正则")
    ap.add_argument("--top-genes", type=int, default=300, help="每数据集 DEG 预过滤基因数(默认 300)")
    ap.add_argument("--top-k", type=int, default=50, help="共识基因 top-K 用于 Jaccard(默认 50)")
    ap.add_argument("--bootstrap", type=int, default=10, help="bootstrap 轮次(默认 10)")
    ap.add_argument("--seed", type=int, default=0, help="随机种子(默认 0)")
    ap.add_argument("--out", help="结果输出目录")
    args = ap.parse_args()

    from core.ml_screening import deg_prefilter

    out_dir = args.out or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "results", "cross_robustness")
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.time()
    consensus_a, deg_a, meta_a = run_one_dataset(
        args.gse_a, args.case_a, args.control_a, args.top_genes, args.bootstrap, args.seed)
    consensus_b, deg_b, meta_b = run_one_dataset(
        args.gse_b, args.case_b, args.control_b, args.top_genes, args.bootstrap, args.seed)

    # Jaccard 重叠
    top_a = set(consensus_a.head(args.top_k).index)
    top_b = set(consensus_b.head(args.top_k).index)
    overlap = top_a & top_b
    jaccard = len(overlap) / len(top_a | top_b) if (top_a | top_b) else 0.0

    print(f"\n{'='*60}")
    print(f"跨数据集稳健性结果")
    print(f"  {args.gse_a} 共识 top-{args.top_k}: {len(top_a)} 个")
    print(f"  {args.gse_b} 共识 top-{args.top_k}: {len(top_b)} 个")
    print(f"  交集: {len(overlap)} 个")
    print(f"  Jaccard 重叠: {jaccard:.3f}")
    if overlap:
        print(f"  共同基因: {sorted(overlap)[:20]}")

    # 落盘
    consensus_a.to_csv(os.path.join(out_dir, f"consensus_{args.gse_a}.csv"),
                       index_label="gene", encoding="utf-8-sig")
    consensus_b.to_csv(os.path.join(out_dir, f"consensus_{args.gse_b}.csv"),
                       index_label="gene", encoding="utf-8-sig")
    overlap_df = pd.DataFrame(sorted(overlap), columns=["gene"])
    overlap_df.to_csv(os.path.join(out_dir, "overlap.csv"), index=False,
                      encoding="utf-8-sig")
    summary = {
        "module": "cross-dataset-robustness",
        "gse_a": args.gse_a, "case_a": args.case_a, "control_a": args.control_a,
        "gse_b": args.gse_b, "case_b": args.case_b, "control_b": args.control_b,
        "top_k": args.top_k, "top_genes": args.top_genes,
        "bootstrap": args.bootstrap, "seed": args.seed,
        "jaccard": round(jaccard, 3),
        "overlap_genes": sorted(overlap),
        "meta_a": meta_a, "meta_b": meta_b,
        "elapsed_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    print(f"\n结果目录: {out_dir}")


if __name__ == "__main__":
    main()
