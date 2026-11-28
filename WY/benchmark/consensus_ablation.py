"""消融对比(二期方案 §4.2,清单 §1.13【本人,产出负责人复核】):
一致子集 vs 单方法模块作为下游种子,量化"地基加固"的实际增益。

实验协议(钉死,可复现):
- 种子变体 S0..S4:S0 = 原始靶点集(基线,一期口径);S1–S3 = 三方法单方法模块;
  S4 = 一致子集(被 ≥k 种方法捞入,Borda 共识)。
- 参考真值 R = DisGeNET curated 高分池 − 当前种子集(留出集:从不进入任何输入
  ——无泄漏;分层复用 1.8 stratify_pool 同一实现,正分数中位数口径)。
- 下游 RWR 腿(α/tol/iter 同一期 config):各变体作种子跑 RWR,统计留出真值在
  top-K 的 Recall(K ∈ 截断集)与命中平均秩——"更好的种子集能否把真关联捞回前排"。
- 邻近度腿(--drug 给出时):proximity(药物靶点, 变体) 的 z 与经验 p——疾病模块
  越扎实,药物-疾病分离度检验越可信;同 ctx 纯计算版,度匹配零模型。
- 与 1.8 的分工:1.8 回测"共识扩张能否捞回隐藏疾病基因"(种子集质量,留一法);
  本脚本量化"更好的种子集对下游检验的实际影响"(下游增益)——前者是根,后者是果。

输出:ablation_table.csv(变体×指标长表)+ summary.json + ablation_chart.png
(Recall@K 分组柱状;matplotlib Figure API,图内英文)。全部确定性:RWR/共识/邻近度
零模型种子均取 config 固定值,无自由随机源。
"""

import argparse
import json
import os
import sys

import pandas as pd


def run_ablation(ctx, seeds, reference, *, n=None, k=None,
                 cutoffs=(50, 100, 200), drug_targets=None, n_perm=None):
    """纯计算核:变体构造 + 下游 RWR Recall@K(+可选邻近度腿)。
    ctx = network_algorithms._lcc_context 注入;reference = 留出真值集合;
    drug_targets 给出时追加邻近度腿(同 ctx 纯计算版,度匹配零模型)。
    返回 {"rows": 变体×指标长表, "reference_size", "consistent_size",
    "cutoffs", "proximity": bool}。"""
    from WY import config
    from WY.core.disease_modules import consensus_targets
    from WY.core.network_algorithms import (_proximity_zscore_from_context,
                                            _rwr_from_context)

    cons = consensus_targets(ctx, seeds, n=n, k=k)
    if not cons.get("ok"):
        raise RuntimeError(f"consensus failed: {cons.get('error')}")
    consistent = sorted(
        cons["consensus"].loc[cons["consensus"]["consistent"].astype(bool),
                              "gene"].astype(str))
    variants: list[tuple[str, list[str]]] = [
        ("S0_seeds", sorted(set(seeds))),
        ("S1_rwr", sorted(cons["modules"]["rwr"])),
        ("S2_diamond", sorted(cons["modules"]["diamond"])),
        ("S3_neighbor", sorted(cons["modules"]["neighbor"])),
        ("S4_consistent", consistent),
    ]

    ref = set(reference)
    alpha, tol, iters = config.RWR_ALPHA, config.RWR_TOL, config.RWR_MAX_ITER
    rows: list[dict] = []
    for name, seedset in variants:
        res = _rwr_from_context(ctx, seedset, alpha, tol, iters, top_nonseed=None)
        if not res.get("ok"):
            raise RuntimeError(f"RWR failed for {name}: {res.get('error')}")
        # 全局排名:种子行与非种子行各按分数降序,合并重排(两者天然有序但互相交错)
        ranked = [r["gene"] for r in sorted(res["seeds"] + res["top_nonseed"],
                                            key=lambda r: -r["score"])]
        row: dict = {"variant": name, "n_variant": len(seedset),
                     "n_variant_in_lcc": int(res.get("n_seeds") or 0)}
        for k_ in cutoffs:
            topk = set(ranked[:k_])
            row[f"recall@{k_}"] = (len(topk & ref) / len(ref)) if ref else None
        hits = [i + 1 for i, g in enumerate(ranked) if g in ref]
        row["ref_hits"] = len(hits)
        row["ref_mean_rank"] = (sum(hits) / len(hits)) if hits else None
        if drug_targets is not None:
            pz = _proximity_zscore_from_context(
                ctx, list(drug_targets), list(seedset),
                n_perm=int(n_perm if n_perm is not None else config.PROXIMITY_N_PERM),
                seed=int(config.PROXIMITY_SEED),
                bin_size=config.PROXIMITY_DEGREE_BIN_SIZE,
                drug_ratio=config.PROXIMITY_MIN_COVERAGE,
                disease_ratio=config.PROXIMITY_MIN_COVERAGE)
            row["prox_ok"] = bool(pz.get("ok"))
            row["prox_z"] = pz.get("z")
            row["prox_p"] = pz.get("p_empirical")
            row["prox_d_c"] = pz.get("d_c")
        rows.append(row)
    return {"rows": rows, "reference_size": len(ref),
            "consistent_size": len(consistent), "cutoffs": list(cutoffs),
            "proximity": drug_targets is not None,
            "k": cons.get("k")}


def render_chart(rows, cutoffs, path: str, title_note: str = "") -> str:
    """Recall@K 分组柱状(变体 × 截断集)。Figure API + Agg,图内英文。"""
    from matplotlib.figure import Figure

    variants = [r["variant"] for r in rows]
    fig = Figure(figsize=(9.6, 5.0))
    ax = fig.add_subplot(111)
    width = 0.8 / max(len(cutoffs), 1)
    for i, k_ in enumerate(cutoffs):
        vals = [float(r[f"recall@{k_}"] or 0.0) for r in rows]
        xs = [j + i * width - 0.4 + width / 2 for j in range(len(variants))]
        ax.bar(xs, vals, width=width, label=f"Recall@{k_}")
    ax.set_xticks(list(range(len(variants))))
    ax.set_xticklabels(variants, fontsize=8)
    ax.set_ylabel("Recall of held-out high-score genes")
    ax.set_title("Consensus ablation: downstream RWR recall by seed-set variant"
                 + (f" ({title_note})" if title_note else ""))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Consensus ablation (plan §4.2): consistent subset vs "
                    "single-method modules as downstream seeds")
    ap.add_argument("--disease", required=True,
                    help="disease name (DisGeNET reference pool + default seed "
                         "source when --genes absent)")
    ap.add_argument("--genes", default="",
                    help="optional explicit seed list (comma-separated or file "
                         "path); default = merged disease target set T_D")
    ap.add_argument("--drug", default="",
                    help="optional drug name for the proximity leg")
    ap.add_argument("--top-n", type=int, default=None,
                    help="module expansion size N (default = config)")
    ap.add_argument("--k", type=int, default=None, help="consistency marker k")
    ap.add_argument("--cutoffs", default="50,100,200",
                    help="comma-separated Recall@K cutoffs")
    ap.add_argument("--out", required=True, help="output directory")
    args = ap.parse_args(argv)

    from WY import config
    from WY.core.disease_modules import _module_size
    from WY.core.interactome import is_built
    from WY.core.network_algorithms import _lcc_context
    from WY.core.targets import (get_disease_targets_merged,
                                 get_drug_targets_merged)
    from WY.core.disgenet import fetch_disease_targets_disgenet
    from .disease_module_benchmark import stratify_pool

    if not is_built():
        print("[ablation] background interactome not built; build it first",
              file=sys.stderr)
        return 1
    cutoffs = tuple(int(c) for c in str(args.cutoffs).split(",") if c.strip())
    n_eff = _module_size(args.top_n)

    seeds: list[str] = []
    if args.genes.strip():
        text = args.genes
        if os.path.isfile(text):
            with open(text, encoding="utf-8-sig") as f:
                text = f.read()
        seeds = sorted({t.strip().upper() for t in
                        text.replace(";", ",").replace("\n", ",").split(",")
                        if t.strip()})
    if not seeds:
        df = get_disease_targets_merged(args.disease)
        seeds = sorted({str(g).strip().upper()
                        for g in df["approvedSymbol"].dropna()}) if (
            df is not None and not df.empty) else []
    if not seeds:
        print("[ablation] no seeds resolved", file=sys.stderr)
        return 1

    assoc = fetch_disease_targets_disgenet(args.disease)
    if assoc is None or assoc.empty:
        print("[ablation] no DisGeNET associations for the reference pool",
              file=sys.stderr)
        return 1
    strata = stratify_pool(assoc)
    reference = sorted(strata["high"] - set(seeds))
    if not reference:
        print("[ablation] reference pool empty after removing seeds; nothing "
              "to measure", file=sys.stderr)
        return 1

    drug_targets = None
    if args.drug.strip():
        ddf = get_drug_targets_merged(args.drug)
        drug_targets = sorted({str(g).strip().upper()
                               for g in ddf["approvedSymbol"].dropna()}) if (
            ddf is not None and not ddf.empty) else []

    ctx = _lcc_context("Human")
    res = run_ablation(ctx, seeds, reference, n=n_eff, k=args.k,
                       cutoffs=cutoffs, drug_targets=drug_targets)

    os.makedirs(args.out, exist_ok=True)
    table_path = os.path.join(args.out, "ablation_table.csv")
    pd.DataFrame(res["rows"]).to_csv(table_path, index=False, encoding="utf-8-sig")
    chart_path = render_chart(res["rows"], cutoffs,
                              os.path.join(args.out, "ablation_chart.png"),
                              title_note=args.disease)
    summary = {
        "schema": "consensus_ablation.v1",
        "disease": args.disease, "drug": args.drug or None,
        "seeds_total": len(seeds),
        "reference_size": res["reference_size"],
        "reference_rule": "DisGeNET curated high pool (>= median score) minus seeds",
        "threshold": float(strata["threshold"]),
        "params": {"n": n_eff, "k": res["k"], "cutoffs": list(cutoffs)},
        "consistent_size": res["consistent_size"],
        "files": ["ablation_table.csv", "ablation_chart.png", "summary.json"],
    }
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(f"[ablation] seeds={len(seeds)} reference(held-out high)={len(reference)} "
          f"N={n_eff} k={res['k']}")
    for r in res["rows"]:
        recalls = " ".join(f"R@{k}={r[f'recall@{k}']:.3f}" for k in cutoffs)
        print(f"  {r['variant']:<14} |{r['n_variant']:>5} | {recalls}")
    print(f"[ablation] written: {table_path} | {chart_path} | "
          f"{os.path.join(args.out, 'summary.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
