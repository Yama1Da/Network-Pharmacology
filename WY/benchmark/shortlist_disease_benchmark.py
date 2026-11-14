#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""疾病相关性验证 benchmark(ML 优化期 B2):回答"ML 选出的基因对疾病真有意义吗"。

验证对象分两层,证据地位不同(工程/生物共同评审定稿,见 文件/ML优化期实施方案):
- **E1 ML 共识臂(主证据)**:共识 top-K(纯数据驱动,从未接触疾病知识)对已知
  疾病基因的富集。OT 与 DisGeNET 对此臂均为有效金标准,双库并报互为复核;
- **E1b 全 shortlist(仅构成说明)**:shortlist 第三臂(平台排序)上游在
  core/targets.py 已聚合 OT+DisGeNET 疾病靶点,故对两库"重叠"不构成独立证据,
  只做每基因证据档案 + 构成统计 + 循环声明,不做任何检验主张;
- **E2 网络邻近**:ML 共识臂相对 OT 疾病模块的邻近度,复用
  network_algorithms.proximity_zscore(度匹配置换检验、LCC、入网率保护、固定
  seed)——与 RWR 模块的 recovery 验收同款指标,数字直接可比。

统计口径(优化期方案 §B2 定稿):
- 背景宇宙 = 该 run 的候选池(ml_consensus_scores.csv 全表,即预过滤 top-N),
  不用全基因组——富集背景错误是本检验最经典的翻车点;
- 经验 p = (1 + #{随机命中 ≥ 观测}) / (1 + N 抽),seed 由 config.ML_SEED
  派生固定偏移,同输入逐值可复现;另报超几何 p 作参数化 sanity 交叉核对;
- 全部格子(K 档 × 库)无条件报告,不挑显著(防多重比较挑樱桃)。

用法(在 WY 目录下):
    python benchmark/shortlist_disease_benchmark.py --list
    python benchmark/shortlist_disease_benchmark.py \
        --run results/mlscreen/GSE281291_20261113_200236 \
        --disease "type 2 diabetes mellitus"

产出 <结果根>/disease_relevance/<tag>_<时间戳>/:
    report.md            人类可读报告(证据独立性声明内置)
    summary.json         全部数字(可复现参数快照)
    per_gene_evidence.csv  最终 shortlist 逐基因证据档案(湿实验交接用)
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

# 经验重抽样轮数(E1):抽 K 个基因成本低,1000 次 p 值分辨粒度足够;
# E2 沿用平台标准 config.PROXIMITY_N_PERM(与 RWR 模块验收数字可比,不另立口径)
_N_EMPIRICAL = 1000
# A3 同款纪律:经验重抽样 seed 固定派生,同输入两次运行逐值一致
_EMP_SEED_OFFSET = 88021


def _config():
    from WY import config
    return config


def _resolve_run(run: str | None, do_list: bool) -> str:
    """定位 run 目录:--run 显式给定,或 --list 打印最近 run 供挑选。"""
    from core.ml_results import list_ml_runs

    rows = list_ml_runs()
    if do_list or not run:
        print("可用 run(新→旧,取前 15):")
        for r in rows[:15]:
            s = r["summary"] or {}
            m = s.get("ml_meta", {})
            print(f"  {r['path']}  | gse={r['gse']} source={r['source']} "
                  f"auc={m.get('cv_auc')} v{s.get('ml_version')}"
                  f"{' 陈旧' if r['stale'] else ''}")
        if not run:
            print("请用 --run <上述路径> 指定要验证的 run")
            sys.exit(0)
    if not os.path.isdir(run):
        sys.exit(f"run 目录不存在: {run}")
    for fn in ("ml_consensus_scores.csv", "summary.json"):
        if not os.path.isfile(os.path.join(run, fn)):
            sys.exit(f"{run} 缺少 {fn}——不是完整的 ML run 目录")
    return run


def _load_universe_and_arms(run: str, k: int, k_scan: list[int]) -> tuple[pd.DataFrame, dict]:
    """读共识全表(=候选宇宙)与各验证臂;k 超出宇宙时截断到宇宙大小。"""
    consensus = pd.read_csv(os.path.join(run, "ml_consensus_scores.csv"),
                            index_col="gene", encoding="utf-8-sig")
    if "consensus_score" not in consensus.columns:
        sys.exit("ml_consensus_scores.csv 缺 consensus_score 列——请用 v3 版 run")
    universe = list(consensus.index.astype(str).str.upper())
    arms = {f"ml_top{k2}": list(consensus.index[:min(k2, len(universe))].astype(str).str.upper())
            for k2 in sorted(set(k_scan + [k]))}
    return consensus, {"universe": universe, "arms": arms,
                       "universe_n": len(universe)}


def _gold_sets(disease: str) -> dict:
    """取两库金标准:OT 沿用上游同款过滤口径(associationScore ≥ OT_MIN_SCORE);
    DisGeNET 学术 key 只含 CURATED 来源,报告须注明。任一库失败/为空不炸整体,
    标注不可用,单库证据链仍可交付。"""
    from core.disgenet import fetch_disease_targets_disgenet
    from core.opentargets import get_disease_targets

    cfg = _config()
    sets: dict[str, dict] = {}
    ot = get_disease_targets(disease)
    if ot is None or ot.empty:
        print(f"  [警告] Open Targets 未返回 '{disease}' 的靶点,该库本轮不可用")
    else:
        sym = ot["approvedSymbol"].astype(str).str.upper().str.strip()
        score = pd.to_numeric(ot["association_score"], errors="coerce").fillna(0)
        mask = score >= cfg.OT_MIN_SCORE
        if "biotype" in ot.columns:
            mask &= ot["biotype"].fillna("").eq("protein_coding")
        keep = sym[mask].drop_duplicates().tolist()
        sets["ot"] = {"genes": set(keep), "n_raw": int(sym.nunique()),
                      "threshold": cfg.OT_MIN_SCORE,
                      "note": f"association_score≥{cfg.OT_MIN_SCORE} 且 protein_coding,口径同上游 targets.py;"
                              "表内已带 genetic_association 等单证据类型分数列,子集金标准玩法可后续启用"}
    dg = fetch_disease_targets_disgenet(disease)
    if dg is None or dg.empty:
        print(f"  [警告] DisGeNET 未返回 '{disease}' 的关联(检查 API key/配额),该库本轮不可用")
    else:
        sym = dg["approvedSymbol"].astype(str).str.upper().str.strip()
        sets["disgenet"] = {"genes": set(sym.dropna().drop_duplicates()),
                            "n_raw": int(sym.nunique()), "threshold": None,
                            "note": "学术 key 仅 CURATED 来源(CLINVAR/UNIPROT/ORPHANET 等)"}
    if not sets:
        sys.exit("两库均不可用,无法做疾病相关性验证——检查网络/API key 后重试")
    return sets


def _empirical_enrichment(arm: list[str], gold: set[str], universe: list[str],
                          n_emp: int, seed: int) -> dict:
    """经验 p + 超几何 sanity:宇宙固定为候选池,同输入逐值可复现(A3 纪律)。"""
    from scipy.stats import hypergeom

    gold_u = gold & set(universe)
    m, n_gold, k = len(universe), len(gold_u), len(arm)
    obs = len(set(arm) & gold_u)
    rng = np.random.default_rng(seed)
    draws_obs = np.empty(n_emp, dtype=np.int32)
    pool_idx = np.arange(m)
    gold_idx = set(np.where([g in gold_u for g in universe])[0])
    for i in range(n_emp):
        draw = rng.choice(pool_idx, size=k, replace=False)
        draws_obs[i] = np.fromiter((ix in gold_idx for ix in draw), dtype=np.int32).sum()
    mean_rand = float(draws_obs.mean())
    p_emp = (1 + int((draws_obs >= obs).sum())) / (1 + n_emp)
    p_hyper = float(hypergeom.sf(obs - 1, m, n_gold, k)) if (n_gold and obs > 0) else 1.0
    return {"k": k, "universe_n": m, "gold_n_in_universe": int(n_gold),
            "hits": obs, "mean_random_hits": round(mean_rand, 3),
            "fold": round(obs / mean_rand, 2) if mean_rand > 0 else None,
            "p_empirical": round(p_emp, 4), "n_empirical": n_emp,
            "p_hypergeometric": round(p_hyper, 6)}


def _shortlist_evidence(run: str, sets: dict) -> tuple[pd.DataFrame, dict]:
    """②最小版:最终 shortlist 逐基因证据档案 + 构成统计(无检验主张)。
    循环性是结构性的:平台排序臂上游(core/targets.py)已聚合两库疾病靶点。"""
    sl_path = os.path.join(run, "shortlist.csv")
    sl = pd.read_csv(sl_path, comment="#", encoding="utf-8-sig")
    genes = sl["gene"].astype(str).str.upper().tolist() if "gene" in sl.columns else []
    dg = sets.get("disgenet", {}).get("genes", set())
    ot = sets.get("ot", {}).get("genes", set())
    rows = []
    for g in genes:
        in_sl = bool(sl.loc[sl["gene"].astype(str).str.upper() == g, "in_shortlist"].iloc[0]) \
            if "in_shortlist" in sl.columns else True
        rows.append({
            "gene": g, "in_shortlist": in_sl,
            "in_ot": g in ot,
            "in_disgenet": g in dg,
        })
    ev = pd.DataFrame(rows)
    final = ev[ev["in_shortlist"]] if "in_shortlist" in ev.columns else ev
    comp = {"shortlist_n": int(len(final)),
            "ot_overlap": int(final["in_ot"].sum()) if len(final) else 0,
            "disgenet_overlap": int(final["in_disgenet"].sum()) if len(final) else 0}
    return ev, comp


def main() -> None:
    ap = argparse.ArgumentParser(description="ML shortlist 疾病相关性验证(B2)")
    ap.add_argument("--run", help="ML run 目录(含 ml_consensus_scores.csv/summary.json)")
    ap.add_argument("--disease", help="疾病名(OT/DisGeNET 检索词,建议标准名)")
    ap.add_argument("--k", type=int, default=None,
                    help="主证据臂 K(默认 config.ML_SHORTLIST_TOP_K)")
    ap.add_argument("--k-scan", default="10,20,50", help="K 敏感性扫描档,逗号分隔")
    ap.add_argument("--list", action="store_true", help="列出最近 run 后退出")
    ap.add_argument("--out", default=None, help="输出根目录(默认 <结果根>/disease_relevance)")
    args = ap.parse_args()

    cfg = _config()
    run = _resolve_run(args.run, args.list)
    if not args.disease:
        sys.exit("需要 --disease(疾病标准名,如 'type 2 diabetes mellitus')")
    k = args.k if args.k is not None else cfg.ML_SHORTLIST_TOP_K
    k_scan = [int(x) for x in str(args.k_scan).split(",") if str(x).strip()]

    tag = os.path.basename(run.rstrip("/\\")).rsplit("_", 1)[0]
    print(f"== B2 疾病相关性验证 ==\nrun: {run}\n疾病: {args.disease}")

    consensus, arms_info = _load_universe_and_arms(run, k, k_scan)
    universe = arms_info["universe"]
    print(f"候选宇宙 {len(universe)} 基因 | 验证臂 {sorted(arms_info['arms'])}")

    sets = _gold_sets(args.disease)
    for name, s in sets.items():
        hit = len(s["genes"] & set(universe))
        print(f"金标准[{name}]: {s['n_raw']} 靶点,其中 {hit} 个落在候选宇宙({s['note']})")

    emp_seed = int(cfg.ML_SEED) + _EMP_SEED_OFFSET
    e1_rows = []
    for arm_name, arm in arms_info["arms"].items():
        for db, s in sets.items():
            r = _empirical_enrichment(arm, s["genes"], universe, _N_EMPIRICAL, emp_seed)
            r.update({"arm": arm_name, "gold_db": db, "gold_threshold": s["threshold"]})
            e1_rows.append(r)
            print(f"E1 {arm_name} × {db}: 命中 {r['hits']}/{r['gold_n_in_universe']} "
                  f"(随机均值 {r['mean_random_hits']}) 倍数 {r['fold']} "
                  f"p_emp={r['p_empirical']} p_hg={r['p_hypergeometric']}")

    e2 = None
    if "ot" in sets:
        from core.network_algorithms import proximity_zscore
        try:
            e2 = proximity_zscore(arms_info["arms"][f"ml_top{k}"], sets["ot"]["genes"])
            print(f"E2 网络邻近: {json.dumps({k2: v for k2, v in e2.items() if k2 != 'details'}, ensure_ascii=False)}")
        except Exception as exc:
            e2 = {"ok": False, "error": str(exc)}
            print(f"E2 网络邻近不可用: {exc}(背景互作组需先构建)")

    ev, comp = _shortlist_evidence(run, sets)
    print(f"② shortlist 构成: {comp['shortlist_n']} 基因,OT 重叠 {comp['ot_overlap']},"
          f"DisGeNET 重叠 {comp['disgenet_overlap']}(仅构成说明——上游 targets.py 已注入两库)")

    out_base = args.out or os.path.join(_default_root(), "disease_relevance")
    out_dir = os.path.join(out_base, f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}")
    os.makedirs(out_dir, exist_ok=True)
    ev.to_csv(os.path.join(out_dir, "per_gene_evidence.csv"), index=False, encoding="utf-8-sig")

    summary = {"module": "shortlist-disease-benchmark", "b2_version": 1,
               "run_dir": run, "disease": args.disease, "universe_n": len(universe),
               "gold_sets": {n: {"n_raw": s["n_raw"], "threshold": s["threshold"],
                                 "note": s["note"]} for n, s in sets.items()},
               "e1_ml_arm": e1_rows, "e2_proximity": e2,
               "e1b_shortlist_composition": comp,
               "independence_tiers": {
                   "ml_arm": "非循环:共识臂仅见表达数据与标签,两库均为有效金标准",
                   "full_shortlist": "循环:平台排序臂上游 core/targets.py 已聚合 OT+DisGeNET,"
                                     "重叠仅作构成说明,不构成独立证据"},
               "args": vars(args), "seed_empirical": emp_seed}
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1, default=str)
    _write_report(out_dir, summary)
    print(f"落盘 → {out_dir}(report.md / summary.json / per_gene_evidence.csv)")


def _default_root() -> str:
    from core.ml_pipeline import default_results_root
    return default_results_root()


def _write_report(out_dir: str, s: dict) -> None:
    """人类可读报告:证据独立性声明写死在首部,防止报告被单独转发后误读。"""
    e1 = s["e1_ml_arm"]
    lines = [
        "# 疾病相关性验证报告(B2)",
        f"- run: `{s['run_dir']}` | 疾病: {s['disease']} | 候选宇宙: {s['universe_n']} 基因",
        f"- 证据独立性:①ML 共识臂=**非循环**(仅见表达数据与标签);"
        f"②全 shortlist=**循环**(平台排序臂上游已聚合 OT+DisGeNET),仅构成说明",
        "",
        "## E1 ML 共识臂富集(主证据)",
        "| 臂 | 金标准库 | K | 命中 | 随机均值 | 倍数 | 经验 p | 超几何 p |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in e1:
        lines.append(f"| {r['arm']} | {r['gold_db']} | {r['k']} | {r['hits']} | "
                     f"{r['mean_random_hits']} | {r['fold']} | {r['p_empirical']} | "
                     f"{r['p_hypergeometric']} |")
    lines += ["", f"## E2 网络邻近(ML top{s['args']['k']} vs OT 疾病模块,度匹配置换)"]
    e2 = s["e2_proximity"]
    if e2 and e2.get("ok"):
        keys = [k2 for k2 in e2 if k2 not in ("ok", "details")]
        lines += [f"- {k2}: {e2[k2]}" for k2 in keys]
    else:
        lines.append(f"- 不可用:{(e2 or {}).get('error', '未知原因')}(背景互作组需先构建)")
    c = s["e1b_shortlist_composition"]
    lines += ["", "## ② 全 shortlist 构成(仅说明,无检验主张)",
              f"- {c['shortlist_n']} 基因:OT 重叠 {c['ot_overlap']},DisGeNET 重叠 {c['disgenet_overlap']}",
              "- 逐基因证据见 per_gene_evidence.csv(湿实验交接用档案)",
              "", f"- 经验重抽样 seed: {s['seed_empirical']}(同输入可复现)",
              f"- 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    with open(os.path.join(out_dir, "report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
