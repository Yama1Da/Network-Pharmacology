"""疾病模块回测(方案 §4.3;基准包第三件,清单 §1.8,队友实现)。

协议(定死,可复现):DisGeNET curated 关联按证据分分层(高分池=留一法对象,
低分池=泛化检验);每轮从 T_D 随机隐藏 20% 已知疾病基因 → 三方法+共识各自
扩张 → Recall@N(被隐藏基因被捞回的比例),重复 100 次取均值±标准差;
分层报告(高分/低分)。预期结论:共识的 Recall@N 优于任何单方法——若不
成立,模块 C 降级为"多算法对照展示"(方案 §8 风险表),本表即裁决依据。

陷阱钉死(清单 §1.8,实现处再各钉一次):
  1) 隐藏池基因绝不能出现在种子集——assert_no_leak 显式断言(防泄漏
     最后一道闸;tests 用违规输入验证闸会响)
  2) 100 次重复的隐藏抽样用 np.random.SeedSequence.spawn 派生独立流——
     全部同种子会产生完全相同的隐藏集(样本量虚而不实)
  3) Recall@N 的 N 与扩张规模 N 是同一个值(下方协议内注释钉死)
  4) DisGeNET 拉取走一期既有缓存与限流(fetch_disease_targets_disgenet,
     不绕开 http_client)

分层口径(负责人复核点):方案未给数值阈值,默认按正分数中位数切分
(高分池 = ≥中位数);分数缺失/≤0 的基因不入池(宁缺毋滥);阈值可显式
传参,敏感性另跑。审计清单(audit)即切分清单,整体进参数快照(§5.1)。
"""

import argparse
import json
import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

BENCHMARK_SCHEMA = "disease_module_benchmark.v1"


def assert_no_leak(hidden: set, seeds, repeat_no: int) -> None:
    """防泄漏最后一道闸(陷阱#1):隐藏池基因绝不能出现在种子集。
    独立成函数以便测试直接验证闸会响(违规输入 → AssertionError)。"""
    assert hidden.isdisjoint(set(seeds)), (
        f"第 {repeat_no} 轮泄漏:隐藏池与种子集相交 "
        f"{sorted(hidden & set(seeds))[:5]}(方案 §4.3 防泄漏,禁止绕过)")


def _split_hidden(pool: list[str], rng, hide_frac: float) -> set:
    """从候选隐藏池随机抽 hide_frac 比例(至少 1 个,至多全池)。"""
    k = min(len(pool), max(1, int(round(len(pool) * float(hide_frac)))))
    idx = rng.choice(len(pool), size=k, replace=False)
    return {pool[i] for i in idx}


def stratify_pool(assoc: pd.DataFrame, *, gene_col: str = "approvedSymbol",
                  score_col: str = "disgenet_score",
                  threshold=None) -> dict:
    """高分/低分池切分:threshold=None 按正分数中位数切(默认口径);
    分数缺失/≤0 不入池(宁缺毋滥)。返回 {"high": set, "low": set,
    "threshold": 实际生效阈值}。"""
    scores = pd.to_numeric(assoc[score_col], errors="coerce")
    df = assoc.assign(_s=scores).dropna(subset=["_s"])
    df = df[df["_s"] > 0]
    if df.empty:
        raise RuntimeError("DisGeNET 关联无有效正分数,无法分层(宁缺毋滥)")
    thr = float(threshold) if threshold is not None else float(df["_s"].median())
    high = set(df.loc[df["_s"] >= thr, gene_col].astype(str))
    low = set(df.loc[df["_s"] < thr, gene_col].astype(str))
    return {"high": high, "low": low, "threshold": thr}


def _consensus_ranking(module_ranks: dict[str, dict[str, int]]) -> pd.DataFrame:
    """三方法秩表 → Borda 排名(与 disease_modules 同一秩语义:秩从 1;
    排序 consensus_score 升序 → n_methods 降序 → 基因序,可复现红线)。"""
    union = sorted(set().union(*[set(r) for r in module_ranks.values()]))
    rows = []
    for gene in union:
        present = [ranks[gene] for ranks in module_ranks.values() if gene in ranks]
        rows.append({"gene": gene, "n_methods": len(present),
                     "consensus_score": float(sum(present) / len(present))})
    return pd.DataFrame(rows).sort_values(
        ["consensus_score", "n_methods", "gene"],
        ascending=[True, False, True]).reset_index(drop=True)


def leave_one_out_protocol(ctx, assoc: pd.DataFrame, t_d, *,
                           n_repeats: int = 100, hide_frac: float = 0.2,
                           n=None, base_seed: int = 0,
                           score_threshold=None) -> dict:
    """留一法协议(§4.3,纯函数,ctx 注入)。返回:
      table: method × stratum(overall/high/low)的 Recall@N 均值±标准差;
      consensus_wins: 共识 overall ≥ 每个单方法 overall(§4.3 预期,裁决列);
      audit: 每轮隐藏清单(切分清单,进参数快照 §5.1)。
    N 语义(陷阱#3):Recall@N 的 N = 扩张规模 N,同一个值,共识模块同取
    前 N 个(Borda 排名),三方法+共识口径完全可比。"""
    from WY.core import disease_modules as dmo

    n_target = dmo._module_size(n)
    _nodes, idx, _Mt, _degree_map = ctx
    pool = sorted({str(g).strip().upper() for g in t_d if str(g).strip().upper() in idx})
    out_of_lcc = sorted({str(g).strip().upper() for g in t_d
                         if str(g).strip().upper() not in idx})
    if len(pool) < 3:
        raise RuntimeError(f"T_D∩LCC 仅 {len(pool)} 基因,留一法无统计意义(拒绝)")
    strata = stratify_pool(assoc, threshold=score_threshold)
    high, low = strata["high"], strata["low"]

    # 陷阱#2:SeedSequence.spawn 派生独立流——每轮独立,绝不同种子
    children = np.random.SeedSequence(base_seed).spawn(int(n_repeats))
    expanders = {
        "rwr": lambda ss: dmo.rwr_module(ctx, ss, n_target),
        "diamond": lambda ss: dmo.diamond_expand(ctx, ss, n_target),
        "neighbor": lambda ss: dmo.neighbor_expansion(ctx, ss, n_target),
    }
    records: list[dict] = []
    audit: list[dict] = []
    for i, child in enumerate(children):
        rng = np.random.default_rng(child)
        hidden = _split_hidden(pool, rng, hide_frac)
        seeds = [g for g in pool if g not in hidden]
        # 防泄漏最后一道闸(陷阱#1):隐藏池与种子集绝不能相交
        assert_no_leak(hidden, seeds, i)
        modules, joins = {}, {}
        for m, expand in expanders.items():
            res = expand(seeds)
            modules[m], joins[m] = res["module"], res["join_order"]
        ranking = _consensus_ranking(
            {m: dmo._method_ranks(modules[m], joins[m]) for m in modules})
        modules["consensus"] = list(ranking["gene"].head(n_target))  # 同一 N(陷阱#3)
        hidden_set = set(hidden)
        rec = {}
        for m, mod in modules.items():
            mod_set = set(mod)
            n_high = len(hidden_set & high)
            n_low = len(hidden_set & low)
            rec[m] = {
                "recall": len(mod_set & hidden_set) / len(hidden_set),
                "recall_high": (len(mod_set & hidden_set & high) / n_high
                                if n_high else None),
                "recall_low": (len(mod_set & hidden_set & low) / n_low
                               if n_low else None),
            }
        records.append(rec)
        audit.append({"repeat": i, "hidden": sorted(hidden_set),
                      "n_seeds": len(seeds)})

    table_rows = []
    for m in records[0]:
        for key, label in (("recall", "overall"), ("recall_high", "high"),
                           ("recall_low", "low")):
            vals = [r[m][key] for r in records if r[m][key] is not None]
            if not vals:
                continue  # 该分层本轮无隐藏基因(如实缺报,不补造)
            table_rows.append({"method": m, "stratum": label,
                               "mean": float(np.mean(vals)),
                               "std": float(np.std(vals)),
                               "n_repeats": len(vals)})
    table = pd.DataFrame(table_rows)
    overall = {r["method"]: r["mean"] for r in table_rows
               if r["stratum"] == "overall"}
    singles = {m: v for m, v in overall.items() if m != "consensus"}
    consensus_wins = (bool(all(overall["consensus"] >= v
                               for v in singles.values()))
                      if "consensus" in overall and singles else None)
    logger_note = (None if consensus_wins is None or consensus_wins
                   else "共识 Recall@N 未优于单方法——按方案 §4.3 回 §8 风险预案,"
                        "模块 C 降级为多算法对照展示")
    return {"schema": BENCHMARK_SCHEMA,
            "n": n_target, "n_repeats": int(n_repeats),
            "hide_frac": float(hide_frac), "base_seed": int(base_seed),
            "pool_size": len(pool), "out_of_lcc": out_of_lcc,
            "strata": {"high": sorted(high), "low": sorted(low),
                       "threshold": strata["threshold"]},
            "table": table, "consensus_wins": consensus_wins,
            "expectation_note": logger_note, "audit": audit}


# ---- 数据适配(一期既有缓存与限流,陷阱#4) ----

def load_associations(disease_name: str) -> pd.DataFrame:
    """DisGeNET curated 关联:直连一期 fetch(缓存 + 429 退避限流),
    不绕开 http_client(陷阱#4)。"""
    from WY.core.disgenet import fetch_disease_targets_disgenet

    df = fetch_disease_targets_disgenet(disease_name)
    if df is None or df.empty:
        raise RuntimeError(f"DisGeNET 无关联数据: {disease_name}"
                           "(检查 DisGeNET_API_KEY 或先在一期抓取缓存)")
    return df


def load_t_d(disease_name: str) -> list[str]:
    """T_D = 一期疾病靶点集(四源整合产物,merged 缓存)。缺失响亮报错——
    回测不得拿 DisGeNET 单源顶替 T_D(四源整合是一期语义)。"""
    from WY import config
    from WY.core.caching import open_named_csv_cache

    _path, cached = open_named_csv_cache(
        config.MERGED_CACHE_DIR, config.MERGED_DISEASE_PREFIX, disease_name,
        label=f"merged disease targets '{disease_name}'")
    if cached is None or cached.empty:
        raise RuntimeError(f"一期疾病靶点集(merged)缺失: {disease_name}"
                           "——先在一期分析中抓取该疾病")
    col = "approvedSymbol" if "approvedSymbol" in cached.columns else cached.columns[0]
    return sorted({str(g).strip().upper() for g in cached[col].dropna()})


def run_benchmark(disease_name: str, *, species: str = "Human",
                  n_repeats: int = 100, hide_frac: float = 0.2, n=None,
                  base_seed: int = 0, score_threshold=None) -> dict:
    """基准入口:取数(一期缓存)→ LCC ctx → 协议 → 结果(附疾病名)。"""
    from WY.core.network_algorithms import _lcc_context

    assoc = load_associations(disease_name)
    t_d = load_t_d(disease_name)
    ctx = _lcc_context(species)
    res = leave_one_out_protocol(ctx, assoc, t_d, n_repeats=n_repeats,
                                 hide_frac=hide_frac, n=n, base_seed=base_seed,
                                 score_threshold=score_threshold)
    res["disease"] = disease_name
    return res


def save_results(res: dict, out_dir: str) -> dict:
    """结果落盘:summary.json(参数快照口径)+ recall_table.csv + split_audit.csv
    (切分清单,§5.1)。返回实际路径 dict。"""
    from WY.core.caching import save_csv_cache, save_json_file

    os.makedirs(out_dir, exist_ok=True)
    summary = {k: v for k, v in res.items()
               if k not in ("table", "audit")}
    summary["table"] = res["table"].to_dict("records")
    save_json_file(os.path.join(out_dir, "summary.json"), summary, indent=1)
    save_csv_cache(res["table"], os.path.join(out_dir, "recall_table.csv"),
                   skip_empty=False)
    audit_df = pd.DataFrame(res["audit"])
    save_csv_cache(audit_df, os.path.join(out_dir, "split_audit.csv"),
                   skip_empty=False)
    return {"summary": os.path.join(out_dir, "summary.json"),
            "table": os.path.join(out_dir, "recall_table.csv"),
            "audit": os.path.join(out_dir, "split_audit.csv")}


def main(argv=None) -> int:
    """CLI 入口(benchmark 脚本,与 lit_pcba 同风格;退出码约定:全成 0/任一败 1)。"""
    ap = argparse.ArgumentParser(description="疾病模块留一法回测(方案 §4.3)")
    ap.add_argument("--disease", required=True, help="疾病名(英文,一期命名)")
    ap.add_argument("--repeats", type=int, default=100)
    ap.add_argument("--hide-frac", type=float, default=0.2)
    ap.add_argument("--n", type=int, default=None, help="扩张规模 N(默认 config)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--score-threshold", type=float, default=None)
    ap.add_argument("--out", required=True, help="结果目录")
    args = ap.parse_args(argv)
    try:
        res = run_benchmark(args.disease, n_repeats=args.repeats,
                            hide_frac=args.hide_frac, n=args.n,
                            base_seed=args.seed,
                            score_threshold=args.score_threshold)
        paths = save_results(res, args.out)
        print(json.dumps({"disease": args.disease,
                          "consensus_wins": res["consensus_wins"],
                          "paths": paths}, ensure_ascii=False))
        return 0
    except Exception as exc:  # noqa: BLE001 —— CLI 边界统一转退出码
        print(f"benchmark failed: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
