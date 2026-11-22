"""协同回收基准:DrugComb 标签聚合 + 防泄漏切分 + 分层报告
(二期方案 §2.3;基准包第二件,清单 §1.6,队友实现)。

定位:benchmark 是尺子——本文件只做"标签聚合 → 药对数据集 → 指标报告",
不改打分算法;预测分数(协同置信)由协同预测管线产出后从外部传入。

定稿口径(全部引方案 §12,不自行拍板):
  - 打分 ZIP,药对×细胞系粒度;ZIP ≥ +10 协同 / ≤ −10 拮抗 / |ZIP|<10 模糊带剔除
    (B6 定稿;±5 平移敏感性已排期,由调用方以 zip_band 传参另跑)
  - 药对级标签 = 多细胞系多数投票;平票或细胞系样本 <3 → 剔除(B6)
  - 投票比例列 = 协同票数/有效票数(12.2-B6:全票协同与多数协同区分可信度)
  - 细胞系→癌种:CCLE 单源主匹配(B5);失败不剔除,挂泛癌对照分层报告
  - 泛癌对照基因集 = 1.9 产物(OT ∪ CancerMine),禁 COSMIC/MSigDB(许可)
  - 防泄漏:调参/测试集按**药对**切分——同一药对的多细胞系记录跨侧即泄漏
    (评测虚高),泄漏闸显式断言(第一审核项)

数据获取适配位:DrugComb 导出格式在 G1 预检核实;DEFAULT_COLUMN_MAP 按门户
常见导出命名给默认映射,可传参覆盖。药名/细胞系名规范化后匹配失败的进
unmatched 报告,不静默丢(陷阱#2/#3)。
"""

import argparse
import json
import logging
import os
import re

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

logger = logging.getLogger(__name__)

BENCHMARK_SCHEMA = "drugcomb_benchmark.v1"

# G1 预检前默认列映射(DrugComb 门户导出常见命名;预检后如不符,调用方传参覆盖)
DEFAULT_COLUMN_MAP = {
    "drug_row": "drug_row",
    "drug_col": "drug_col",
    "cell_line": "cell_line_name",
    "zip": "synergy_zip",
}

_vote_re = re.compile(r"[\s\-_]+")


def normalize_drug_name(name) -> str:
    """药名规范化(陷阱#2):小写化 + 连字符/空格收敛为空 —— aspirin-ish 变体
    归一后再与 ChEMBL/预测表对齐;失败进 unmatched 报告,不静默丢。"""
    s = str(name or "").strip().lower()
    return _vote_re.sub("", s)


def normalize_cell_line(name) -> str:
    """细胞系规范化(陷阱#3):MCF-7 vs MCF7 → 大写去连字符/空格。"""
    s = str(name or "").strip().upper()
    return _vote_re.sub("", s)


def aggregate_pair_labels(df: pd.DataFrame, *, col_map: dict | None = None,
                          zip_band: float = 10.0,
                          min_cell_lines: int = 3,
                          block_agg: str = "median") -> dict:
    """DrugComb 记录 → 药对级标签表(§2.3 定稿口径)。
    步骤:规范化 → (药对,细胞系)粒度聚合(block 多记录取中位数,复核点)→
    模糊带剔除 → 多数投票(平票或 <min_cell_lines 剔除)→ 投票比例列。
    返回 {"labels": DataFrame, "stats": dict(各剔除环节计数,不静默丢)}。
    labels 列:drug_row/drug_col(规范化)/n_cell_lines/n_votes/vote/
    vote_label(协同=1,拮抗=0)/vote_fraction。"""
    cm = dict(DEFAULT_COLUMN_MAP)
    cm.update(col_map or {})
    need = ["drug_row", "drug_col", "cell_line", "zip"]
    # 检查源列名(cm[c])是否在原始表里——归一名是本函数产出的列,拿它查原始表
    # 永远缺失(首跑即被手算测试抓出的检查键错位)
    missing = [c for c in need if cm[c] not in df.columns]
    if missing:
        raise RuntimeError(f"列映射后缺失必需列: {missing}(G1 预检核实导出格式)")
    d = df.copy()
    d["drug_row"] = d[cm["drug_row"]].map(normalize_drug_name)
    d["drug_col"] = d[cm["drug_col"]].map(normalize_drug_name)
    d["cell_line"] = d[cm["cell_line"]].map(normalize_cell_line)
    d["zip"] = pd.to_numeric(d[cm["zip"]], errors="coerce")
    stats = {"rows_in": int(len(d)),
             "rows_invalid_zip": int(d["zip"].isna().sum()),
             "rows_fuzzy_band": 0, "pairs_dropped_small": 0,
             "pairs_dropped_tie": 0, "pairs_kept": 0}
    d = d.dropna(subset=["zip"])
    # 药对 × 细胞系粒度:同格多 block 取中位数(稳健;复核点)
    grain = (d.groupby(["drug_row", "drug_col", "cell_line"], as_index=False)
              .agg(zip=("zip", block_agg)))
    # 模糊带剔除(|ZIP|<zip_band)——投票前剔除,不进票
    fuzzy = grain[grain["zip"].abs() < float(zip_band)]
    stats["rows_fuzzy_band"] = int(len(fuzzy))
    grain = grain[grain["zip"].abs() >= float(zip_band)]
    grain["vote"] = np.where(grain["zip"] >= float(zip_band), 1, 0)
    rows = []
    for (row, col), g in grain.groupby(["drug_row", "drug_col"]):
        n = int(len(g))
        n_syn = int(g["vote"].sum())
        if n < int(min_cell_lines):
            stats["pairs_dropped_small"] += 1
            continue
        label = 1 if n_syn * 2 > n else 0
        if n_syn * 2 == n:
            stats["pairs_dropped_tie"] += 1
            continue  # 平票剔除(B6)
        rows.append({"drug_row": row, "drug_col": col, "n_cell_lines": n,
                     "n_votes": n, "n_syn": n_syn, "vote": label,
                     "vote_label": label,
                     "vote_fraction": n_syn / n})
    labels = pd.DataFrame(rows, columns=[
        "drug_row", "drug_col", "n_cell_lines", "n_votes", "n_syn",
        "vote", "vote_label", "vote_fraction"])
    stats["pairs_kept"] = int(len(labels))
    logger.info("DrugComb 标签聚合: %s", stats)
    return {"labels": labels, "stats": stats}


def assert_pair_split_clean(train_pairs: set, test_pairs: set) -> None:
    """防泄漏闸(陷阱#1,第一审核项):同一药对绝不允许跨切分两侧。
    独立成函数以便测试验证闸会响(违规输入 → AssertionError)。"""
    leak = set(train_pairs) & set(test_pairs)
    assert not leak, f"药对切分泄漏: {sorted(leak)[:5]} 同药对跨调参/测试两侧"


def split_pairs(labels: pd.DataFrame, *, test_frac: float = 0.3,
                seed: int = 0) -> dict:
    """药对粒度切分(§2.3:固定种子,切分清单整体进参数快照)。
    返回 {"train": labels 子集, "test": labels 子集,
          "pairs": {"train": [...], "test": [...]}};泄漏闸显式断言。"""
    pairs = (labels[["drug_row", "drug_col"]]
             .drop_duplicates().apply(tuple, axis=1).tolist())
    pairs = sorted(pairs)
    rng = np.random.default_rng(seed)
    rng.shuffle(pairs)
    k_test = max(1, int(round(len(pairs) * float(test_frac))))
    test_pairs = set(pairs[:k_test])
    train_pairs = set(pairs) - test_pairs
    assert_pair_split_clean(train_pairs, test_pairs)  # 泄漏闸(陷阱#1)
    key = labels.apply(lambda r: (r["drug_row"], r["drug_col"]), axis=1)
    train = labels[key.isin(train_pairs)].reset_index(drop=True)
    test = labels[key.isin(test_pairs)].reset_index(drop=True)
    return {"train": train, "test": test,
            "pairs": {"train": sorted(train_pairs), "test": sorted(test_pairs)}}


def map_lineage(cell_lines, ccle: pd.DataFrame, *,
                ccle_cell_col: str = "cell_line",
                ccle_lineage_col: str = "lineage",
                pancancer_label: str = "pancancer") -> dict:
    """CCLE 谱系映射(B5:单源主匹配;失败挂泛癌对照,不剔除)。
    ccle 表规范化后精确匹配;返回 {"lineage": {cell_line: 谱系或 pancancer},
    "unmatched": [...]}——unmatched 由调用方写报告(不静默丢,陷阱#3)。"""
    table = {normalize_cell_line(c): str(l)
             for c, l in zip(ccle[ccle_cell_col], ccle[ccle_lineage_col])}
    lineage, unmatched = {}, []
    for c in cell_lines:
        norm = normalize_cell_line(c)
        if norm in table:
            lineage[str(c)] = table[norm]
        else:
            lineage[str(c)] = pancancer_label  # B5:失败挂泛癌,不剔除
            unmatched.append(str(c))
    return {"lineage": lineage, "unmatched": unmatched}


def topk_recall(scores: pd.DataFrame, *, score_col: str, label_col: str,
                k: int) -> float:
    """Top-K 回收率:按预测分降序取前 K,其中协同(标签 1)占比
    / 全体协同数。K 超过行数时按行数计。"""
    d = scores.sort_values(score_col, ascending=False).head(int(k))
    positives = int((scores[label_col] == 1).sum())
    if positives == 0:
        raise RuntimeError("标签无协同正例,Top-K 回收无定义(拒绝)")
    return float((d[label_col] == 1).sum() / positives)


def stratified_report(labels_with_pred: pd.DataFrame, *,
                      score_col: str, label_col: str = "vote_label",
                      geometry_col: str | None = None,
                      lineage_col: str | None = None,
                      top_k: int | None = None) -> pd.DataFrame:
    """分层报告(§2.3):AUROC(sklearn,不手写)+ Top-K 回收;分层维度 =
    (重叠/互补,若提供 geometry 列)×(映射成功/泛癌,若提供 lineage 列)。
    返回 DataFrame[stratum, n_pairs, auroc, topk_recall(可空)]。"""
    rows = []
    df = labels_with_pred.copy()
    strata = {"overall": np.ones(len(df), dtype=bool)}
    if geometry_col and geometry_col in df.columns:
        strata["overlap"] = df[geometry_col] < 0
        strata["separated"] = df[geometry_col] >= 0
    if lineage_col and lineage_col in df.columns:
        strata["lineage_mapped"] = df[lineage_col] != "pancancer"
        strata["pancancer"] = df[lineage_col] == "pancancer"
    for name, mask in strata.items():
        sub = df[mask]
        y = sub[label_col]
        if y.nunique() < 2:
            rows.append({"stratum": name, "n_pairs": int(len(sub)),
                         "auroc": None, "topk_recall": None})
            continue
        auc = float(roc_auc_score(y, sub[score_col]))
        tk = (topk_recall(sub, score_col=score_col, label_col=label_col,
                          k=top_k) if top_k else None)
        rows.append({"stratum": name, "n_pairs": int(len(sub)),
                     "auroc": auc, "topk_recall": tk})
    return pd.DataFrame(rows)


def load_table(path: str, col_map: dict | None = None) -> pd.DataFrame:
    """DrugComb 导出文件加载(G1 预检后回填格式;CSV/TSV 按扩展名)。"""
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, encoding="utf-8-sig")
    if col_map:
        df = df.rename(columns={v: k for k, v in col_map.items()})
    return df


def save_results(res: dict, out_dir: str) -> dict:
    """结果落盘:summary.json(切分清单/统计进参数快照)+ labels.csv +
    unmatched.csv(药名/细胞系匹配失败报告,不静默丢)。"""
    from WY.core.caching import save_csv_cache, save_json_file

    os.makedirs(out_dir, exist_ok=True)
    summary = {k: v for k, v in res.items() if not hasattr(v, "to_dict")}
    save_json_file(os.path.join(out_dir, "summary.json"), summary, indent=1)
    save_csv_cache(res["labels"], os.path.join(out_dir, "labels.csv"),
                   skip_empty=False)
    unmatched = res.get("unmatched")
    if unmatched is not None and not unmatched.empty:
        save_csv_cache(unmatched, os.path.join(out_dir, "unmatched.csv"),
                       skip_empty=False)
    return {"summary": os.path.join(out_dir, "summary.json"),
            "labels": os.path.join(out_dir, "labels.csv"),
            "unmatched": os.path.join(out_dir, "unmatched.csv")}


def main(argv=None) -> int:
    """CLI 入口(退出码约定:全成 0/任一败 1)。输入 = DrugComb 导出文件 +
    预测分数表(协同预测管线产出);报告 = 分层 AUROC/Top-K。"""
    ap = argparse.ArgumentParser(description="DrugComb 协同回收基准(方案 §2.3)")
    ap.add_argument("--labels", required=True, help="DrugComb 导出文件(CSV/TSV)")
    ap.add_argument("--predictions", required=True,
                    help="药对预测分数表:drug_row,drug_col,<score_col>[,geometry]")
    ap.add_argument("--score-col", default="confidence")
    ap.add_argument("--zip-band", type=float, default=10.0)
    ap.add_argument("--min-cell-lines", type=int, default=3)
    ap.add_argument("--test-frac", type=float, default=0.3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--top-k", type=int, default=None)
    ap.add_argument("--out", required=True, help="结果目录")
    args = ap.parse_args(argv)
    try:
        raw = load_table(args.labels)
        agg = aggregate_pair_labels(raw, zip_band=args.zip_band,
                                    min_cell_lines=args.min_cell_lines)
        split = split_pairs(agg["labels"], test_frac=args.test_frac,
                            seed=args.seed)
        preds = load_table(args.predictions)
        preds["drug_row"] = preds["drug_row"].map(normalize_drug_name)
        preds["drug_col"] = preds["drug_col"].map(normalize_drug_name)
        merged = split["test"].merge(preds, on=["drug_row", "drug_col"],
                                     how="inner")
        if merged.empty:
            raise RuntimeError("预测表与测试药对无交集(检查药名规范化)")
        report = stratified_report(
            merged, score_col=args.score_col,
            geometry_col="geometry" if "geometry" in merged.columns else None,
            top_k=args.top_k)
        res = {"schema": BENCHMARK_SCHEMA, "stats": agg["stats"],
               "pairs": split["pairs"], "report": report,
               "labels": agg["labels"]}
        paths = save_results(res, args.out)
        print(json.dumps({"paths": paths,
                          "report": report.to_dict("records")},
                         ensure_ascii=False, default=str))
        return 0
    except Exception as exc:  # noqa: BLE001 —— CLI 边界统一转退出码
        print(f"benchmark failed: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
