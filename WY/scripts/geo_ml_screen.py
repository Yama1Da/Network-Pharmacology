"""GEO 表达谱与 ML 共识筛选的命令行入口(实施方案 §6,仿 lit_pcba_benchmark.py 先例)。

M1:geo(数据集检索)/ deg(下载→解析→分组→预处理→DEG)两个子命令;
M2:screen 子命令(在 deg 产物之上叠加六算法共识与 shortlist)。

用法(在 WY 目录下):
    python scripts/geo_ml_screen.py geo --disease "type 2 diabetes" --limit 20
    python scripts/geo_ml_screen.py deg --gse GSE63067 --case-pattern "steatosis" --control-pattern "healthy"
    python scripts/geo_ml_screen.py deg --gse GSE63067            # 不带分组:仅打印分组字段总览
    python scripts/geo_ml_screen.py screen --gse GSE63067 --case-pattern "steatosis" --control-pattern "healthy"
        --ranked-csv "path/to/rwr_output.csv" --top-k 20
    python scripts/geo_ml_screen.py screen --expr expr.csv --labels labels.csv --case-label disease
        # 自备数据直跑(方案 §6,不经 GEO):expr.csv 首列为基因符号,样本轴自动识别;
        # labels 两列(样本ID, 分组),分组写法宽松(disease/control、case/ctrl、1/0)
"""

import argparse
import json
import os
import sys
import time
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

import config
from core.expression import differential_expression, map_probe_to_symbol, preprocess
from core.geo import (
    assign_groups,
    fetch_geo_dataset,
    group_fields_overview,
    search_geo_datasets,
)


def _results_dir(args_out: str | None, tag: str) -> str:
    """结果目录优先级:--out > API.env 的 CLI_RESULTS_DIR > 当前目录 results。"""
    base = args_out or config.CLI_RESULTS_DIR or os.path.join(os.getcwd(), "results")
    out = os.path.join(base, "mlscreen", f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}")
    os.makedirs(out, exist_ok=True)
    return out


def cmd_geo(args) -> int:
    df = search_geo_datasets(args.disease, organism=args.organism, limit=args.limit)
    if df.empty:
        print(f"未检索到匹配数据集(疾病={args.disease!r},物种={args.organism!r})。"
              "可尝试换同义词/英文全名,或放宽检索词。")
        return 0
    with pd.option_context("display.max_colwidth", 60, "display.width", 200):
        print(df.to_string(index=False))
    print("\n提示:请核对「标题」列确认数据集内容;possible_cellline=表示可能为细胞系实验"
          "(组织与细胞模型的表达证据迁移性不同)。选定后用 deg 子命令跑该 GSE。")
    return 0


def cmd_deg(args) -> int:
    print(f"[1/4] 获取 {args.gse} …")
    ds = fetch_geo_dataset(args.gse, force_refresh=args.force_refresh)
    if ds.issues:
        raise RuntimeError(f"{args.gse}: " + ";".join(ds.issues)
                           + " —— 请换数据集(元数据型提交无表达矩阵)")
    print(f"    样本 {len(ds.samples)} | 矩阵 {ds.table.shape[0]}x{ds.table.shape[1]} "
          f"| 平台 {ds.gpl}")

    overview = group_fields_overview(ds)
    if not args.case_pattern or not args.control_pattern:
        print("\n[分组总览] 未提供 --case-pattern/--control-pattern,请按以下字段取值确定分组:")
        with pd.option_context("display.max_colwidth", 90, "display.width", 200):
            print(overview.to_string(index=False))
        print("\n正则匹配对象 = 样本注释全文(source_name + title + characteristics 各字段拼接)。")
        return 0

    print("[2/4] 分组 …")
    labels = assign_groups(ds, args.case_pattern, args.control_pattern)
    print(f"    case={int((labels == 1).sum())} control={int((labels == 0).sum())} "
          f"丢弃={len(ds.samples) - len(labels)}")

    print("[3/4] 预处理与探针映射 …")
    matrix, log_branch = preprocess(ds.table, force_log=args.force_log)
    matrix = map_probe_to_symbol(matrix, ds.gpl)
    print(f"    log2 决策: {log_branch} | 映射后基因数 {len(matrix)}")

    print("[4/4] 差异表达(Welch t + BH)…")
    deg = differential_expression(matrix, labels)
    out_dir = _results_dir(args.out, args.gse)
    deg_path = os.path.join(out_dir, "deg_results.csv")
    deg.to_csv(deg_path, index=False, encoding="utf-8-sig")
    sig = int((deg["FDR"] < 0.05).sum())
    summary = {
        "module": "geo-deg", "gse": ds.gse, "gpl": ds.gpl,
        "n_samples": int(len(ds.samples)), "n_genes": int(len(matrix)),
        "groups": {"case": int((labels == 1).sum()), "control": int((labels == 0).sum()),
                   "case_pattern": args.case_pattern, "control_pattern": args.control_pattern},
        "preprocess": log_branch,
        "deg_significant_fdr05": sig,
        "args": {k: v for k, v in vars(args).items() if k != "func"},
    }
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)
    print(f"    显著基因(FDR<0.05): {sig} | 结果: {deg_path}")
    print(f"    summary: {os.path.join(out_dir, 'summary.json')}")
    return 0


def cmd_screen(args) -> int:
    """M2 全链:获取→分组→预处理→(映射)→DEG→预过滤→六算法共识→shortlist→落盘。

    流程本体在 core/ml_pipeline.run_screen_chain(CLI 与 M4 页面共用,方案 §M4
    "页面不新增能力、只换交互形态");本函数只做参数收集与进度打印。
    """
    from core.ml_pipeline import STAGES, run_screen_chain

    # 抑制小样本 LinearSVC 偶发收敛警告(不影响结果,代码内已处理)
    warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

    ranked_csvs = {}
    for rc in (args.ranked_csv or []):
        name = os.path.splitext(os.path.basename(rc))[0]
        ranked_csvs[name] = rc

    base = args.out or config.CLI_RESULTS_DIR or os.path.join(os.getcwd(), "results")
    step = {"i": 0}

    def _progress(stage, message):
        if stage is None:
            print(f"    {message}")
        else:
            step["i"] += 1
            print(f"[{step['i']}/{len(STAGES)}] {message}")

    res = run_screen_chain(
        gse=args.gse, case_pattern=args.case_pattern,
        control_pattern=args.control_pattern,
        expr=args.expr, labels=args.labels, case_label=args.case_label,
        force_refresh=args.force_refresh, force_log=args.force_log,
        prefilter_top=args.prefilter_top, bootstrap_rounds=args.bootstrap_rounds,
        ranked_csvs=ranked_csvs, top_k=args.top_k, out_base=base,
        progress=_progress)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="GEO 表达谱与 ML 共识筛选(M1:geo/deg;M2 增加 screen)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_geo = sub.add_parser("geo", help="按疾病检索 GEO 表达谱数据集")
    p_geo.add_argument("--disease", required=True, help="疾病名(英文检索效果更好)")
    p_geo.add_argument("--organism", default="Homo sapiens")
    p_geo.add_argument("--limit", type=int, default=20)
    p_geo.set_defaults(func=cmd_geo)

    p_deg = sub.add_parser("deg", help="下载→解析→分组→预处理→DEG")
    p_deg.add_argument("--gse", required=True, help="GSE 编号")
    p_deg.add_argument("--case-pattern", help="病例组正则(不填则仅打印分组字段总览)")
    p_deg.add_argument("--control-pattern", help="对照组正则")
    p_deg.add_argument("--force-log", dest="force_log", default=None, type=str,
                       choices=["true", "false"], help="覆盖 log2 启发式(true/false)")
    p_deg.add_argument("--force-refresh", action="store_true", help="跳过缓存重新下载")
    p_deg.add_argument("--out", help="结果根目录(默认 CLI_RESULTS_DIR 或 ./results)")
    p_deg.set_defaults(func=cmd_deg)

    p_scr = sub.add_parser("screen", help="全链:DEG + 六算法共识 + shortlist(M2;--gse 或 --expr 二选一)")
    p_scr.add_argument("--gse", default=None, help="GSE 编号(与 --expr 二选一)")
    p_scr.add_argument("--expr", help="本地表达矩阵 CSV(首列=基因符号,样本轴自动识别;与 --gse 二选一)")
    p_scr.add_argument("--labels", help="两列 CSV(样本ID, 分组);--expr 模式必填")
    p_scr.add_argument("--case-label", default=None,
                       help="--expr 模式:分组列中代表病例组的值"
                            "(缺省按 1/case/disease/treated/patient 推断,推断不了报错)")
    p_scr.add_argument("--case-pattern", default=None, help="病例组正则(--gse 模式必填)")
    p_scr.add_argument("--control-pattern", default=None, help="对照组正则(--gse 模式必填)")
    p_scr.add_argument("--ranked-csv", action="append",
                       help="平台排序 CSV(可多次传入,如 --ranked-csv rwr.csv)")
    p_scr.add_argument("--top-k", type=int, default=None,
                       help="shortlist 截断(默认 config.ML_SHORTLIST_TOP_K=20)")
    p_scr.add_argument("--prefilter-top", type=int, default=None,
                       help="DEG 预过滤基因数(默认 config.ML_DEG_PREFILTER_TOP_N=2000)")
    p_scr.add_argument("--bootstrap-rounds", type=int, default=None,
                       help="bootstrap 轮次(默认 config.ML_BOOTSTRAP_ROUNDS=100;冒烟可 10)")
    p_scr.add_argument("--force-log", dest="force_log", default=None, type=str,
                       choices=["true", "false"], help="覆盖 log2 启发式")
    p_scr.add_argument("--force-refresh", action="store_true")
    p_scr.add_argument("--out", help="结果根目录")
    p_scr.set_defaults(func=cmd_screen)

    args = ap.parse_args()
    if getattr(args, "force_log", None) is not None:
        args.force_log = {"true": True, "false": False}[args.force_log]
    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("\n[中断] 本次运行未写入结果目录(产物只在全链成功后一次性落盘);"
              "已下载的 GEO 缓存文件不受影响,重跑同一命令自动复用")
        return 130
    except Exception as exc:
        print(f"[错误] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
