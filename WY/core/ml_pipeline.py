"""GEO/ML screen 全链编排:CLI(geo_ml_screen.py)与 M4 页面共用的同一套流程。

方案 §M4「本页面不新增能力、只换交互形态」的落地方式:把 获取→分组→预处理→
(探针映射)→DEG→预过滤→六算法共识→shortlist→落盘 整条链收拢在此处,
CLI 与页面各自只做参数收集与结果呈现,防止两处各改各的造成漂移。

- 数据源二选一:gse(GEO 下载链,需 case/control 正则)或 expr+labels(本地矩阵,
  方案 §6,可选 case_label);硬边界与 core/geo.assign_groups 同款;
- 产物一次性落盘 <结果根>/mlscreen/<tag>_<时间戳>/,summary.json 含算法版本与
  排序 CSV 指纹(mtime+size),供历史结果列表做陈旧检测(方案失效场景③);
- progress(stage, message) 回调:stage 为七阶段键之一(见 _STAGES),传 None
  表示不计编号的补充信息行;CLI 打印全部行,页面按 stage 映射 i18n 文案、
  忽略补充行(message 为中文,页面直显会与界面语言冲突)。
"""

import json
import os
import time

import pandas as pd

# 全链七阶段键(progress 回调的 stage 取值;CLI 进度编号分母也取自它)
STAGES = ("fetch", "group", "preprocess", "deg", "ml", "shortlist", "save")

# shortlist.csv 首行定位注释(方案 §11"表头带定位注释行"):CSV 会单独流传,
# 离开 docstring 与 summary.json 后须自带定位声明,防被误读为交叉验证/预测结论。
# 行首 "#" 是刻意选择:pd.read_csv(comment="#") 可干净跳过(见 selfcheck 回归),
# Excel 打开仅多一行说明,不会混入数据
SHORTLIST_POSITIONING_NOTE = (
    "定位:湿实验优先级证据(多算法共识特征筛选),非交叉验证结论、非预测金标准;"
    "生成参数与算法版本见同目录 summary.json;机器读取本表请用 skiprows=1 或 comment='#'"
)


def default_results_root() -> str:
    """结果根目录:CLI_RESULTS_DIR > <项目 WY>/results。

    用 BASE_DIR 而非 cwd 兜底:Streamlit 服务的 cwd 不保证是 WY,锚到项目目录
    后,CLI(文档约定在 WY 下运行)与页面产物才落在同一处,历史列表才能同源。
    """
    import config
    return config.CLI_RESULTS_DIR or os.path.join(config.BASE_DIR, "results")


# ---------------- 本地矩阵输入(--expr 模式,方案 §6) ----------------


def load_local_expr(expr_path: str, labels_path: str,
                    case_label: str | None = None,
                    progress=None) -> tuple[pd.DataFrame, pd.Series, str, dict]:
    """读本地表达矩阵与两列标签,产出与 GEO 链同构的(矩阵, 标签)。

    labels 为两列 CSV(sample_id, group),分组写法宽松(disease/control、
    case/ctrl、1/0),case_label 指定病例组取值(缺省按常见写法推断,推断不了
    报错);case 匹配 0 个或全部样本硬报错、每组样本量下限,与 assign_groups
    同款边界。表达矩阵样本轴自动识别:列名或行索引与标签样本ID重合多者视为
    样本轴,「行=基因」与「行=样本」两种布局都接受;行标识须为基因符号级
    (本地模式无 GPL 注释表,探针级无法映射 —— 方案 §2"映射失败明确报错")。

    返回 (矩阵, 标签, tag, info);info 记录 case/control 取值与样本数,供 summary。
    """
    from .expression import looks_symbol_level

    def _note(message: str) -> None:
        if progress:
            progress(None, message)

    expr_df = pd.read_csv(expr_path, index_col=0, encoding="utf-8-sig")
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    if expr_df.empty:
        raise RuntimeError(f"{expr_path}: 未读到数值数据(首列应为基因符号,其余列为样本数值)")

    lab = pd.read_csv(labels_path, encoding="utf-8-sig")
    if lab.shape[1] < 2:
        raise RuntimeError(f"{labels_path}: 标签文件需两列(样本ID, 分组),实际 {lab.shape[1]} 列")
    sids = lab.iloc[:, 0].astype(str).str.strip()
    if sids.duplicated().any():
        dup = sorted(set(sids[sids.duplicated()]))
        raise RuntimeError(f"{labels_path}: 样本ID重复: {dup[:5]}")
    labels_raw = pd.Series(lab.iloc[:, 1].astype(str).str.strip().values, index=sids.values)

    id_set = set(labels_raw.index)
    col_hits = sum(str(c).strip() in id_set for c in expr_df.columns)
    row_hits = sum(str(i).strip() in id_set for i in expr_df.index)
    if col_hits == 0 and row_hits == 0:
        raise RuntimeError(
            f"表达矩阵与标签的样本ID无交集(矩阵列名如 {list(expr_df.columns[:3])},"
            f"标签ID如 {list(labels_raw.index[:3])}),请核对两文件的样本命名是否一致")
    if row_hits > col_hits:
        _note(f"样本方向:行=样本(行索引命中 {row_hits} > 列名命中 {col_hits}),已转置")
        expr_df = expr_df.T
    expr_df.columns = [str(c).strip() for c in expr_df.columns]

    vals = sorted(set(labels_raw))
    if len(vals) != 2:
        raise RuntimeError(
            f"分组列应恰好 2 类(病/对照),实际 {len(vals)} 类: {vals}——"
            "如数据含中间分级,请先在 labels 文件里只保留目标两类(其余样本不会自动丢弃)")
    case = (case_label or "").strip()
    if case:
        hit = [v for v in vals if v.lower() == case.lower()]
        if not hit:
            raise RuntimeError(f"case_label {case_label!r} 不在分组取值 {vals} 中")
        case = hit[0]
    else:
        prefer = ("1", "case", "disease", "treated", "patient")
        lower = {v.lower(): v for v in vals}
        case = next((lower[p] for p in prefer if p in lower), "")
        if not case:
            raise RuntimeError(
                f"无法从分组取值 {vals} 推断病例组,请显式指定 case_label"
                f"(如 case_label={vals[0]!r})")
    control = vals[0] if vals[1] == case else vals[1]
    labels = labels_raw.map({case: 1, control: 0}).astype(int)
    n_case, n_ctrl = int((labels == 1).sum()), int((labels == 0).sum())
    if n_case == 0 or n_ctrl == 0:
        raise RuntimeError(f"分组结果为空(case={n_case}, control={n_ctrl}),标签方向无效")
    import config
    for name, n in (("case", n_case), ("control", n_ctrl)):
        if n < config.ML_MIN_SAMPLES_PER_GROUP:
            raise RuntimeError(
                f"{name} 组仅 {n} 个样本(低于 {config.ML_MIN_SAMPLES_PER_GROUP}):"
                "小样本下统计与 ML 均无意义,请补充样本或换数据")
    _note(f"case({case})={n_case} control({control})={n_ctrl}")

    if not looks_symbol_level(expr_df.index):
        raise RuntimeError(
            "本地矩阵行标识疑似探针级(_at/ENSG/ILMN_/A_ 等),本地模式无 GPL 注释表无法映射,"
            "请先转为基因符号矩阵(方案边界:映射失败明确报错,不做激进猜测)")
    n_dup = int(expr_df.index.duplicated().sum())
    if n_dup:
        _note(f"注意:重复基因符号 {n_dup} 个,各保留第一条")
        expr_df = expr_df[~expr_df.index.duplicated(keep="first")]
    missing = [s for s in labels.index if s not in set(expr_df.columns)]
    if missing:
        raise RuntimeError(f"标签中的样本在表达矩阵中缺失: {missing[:5]}"
                           f"{' …' if len(missing) > 5 else ''}")
    tag = os.path.splitext(os.path.basename(expr_path))[0]
    info = {"case_value": case, "control_value": control,
            "n_case": n_case, "n_control": n_ctrl}
    return expr_df, labels, tag, info


# ---------------- 全链主入口 ----------------


def run_screen_chain(*, gse: str | None = None,
                     case_pattern: str | None = None,
                     control_pattern: str | None = None,
                     expr: str | None = None,
                     labels: str | None = None,
                     case_label: str | None = None,
                     force_refresh: bool = False,
                     force_log: bool | None = None,
                     prefilter_top: int | None = None,
                     bootstrap_rounds: int | None = None,
                     ranked_csvs: dict[str, str] | None = None,
                     top_k: int | None = None,
                     out_base: str | None = None,
                     out_subdir: str | None = "mlscreen",
                     progress=None) -> dict:
    """执行 screen 全链并落盘,返回内存产物供 CLI/页面呈现。

    返回 dict:consensus/shortlist/funnel/deg(DataFrame)、meta/roc(ML 产物)、
    gse/gpl/tag/out_dir/summary(已写入的 summary.json 内容)、log_branch/leak_note。

    out_base:结果根目录(None → default_results_root);out_subdir:CLI 用默认
    "mlscreen"(方案 §6 的 results/mlscreen 布局),页面传 None 使产物直接落
    <out_base>(GEOresult)根下,与 MDresult 平级。产物目录为
    <out_base>[/<out_subdir>]/<tag>_<时间戳>/,在全部计算成功后才创建
    (失败语义:不部分落盘,summary.json 最后写入)。
    """
    import config
    from .expression import differential_expression, map_probe_to_symbol, preprocess
    from .geo import assign_groups, fetch_geo_dataset, sample_source_label
    from .ml_screening import build_shortlist, deg_prefilter, run_ml_screening

    def _report(stage: str | None, message: str) -> None:
        if progress:
            progress(stage, message)

    t0 = time.time()

    # 参数边界(CLI 与页面共用同一套防呆)
    if bool(expr) == bool(gse):
        raise RuntimeError("数据源必须二选一:gse(GEO 数据集)或 expr(本地矩阵 + labels)")
    if expr:
        if not labels:
            raise RuntimeError("本地矩阵模式需同时提供 labels(两列 CSV:样本ID, 分组)")
        if case_pattern or control_pattern:
            raise RuntimeError("本地矩阵模式不用 case_pattern/control_pattern:"
                               "分组由 labels 文件与 case_label 确定")
    elif not case_pattern or not control_pattern:
        raise RuntimeError("GSE 模式需提供 case_pattern/control_pattern"
                           "(可先跑 deg 子命令不带分组,打印分组字段总览再定正则)")

    ds = None
    local_info: dict = {}
    gse_id = gpl_id = None
    if expr:
        _report("fetch", f"读取本地表达矩阵 {expr}")
        matrix_raw, label_series, tag, local_info = load_local_expr(
            expr, labels, case_label, progress=progress)
        _report(None, f"矩阵 {matrix_raw.shape[0]} 基因 × {matrix_raw.shape[1]} 样本")
    else:
        gse_id = gse.strip().upper()
        _report("fetch", f"获取 {gse_id}")
        ds = fetch_geo_dataset(gse_id, force_refresh=force_refresh)
        if ds.issues:
            raise RuntimeError(f"{gse_id}: " + ";".join(ds.issues))
        gpl_id = ds.gpl
        tag = gse_id
        matrix_raw = ds.table
        _report(None, f"样本 {len(ds.samples)} | 矩阵 {ds.table.shape} | 平台 {ds.gpl}")
        if ds.cache_fallback:
            _report(None, f"警告: GEO 文件下载失败,已沿用本地旧缓存(缓存日期 "
                          f"{ds.cache_date});CLI 可加 --force-refresh 强制重下")

    # 样本来源标注(承检索启发式,方案 §4.3 shortlist 的 evidence_scope 列;
    # geo 分支取数据集标题/摘要判定,本地矩阵无数据集元数据 → 未标注)
    if ds is not None:
        evidence_scope = sample_source_label(ds.meta.get("Series_title", ""),
                                             ds.meta.get("Series_summary", ""))
        _report(None, f"样本来源(启发式): {evidence_scope}")
    else:
        evidence_scope = "未标注(本地矩阵)"

    _report("group", "分组" if ds is not None else "分组(labels 文件,硬边界已校验)")
    if ds is not None:
        label_series = assign_groups(ds, case_pattern, control_pattern)
        n_case, n_ctrl = int((label_series == 1).sum()), int((label_series == 0).sum())
        _report(None, f"case={n_case} control={n_ctrl} 丢弃={len(ds.samples) - len(label_series)}")
    else:
        n_case, n_ctrl = local_info["n_case"], local_info["n_control"]

    _report("preprocess", "预处理与映射")
    matrix, log_branch = preprocess(matrix_raw, force_log=force_log)
    if ds is not None:
        matrix = map_probe_to_symbol(matrix, ds.gpl)
    _report(None, f"log2 决策: {log_branch} | 基因数 {len(matrix)}")

    _report("deg", "差异表达")
    deg = differential_expression(matrix, label_series)
    sig_n = int((deg["FDR"] < 0.05).sum())
    _report(None, f"FDR<0.05: {sig_n} | 全部基因 {len(deg)}")

    _report("ml", "DEG 预过滤 + ML 共识(六算法 × 三验证)")
    deg_top, leak_note = deg_prefilter(deg, top_n=prefilter_top)
    genes = deg_top["gene"].tolist()
    # 只取有标签的样本(丢弃的中间分级样本不在标签中,全样本矩阵会 KeyError,审核 BUG-7)
    labeled = list(label_series[label_series == 1].index) + list(label_series[label_series == 0].index)
    X = matrix[labeled].T[genes]  # 转置为样本×基因(仅已标记样本)
    res = run_ml_screening(X, label_series.loc[X.index], bootstrap_rounds=bootstrap_rounds)
    m = res["meta"]
    _report(None, f"ML 输入: {m['n_samples']}样本 × {m['n_genes']}基因 | "
                  f"{m['cv_mode']} AUC={m['cv_auc']} | 算法一致性={m['algo_spearman_mean']} | "
                  f"耗时 {m['elapsed_s']}s")

    _report("shortlist", "shortlist 与漏斗")
    ranked = ranked_csvs or {}
    resolved_top_k = top_k if top_k is not None else config.ML_SHORTLIST_TOP_K
    sl, funnel = build_shortlist(res["consensus"], deg, ranked, top_k=resolved_top_k,
                                 evidence_scope=evidence_scope)
    n_hit = int(sl["in_shortlist"].sum())
    _report(None, f"shortlist: {n_hit} 个全证据交集基因")

    base = out_base or default_results_root()
    # out_subdir 语义见 docstring:CLI 落 results/mlscreen,页面直落 GEOresult 根
    out_parent = os.path.join(base, out_subdir) if out_subdir else base
    out_dir = os.path.join(out_parent, f"{tag}_{time.strftime('%Y%m%d_%H%M%S')}")
    _report("save", f"落盘 → {out_dir}")
    os.makedirs(out_dir, exist_ok=True)
    res["consensus"].to_csv(os.path.join(out_dir, "ml_consensus_scores.csv"),
                            index_label="gene", encoding="utf-8-sig")
    # shortlist.csv 带首行定位注释(方案 §11;读取方式见 SHORTLIST_POSITIONING_NOTE)
    with open(os.path.join(out_dir, "shortlist.csv"), "w",
              encoding="utf-8-sig", newline="") as f:
        f.write(f"# {SHORTLIST_POSITIONING_NOTE}\n")
        sl.to_csv(f, index=False)
    funnel.to_csv(os.path.join(out_dir, "funnel.csv"), index=False, encoding="utf-8-sig")
    deg.to_csv(os.path.join(out_dir, "deg_results.csv"), index=False, encoding="utf-8-sig")
    from .ml_screening import write_consensus_curve_png, write_roc_png
    write_roc_png(os.path.join(out_dir, f"ml_roc_{tag}_{m['cv_mode']}.png"),
                  res["roc"]["y_true"], res["roc"]["scores"], str(tag))
    write_consensus_curve_png(os.path.join(out_dir, "consensus_curve.png"),
                              res["consensus"], resolved_top_k)
    summary = {
        "module": "geo-ml-screen", "ml_version": config.ML_SCREENING_VERSION,
        "source": "local-expr" if expr else "geo",
        "gse": gse_id, "gpl": gpl_id, "expr": expr, "labels": labels,
        "n_samples": len(ds.samples) if ds is not None else n_case + n_ctrl,
        "groups": {"case": n_case, "control": n_ctrl,
                   "case_pattern": case_pattern, "control_pattern": control_pattern,
                   "case_label": (local_info or {}).get("case_value", case_label)},
        "preprocess": log_branch,
        "deg_significant_fdr05": sig_n,
        "deg_prefilter": {"top_n": len(deg_top), "leak_note": leak_note},
        "ml_meta": m,
        "ranked_csvs": {n: {"path": p, "mtime": os.path.getmtime(p),
                            "size": os.path.getsize(p)}
                        for n, p in ranked.items()},
        "shortlist_count": n_hit,
        "evidence_scope": evidence_scope,
        "geo_cache_fallback": bool(getattr(ds, "cache_fallback", False)),
        "geo_cache_date": getattr(ds, "cache_date", "") or "",
        "elapsed_total_s": round(time.time() - t0, 1),
        "args": {"gse": gse_id, "expr": expr, "labels": labels,
                 "case_label": case_label, "case_pattern": case_pattern,
                 "control_pattern": control_pattern,
                 "force_refresh": force_refresh, "force_log": force_log,
                 "prefilter_top": prefilter_top,
                 "bootstrap_rounds": bootstrap_rounds, "top_k": top_k,
                 "ranked_csvs": ranked, "out_base": out_base},
    }
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1)

    return {"consensus": res["consensus"], "shortlist": sl, "funnel": funnel,
            "deg": deg, "meta": m, "roc": res["roc"], "roc_png":
            os.path.join(out_dir, f"ml_roc_{tag}_{m['cv_mode']}.png"),
            "gse": gse_id, "gpl": gpl_id, "tag": tag, "out_dir": out_dir,
            "summary": summary, "log_branch": log_branch, "leak_note": leak_note}


# ---------------- 结果陈旧判定(方案 §M4 失效场景③) ----------------


def stale_reasons(summary: dict) -> list[dict]:
    """判定一次历史结果是否已陈旧,返回机器可读原因列表(空列表 = 未陈旧)。

    两种触发(方案 §M4):算法版本升级(summary.ml_version != 当前
    ML_SCREENING_VERSION);平台排序 CSV 在结果生成后被重跑(mtime/size 指纹
    对不上)或已删除。页面把原因映射为 i18n 文案,自检直接断言原因类型。
    """
    import config
    reasons: list[dict] = []
    if summary.get("ml_version") != config.ML_SCREENING_VERSION:
        reasons.append({"type": "version", "old": summary.get("ml_version"),
                        "new": config.ML_SCREENING_VERSION})
    for name, fp in (summary.get("ranked_csvs") or {}).items():
        path = fp.get("path") if isinstance(fp, dict) else None
        if not path or not os.path.isfile(path):
            reasons.append({"type": "ranked_missing", "name": name})
        elif (os.path.getmtime(path) != fp.get("mtime")
              or os.path.getsize(path) != fp.get("size")):
            reasons.append({"type": "ranked_changed", "name": name})
    return reasons
