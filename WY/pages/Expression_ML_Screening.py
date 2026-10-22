"""GEO 表达谱与 ML 共识筛选页面(M4,模块收官 GUI)。

四区布局(方案 §M4):
① 输入区:疾病名检索 GEO 数据集 / 直接输入 GSE 编号 / 本地矩阵直跑(方案 §6)
② 配置区:分组正则 + 平台排序 CSV 多选(按 analysis_meta.json 元数据展示与过滤,
   只列排序表 *_query_genes.csv)+ top-K/预过滤/bootstrap 参数
③ 运行区:一键执行(DEG → ML 共识 → shortlist),产物经共用管线落盘
   Static/data/downloads/GEOresult/<GSE>_<时间戳>/,与 PPI(_data)/对接(MDresult)
   同属 Web 端产物区;CLI 默认仍落 results/mlscreen(方案 §6),历史列表双根同扫;
   GEO 下载失败降级用旧缓存时,显示黄色警告条 + 强制更新按钮(失效场景②)
④ 结果区:共识表/ROC/shortlist 下载/历史结果列表(读各次 summary.json,
   含陈旧检测:算法版本升级或排序 CSV 重跑即提示建议重算,方案失效场景③)

本页面不新增能力、只换交互形态:全链流程在 core/ml_pipeline.run_screen_chain,
CLI 与页面共用同一套编排,防止两处漂移。
"""
import json
import os
import warnings

import streamlit as st

import config
from core.geo import fetch_geo_dataset
from core.ml_pipeline import default_results_root, run_screen_chain, stale_reasons
from core.ppi import base_disease_of, list_analyzed_pairs
from i18n import t
from ui.widgets import back_button, hide_sidebar

st.set_page_config(page_title=t("eml.page_title"), layout="wide")

hide_sidebar()

st.title(t("eml.title"))
st.caption(t("eml.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    back_button()

st.divider()

# ---------- 会话状态 ----------
ss = st.session_state
ss.setdefault("eml", {})
ss.setdefault("eml_busy", False)


def _fmt_size(n: float) -> str:
    """字节数 → KB/MB 显示(下载耗时预估用)。"""
    return f"{n / 1024:.0f} KB" if n < 1048576 else f"{n / 1048576:.1f} MB"


def _fmt_eta(seconds: float) -> str:
    """秒数 → 秒 / 分钟 / 小时+分钟 显示。"""
    if seconds < 90:
        return t("eml.eta_sec", n=max(1, round(seconds)))
    if seconds < 3600:
        return t("eml.eta_min", n=max(1, round(seconds / 60)))
    return t("eml.eta_hr_min", h=int(seconds // 3600), m=round(seconds % 3600 / 60))

# ---------- ① 输入区 ----------
st.subheader(t("eml.step1"))

tab_search, tab_direct, tab_local = st.tabs(
    [t("eml.tab_search"), t("eml.tab_direct"), t("eml.tab_local")])

with tab_search:
    # 疾病名优先从本机 App 分析记录取(用户跑过靶点的疾病):下拉列出干净的
    # 基准疾病名(剥掉 rwrNN 变体后缀),accept_new_options 仍允许手输新名称;
    # 无记录时退回文本框直接输入
    _pairs = list_analyzed_pairs()
    _seen: set = set()
    disease_options = []
    for _p in _pairs:
        _key = _p["base_disease"].strip().lower()
        if _key and _key not in _seen:
            _seen.add(_key)
            disease_options.append(_p["base_disease"].strip())
    col_d, col_l = st.columns([3, 1])
    with col_d:
        if disease_options:
            disease = st.selectbox(t("eml.disease_label"), options=disease_options,
                                   index=None, accept_new_options=True,
                                   key="eml_disease_sel") or ""
            st.caption(t("eml.analyzed_hint", n=len(disease_options)))
        else:
            disease = st.text_input(t("eml.disease_label"), key="eml_disease_txt",
                                    placeholder=t("eml.disease_ph"))
            st.caption(t("eml.analyzed_empty"))
    with col_l:
        limit = st.number_input(t("eml.limit_label"), min_value=5, max_value=50,
                                value=15, key="eml_limit")
    if st.button(t("eml.search_btn"), key="eml_search", type="primary"):
        if not disease.strip():
            st.warning(t("eml.disease_empty"))
        else:
            with st.spinner(t("eml.searching")):
                try:
                    from core.geo import search_geo_datasets
                    results = search_geo_datasets(disease.strip(), limit=int(limit))
                    ss.eml["search_results"] = results
                except Exception as exc:
                    st.error(t("eml.search_failed", exc=exc))
    if "search_results" in ss.eml and not ss.eml["search_results"].empty:
        st.dataframe(ss.eml["search_results"], use_container_width=True,
                     hide_index=True)
        st.caption(t("eml.search_hint"))
    elif "search_results" in ss.eml:
        st.info(t("eml.no_results"))

with tab_direct:
    st.caption(t("eml.direct_hint"))

with tab_local:
    st.caption(t("eml.local_hint"))
    col_le, col_ll = st.columns(2)
    with col_le:
        expr_path = st.text_input(t("eml.local_expr"), key="eml_local_expr",
                                  placeholder=r"D:\data\GSE12345_expr.csv")
    with col_ll:
        labels_path = st.text_input(t("eml.local_labels"), key="eml_local_labels",
                                    placeholder=r"D:\data\labels.csv")
    case_label = st.text_input(t("eml.local_case"), key="eml_local_case",
                               placeholder=t("eml.local_case_ph"))

# 数据集选择:检索后有结果 → 下拉点选(accept_new_options 仍可直填其他编号);
# 未检索 → 文本框直填(tab_direct 路径)
_srch = ss.eml.get("search_results")
_gse_opts = []
if _srch is not None and not _srch.empty and "gse" in _srch.columns:
    _gse_opts = [str(g).strip() for g in _srch["gse"].tolist() if str(g).strip()]
if _gse_opts:
    _sel = st.selectbox(t("eml.select_gse"), options=_gse_opts, index=None,
                        accept_new_options=True, key="eml_selected_gse_sel",
                        placeholder="GSE63067")
    selected_gse = (_sel or "").strip().upper()
else:
    selected_gse = st.text_input(t("eml.select_gse"), key="eml_selected_gse",
                                 placeholder="GSE63067").strip().upper()
if selected_gse:
    ss.eml["selected_gse"] = selected_gse
    st.info(t("eml.gse_confirmed", gse=selected_gse))
    # 下载耗时预估(方案 §4.1:跨境带宽波动大,选型阶段先知道要等多久)。
    # 每个 GSE 每会话只探测一次(存 session),网络变化后换 GSE 或刷新页面重测;
    # 本地矩阵模式不探测。单个小 Range 请求,不落盘、不写缓存、不影响正式下载。
    if not expr_path.strip() and ss.eml.get("probe", {}).get("gse") != selected_gse:
        with st.spinner(t("eml.est_probing")):
            try:
                from core.geo import probe_series_matrix
                ss.eml["probe"] = {"gse": selected_gse,
                                   "r": probe_series_matrix(selected_gse), "err": ""}
            except Exception as _exc:
                ss.eml["probe"] = {"gse": selected_gse, "r": None, "err": str(_exc)}
    _probe = ss.eml.get("probe", {})
    if not expr_path.strip() and _probe.get("gse") == selected_gse:
        _cache_dir = os.path.join(config.CACHE_DIR, config.GEO_SERIES_CACHE_SUBDIR)
        if _probe.get("r"):
            _pr = _probe["r"]
            _eta = _pr["size_bytes"] / _pr["speed_bps"] if _pr["speed_bps"] > 0 else 10 ** 9
            st.caption(t("eml.est_line",
                         size=_fmt_size(_pr["size_bytes"]),
                         speed=f"{_pr['speed_bps'] / 1024:.0f} KB/s",
                         eta=_fmt_eta(_eta)))
            if not _pr["reliable"] or _eta > 300:  # 预计超 5 分钟或测速被时限截断
                st.caption(t("eml.est_slow", path=_cache_dir))
        elif _probe.get("err"):
            st.caption(t("eml.est_failed", path=_cache_dir))

# 数据源判定:本地矩阵(两路都填)与 GSE 二选一(方案 §6,与管线防呆一致)
use_local = bool(expr_path.strip()) and bool(labels_path.strip())
local_partial = bool(expr_path.strip()) != bool(labels_path.strip())
has_gse = bool(selected_gse.strip())
source_conflict = use_local and has_gse

st.divider()

# ---------- ② 配置区 ----------
st.subheader(t("eml.step2"))

col_case, col_ctrl = st.columns(2)
with col_case:
    case_pat = st.text_input(t("eml.case_pattern"), key="eml_case",
                             placeholder="steatosis|non-alcoholic")
with col_ctrl:
    ctrl_pat = st.text_input(t("eml.ctrl_pattern"), key="eml_ctrl",
                             placeholder="healthy")

# 分组总览按钮(运行中禁用:重复触发会对同一数据集并发下载,文件锁只能保正确、不能省时间)
if selected_gse.strip() and st.button(t("eml.overview_btn"), key="eml_overview",
                                      disabled=ss.eml_busy):
    ss.eml_busy = True
    try:
        with st.spinner(t("eml.loading")):
            try:
                from core.geo import fetch_geo_dataset, group_fields_overview
                ds = fetch_geo_dataset(selected_gse.strip())
                ss.eml["dataset"] = ds
                if ds.issues:
                    st.warning(t("eml.dataset_issues", issues="; ".join(ds.issues)))
                else:
                    overview = group_fields_overview(ds)
                    st.markdown(t("eml.overview_title"))
                    st.dataframe(overview, use_container_width=True, hide_index=True)
                    st.caption(t("eml.overview_hint"))
            except Exception as exc:
                st.error(t("eml.overview_failed", exc=exc))
    finally:
        ss.eml_busy = False

# 平台排序 CSV 自动发现
col_rank, col_param = st.columns([2, 1])
with col_rank:
    st.markdown(t("eml.ranked_csv_title"))
    # 自动发现:与疾病下拉同源(list_analyzed_pairs),只收排序表 *_query_genes.csv
    # (nodes/edges 是网络结构表,混入会被当排序用,语义错误)。列表跟随上方选定
    # 的疾病自动限定 —— rwr 等参数变体的排序以文件维度并列展示,跨疾病的排序
    # 不再出现,从源头排除"把别的疾病的排序勾进交集"的可能(方案 §M4 按元数据
    # 匹配,不用文件名猜)
    ranked_files = []
    for _rec in list_analyzed_pairs():
        try:
            _dir_files = os.listdir(_rec["directory"])
        except OSError:
            continue
        for f in sorted(_dir_files):
            if not f.endswith("_query_genes.csv"):
                continue
            title = " · ".join(x for x in (
                os.path.basename(_rec["directory"])[:-len("_data")],
                _rec["analyzed_at"][:19],
                f"STRING {_rec['string_version']}" if _rec["string_version"] else "",
            ) if x)
            ranked_files.append({"label": title,
                                 "path": os.path.join(_rec["directory"], f),
                                 "drug": _rec["drug"],
                                 "base_disease": _rec["base_disease"]})
    # 当前疾病上下文(下拉或手输,同样剥掉可能手输的 rwrNN 后缀);为空则展示全部。
    # 匹配只允许"用户输的更泛 → 命中更具体的记录"(如输 cancer 带出 gastric cancer),
    # 不做反向:选 gastric cancer 不能把泛癌 cancer 的排序带进交集(疾病混淆)
    _ctx = base_disease_of(disease).lower() if disease.strip() else ""
    if ranked_files:
        f_drug = st.text_input(t("eml.filter_drug"), key="eml_ranked_drug").strip().lower()
        shown = [rf for rf in ranked_files
                 if (not _ctx
                     or _ctx == rf["base_disease"].lower()
                     or _ctx in rf["base_disease"].lower())
                 and (not f_drug or f_drug in rf["drug"].lower())]
        if shown:
            selected_files = st.multiselect(
                t("eml.ranked_csv_select"),
                options=[rf["label"] for rf in shown],
                key="eml_ranked")
            ss.eml["ranked_paths"] = [
                rf["path"] for rf in shown if rf["label"] in selected_files]
        else:
            st.caption(t("eml.no_ranked_for_disease"))
    else:
        st.caption(t("eml.no_ranked_csv"))

with col_param:
    top_k = st.number_input(t("eml.topk_label"), min_value=5, max_value=100,
                            value=config.ML_SHORTLIST_TOP_K, key="eml_topk")
    # 默认随 config(500 基因 × 30 轮,页面/CLI 统一口径):本机实测 500×30 ML 段约 2 分钟
    # (LOOCV 小样本数据集约 2.5 分钟);SVM-RFE 对特征数近平方级,调大候选池会显著变慢
    prefilter_n = st.number_input(
        t("eml.prefilter_label"), min_value=50, max_value=5000,
        value=config.ML_DEG_PREFILTER_TOP_N, key="eml_prefilter")
    bootstrap_n = st.number_input(
        t("eml.bootstrap_label"), min_value=5, max_value=200,
        value=config.ML_BOOTSTRAP_ROUNDS, key="eml_bootstrap")
    st.caption(t("eml.perf_hint"))

st.divider()

# ---------- ③ 运行区 ----------
st.subheader(t("eml.step3"))
can_run = ((has_gse and case_pat.strip() and ctrl_pat.strip()) or use_local) \
    and not source_conflict
if source_conflict:
    st.error(t("eml.source_conflict"))
elif local_partial:
    st.warning(t("eml.local_missing"))
elif not can_run:
    st.info(t("eml.fill_required"))
if ss.eml_busy:
    st.caption(t("eml.busy_hint"))

if st.button(t("eml.run_btn"), key="eml_run", type="primary",
             disabled=(not can_run) or ss.eml_busy):
    ss.eml_busy = True
    # 进度回调:管线发七阶段键,页面按 i18n 映射展示;补充信息行(message 为中文)
    # 不直显,避免与界面语言冲突
    _progress_map = {
        "fetch": t("eml.progress_fetch"),
        "group": t("eml.progress_group"),
        "preprocess": t("eml.progress_preprocess"),
        "deg": t("eml.progress_deg"),
        "ml": t("eml.progress_ml"),
        "shortlist": t("eml.progress_shortlist"),
        "save": t("eml.progress_save"),
    }

    def _on_progress(stage, _message):
        if stage in _progress_map:
            st.markdown(_progress_map[stage])

    try:
        with st.spinner(t("eml.running")):
            try:
                ranked_csvs = {os.path.splitext(os.path.basename(p))[0]: p
                               for p in ss.eml.get("ranked_paths", [])}
                common = dict(prefilter_top=int(prefilter_n),
                              bootstrap_rounds=int(bootstrap_n),
                              ranked_csvs=ranked_csvs, top_k=int(top_k),
                              out_base=config.GEO_RESULT_DIR,
                              out_subdir=None,
                              progress=_on_progress)
                if use_local:
                    res = run_screen_chain(
                        expr=expr_path.strip(), labels=labels_path.strip(),
                        case_label=case_label.strip() or None, **common)
                else:
                    res = run_screen_chain(
                        gse=selected_gse.strip().upper(),
                        case_pattern=case_pat.strip(),
                        control_pattern=ctrl_pat.strip(), **common)
                ss.eml["result"] = {
                    "consensus": res["consensus"],
                    "shortlist": res["shortlist"],
                    "funnel": res["funnel"],
                    "deg": res["deg"],
                    "meta": res["meta"],
                    "roc": res["roc"],
                    "gse": res["gse"] or res["tag"],
                    "out_dir": res["out_dir"],
                    "elapsed": res["summary"]["elapsed_total_s"],
                }
                st.success(t("eml.run_success",
                             gse=res["gse"] or res["tag"],
                             elapsed=res["summary"]["elapsed_total_s"]))
                st.caption(t("eml.saved_to", dir=res["out_dir"]))
                # 降级标记入会话(方案失效场景②):强制更新按钮必须渲染在运行按钮
                # 块之外 —— 块内按钮在本次 run 结束后不再重建,点击会被 Streamlit
                # 静默丢弃(审查 BUG-12)
                if res["summary"].get("geo_cache_fallback"):
                    ss.eml["geo_fallback"] = {
                        "gse": selected_gse.strip().upper(),
                        "date": res["summary"].get("geo_cache_date", "")}
                else:
                    ss.eml.pop("geo_fallback", None)
            except Exception as exc:
                st.error(t("eml.run_failed", exc=exc))
    finally:
        ss.eml_busy = False

# ---------- GEO 下载降级提示(方案失效场景②) ----------
# 常驻渲染:强制更新按钮任何时候可点(点击触发重跑,本块会重新执行,点击不会丢)
if ss.eml.get("geo_fallback"):
    _fb = ss.eml["geo_fallback"]
    st.warning(t("eml.cache_fallback", date=_fb["date"]))
    if st.button(t("eml.force_redownload"), key="eml_force_dl",
                 disabled=ss.eml_busy):
        with st.spinner(t("eml.loading")):
            try:
                fetch_geo_dataset(_fb["gse"], force_refresh=True)
                ss.eml.pop("geo_fallback", None)
                st.success(t("eml.force_dl_ok"))
            except Exception as exc:
                st.error(t("eml.run_failed", exc=exc))

st.divider()

# ---------- ④ 结果区 ----------
st.subheader(t("eml.step4"))

if "result" not in ss.eml:
    st.info(t("eml.no_result"))
else:
    r = ss.eml["result"]
    m = r["meta"]

    # 指标卡
    col1, col2, col3, col4 = st.columns(4)
    col1.metric(t("eml.metric_samples"), m["n_samples"])
    col2.metric(t("eml.metric_genes"), m["n_genes"])
    col3.metric(t("eml.metric_auc"), m["cv_auc"] if m["cv_auc"] is not None else "N/A")
    col4.metric(t("eml.metric_spearman"), m["algo_spearman_mean"])

    # ROC 曲线
    st.markdown(t("eml.roc_title"))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.metrics import roc_curve, roc_auc_score
        fpr, tpr, _ = roc_curve(r["roc"]["y_true"], r["roc"]["scores"])
        auc_val = roc_auc_score(r["roc"]["y_true"], r["roc"]["scores"])
        fig, ax = plt.subplots(figsize=(4.2, 4.0), dpi=120)
        ax.plot(fpr, tpr, lw=1.6, color="#1f6fb4")
        ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color="#999")
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_title(f"CV-ROC (AUC={auc_val:.3f})", fontsize=10)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
    except Exception:
        st.warning(t("eml.roc_failed"))

    # 共识表(top 50,可排序)
    st.markdown(t("eml.consensus_title"))
    top50 = r["consensus"].head(50).copy()
    top50.insert(0, "gene", top50.index)
    st.dataframe(top50, use_container_width=True, hide_index=True)

    # 漏斗
    st.markdown(t("eml.funnel_title"))
    st.dataframe(r["funnel"], use_container_width=True, hide_index=True)

    # shortlist 下载
    st.markdown(t("eml.shortlist_title",
                  n=int(r["shortlist"]["in_shortlist"].sum())))
    st.dataframe(r["shortlist"], use_container_width=True, hide_index=True)

    # CSV 下载按钮
    col_dl1, col_dl2, col_dl3 = st.columns(3)
    with col_dl1:
        st.download_button(
            t("eml.dl_consensus"),
            data=r["consensus"].to_csv(index_label="gene").encode("utf-8-sig"),
            file_name=f"ml_consensus_{r['gse']}.csv",
            mime="text/csv", key="eml_dl_consensus")
    with col_dl2:
        st.download_button(
            t("eml.dl_shortlist"),
            data=r["shortlist"].to_csv(index=False).encode("utf-8-sig"),
            file_name=f"shortlist_{r['gse']}.csv",
            mime="text/csv", key="eml_dl_shortlist")
    with col_dl3:
        st.download_button(
            t("eml.dl_deg"),
            data=r["deg"].to_csv(index=False).encode("utf-8-sig"),
            file_name=f"deg_{r['gse']}.csv",
            mime="text/csv", key="eml_dl_deg")

st.divider()

# ---------- 历史结果列表(方案 §M4 ④:读各次 summary.json;陈旧检测=失效场景③) ----------
st.subheader(t("eml.history_title"))

# 双根同扫:页面新产物在 downloads/GEOresult,CLI 产物在 results/mlscreen,都不丢
_history_roots = [
    ("GEOresult", config.GEO_RESULT_DIR),
    ("mlscreen", os.path.join(default_results_root(), "mlscreen")),
]
_entries = []
_seen_labels = set()
for _root_tag, _root in _history_roots:
    if not os.path.isdir(_root):
        continue
    for _d in sorted(os.listdir(_root), reverse=True):
        _spath = os.path.join(_root, _d, "summary.json")
        if not os.path.isfile(_spath):
            continue
        try:
            with open(_spath, encoding="utf-8") as _f:
                _summary = json.load(_f)
        except Exception:
            continue  # 损坏的 summary 不挡其余历史
        _label = _d if _d not in _seen_labels else f"{_d}·{_root_tag}"
        _seen_labels.add(_label)
        _entries.append({"label": _label, "dir": _d, "root": _root,
                         "summary": _summary})

if not _entries:
    st.caption(t("eml.history_empty"))
else:
    _sel = st.selectbox(t("eml.history_select"),
                        options=[e["label"] for e in _entries], key="eml_history")
    _entry = next(e for e in _entries if e["label"] == _sel)
    _s = _entry["summary"]
    _m = _s.get("ml_meta") or {}
    _reasons = stale_reasons(_s)
    if _reasons:
        _msgs = []
        for _r in _reasons:
            if _r["type"] == "version":
                _msgs.append(t("eml.stale_version", old=_r.get("old"), new=_r.get("new")))
            elif _r["type"] == "ranked_changed":
                _msgs.append(t("eml.stale_ranked", name=_r["name"]))
            else:
                _msgs.append(t("eml.stale_ranked_missing", name=_r["name"]))
        st.warning(f"{t('eml.stale_title')}: " + "; ".join(_msgs))
    _c1, _c2, _c3, _c4, _c5 = st.columns(5)
    _c1.metric(t("eml.metric_run"), _s.get("gse")
               or os.path.basename(str(_s.get("expr") or "?")))
    _c2.metric(t("eml.metric_samples"), _s.get("n_samples", "?"))
    _c3.metric(t("eml.metric_auc"),
               _m.get("cv_auc") if _m.get("cv_auc") is not None else "N/A")
    _c4.metric(t("eml.metric_shortlist"), _s.get("shortlist_count", "?"))
    _c5.metric(t("eml.metric_version"), f"v{_s.get('ml_version', '?')}")
    st.caption(f"{t('eml.history_files')}: "
               f"{os.path.join(_entry['root'], _entry['dir'])}")

st.divider()
st.caption(t("eml.footer"))
