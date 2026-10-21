"""GEO 表达谱与 ML 共识筛选页面(M4,模块收官 GUI)。

四区布局(方案 §M4):
① 输入区:疾病名检索 GEO 数据集(或直接输入 GSE 编号)
② 配置区:分组正则 + 平台排序 CSV 多选 + top-K/预过滤/bootstrap 参数
③ 运行区:一键执行(DEG → ML 共识 → shortlist)
④ 结果区:共识表/ROC/shortlist 下载/历史结果

本页面不新增能力、只换交互形态,底层调用与 CLI(geo_ml_screen.py)完全同一套 core 模块。
"""
import os
import warnings

import streamlit as st

import config
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

# ---------- ① 输入区 ----------
st.subheader(t("eml.step1"))

tab_search, tab_direct = st.tabs([t("eml.tab_search"), t("eml.tab_direct")])

with tab_search:
    col_d, col_l = st.columns([3, 1])
    with col_d:
        disease = st.text_input(t("eml.disease_label"), key="eml_disease",
                                placeholder=t("eml.disease_ph"))
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

# 数据集选择(检索结果点选或直填)
selected_gse = st.text_input(
    t("eml.select_gse"),
    value=ss.eml.get("selected_gse", ""),
    key="eml_selected_gse",
    placeholder="GSE63067")
if selected_gse.strip():
    ss.eml["selected_gse"] = selected_gse.strip().upper()
    st.info(t("eml.gse_confirmed", gse=selected_gse.strip().upper()))

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

# 分组总览按钮
if selected_gse.strip() and st.button(t("eml.overview_btn"), key="eml_overview"):
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

# 平台排序 CSV 自动发现
col_rank, col_param = st.columns([2, 1])
with col_rank:
    st.markdown(t("eml.ranked_csv_title"))
    # 自动发现 downloads 下的 *_query_genes.csv 和 *_nodes.csv
    ranked_files = []
    dl_dir = config.DOWNLOADS_DIR
    if os.path.isdir(dl_dir):
        for d in sorted(os.listdir(dl_dir), reverse=True):
            dpath = os.path.join(dl_dir, d)
            if not os.path.isdir(dpath) or not d.endswith("_data"):
                continue
            for f in sorted(os.listdir(dpath)):
                if f.endswith("_query_genes.csv") or f.endswith("_nodes.csv"):
                    ranked_files.append({
                        "label": f"{d}/{f}",
                        "path": os.path.join(dpath, f),
                    })
    if ranked_files:
        selected_files = st.multiselect(
            t("eml.ranked_csv_select"),
            options=[rf["label"] for rf in ranked_files],
            key="eml_ranked")
        ss.eml["ranked_paths"] = [
            rf["path"] for rf in ranked_files if rf["label"] in selected_files]
    else:
        st.caption(t("eml.no_ranked_csv"))

with col_param:
    top_k = st.number_input(t("eml.topk_label"), min_value=5, max_value=100,
                            value=config.ML_SHORTLIST_TOP_K, key="eml_topk")
    # GUI 默认值小于 CLI:网页交互场景下 2000 基因 + 20 轮会让 ML 段跑数十分钟
    # (SVM-RFE 对特征数近平方级;实测 2000×20 轮 >3 分钟未完成,500×10 轮约 1-2 分钟)
    prefilter_n = st.number_input(
        t("eml.prefilter_label"), min_value=50, max_value=5000,
        value=500, key="eml_prefilter")
    bootstrap_n = st.number_input(
        t("eml.bootstrap_label"), min_value=5, max_value=200,
        value=10, key="eml_bootstrap")
    st.caption(t("eml.perf_hint"))

st.divider()

# ---------- ③ 运行区 ----------
st.subheader(t("eml.step3"))
can_run = bool(selected_gse.strip() and case_pat.strip() and ctrl_pat.strip())
if not can_run:
    st.info(t("eml.fill_required"))

if st.button(t("eml.run_btn"), key="eml_run", type="primary",
             disabled=not can_run):
    from core.geo import fetch_geo_dataset, assign_groups
    from core.expression import preprocess, map_probe_to_symbol, differential_expression
    from core.ml_screening import deg_prefilter, run_ml_screening, build_shortlist
    import time

    warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")
    t_start = time.time()

    with st.spinner(t("eml.running")):
        try:
            # Step 1: 获取数据集
            st.markdown(t("eml.progress_fetch"))
            ds = fetch_geo_dataset(selected_gse.strip())
            if ds.issues:
                raise RuntimeError("; ".join(ds.issues))

            # Step 2: 分组
            st.markdown(t("eml.progress_group"))
            labels = assign_groups(ds, case_pat.strip(), ctrl_pat.strip())

            # Step 3: 预处理+映射+DEG
            st.markdown(t("eml.progress_deg"))
            matrix, log_branch = preprocess(ds.table)
            matrix = map_probe_to_symbol(matrix, ds.gpl)
            deg = differential_expression(matrix, labels)

            # Step 4: 预过滤+ML
            st.markdown(t("eml.progress_ml"))
            deg_top, leak_note = deg_prefilter(deg, top_n=int(prefilter_n))
            genes = deg_top["gene"].tolist()
            labeled = list(labels[labels == 1].index) + list(labels[labels == 0].index)
            X = matrix[labeled].T[genes]
            res = run_ml_screening(X, labels.loc[X.index],
                                   bootstrap_rounds=int(bootstrap_n))

            # Step 5: shortlist
            st.markdown(t("eml.progress_shortlist"))
            ranked_csvs = {}
            for path in ss.eml.get("ranked_paths", []):
                name = os.path.splitext(os.path.basename(path))[0]  # 用文件名(去后缀)作列名前缀
                ranked_csvs[name] = path
            sl, funnel = build_shortlist(res["consensus"], deg, ranked_csvs,
                                         top_k=int(top_k))

            # 保存结果到会话
            ss.eml["result"] = {
                "consensus": res["consensus"],
                "shortlist": sl,
                "funnel": funnel,
                "deg": deg,
                "meta": res["meta"],
                "roc": res["roc"],
                "gse": selected_gse.strip(),
                "elapsed": round(time.time() - t_start, 1),
                "log_branch": log_branch,
            }
            st.success(t("eml.run_success",
                         gse=selected_gse.strip(),
                         elapsed=round(time.time() - t_start, 1)))
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
st.caption(t("eml.footer"))
