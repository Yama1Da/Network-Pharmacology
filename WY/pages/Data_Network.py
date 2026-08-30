import os
import re

import pandas as pd
import streamlit as st

import config
import Utils
from i18n import t, render_lang_toggle

st.set_page_config(page_title=t("net.page_title"), layout="wide")
render_lang_toggle()

# ---------- Hide auto navigation sidebar for full-screen display ----------
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        display: none;
    }
    /* Run button colors (#167A9C bg, #C6E7F5 text) */
    .st-key-run_btn button {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
    }
    /* Data Management button colors (#6347A1 bg, #E0D8F2 text) */
    .st-key-data_mgmt_btn button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }
</style>
""", unsafe_allow_html=True)

st.title(t("net.title"))
st.caption(t("net.caption"))

col_back, col_mgmt = st.columns([1, 1])
with col_back:
    if st.button(t("common.back_main"), key="back_btn"):
        st.switch_page("App.py")
with col_mgmt:
    if st.button(t("net.data_mgmt"), key="data_mgmt_btn"):
        st.switch_page("pages/Data_Management.py")

st.divider()

# ---------- 输入区 ----------
input_mode = st.radio(
    t("net.input_mode"),
    ["names", "genes"],
    format_func=lambda v: t("net.mode_names") if v == "names" else t("net.mode_genes"),
    horizontal=True,
)

col_input, col_params = st.columns([2, 1])
with col_input:
    if input_mode == "names":
        col_drug, col_dis = st.columns(2)
        with col_drug:
            drug_name = st.text_input(t("net.drug_label"), value="aspirin")
        with col_dis:
            disease_name = st.text_input(t("net.disease_label"), value="cancer")
        gene_text = ""
    else:
        drug_name, disease_name = "", ""
        gene_text = st.text_area(
            t("net.gene_text_label"),
            placeholder="TP53\nBRCA1\nEGFR",
            height=140,
        )

with col_params:
    # 物种选择：当前仅开放“人类（Human）”。
    # 多物种支持待后续优化后再开放（届时把下面选项恢复为
    # list(config.SPECIES_MAP.keys())，并同步处理疾病靶点数据源
    # Open Targets 只有人类数据的问题）。
    # 注意：选项值 "Human" 同时是 config.SPECIES_MAP 的键，不参与翻译。
    species_label = st.selectbox(t("net.species"), ["Human"])
    score = st.slider(t("net.confidence"), 0.1, 0.9, config.PPI_CONFIDENCE, 0.05)
    min_degree = st.number_input(t("net.min_degree"), 1, 20, config.MIN_DEGREE_FILTER, 1)

run_btn = st.button(t("net.build_report"), key="run_btn", type="primary", width="stretch")


# ---------- 分析执行 ----------
def _run_ppi_analysis(mode, drug_name, disease_name, gene_text, species_id, organism, score, min_degree) -> dict:
    """执行完整分析（交集→STRING→建图→指标→Hub→富集→可视化→保存），返回结果字典。"""
    # ============ 物种策略说明（待优化时确认） ============
    # 当前疾病靶点仅来自 Open Targets（只有人类数据），因此整个分析
    # 链路按【人类】处理；界面侧也暂时只开放 Human 选项。
    # 后续优化时再决定：是否定死人类，或放开多物种（Mouse/Rat）。
    # 若放开多物种，需同步解决疾病靶点数据源的多物种支持
    # （如接入 DisGeNET / MGI / RGD 等，并处理人-鼠同源映射）。
    # =====================================================
    # Step 0: 确定输入基因列表
    if mode == "names":
        genes = Utils.get_common_targets(drug_name, disease_name, top_n=config.STRING_QUERY_LIMIT)
        # 结果目录/文件前缀：药物名_疾病名（save_ppi_results 会生成 <prefix>_data 目录）
        prefix = re.sub(r"[^\w\-]+", "_", f"{drug_name}_{disease_name}")
    else:
        tokens = re.split(r"[\s,;]+", gene_text)
        genes = sorted({t.strip().upper() for t in tokens if t.strip()})
        prefix = "manual"

    # Step 1: STRING 获取相互作用（缓存命名含参数，如 StringPPI_aspirin_cancer_0.7_3.csv）
    edges_df = Utils.fetch_string_ppi(genes, species=species_id, score=score,
                                      label=prefix, min_degree=int(min_degree))

    # Step 2-4: 建图 + 拓扑指标 + Hub 筛选 + 可视化
    graph = Utils.build_ppi_graph(edges_df, genes)
    metrics_df = Utils.compute_ppi_metrics(graph)
    metrics_df = Utils.mark_hub_genes(metrics_df, min_degree=min_degree)
    hub_df = metrics_df[metrics_df["is_hub"]]
    html_content = Utils.render_pyvis_html(graph, metrics_df)
    static_fig = Utils.render_static_plot(graph, metrics_df)

    # Step 4.5: 富集分析（g:Profiler 主 + Enrichr 备，输入为交集靶点）
    enrichment_df, enrichment_source = Utils.run_enrichment(
        genes, organism=organism, pvalue=config.ENRICHMENT_PVALUE, label=prefix
    )
    enrichment_bar_fig = enrichment_dot_fig = None
    if not enrichment_df.empty:
        enrichment_bar_fig, enrichment_dot_fig = Utils.render_enrichment_plots(enrichment_df)

    # Step 5: 保存到 DOWNLOADS_DIR（替代式：同 prefix 旧目录会被替换）
    saved = Utils.save_ppi_results(
        graph, static_fig, metrics_df, edges_df, html_content, prefix=prefix,
        enrichment_df=enrichment_df if not enrichment_df.empty else None,
        enrichment_bar_fig=enrichment_bar_fig,
        enrichment_dot_fig=enrichment_dot_fig,
    )

    return {
        "genes": genes,
        "edges_df": edges_df,
        "graph": graph,
        "metrics_df": metrics_df,
        "hub_df": hub_df,
        "html_content": html_content,
        "static_fig": static_fig,
        "enrichment_df": enrichment_df,
        "enrichment_source": enrichment_source,
        "enrichment_bar_fig": enrichment_bar_fig,
        "enrichment_dot_fig": enrichment_dot_fig,
        "saved": saved,
        "prefix": prefix,
        "params": (mode, drug_name, disease_name, gene_text, species_id, organism, score, min_degree),
    }


# ---------- 结果持久化：存入 session_state，点击下载等触发的 rerun 不再重跑分析 ----------
if "ppi" not in st.session_state:
    st.session_state.ppi = None

current_params = (
    input_mode,
    drug_name.strip(),
    disease_name.strip(),
    gene_text.strip(),
    config.SPECIES_MAP[species_label]["ncbi_taxon"],
    config.SPECIES_MAP[species_label]["gprofiler"],
    score,
    int(min_degree),
)

if run_btn:
    if input_mode == "names" and (not drug_name.strip() or not disease_name.strip()):
        st.warning(t("net.warn_missing"))
    else:
        with st.spinner(t("net.spinner_build")):
            st.session_state.ppi = _run_ppi_analysis(
                input_mode,
                drug_name.strip(),
                disease_name.strip(),
                gene_text.strip(),
                config.SPECIES_MAP[species_label]["ncbi_taxon"],
                config.SPECIES_MAP[species_label]["gprofiler"],
                score,
                int(min_degree),
            )
        st.rerun()

# ---------- 结果展示（读取 session_state，不重新分析） ----------
ppi = st.session_state.ppi
if ppi is not None:
    if current_params != ppi["params"]:
        st.info(t("net.params_changed"))

    genes = ppi["genes"]
    edges_df = ppi["edges_df"]
    metrics_df = ppi["metrics_df"]
    hub_df = ppi["hub_df"]
    html_content = ppi["html_content"]
    static_fig = ppi["static_fig"]
    saved = ppi["saved"]
    prefix = ppi["prefix"]

    if not genes:
        st.warning(t("net.no_common"))
    elif edges_df.empty:
        st.error(t("net.no_interactions"))
    else:
        st.info(t("net.input_genes", n=len(genes)))

        col_m1, col_m2, col_m3, col_m4 = st.columns(4)
        col_m1.metric(t("net.metric_input"), len(genes))
        col_m2.metric(t("net.metric_nodes"), len(metrics_df))
        col_m3.metric(t("net.metric_edges"), len(edges_df))
        col_m4.metric(t("net.metric_hubs"), len(hub_df))

        st.divider()

        st.subheader(t("net.hub_subheader", d=int(min_degree)))
        if hub_df.empty:
            st.info(t("net.no_hubs"))
        else:
            st.dataframe(hub_df, width="stretch", hide_index=True)

        st.subheader(t("net.interactive"))
        st.iframe(html_content, height=640)

        st.subheader(t("net.static"))
        st.pyplot(static_fig)

        # ---------- 富集分析展示 ----------
        st.divider()
        st.subheader(t("net.enrich_subheader"))
        enrichment_df = ppi["enrichment_df"]
        enrichment_source = ppi["enrichment_source"]
        enrichment_bar_fig = ppi["enrichment_bar_fig"]
        enrichment_dot_fig = ppi["enrichment_dot_fig"]
        if enrichment_df is None or enrichment_df.empty:
            st.info(t("net.no_enrich"))
        else:
            st.caption(
                t("net.enrich_caption", src=enrichment_source, n=len(enrichment_df), m=len(genes))
            )
            display_cols = ["source", "term_id", "term_name", "adjusted_p_value",
                            "intersection_size", "query_size"]
            st.dataframe(enrichment_df[display_cols], width="stretch", hide_index=True)
            if enrichment_bar_fig is not None:
                st.caption(t("net.bar_caption"))
                st.pyplot(enrichment_bar_fig)
            if enrichment_dot_fig is not None:
                st.caption(t("net.bubble_caption"))
                st.pyplot(enrichment_dot_fig)

        # ---------- 下载区 ----------
        with st.expander(t("net.download_expander")):
            col_d1, col_d2, col_d3 = st.columns(3)
            with col_d1:
                st.download_button(
                    t("net.dl_nodes"),
                    metrics_df.to_csv(index=False).encode("utf-8-sig"),
                    f"{prefix}_nodes.csv", "text/csv", width="stretch",
                )
            with col_d2:
                st.download_button(
                    t("net.dl_edges"),
                    edges_df.to_csv(index=False).encode("utf-8-sig"),
                    f"{prefix}_edges.csv", "text/csv", width="stretch",
                )
            with col_d3:
                st.download_button(
                    t("net.dl_html"),
                    html_content,
                    f"{prefix}_network.html", "text/html", width="stretch",
                )
            col_d4, col_d5, col_d6 = st.columns(3)
            with col_d4:
                st.download_button(
                    t("net.dl_graphml"),
                    open(saved["graphml"], "rb").read(),
                    f"{prefix}_graph.graphml", "application/xml", width="stretch",
                )
            with col_d5:
                st.download_button(
                    t("net.dl_svg"),
                    open(saved["svg"], "rb").read(),
                    f"{prefix}_network.svg", "image/svg+xml", width="stretch",
                )
            with col_d6:
                # 富集气泡图（与 ZIP 换位到第二行第三列；有富集结果时显示）
                if "enrichment_csv" in saved:
                    st.download_button(
                        t("net.dl_bubble"),
                        open(saved["enrichment_dot_png"], "rb").read(),
                        f"{prefix}_enrichment_dot.png", "image/png", width="stretch",
                    )
            # 第三行：富集 CSV / Bar（有富集结果时）+ 右下角 ZIP（始终显示）
            col_e1, col_e2, col_e3 = st.columns(3)
            if "enrichment_csv" in saved:
                with col_e1:
                    st.download_button(
                        t("net.dl_enrich_csv"),
                        open(saved["enrichment_csv"], "rb").read(),
                        f"{prefix}_enrichment.csv", "text/csv", width="stretch",
                    )
                with col_e2:
                    st.download_button(
                        t("net.dl_enrich_bar"),
                        open(saved["enrichment_bar_png"], "rb").read(),
                        f"{prefix}_enrichment_bar.png", "image/png", width="stretch",
                    )
            with col_e3:
                st.download_button(
                    t("net.dl_zip"),
                    Utils.make_ppi_zip(saved),
                    f"{prefix}_results.zip", "application/zip", width="stretch",
                )
            st.caption(
                t("net.files_grouped", dir=saved["directory"])
            )

        st.success(t("net.success"))
