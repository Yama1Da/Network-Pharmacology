import os
import re

import pandas as pd
import streamlit as st

import config
import Utils

st.set_page_config(page_title="PPI Network Analysis", layout="wide")

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
</style>
""", unsafe_allow_html=True)

st.title("🧬 PPI Network Analysis")
st.caption("Drug ∩ Disease targets → STRING interactions → topology analysis → network visualization")

col_back, _ = st.columns([1, 4])
with col_back:
    if st.button("← Back to Main", key="back_btn"):
        st.switch_page("App.py")

st.divider()

# ---------- 输入区 ----------
input_mode = st.radio(
    "Input mode",
    ["Drug + Disease names", "Paste gene list"],
    horizontal=True,
)

col_input, col_params = st.columns([2, 1])
with col_input:
    if input_mode == "Drug + Disease names":
        col_drug, col_dis = st.columns(2)
        with col_drug:
            drug_name = st.text_input("Drug name", value="aspirin")
        with col_dis:
            disease_name = st.text_input("Disease name", value="cancer")
        gene_text = ""
    else:
        drug_name, disease_name = "", ""
        gene_text = st.text_area(
            "Gene symbols (one per line or comma-separated)",
            placeholder="TP53\nBRCA1\nEGFR",
            height=140,
        )

with col_params:
    species_label = st.selectbox("Species", list(config.SPECIES_MAP.keys()))
    score = st.slider("STRING confidence", 0.1, 0.9, config.PPI_CONFIDENCE, 0.05)
    min_degree = st.number_input("Min degree (hub)", 1, 20, config.MIN_DEGREE_FILTER, 1)

run_btn = st.button("Build PPI Network", key="run_btn", type="primary", width="stretch")


# ---------- 分析执行 ----------
def _run_ppi_analysis(mode, drug_name, disease_name, gene_text, species_id, score, min_degree) -> dict:
    """执行完整 PPI 分析（交集→STRING→建图→指标→Hub→可视化→保存），返回结果字典。"""
    # Step 0: 确定输入基因列表
    if mode == "Drug + Disease names":
        genes = Utils.get_common_targets(drug_name, disease_name, top_n=config.STRING_QUERY_LIMIT)
        prefix = re.sub(r"[^\w\-]+", "_", f"PPI_{drug_name}_{disease_name}")
    else:
        tokens = re.split(r"[\s,;]+", gene_text)
        genes = sorted({t.strip().upper() for t in tokens if t.strip()})
        prefix = "PPI_manual"

    # Step 1: STRING 获取相互作用
    edges_df = Utils.fetch_string_ppi(genes, species=species_id, score=score)

    # Step 2-4: 建图 + 拓扑指标 + Hub 筛选 + 可视化
    graph = Utils.build_ppi_graph(edges_df, genes)
    metrics_df = Utils.compute_ppi_metrics(graph)
    metrics_df = Utils.mark_hub_genes(metrics_df, min_degree=min_degree)
    hub_df = metrics_df[metrics_df["is_hub"]]
    html_content = Utils.render_pyvis_html(graph, metrics_df)
    static_fig = Utils.render_static_plot(graph, metrics_df)

    # Step 5: 保存到 DOWNLOADS_DIR（替代式：同 prefix 旧目录会被替换）
    saved = Utils.save_ppi_results(graph, static_fig, metrics_df, edges_df, html_content, prefix=prefix)

    return {
        "genes": genes,
        "edges_df": edges_df,
        "graph": graph,
        "metrics_df": metrics_df,
        "hub_df": hub_df,
        "html_content": html_content,
        "static_fig": static_fig,
        "saved": saved,
        "prefix": prefix,
        "params": (mode, drug_name, disease_name, gene_text, species_id, score, min_degree),
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
    score,
    int(min_degree),
)

if run_btn:
    if input_mode == "Drug + Disease names" and (not drug_name.strip() or not disease_name.strip()):
        st.warning("Please enter both drug and disease names")
    else:
        with st.spinner("Building PPI network…"):
            st.session_state.ppi = _run_ppi_analysis(
                input_mode,
                drug_name.strip(),
                disease_name.strip(),
                gene_text.strip(),
                config.SPECIES_MAP[species_label]["ncbi_taxon"],
                score,
                int(min_degree),
            )
        st.rerun()

# ---------- 结果展示（读取 session_state，不重新分析） ----------
ppi = st.session_state.ppi
if ppi is not None:
    if current_params != ppi["params"]:
        st.info("Parameters have changed — click Build PPI Network to re-analyze.")

    genes = ppi["genes"]
    edges_df = ppi["edges_df"]
    metrics_df = ppi["metrics_df"]
    hub_df = ppi["hub_df"]
    html_content = ppi["html_content"]
    static_fig = ppi["static_fig"]
    saved = ppi["saved"]
    prefix = ppi["prefix"]

    if not genes:
        st.warning("No common targets found between the two datasets.")
    elif edges_df.empty:
        st.error("STRING returned no interactions above the confidence threshold.")
    else:
        st.info(f"Input: {len(genes)} target genes")

        col_m1, col_m2, col_m3, col_m4 = st.columns(4)
        col_m1.metric("Input targets", len(genes))
        col_m2.metric("Network nodes", len(metrics_df))
        col_m3.metric("Edges", len(edges_df))
        col_m4.metric("Hub genes", len(hub_df))

        st.divider()

        st.subheader(f"⭐ Hub genes (degree > median, degree ≥ {int(min_degree)})")
        if hub_df.empty:
            st.info("No hub genes met the threshold.")
        else:
            st.dataframe(hub_df, width="stretch", hide_index=True)

        st.subheader("🕸️ Interactive network (drag / zoom)")
        st.iframe(html_content, height=640)

        st.subheader("📊 Static network")
        st.pyplot(static_fig)

        # ---------- 下载区 ----------
        with st.expander("Download results"):
            col_d1, col_d2, col_d3 = st.columns(3)
            with col_d1:
                st.download_button(
                    "⬇️ Nodes CSV",
                    metrics_df.to_csv(index=False).encode("utf-8-sig"),
                    f"{prefix}_nodes.csv", "text/csv", width="stretch",
                )
            with col_d2:
                st.download_button(
                    "⬇️ Edges CSV",
                    edges_df.to_csv(index=False).encode("utf-8-sig"),
                    f"{prefix}_edges.csv", "text/csv", width="stretch",
                )
            with col_d3:
                st.download_button(
                    "⬇️ Network HTML",
                    html_content,
                    f"{prefix}_network.html", "text/html", width="stretch",
                )
            col_d4, col_d5, col_d6 = st.columns(3)
            with col_d4:
                st.download_button(
                    "⬇️ GraphML (Cytoscape)",
                    open(saved["graphml"], "rb").read(),
                    f"{prefix}_graph.graphml", "application/xml", width="stretch",
                )
            with col_d5:
                st.download_button(
                    "⬇️ SVG (vector)",
                    open(saved["svg"], "rb").read(),
                    f"{prefix}_network.svg", "image/svg+xml", width="stretch",
                )
            with col_d6:
                st.download_button(
                    "⬇️ All results (ZIP)",
                    Utils.make_ppi_zip(saved),
                    f"{prefix}_results.zip", "application/zip", width="stretch",
                )
            st.caption(
                f"All files of this analysis are grouped in: `{saved['directory']}`"
            )

        st.success("✅ PPI network analysis complete")
