import logging
import os
import re

import pandas as pd
import streamlit as st

import config
import Utils
from i18n import t
from core.network_algorithms import (
    detect_louvain_modules,
    enrich_modules,
    run_proximity_cached,
    run_rwr_cached,
)
from core.enrichment import (
    build_ranked_list,
    fetch_enrichment_gsea_prerank,
    gsea_gene_sets,
    render_gsea_plot,
)
from core.ppi import render_pyvis_html, save_gsea_results, save_module_results, save_rwr_results

st.set_page_config(page_title=t("net.page_title"), layout="wide")

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

    # Step 4.6: 排序富集（GSEA prerank，仅 names 模式）——Build Report 自动跑
    # 药物侧（pChEMBL 排秩）与疾病侧（OT 关联分排秩），CSV 由 save_gsea_results
    # 落入 <prefix>_data 目录并随后纳入 ZIP。整段 try/except 兜底：GSEA 出任何
    # 意外都不能拖垮主报告链路（主分析已有缓存，此时只会少两份 CSV）。
    gsea_results: dict = {}
    gsea_drug_df = gsea_disease_df = None
    if mode == "names":
        try:
            libs = gsea_gene_sets()
            for side_label, src_df, score_col in (
                ("drug", Utils.get_drug_targets_merged(drug_name.strip()), "score_pchembl"),
                ("disease", Utils.get_disease_targets_merged(disease_name.strip()), "score_ot"),
            ):
                rnk, rstats = build_ranked_list(src_df, "approvedSymbol", score_col)
                gsea_status: dict = {}
                gdf = fetch_enrichment_gsea_prerank(
                    rnk, libs=libs, label=f"{prefix}_{side_label}", status=gsea_status)
                gsea_results[side_label] = {
                    "df": gdf,
                    "stats": rstats,
                    "reason": gsea_status.get("reason", ""),
                }
                if side_label == "drug":
                    gsea_drug_df = gdf if not gdf.empty else None
                else:
                    gsea_disease_df = gdf if not gdf.empty else None
        except Exception:
            logging.exception("GSEA step failed during report build — skipped gracefully")
            gsea_results = {}
            gsea_drug_df = gsea_disease_df = None

    # Step 4.7: 网络邻近度检验（M3，仅 names 模式）——依赖背景互作组本地缓存
    # （未构建时返回 not_built 提示，不阻断）；结果带文件缓存（同参数秒级复用）。
    # 整段 try/except 兜底：邻近度出任何意外都不能拖垮主报告链路。
    proximity_result: dict | None = None
    if mode == "names":
        try:
            drug_df_full = Utils.get_drug_targets_merged(drug_name.strip())
            disease_df_full = Utils.get_disease_targets_merged(disease_name.strip())
            proximity_result = run_proximity_cached(
                drug_df_full["approvedSymbol"].tolist(),
                disease_df_full["approvedSymbol"].tolist(),
                label=prefix,
            )
        except Exception:
            logging.exception("Proximity step failed during report build — skipped gracefully")
            proximity_result = {"ok": False, "error": "exception", "coverage": {}, "cached": False}

    # Step 4.8: RWR 靶点优先级排序（M4，仅 names 模式）——双侧种子各跑一次
    # （幂迭代毫秒级，计算近零成本）。标注集（is_drug/is_disease）与缓存解耦，
    # 由 run_rwr_cached 在缓存命中后按当前对侧靶点集附加。
    # CSV 落盘放在 Step 5 saved 生成之后（与 GSEA 落盘并排），避免前向引用。
    rwr_results: dict = {}
    if mode == "names":
        try:
            drug_syms = {str(g).strip().upper() for g in drug_df_full["approvedSymbol"]}
            disease_syms = {str(g).strip().upper() for g in disease_df_full["approvedSymbol"]}
            rwr_results["disease_seeds"] = run_rwr_cached(
                list(disease_syms), label=f"{prefix}_disease_seeds", annotate=drug_syms)
            rwr_results["drug_seeds"] = run_rwr_cached(
                list(drug_syms), label=f"{prefix}_drug_seeds", annotate=disease_syms)
        except Exception:
            logging.exception("RWR step failed during report build — skipped gracefully")
            rwr_results = {}

    # Step 5: 保存到 DOWNLOADS_DIR（替代式：同 prefix 旧目录会被替换）
    saved = Utils.save_ppi_results(
        graph, static_fig, metrics_df, edges_df, html_content, prefix=prefix,
        enrichment_df=enrichment_df if not enrichment_df.empty else None,
        enrichment_bar_fig=enrichment_bar_fig,
        enrichment_dot_fig=enrichment_dot_fig,
    )
    # GSEA 结果 CSV 落盘进 <prefix>_data（空表不落盘），并入 saved 使 ZIP 自动包含
    gsea_saved = save_gsea_results(saved["directory"], prefix, gsea_drug_df, gsea_disease_df)
    saved.update(gsea_saved)
    # RWR 结果 CSV 落盘（空表/失败不落盘），并入 saved 使 ZIP 自动包含。
    # 注意 result 结构为 seeds（全量）+ top_nonseed（截断）两个分列表——
    # 两者合并为一张表落盘，is_seed 列区分；annotated 仅非种子行携带。
    if rwr_results:

        def _rwr_df(r: dict):
            if not r.get("ok"):
                return None
            frames = []
            if r.get("seeds"):
                frames.append(pd.DataFrame(r["seeds"]).assign(is_seed=True))
            if r.get("top_nonseed"):
                frames.append(pd.DataFrame(r["top_nonseed"]).assign(is_seed=False))
            if not frames:
                return None
            return pd.concat(frames, ignore_index=True)

        rwr_saved = save_rwr_results(
            saved["directory"], prefix,
            disease_df=_rwr_df(rwr_results.get("disease_seeds", {})),
            drug_df=_rwr_df(rwr_results.get("drug_seeds", {})),
        )
        saved.update(rwr_saved)

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
        "gsea": gsea_results,
        "proximity": proximity_result,
        "rwr": rwr_results,
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
            ppi_result = _run_ppi_analysis(
                input_mode,
                drug_name.strip(),
                disease_name.strip(),
                gene_text.strip(),
                config.SPECIES_MAP[species_label]["ncbi_taxon"],
                config.SPECIES_MAP[species_label]["gprofiler"],
                score,
                int(min_degree),
            )
            st.session_state.ppi = ppi_result
            # Build Report 已自动跑完 GSEA 双侧 → 直接作为交互区的初始结果
            # （按侧存储；缓存键相同，交互区重跑同侧也是秒级命中）。
            st.session_state.ppi_gsea = {
                side: {**info, "ppi_params": ppi_result["params"]}
                for side, info in ppi_result.get("gsea", {}).items()
            }
        st.rerun()

# ---------- 结果展示（读取 session_state，不重新分析） ----------
ppi = st.session_state.ppi
if ppi is not None:
    if current_params != ppi["params"]:
        st.info(t("net.params_changed"))

    genes = ppi["genes"]
    edges_df = ppi["edges_df"]
    graph = ppi["graph"]  # 模块检测与按模块着色渲染需要原图对象
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

        # ---------- 网络邻近度检验（M3，仅 names 模式） ----------
        # 随主分析自动计算（文件缓存后秒级复用）；genes 模式无药物/疾病分侧，不展示。
        proximity_result = ppi.get("proximity")
        if proximity_result is not None and input_mode == "names":
            st.divider()
            st.subheader(t("net.prox_subheader"))
            st.caption(t("net.prox_caption", n_perm=config.PROXIMITY_N_PERM,
                         seed=config.PROXIMITY_SEED,
                         bin_size=config.PROXIMITY_DEGREE_BIN_SIZE))
            prox_err = proximity_result.get("error")
            if prox_err == "not_built":
                st.warning(t("net.prox_not_built"))
            elif prox_err in ("low_coverage_drug", "low_coverage_disease", "empty_after_filter"):
                cov = proximity_result.get("coverage", {})
                st.warning(t("net.prox_low_coverage",
                             d_in=cov.get("drug_in_lcc", 0), d_tot=cov.get("drug_total", 0),
                             s_in=cov.get("disease_in_lcc", 0), s_tot=cov.get("disease_total", 0),
                             thr=config.PROXIMITY_MIN_COVERAGE))
            elif prox_err:
                st.error(t("net.prox_fail", reason=prox_err))
            elif proximity_result.get("ok"):
                col_p1, col_p2, col_p3 = st.columns(3)
                col_p1.metric(t("net.prox_dc"), f"{proximity_result['d_c']:.3f}")
                col_p2.metric(t("net.prox_z"), f"{proximity_result['z']:.3f}")
                col_p3.metric(t("net.prox_p"), f"{proximity_result['p_empirical']:.4f}")
                if proximity_result["z"] < 0:
                    st.success(t("net.prox_support"))
                else:
                    st.info(t("net.prox_nosupport"))
                cov = proximity_result.get("coverage", {})
                st.caption(t("net.prox_coverage",
                             ver=config.STRING_VERSION, conf=config.INTERACTOME_CONFIDENCE,
                             lcc=proximity_result.get("lcc_nodes"),
                             d_in=cov.get("drug_in_lcc", 0), d_tot=cov.get("drug_total", 0),
                             s_in=cov.get("disease_in_lcc", 0), s_tot=cov.get("disease_total", 0),
                             n_perm=proximity_result.get("n_perm"),
                             cached=t("net.prox_cached") if proximity_result.get("cached") else ""))
                # 去交集版本：剔除共同靶点后仍显著，才能排除"纯属直接作用"的解释
                if proximity_result.get("d_c_excl") is not None:
                    st.caption(t("net.prox_excl_caption", n=proximity_result.get("n_excl") or 0))
                    col_q1, col_q2, col_q3 = st.columns(3)
                    col_q1.metric(t("net.prox_dc_excl"), f"{proximity_result['d_c_excl']:.3f}")
                    col_q2.metric(t("net.prox_z_excl"),
                                  f"{proximity_result['z_excl']:.3f}"
                                  if proximity_result.get("z_excl") is not None else "—")
                    col_q3.metric(t("net.prox_p_excl"),
                                  f"{proximity_result['p_excl']:.4f}"
                                  if proximity_result.get("p_excl") is not None else "—")
                elif proximity_result.get("ok"):
                    st.caption(t("net.prox_excl_unavailable"))

        # ---------- 靶点优先级排序（RWR，M4，仅 names 模式） ----------
        # 随主分析自动计算双侧（种子=疾病 / 种子=药物）；radio 仅切换展示，无需重算。
        rwr_results = ppi.get("rwr") or {}
        if rwr_results and input_mode == "names":
            st.divider()
            with st.expander(t("net.rwr_subheader")):
                st.caption(t("net.rwr_caption", alpha=config.RWR_ALPHA))
                rwr_side = st.radio(
                    t("net.rwr_side"), ["disease_seeds", "drug_seeds"],
                    format_func=lambda v: (t("net.rwr_side_disease") if v == "disease_seeds"
                                           else t("net.rwr_side_drug")),
                    horizontal=True, key="rwr_side",
                )
                rwr_res = rwr_results.get(rwr_side) or {}
                rwr_err = rwr_res.get("error")
                if rwr_err == "not_built":
                    st.warning(t("net.rwr_not_built"))
                elif rwr_err:
                    st.error(t("net.rwr_fail", reason=rwr_err))
                elif rwr_res.get("ok"):
                    cov = rwr_res.get("coverage", {})
                    conv_txt = t("net.rwr_conv_yes") if rwr_res.get("converged") else t("net.rwr_conv_no")
                    anno_label = (t("net.rwr_col_is_drug") if rwr_side == "disease_seeds"
                                  else t("net.rwr_col_is_disease"))
                    st.caption(t("net.rwr_stats",
                                 alpha=rwr_res.get("alpha"), iters=rwr_res.get("n_iter"),
                                 conv=conv_txt,
                                 n_seed=cov.get("seed_in_lcc", 0), n_tot=cov.get("seed_total", 0),
                                 lcc=rwr_res.get("lcc_nodes"), mass=rwr_res.get("seed_mass", 0.0),
                                 cached=t("net.rwr_cached") if rwr_res.get("cached") else ""))
                    # 度匹配校正：标注集平均得分 vs 度匹配随机集（排除枢纽度数假象）
                    enr = rwr_res.get("enrichment")
                    if enr and enr.get("z") is not None:
                        if enr.get("p") is not None and enr["p"] < 0.05:
                            st.success(t("net.rwr_enrich_sig",
                                         n=enr.get("n_annotated", 0), z=enr["z"], p=enr["p"]))
                        else:
                            st.info(t("net.rwr_enrich_nosig",
                                      n=enr.get("n_annotated", 0), z=enr["z"], p=enr["p"]))
                    top_df = pd.concat([
                        pd.DataFrame(rwr_res.get("seeds", [])).assign(is_seed=True),
                        pd.DataFrame(rwr_res.get("top_nonseed", [])).assign(is_seed=False),
                    ], ignore_index=True) if (rwr_res.get("seeds") or rwr_res.get("top_nonseed")) else pd.DataFrame()
                    if top_df.empty:
                        st.info(t("net.rwr_none"))
                    else:
                        top_df = top_df.sort_values("score", ascending=False).reset_index(drop=True)
                        seeds_df = top_df[top_df["is_seed"]]
                        nonseed_df = top_df[~top_df["is_seed"]].head(config.RWR_TOP_N)
                        st.markdown(t("net.rwr_top_header", n=len(nonseed_df)))
                        show_df = nonseed_df[["gene", "score"]].copy()
                        if "annotated" in nonseed_df.columns:
                            show_df[anno_label] = nonseed_df["annotated"]
                        st.dataframe(show_df, width="stretch", hide_index=True)
                        with st.expander(t("net.rwr_seeds_header", n=len(seeds_df))):
                            st.dataframe(seeds_df[["gene", "score"]],
                                         width="stretch", hide_index=True)
                        st.download_button(
                            t("net.dl_rwr_csv"),
                            top_df.to_csv(index=False).encode("utf-8-sig"),
                            f"{prefix}_rwr_{rwr_side}.csv", "text/csv", width="stretch",
                        )
                        if "rwr_disease_csv" in saved or "rwr_drug_csv" in saved:
                            st.caption(t("net.rwr_in_report", dir=saved["directory"]))

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

        # ---------- 排序富集（GSEA prerank，仅智人） ----------
        # Build Report 已自动跑双侧并落盘 ZIP；此交互区可换参数/换侧重跑。
        # 会话状态按"侧"存储（ppi_gsea[side]）：切换药物/疾病单选只会显示
        # 对应侧的结果，不会误展示另一侧的旧结果。
        st.divider()
        with st.expander(t("net.gsea_subheader")):
            st.caption(t("net.gsea_caption", libs="、".join(config.GSEA_ENRICHR_LIBS)))
            # 按侧存储；兼容清理旧版（单结果字典含 "side" 键）的会话残留
            if ("ppi_gsea" not in st.session_state
                    or not isinstance(st.session_state.ppi_gsea, dict)
                    or "side" in st.session_state.ppi_gsea):
                st.session_state.ppi_gsea = {}
            if "gsea_drug_csv" in saved or "gsea_disease_csv" in saved:
                st.caption(t("net.gsea_in_report", dir=saved["directory"]))

            if input_mode == "genes":
                st.info(t("net.gsea_manual_unsupported"))
            else:
                gsea_side = st.radio(
                    t("net.gsea_side"), ["drug", "disease"],
                    format_func=lambda v: (t("net.gsea_side_drug") if v == "drug"
                                           else t("net.gsea_side_disease")),
                    horizontal=True, key="gsea_side",
                )
                if st.button(t("net.gsea_run"), key="gsea_run_btn", width="stretch"):
                    if current_params != ppi["params"]:
                        # 防呆：主分析参数已变，先重跑主分析再 GSEA
                        st.warning(t("net.params_changed"))
                    else:
                        with st.spinner(t("net.gsea_spinner")):
                            if gsea_side == "drug":
                                src_df = Utils.get_drug_targets_merged(drug_name.strip())
                                score_col = "score_pchembl"
                            else:
                                src_df = Utils.get_disease_targets_merged(disease_name.strip())
                                score_col = "score_ot"
                            rnk, rstats = build_ranked_list(src_df, "approvedSymbol", score_col)
                            gsea_status: dict = {}
                            gsea_df = fetch_enrichment_gsea_prerank(
                                rnk,
                                libs=gsea_gene_sets(),
                                label=f"{prefix}_{gsea_side}",
                                status=gsea_status,
                            )
                            st.session_state.ppi_gsea[gsea_side] = {
                                "df": gsea_df,
                                "stats": rstats,
                                "reason": gsea_status.get("reason", ""),
                                "ppi_params": ppi["params"],
                            }
                        st.rerun()

                gsea_res = st.session_state.ppi_gsea.get(gsea_side)
                if gsea_res is not None:
                    if gsea_res["ppi_params"] != ppi["params"]:
                        st.info(t("net.params_changed"))
                    else:
                        stats = gsea_res["stats"]
                        gsea_df = gsea_res["df"]
                        reason = gsea_res.get("reason", "")
                        side_txt = (t("net.gsea_side_drug") if gsea_side == "drug"
                                    else t("net.gsea_side_disease"))
                        # 打分缺失去除是多源合并的预期行为（如仅 DGIdb 收录的靶点无
                        # pChEMBL），作中性说明而非预警
                        if stats.get("n_dropped", 0) > 0:
                            st.caption(t("net.gsea_drop_note",
                                         n=stats.get("n_dropped", 0),
                                         m=stats.get("n_input", 0)))
                        if reason == "insufficient_input":
                            st.info(t("net.gsea_insufficient",
                                      n=stats.get("n_kept", 0), min_n=config.GSEA_MIN_RANKED_GENES))
                        elif reason == "request_failed":
                            st.error(t("net.gsea_fail", reason=reason))
                        elif gsea_df is None or gsea_df.empty:
                            st.info(t("net.gsea_none", thr=config.GSEA_FDR_SCREEN))
                        else:
                            st.caption(t("net.gsea_result_side", side=side_txt,
                                         n_in=stats.get("n_input", 0), n_kept=stats.get("n_kept", 0)))
                            screened = gsea_df[gsea_df["fdr_q_value"] < config.GSEA_FDR_SCREEN]
                            if screened.empty:
                                st.info(t("net.gsea_none", thr=config.GSEA_FDR_SCREEN))
                            else:
                                st.caption(t("net.gsea_table_caption", n=len(gsea_df),
                                             m=len(screened), thr=config.GSEA_FDR_SCREEN))
                                show_cols = ["source", "term_name", "nes", "fdr_q_value",
                                             "nominal_p_value", "tag_pct", "gene_pct"]
                                st.dataframe(screened[show_cols], width="stretch", hide_index=True)
                                gsea_fig = render_gsea_plot(screened)
                                if gsea_fig is not None:
                                    st.caption(t("net.gsea_bar_caption"))
                                    st.pyplot(gsea_fig)
                            st.download_button(
                                t("net.dl_gsea_csv"),
                                gsea_df.to_csv(index=False).encode("utf-8-sig"),
                                f"{prefix}_gsea_{gsea_side}.csv", "text/csv",
                                width="stretch",
                            )

        # ---------- 蛋白模块检测（Louvain + 模块富集） ----------
        # 基于 session_state 中已构建的 PPI 图，独立按钮触发；结果连同其依赖的
        # 主分析参数指纹一起存 session_state——主分析参数变更后模块结果自动标记过期。
        st.divider()
        with st.expander(t("net.modules_subheader")):
            st.caption(t("net.modules_caption",
                         min_size=config.LOUVAIN_MIN_MODULE_SIZE,
                         top_n=config.LOUVAIN_MAX_ENRICH_MODULES))
            if "ppi_modules" not in st.session_state:
                st.session_state.ppi_modules = None

            if st.button(t("net.modules_run"), key="modules_run_btn", width="stretch"):
                if current_params != ppi["params"]:
                    # 防呆：主分析参数已变，图与旧结果不对应，先重跑主分析
                    st.warning(t("net.params_changed"))
                else:
                    with st.spinner(t("net.modules_spinner_detect")):
                        modules = detect_louvain_modules(
                            graph,
                            resolution=config.LOUVAIN_RESOLUTION,
                            seed=config.LOUVAIN_SEED,
                        )
                        summary_df = pd.DataFrame([
                            {
                                "module_id": m["module_id"],
                                "size": m["size"],
                                "avg_degree": m["avg_degree"],
                                "genes_preview": ", ".join(m["genes"][:8]) + ("…" if m["size"] > 8 else ""),
                            }
                            for m in modules if m["size"] >= 2
                        ])
                        # 仅多基因模块参与着色；单点模块节点渲染为浅灰
                        module_map = {
                            gene: m["module_id"]
                            for m in modules if m["size"] >= 2
                            for gene in m["genes"]
                        }
                        # params[5] = organism（见上方 current_params 的元组顺序）
                        long_df, failures, no_sig = enrich_modules(
                            modules,
                            organism=current_params[5],
                            min_module_size=config.LOUVAIN_MIN_MODULE_SIZE,
                            top_n_modules=config.LOUVAIN_MAX_ENRICH_MODULES,
                        )
                        # 结果落盘到本次分析目录（<prefix>_data），空表不落盘；
                        # saved 即 session_state 中的字典，追加后页面 ZIP 自动包含模块文件，
                        # Data Management 的整目录 zip 也会一并打包。
                        mod_saved = save_module_results(
                            saved["directory"], prefix, summary_df, long_df)
                        saved.update(mod_saved)
                        st.session_state.ppi_modules = {
                            "modules": modules,
                            "summary_df": summary_df,
                            "module_map": module_map,
                            "long_df": long_df,
                            "failures": failures,
                            "no_sig": no_sig,
                            "saved_paths": mod_saved,
                            "ppi_params": ppi["params"],
                        }
                    st.rerun()

            mods_res = st.session_state.ppi_modules
            if mods_res is not None:
                if mods_res["ppi_params"] != ppi["params"]:
                    st.info(t("net.params_changed"))
                else:
                    summary_df = mods_res["summary_df"]
                    if summary_df.empty:
                        st.info(t("net.modules_none"))
                    else:
                        # 退化防护：最大模块吞掉大多数节点时提示调参
                        total_nodes = len(metrics_df)
                        if total_nodes and mods_res["modules"]:
                            ratio = mods_res["modules"][0]["size"] / total_nodes
                            if ratio > config.LOUVAIN_MAX_MODULE_RATIO:
                                st.warning(t("net.modules_degenerate", ratio=ratio))
                        n_singles = sum(1 for m in mods_res["modules"] if m["size"] < 2)
                        st.markdown(t("net.modules_summary", n_singles=n_singles))
                        st.dataframe(summary_df, width="stretch", hide_index=True)

                        # 按模块着色的网络视图（默认关，避免重复渲染大图）
                        if st.toggle(t("net.modules_color_toggle"), key="modules_color_toggle"):
                            colored_html = render_pyvis_html(
                                graph, metrics_df, module_map=mods_res["module_map"])
                            st.iframe(colored_html, height=620)

                    # ---- 模块富集结果（失败→错误提示；无显著→一般提示） ----
                    long_df = mods_res["long_df"]
                    failures = mods_res["failures"]
                    no_sig = mods_res["no_sig"]
                    if failures:
                        ids = ", ".join(f"#{f['module_id']}" for f in failures)
                        st.error(t("net.modules_fail", ids=ids))
                    if no_sig:
                        ids = ", ".join(f"#{mid}" for mid in no_sig)
                        st.info(t("net.modules_nosig", ids=ids))
                    if long_df is not None and not long_df.empty:
                        k = config.LOUVAIN_TOP_TERMS_PER_MODULE
                        show_df = long_df.groupby("module_id", group_keys=False).head(k)
                        st.markdown(t("net.modules_enrich_subheader"))
                        st.caption(t("net.modules_enrich_caption",
                                     n=len(long_df), m=long_df["module_id"].nunique(), k=k))
                        cols = ["module_id", "source", "term_id", "term_name",
                                "adjusted_p_value", "intersection_size", "query_size"]
                        st.dataframe(show_df[cols], width="stretch", hide_index=True)
                        col_ms, col_me = st.columns(2)
                        with col_ms:
                            st.download_button(
                                t("net.dl_module_summary_csv"),
                                summary_df.to_csv(index=False).encode("utf-8-sig"),
                                f"{prefix}_modules.csv", "text/csv", width="stretch",
                            )
                        with col_me:
                            st.download_button(
                                t("net.dl_module_enrich_csv"),
                                long_df.to_csv(index=False).encode("utf-8-sig"),
                                f"{prefix}_module_enrichment.csv", "text/csv", width="stretch",
                            )
                        if mods_res.get("saved_paths"):
                            st.caption(t("net.modules_saved", dir=saved["directory"]))
                    elif not failures and not no_sig and not summary_df.empty:
                        st.info(t("net.modules_enrich_empty"))

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
