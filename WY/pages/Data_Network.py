import logging
import re

import pandas as pd
import streamlit as st

from WY import config
from WY.core.files import make_fig_zip, make_ppi_zip
from WY.ui.widgets import (back_button, download_button,
                        download_file_button, enter_busy, finish_busy, hide_sidebar,
                        is_busy, show_busy_notice)
from WY.i18n import t
from WY.core.enrichment import (
    enrichment_retry_hint_minutes,
    render_gsea_plot,
)
from WY.core.ppi import render_pyvis_html
from WY.core.analysis import (compute_ppi_modules, gsea_prerank_side,
                           run_ppi_analysis)

st.set_page_config(page_title=t("net.page_title"), layout="wide")

# RWR 扩展靶点数上限 100;组件不设 max,任何超限输入都落到下方中文警告并钳制,
# 否则 Streamlit 前端会对超出 max 的输入弹内置英文报错(经 i18n 无法翻译)。
RWR_EXPAND_MAX = 100
# Hub 数量合法区间;控件同样不设 max,超限输入由下方中文警告钳制(组件 min=1 仅挡负数/0)。
HUB_TOP_N_MIN, HUB_TOP_N_MAX = 5, 20

# ---------- 隐藏侧边栏全屏显示;Run 按钮(#167A9C)、Data Management 按钮(#6347A1)配色 ----------
hide_sidebar("""
    .st-key-run_btn button {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
    }
    .st-key-data_mgmt_btn button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }
""")

st.title(t("net.title"))
st.caption(t("net.caption"))

# Build Report(busy)运行期间禁用页面交互:防止新点击/参数变动把分钟级分析中断丢弃
# (两段式 busy 状态机见 ui.widgets;disabled 不影响控件取值,翻转安全)
net_busy = is_busy("net_build_busy")

col_back, col_mgmt = st.columns([1, 1])
with col_back:
    back_button(disabled=net_busy)
with col_mgmt:
    if st.button(t("net.data_mgmt"), key="data_mgmt_btn", disabled=net_busy):
        st.switch_page("pages/Data_Management.py")

# 上一次分析被在途点击打断时给出提示(busy 复位后的下一轮才能看到,不能门在 busy 上)
show_busy_notice("net_build_busy")

st.divider()

# ---------- 输入区 ----------
input_mode = st.radio(
    t("net.input_mode"),
    ["names", "genes"],
    format_func=lambda v: t("net.mode_names") if v == "names" else t("net.mode_genes"),
    horizontal=True,
    disabled=net_busy,
)

col_input, col_params = st.columns([2, 1])
with col_input:
    if input_mode == "names":
        col_drug, col_dis = st.columns(2)
        with col_drug:
            drug_name = st.text_input(t("net.drug_label"), value="metformin",
                                      disabled=net_busy)
        with col_dis:
            disease_name = st.text_input(t("net.disease_label"), value="Type 2 diabetes mellitus",
                                         disabled=net_busy)
        gene_text = ""
    else:
        drug_name, disease_name = "", ""
        gene_text = st.text_area(
            t("net.gene_text_label"),
            placeholder="TP53\nBRCA1\nEGFR",
            height=140,
            disabled=net_busy,
        )

with col_params:
    # 物种固定为人类:疾病靶点数据源(Open Targets/DisGeNET)只有人类数据,
    # 药物靶点与富集库亦以人类为准;多物种支持待接入相应疾病数据源后再议。
    score = st.slider(t("net.confidence"), 0.1, 0.9, config.PPI_CONFIDENCE, 0.05,
                      disabled=net_busy)
    # 发起分析的位置常驻显示数据版本:切版本后不忘当前用的是哪版 STRING 数据
    st.caption(t("net.string_ver_caption", ver=config.STRING_VERSION))
    network_type = st.selectbox(
        t("net.network_type"), ("full", "physical"),
        format_func=lambda v: (t("net.network_type_full") if v == "full"
                               else t("net.network_type_physical")),
        help=t("net.network_type_help"),
        disabled=net_busy)
    rwr_expand = st.number_input(
        t("net.rwr_expand"), 0, None, config.RWR_EXPAND_TOP_N, 1,
        help=t("net.rwr_expand_help"), disabled=net_busy)
    if rwr_expand > RWR_EXPAND_MAX:
        st.warning(t("net.rwr_expand_clamp", n=rwr_expand))
        rwr_expand = RWR_EXPAND_MAX
    hub_top_n = st.number_input(
        t("net.hub_top_n"), 1, None, config.HUB_TOP_N, 1,
        help=t("net.hub_top_n_help", factor=config.HUB_SMALL_NET_FACTOR),
        disabled=net_busy)
    if not (HUB_TOP_N_MIN <= hub_top_n <= HUB_TOP_N_MAX):
        clamped = min(max(hub_top_n, HUB_TOP_N_MIN), HUB_TOP_N_MAX)
        st.warning(t("net.hub_top_n_clamp", n=hub_top_n, m=clamped))
        hub_top_n = clamped
    hub_method = st.selectbox(
        t("net.hub_method"), ("degree", "mcc"),
        index=0 if config.HUB_METHOD == "degree" else 1,
        format_func=lambda v: t("net.hub_method_degree") if v == "degree" else t("net.hub_method_mcc"),
        help=t("net.hub_method_help"),
        disabled=net_busy)

run_btn = st.button(t("net.build_report"), key="run_btn", type="primary",
                    width="stretch", disabled=net_busy)

# ---------- 点击轮:校验输入后置位 busy,任务轮再执行(两段式,见 ui.widgets) ----------
if run_btn and not net_busy:
    if input_mode == "names" and (not drug_name.strip() or not disease_name.strip()):
        st.warning(t("net.warn_missing"))
    else:
        enter_busy("net_build_busy")


# ---------- 分析执行 ----------
class _ReportProgress:
    """Build Report 进度面板:总进度条 + 各阶段状态行(⏳ 待执行 / 🔄 进行中 / ✅ 完成)。

    Streamlit 的 ScriptRunContext 绑定在脚本线程上,工作线程里更新组件会被静默忽略,
    因此所有 st 调用都只发生在主线程:线程任务只管计算,完成情况由主线程轮询获取后刷新。
    """

    _ICON = {"pending": "⏳", "running": "🔄", "done": "✅"}

    def __init__(self, steps: list[str]):
        self._states = {k: "pending" for k in steps}
        with st.container(border=True):
            self._caption = st.caption(t("net.prog.caption"))
            self._bar = st.progress(0.0)
            self._rows = {k: st.empty() for k in steps}
        for k in steps:
            self._redraw(k)

    def _redraw(self, key: str) -> None:
        self._rows[key].markdown(
            f"{self._ICON[self._states[key]]} {t(f'net.prog.{key}')}")

    def set(self, key: str, state: str) -> None:
        """更新某阶段状态(pending/running/done);键不在面板里时忽略。"""
        if key not in self._states:
            return
        self._states[key] = state
        self._redraw(key)

    def advance(self, fraction: float) -> None:
        self._bar.progress(min(max(fraction, 0.0), 1.0))


def _enrich_retry_suffix() -> str:
    """富集服务失败提示的重试后缀:有熔断冷却信息时给"约 X 分钟后自动重试"的大致时长,
    否则退回通用文案。"""
    mins = enrichment_retry_hint_minutes()
    return (t("net.enrich_retry_hint", mins=mins) if mins is not None
            else t("net.enrich_retry_later"))


# ---------- 结果暂存：页面上点下载等操作引起刷新时，不会重复跑分析 ----------
st.session_state.setdefault("ppi", None)

current_params = (
    input_mode,
    drug_name.strip(),
    disease_name.strip(),
    gene_text.strip(),
    score,
    network_type,
    int(rwr_expand),
    int(hub_top_n),
    hub_method,
)

# ---------- 任务轮(busy):跑完入库后复位 busy 重渲染,结果展示区只读会话 ----------
if net_busy:
    completed = False
    try:
        # 进度面板:names 模式 8 个阶段,genes 模式无 GSEA 与网络检验
        steps = ["targets", "string", "graph", "enrichment",
                 *(["gsea", "network"] if input_mode == "names" else []),
                 "modules", "save"]
        ppi_result = run_ppi_analysis(
            input_mode,
            drug_name.strip(),
            disease_name.strip(),
            gene_text.strip(),
            score,
            network_type=network_type,
            rwr_expand=int(rwr_expand),
            hub_top_n=int(hub_top_n),
            hub_method=hub_method,
            progress=_ReportProgress(steps),
        )
    except Exception:
        # 面向用户只提示友好信息;完整堆栈经 logging.exception 写入日志文件供排查
        logging.exception("PPI analysis failed")
        st.session_state["net_build_error"] = t("net.analysis_fail")
        completed = True  # 失败信息已入库,按已完成收尾(不显示"被中断")
    else:
        st.session_state.ppi = ppi_result
        # Build Report 已自动跑完药物/疾病两侧 GSEA，直接作为下方交互区的初始结果
        st.session_state.ppi_gsea = {
            side: {**info, "ppi_params": ppi_result["params"]}
            for side, info in ppi_result.get("gsea", {}).items()
        }
        # Build Report 已自动跑完模块检测，作为交互区初始结果；失败时可在交互区手动重跑
        mods_auto = ppi_result.get("modules")
        st.session_state.ppi_modules = (
            {**mods_auto, "ppi_params": ppi_result["params"]}
            if mods_auto is not None else None
        )
        # 新报告生成后,网络图视图选择回到默认交互式:视图选择框的 key 绑定报告实例
        # (见下方显示段),这里清掉旧报告的键,避免会话内键随分析次数无限累积
        for _k in [k for k in st.session_state if str(k).startswith("ppi_view::")]:
            st.session_state.pop(_k, None)
        st.session_state.pop("net_build_error", None)
        completed = True
    finally:
        finish_busy("net_build_busy", completed=completed)
    st.rerun()

# 上一次分析失败的友好提示(成功轮已清掉该键)
net_error = st.session_state.pop("net_build_error", None)
if net_error:
    st.error(net_error)

# ---------- 结果展示（直接读取暂存结果，不重新分析） ----------
ppi = st.session_state.ppi


def _params_stale() -> bool:
    """当前页面参数与暂存分析结果的参数是否不一致（结果已过期）。

    结果展示区各小节与 GSEA/模块手动重跑处共用此防呆判断。
    """
    return current_params != ppi["params"]


def _render_proximity(ppi: dict, input_mode: str) -> None:
    """网络邻近度检验展示(仅 names 模式;手动输入基因模式无药物/疾病之分,不展示)。"""
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
            prox_z = proximity_result.get("z")
            col_p2.metric(t("net.prox_z"),
                          f"{prox_z:.3f}" if prox_z is not None else "—")
            col_p3.metric(t("net.prox_p"), f"{proximity_result['p_empirical']:.4f}")
            if prox_z is None:
                # 零分布退化时 z 无意义,不给出支持/不支持的结论
                st.info(t("net.prox_degenerate"))
            elif prox_z < 0:
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
                # 外层分支已保证 ok=True,走到这里说明 d_c_excl 缺失
                st.caption(t("net.prox_excl_unavailable"))


def _render_rwr(ppi: dict, input_mode: str, prefix: str, saved: dict) -> None:
    """靶点优先级排序展示(仅 names 模式):两个方向随主分析自动算完,单选只切换查看,不重算。"""
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
                    download_button(
                        "net.dl_rwr_csv",
                        top_df.to_csv(index=False).encode("utf-8-sig"),
                        f"{prefix}_rwr_{rwr_side}.csv", "text/csv",
                    )
                    if "rwr_disease_csv" in saved or "rwr_drug_csv" in saved:
                        st.caption(t("net.rwr_in_report", dir=saved["directory"]))


def _render_enrichment(ppi: dict, genes: list[str]) -> None:
    """富集分析展示:结果表 + 柱状/气泡图 + 药物-靶点-通路网络图。"""
    st.divider()
    st.subheader(t("net.enrich_subheader"))
    enrichment_df = ppi["enrichment_df"]
    enrichment_source = ppi["enrichment_source"]
    enrichment_bar_fig = ppi["enrichment_bar_fig"]
    enrichment_dot_fig = ppi["enrichment_dot_fig"]
    if enrichment_df is None or enrichment_df.empty:
        # 双源都不可用(或无兜底物种的主源失败)时显著告警,与"真无显著项"的普通提示区分开;
        # 并附上大致的自动重试时间,避免用户反复手动重跑
        if ppi.get("enrichment_reason") in ("both_sources_down", "gprofiler_down_no_fallback"):
            st.error(t("net.enrich_both_down") + _enrich_retry_suffix())
        else:
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
        # 药物-靶点-通路网络图（基于 KEGG 富集结果构建）
        dtp_fig = ppi.get("dtp_fig")
        if dtp_fig is not None:
            st.caption(t("net.dtp_caption"))
            st.pyplot(dtp_fig)


def _render_gsea(ppi: dict, input_mode: str, drug_name: str, disease_name: str,
                 prefix: str, saved: dict, params_stale: bool) -> None:
    """排序富集(GSEA prerank)展示与手动重跑:换侧只切换已存结果,不触发重新计算。"""
    st.divider()
    with st.expander(t("net.gsea_subheader")):
        st.caption(t("net.gsea_caption", libs="、".join(config.GSEA_ENRICHR_LIBS)))
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
            # 该侧已有有效结果(参数一致且非请求失败)时禁用重跑按钮,避免误点重复
            # 计算;自动步骤失败时按钮保持可用,作为手动兜底入口。
            gsea_res = st.session_state.ppi_gsea.get(gsea_side)
            gsea_reason = (gsea_res or {}).get("reason", "")
            gsea_done = (gsea_res is not None and not params_stale
                         and gsea_reason not in ("request_failed", "circuit_open",
                                                 "normalize_empty"))
            run_clicked = st.button(
                t("net.gsea_run"), key="gsea_run_btn", width="stretch",
                disabled=gsea_done or params_stale or net_busy,
                help=(t("net.params_changed") if params_stale
                      else t("net.gsea_done_hint") if gsea_done else None),
            )
            if run_clicked:
                if params_stale:
                    # 防呆：主分析参数已变，先重跑主分析再 GSEA
                    st.warning(t("net.params_changed"))
                else:
                    with st.spinner(t("net.gsea_spinner")):
                        try:
                            gsea_df, rstats, reason = gsea_prerank_side(
                                gsea_side, drug_name, disease_name, prefix)
                        except Exception:
                            # 手动入口与自动流程同等兜底:网络/文件链路异常不白屏,
                            # 按请求失败记录,由下方 reason 分支提示并允许重试
                            logging.exception("Manual GSEA rerun failed (%s side)", gsea_side)
                            gsea_df, rstats, reason = None, {}, "request_failed"
                        st.session_state.ppi_gsea[gsea_side] = {
                            "df": gsea_df,
                            "stats": rstats,
                            "reason": reason,
                            "ppi_params": ppi["params"],
                        }
                    st.rerun()

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
                            # 展示阈值只作用于页面表格；导出始终为完整计算结果，
                            # 在此说明以免"提示无显著却仍能下载"造成误解
                            st.caption(t("net.gsea_full_export_note", thr=config.GSEA_FDR_SCREEN))
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
                        download_button(
                            "net.dl_gsea_csv",
                            gsea_df.to_csv(index=False).encode("utf-8-sig"),
                            f"{prefix}_gsea_{gsea_side}.csv", "text/csv",
                        )
                        # 该侧显著通路条形图已随报告构建落盘时,提供 PNG+SVG 小包下载;
                        # 筛后无显著(未落盘)则按钮不出现
                        fig_png = saved.get(f"gsea_{gsea_side}_png")
                        fig_svg = saved.get(f"gsea_{gsea_side}_svg")
                        if fig_png and fig_svg:
                            # key 与下载区同名按钮区分:两侧图包参数完全相同,
                            # 不加 key 会因元素 ID 冲突截断下载区的渲染(切侧必现报错)
                            download_button(
                                f"net.dl_gsea_fig_{gsea_side}",
                                make_fig_zip(fig_png, fig_svg),
                                f"{prefix}_gsea_{gsea_side}.zip", "application/zip",
                                key=f"gsea_fig_btn_{gsea_side}",
                            )


def _render_modules(ppi: dict, graph, metrics_df, prefix: str, saved: dict,
                    params_stale: bool) -> None:
    """蛋白模块检测展示与手动重跑(Louvain + 模块富集):重跑仅作自动步骤失败时的兜底。"""
    st.divider()
    with st.expander(t("net.modules_subheader")):
        st.caption(t("net.modules_caption",
                     min_size=config.LOUVAIN_MIN_MODULE_SIZE,
                     top_n=config.LOUVAIN_MAX_ENRICH_MODULES))
        st.session_state.setdefault("ppi_modules", None)

        if st.button(t("net.modules_run"), key="modules_run_btn", width="stretch",
                     disabled=net_busy):
            if params_stale:
                # 防呆：主分析参数已变，图与旧结果不对应，先重跑主分析
                st.warning(t("net.params_changed"))
            else:
                with st.spinner(t("net.modules_spinner_detect")):
                    # 物种固定人类;organism 取 Human 的 gProfiler 代码
                    try:
                        mods = compute_ppi_modules(
                            graph, config.SPECIES_MAP["Human"]["gprofiler"],
                            saved["directory"], prefix)
                    except Exception:
                        logging.exception("Module analysis failed")
                        st.error(t("net.analysis_fail"))
                    else:
                        saved.update(mods["saved_paths"])
                        st.session_state.ppi_modules = {
                            **mods,
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
                    st.error(t("net.modules_fail", ids=ids) + _enrich_retry_suffix())
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
                        download_button(
                            "net.dl_module_summary_csv",
                            summary_df.to_csv(index=False).encode("utf-8-sig"),
                            f"{prefix}_modules.csv", "text/csv",
                        )
                    with col_me:
                        download_button(
                            "net.dl_module_enrich_csv",
                            long_df.to_csv(index=False).encode("utf-8-sig"),
                            f"{prefix}_module_enrichment.csv", "text/csv",
                        )
                    if mods_res.get("saved_paths"):
                        st.caption(t("net.modules_saved", dir=saved["directory"]))
                elif not failures and not no_sig and not summary_df.empty:
                    st.info(t("net.modules_enrich_empty"))


def _render_downloads(metrics_df, edges_df, html_content: str, saved: dict, prefix: str) -> None:
    """下载区:先收集全部可用下载项,再按 3 列网格逐行铺满——
    条件产物(各风格图/富集/GSEA 图)缺失时不会留下空洞。
    图类一图一钮:PNG(300 dpi)+SVG 打成小包一次带走(总 ZIP 仍是散文件)。"""
    with st.expander(t("net.download_expander")):
        items: list[tuple[str, object, str, str]] = [
            ("net.dl_nodes", metrics_df.to_csv(index=False).encode("utf-8-sig"),
             f"{prefix}_nodes.csv", "text/csv"),
            ("net.dl_edges", edges_df.to_csv(index=False).encode("utf-8-sig"),
             f"{prefix}_edges.csv", "text/csv"),
            ("net.dl_html", html_content.encode("utf-8"),
             f"{prefix}_network.html", "text/html"),
            ("net.dl_graphml", saved["graphml"], f"{prefix}_graph.graphml",
             "application/xml"),
        ]
        if "enrichment_csv" in saved:
            items.append(("net.dl_enrich_csv", saved["enrichment_csv"],
                          f"{prefix}_enrichment.csv", "text/csv"))
        if "png" in saved and "svg" in saved:
            items.append(("net.dl_static_fig", make_fig_zip(saved["png"], saved["svg"]),
                          f"{prefix}_network.zip", "application/zip"))
        if "string_png" in saved:
            items.append(("net.dl_string_fig", make_fig_zip(saved["string_png"], saved["string_svg"]),
                          f"{prefix}_network_string.zip", "application/zip"))
        if "cyto_png" in saved:
            items.append(("net.dl_cyto_fig", make_fig_zip(saved["cyto_png"], saved["cyto_svg"]),
                          f"{prefix}_network_cytoscape.zip", "application/zip"))
        if "core_png" in saved:
            items.append(("net.dl_core_fig", make_fig_zip(saved["core_png"], saved["core_svg"]),
                          f"{prefix}_core_targets.zip", "application/zip"))
        if "dtp_png" in saved:
            items.append(("net.dl_dtp_fig", make_fig_zip(saved["dtp_png"], saved["dtp_svg"]),
                          f"{prefix}_enrichment_network.zip", "application/zip"))
        if "enrichment_bar_png" in saved:
            items.append(("net.dl_enrich_bar",
                          make_fig_zip(saved["enrichment_bar_png"], saved["enrichment_bar_svg"]),
                          f"{prefix}_enrichment_bar.zip", "application/zip"))
        if "enrichment_dot_png" in saved:
            items.append(("net.dl_bubble",
                          make_fig_zip(saved["enrichment_dot_png"], saved["enrichment_dot_svg"]),
                          f"{prefix}_enrichment_dot.zip", "application/zip"))
        for side in ("drug", "disease"):
            fig_png, fig_svg = saved.get(f"gsea_{side}_png"), saved.get(f"gsea_{side}_svg")
            if fig_png and fig_svg:
                items.append((f"net.dl_gsea_fig_{side}", make_fig_zip(fig_png, fig_svg),
                              f"{prefix}_gsea_{side}.zip", "application/zip"))
        items.append(("net.dl_zip", make_ppi_zip(saved),
                      f"{prefix}_results.zip", "application/zip"))

        for start in range(0, len(items), 3):
            chunk = items[start:start + 3]
            # 末行不足 3 个时按实际数量分列:按钮拉宽占满整行,不留缺角
            cols = st.columns(len(chunk))
            for col, (label_key, data, file_name, mime) in zip(cols, chunk):
                with col:
                    if isinstance(data, (bytes, bytearray)):
                        download_button(label_key, data, file_name, mime)
                    else:
                        download_file_button(label_key, data, file_name, mime)
        st.caption(
            t("net.files_grouped", dir=saved["directory"])
        )


if ppi is not None:
    if _params_stale():
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
    else:
        # 交集过小(如降噪后)可能没有任何达到置信度的互作边:此时 PPI 网络为空,
        # 但富集/GSEA/网络检验/下载等结果仍然有效,照常展示而非整页拦截。
        # physical 网络只保留物理结合证据,空网概率远高于 full,提示要更有针对性
        if edges_df.empty:
            # params 元组第 6 位是 network_type(旧版暂存结果可能更短,按 full 兜底)
            net_type = ppi["params"][5] if len(ppi["params"]) > 5 else "full"
            st.warning(t("net.no_interactions_physical" if net_type == "physical"
                         else "net.no_interactions"))
        st.info(t("net.input_genes", n=len(genes)))

        col_m1, col_m2, col_m3, col_m4 = st.columns(4)
        col_m1.metric(t("net.metric_input"), len(genes))
        col_m2.metric(t("net.metric_nodes"), len(metrics_df))
        col_m3.metric(t("net.metric_edges"), len(edges_df))
        col_m4.metric(t("net.metric_hubs"), len(hub_df))

        expansion_df = ppi.get("expansion_df")
        if expansion_df is not None:
            n_seed = int((expansion_df["source"] == "seed").sum())
            n_exp = int((expansion_df["source"] == "expanded").sum())
            st.caption(t("net.expansion_caption", s=n_seed, e=n_exp, k=n_seed + n_exp,
                         alpha=config.RWR_ALPHA, f=ppi.get("prefix", "")))

        st.divider()

        hub_rule = ppi.get("hub_rule", "")
        if hub_rule.startswith("top"):
            st.subheader(t("net.hub_top_subheader", n=int(hub_top_n),
                           method=ppi.get("hub_method_eff", "degree"),
                           t=config.HUB_SMALL_NET_FACTOR * int(hub_top_n)))
        else:
            st.subheader(t("net.hub_subheader"))
        if hub_df.empty:
            st.info(t("net.no_hubs"))
        else:
            st.dataframe(hub_df, width="stretch", hide_index=True)

        # ---------- 网络图展示:下拉切换,默认 Cytoscape 风格静态图(交互式 HTML 加载慢) ----------
        # 四种图(交互/静态/STRING/Cytoscape)无论当前选中哪个都已在 save_ppi_results 阶段落盘,
        # 这里只控制"看哪张";渲染失败的风格图不出现在选项里,不阻塞其余选项。
        views: dict[str, tuple[str, str | None]] = {
            "interactive": (t("net.view_interactive"), None),
            "static": (t("net.view_static"), None),
            "string": (t("net.view_string"), t("net.string_style_caption")),
            "cytoscape": (t("net.view_cyto"), t("net.cyto_style_caption")),
        }
        available = ["interactive", "static"]
        if ppi.get("string_fig") is not None:
            available.append("string")
        if ppi.get("cyto_fig") is not None:
            available.append("cytoscape")
        # 上次选中的图这次不存在(旧结果/渲染失败)时清掉,让选择框回到默认风格
        # key 绑定报告实例(id 每次新分析都不同;旧 ppi 被替换前一直存活,id 不会复用):
        # 直接复用固定 key 时,streamlit 前端控件不会跟随服务端 pop 重置,仍显示用户
        # 上次选的旧选项,与服务端实际渲染的默认风格图脱节;换 key 即换新控件,彻底同步。
        view_key = f"ppi_view::{id(ppi)}"
        if st.session_state.get(view_key) not in available:
            st.session_state.pop(view_key, None)
        st.subheader(t("net.view_subheader"))
        default_view = "cytoscape" if "cytoscape" in available else "interactive"
        view = st.selectbox(t("net.view_label"), available,
                            index=available.index(default_view),
                            format_func=lambda v: views[v][0], key=view_key)
        view_caption = views[view][1]
        if view_caption:
            st.caption(view_caption)
        if view == "interactive":
            st.iframe(html_content, height=640)
        elif view == "static":
            st.pyplot(static_fig)
        elif view == "string":
            st.pyplot(ppi["string_fig"])
        else:
            st.pyplot(ppi["cyto_fig"])
        st.caption(t("net.view_saved_note"))

        # ---------- 核心靶点子网络（CytoHubba 式 MCC Top N） ----------
        core_fig = ppi.get("core_fig")
        if core_fig is not None:
            st.subheader(t("net.core_subheader"))
            st.caption(t("net.core_caption", n=config.MCC_TOP_N))
            st.pyplot(core_fig)

        # ---------- 网络邻近度 / RWR / 富集 / GSEA / 模块 / 下载区 ----------
        params_stale = _params_stale()
        _render_proximity(ppi, input_mode)

        _render_rwr(ppi, input_mode, prefix, saved)
        _render_enrichment(ppi, genes)

        _render_gsea(ppi, input_mode, drug_name, disease_name, prefix, saved, params_stale)

        _render_modules(ppi, graph, metrics_df, prefix, saved, params_stale)

        _render_downloads(metrics_df, edges_df, html_content, saved, prefix)

        st.success(t("net.success"))
