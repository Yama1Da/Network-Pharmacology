import pandas as pd
import streamlit as st
import Utils
from i18n import t, set_lang, get_lang
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# 初始化日志:每天一个日志文件,保留两年。
Utils.setup_logging()

#页面布局
st.set_page_config(page_title=t("app.page_title"), layout="centered")

is_en = get_lang() == "en"

# ---------- Custom Minimal CSS ----------
st.markdown("""
<style>
    /* Center container */
    .main > div {
        max-width: 720px;
        margin: 0 auto;
    }
    /* Title area */
    .title-area {
        text-align: center;
        padding: 2.5rem 0 1rem 0;
    }
    .title-area h1 {
        font-size: 2rem;
        font-weight: 600;
        letter-spacing: 1px;
    }
    .title-area p {
        color: #6b7280;
        font-size: 0.95rem;
        margin-top: 0.3rem;
    }
    /* Enlarge input */
    .stTextInput input {
        font-size: 1.2rem;
        padding: 0.7rem 1rem;
    }
    /* Enlarge button */
    .stButton button {
        font-size: 1.1rem;
        padding: 0.5rem 1rem;
    }
    /* Result card */
    .result-card {
        background: #ffffff;
        border-radius: 12px;
        padding: 0.3rem 0.3rem;
    }
    /* Alert styling */
    .stAlert {
        border-radius: 8px;
    }
    /* Hide auto multipage navigation sidebar */
    [data-testid="stSidebar"] {
        display: none;
    }
    /* Top-right corner language + Logs buttons (fixed overlay) */
    div.st-key-top_right_logs {
        position: fixed;
        top: 3.5rem;
        right: 1rem;
        left: 0;
        z-index: 9999;
        display: flex;
        flex-direction: row;
        justify-content: flex-end;
        align-items: center;
        gap: 0.5rem;
        pointer-events: none;
    }
    /* 容器不拦截点击，按钮自身恢复可点击 */
    div.st-key-top_right_logs .stButton {
        pointer-events: auto;
    }
    /* Logs button: #909399 background, white text；高度与蓝色按钮一致（40px、单行） */
    div.st-key-top_right_logs .stButton button {
        background-color: #909399;
        color: #FFFFFF;
        border: 1px solid #909399;
        box-shadow: none;
        height: 40px;
        padding: 4px 12px;
        font-size: 16px;
        line-height: 1;
        white-space: nowrap;
    }
    div.st-key-top_right_logs .stButton button:hover,
    div.st-key-top_right_logs .stButton button:active {
        background-color: #909399;
        color: #FFFFFF;
        border: 1px solid #909399;
        box-shadow: none;
    }
    /* 语言切换：iOS 风格滑动开关，文本直接嵌在开关轨道内
       （左=中文红 / 右=English蓝，文字始终位于滑块对侧） */
    .st-key-lang_switch_wrap {
        display: flex;
        flex-direction: row;
        align-items: center;
        justify-content: flex-end;
        gap: 0.5rem !important;
        /* 在上层纵向布局中靠右对齐，且宽度收缩为内容宽 */
        align-self: flex-end;
    }
    /* fixed 容器默认不拦截点击，语言开关整体恢复可点击 */
    div.st-key-top_right_logs .st-key-lang_switch_wrap {
        pointer-events: auto;
    }
    /* 轨道：显式放大为带文字的宽开关。Streamlit 1.62 开关结构：
       轨道是 label 内第一个 div，选中态在 label[data-selected] 上 */
    .st-key-lang_switch_wrap label > div:first-of-type {
        position: relative;
        width: 80px !important;
        height: 26px !important;
        padding: 0 !important;
    }
    /* 轨道配色：未选中（中文模式）=红；选中（English模式）=蓝 */
    .st-key-lang_switch_wrap label:not([data-selected="true"]) > div:first-of-type {
        background-color: #E53E3E !important;
        border-color: #E53E3E !important;
    }
    .st-key-lang_switch_wrap label[data-selected="true"] > div:first-of-type {
        background-color: #3B82F6 !important;
        border-color: #3B82F6 !important;
    }
    /* 滑块：绝对定位（20x20，3px 内边距），选中时滑到右端 */
    .st-key-lang_switch_wrap label > div:first-of-type > div {
        position: absolute !important;
        top: 3px !important;
        left: 3px !important;
        width: 20px !important;
        height: 20px !important;
        transform: translateX(0) !important;
    }
    .st-key-lang_switch_wrap label[data-selected="true"] > div:first-of-type > div {
        transform: translateX(54px) !important;
    }
    /* 轨道内文字：显示当前语言，位于滑块对侧 */
    .st-key-lang_switch_wrap label > div:first-of-type::after {
        position: absolute;
        top: 50%;
        transform: translateY(-50%);
        font-size: 12px;
        font-weight: 700;
        line-height: 1;
        letter-spacing: 0.3px;
        color: #FFFFFF;
        pointer-events: none;
    }
    .st-key-lang_switch_wrap label:not([data-selected="true"]) > div:first-of-type::after {
        content: "中文";
        right: 8px;
    }
    .st-key-lang_switch_wrap label[data-selected="true"] > div:first-of-type::after {
        content: "English";
        left: 8px;
    }
    /* Top-left corner management buttons (fixed overlay, no card background) */
    div.st-key-top_left_buttons {
        position: fixed;
        top: 3.5rem;
        left: 1rem;
        z-index: 9999;
        display: flex;
        flex-direction: row;
        gap: 0.6rem;
        width: fit-content;
    }
    /* Management buttons: #167A9C background, #C6E7F5 text, gently float up/down */
    div.st-key-top_left_buttons .stButton button {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
        box-shadow: none;
        animation: gentle-float 2s ease-in-out infinite;
    }
    /* On hover: stop floating, keep colors unchanged */
    div.st-key-top_left_buttons .stButton button:hover,
    div.st-key-top_left_buttons .stButton button:active {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
        box-shadow: none;
        animation: none;
    }
    /* Popover 触发按钮（Source Information Management / PPI-Enrichment & Molecular Docking）同色 */
    .st-key-source_mgmt_popover button,
    .st-key-docking_mgmt_popover button {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
        box-shadow: none;
    }
    .st-key-source_mgmt_popover button:hover,
    .st-key-source_mgmt_popover button:active,
    .st-key-docking_mgmt_popover button:hover,
    .st-key-docking_mgmt_popover button:active {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
        box-shadow: none;
    }
    /* 下拉面板选项：浅灰底、深灰 hover 变色，醒目可辨 */
    .st-key-drug_mgmt_btn button,
    .st-key-disease_mgmt_btn button,
    .st-key-ppi_info_btn button,
    .st-key-enrichment_info_btn button,
    .st-key-data_network_btn button,
    .st-key-molecular_docking_btn button {
        background-color: #F0F0F0;
        color: #333333;
        border: 1px solid #D0D0D0;
        width: 100%;
    }
    .st-key-drug_mgmt_btn button:hover,
    .st-key-disease_mgmt_btn button:hover,
    .st-key-ppi_info_btn button:hover,
    .st-key-enrichment_info_btn button:hover,
    .st-key-data_network_btn button:hover,
    .st-key-molecular_docking_btn button:hover {
        background-color: #B8B8B8;
        color: #000000;
        border-color: #888888;
    }
    @keyframes gentle-float {
        0%   { transform: translateY(0); }
        50%  { transform: translateY(-4px); }
        100% { transform: translateY(0); }
    }
</style>
""", unsafe_allow_html=True)

# ---------- Title ----------
st.markdown(f'<div class="title-area"><h1>{t("app.title")}</h1><p>{t("app.subtitle")}</p></div>',
            unsafe_allow_html=True)

# ---------- Language sliding switch + Logs button (fixed top-right) ----------
with st.container(key="top_right_logs"):
    # 语言切换开关:左为中文(红色),右为英文(蓝色),文字显示在开关轨道上。
    with st.container(key="lang_switch_wrap"):
        # 若语言在其它页面被切换过,让开关跟随当前语言。
        if st.session_state.get("_lang_switch_synced") != get_lang():
            st.session_state["_lang_switch"] = is_en
            st.session_state["_lang_switch_synced"] = get_lang()
        lang_is_en = st.toggle("Language", value=is_en, key="_lang_switch",
                               label_visibility="collapsed")
        if lang_is_en != is_en:
            set_lang("en" if lang_is_en else "zh")
            st.session_state["_lang_switch_synced"] = get_lang()
            st.rerun()
    if st.button(t("app.logs_btn"), key="logs_btn"):
        st.switch_page("pages/Logs.py")

# ---------- Management buttons (fixed top-left, open full-screen pages) ----------
with st.container(key="top_left_buttons"):
    # 下拉菜单:跳转到四个信息管理页面。
    with st.popover(t("app.popover_source"), key="source_mgmt_popover"):
        if st.button(t("app.btn_drug_mgmt"), key="drug_mgmt_btn", width="stretch"):
            st.switch_page("pages/Drug_Target_Information_Management.py")
        if st.button(t("app.btn_disease_mgmt"), key="disease_mgmt_btn", width="stretch"):
            st.switch_page("pages/Disease_Target_Information_Management.py")
        if st.button(t("app.btn_ppi_info"), key="ppi_info_btn", width="stretch"):
            st.switch_page("pages/PPI_Information.py")
        if st.button(t("app.btn_enrichment_info"), key="enrichment_info_btn", width="stretch"):
            st.switch_page("pages/Enrichment_Information.py")
    # 下拉菜单:主分析(PPI/富集)与分子对接入口。
    with st.popover(t("app.popover_docking"), key="docking_mgmt_popover"):
        if st.button(t("app.btn_data_network"), key="data_network_btn", width="stretch"):
            st.switch_page("pages/Data_Network.py")
        if st.button(t("app.btn_molecular_docking"), key="molecular_docking_btn", width="stretch"):
            st.session_state["_docking_entry"] = True  # 通知对接页:从入口进入并重置状态
            st.switch_page("pages/Molecular_Docking.py")

st.divider()# 分隔管理按钮和搜索栏

# ---------- Search bar ----------
col_drug, col_disease = st.columns(2)# 两列分别输入药物名称和疾病名称
with col_drug:
    drug_name = st.text_input(t("app.drug_label"), value="aspirin",
                               label_visibility="collapsed",
                               placeholder=t("app.drug_placeholder"))
with col_disease:
    disease_name = st.text_input(t("app.disease_label"), value="cancer",
                                  label_visibility="collapsed",
                                  placeholder=t("app.disease_placeholder"))

analyze_btn = st.button(t("app.analyze_btn"), type="primary", width="stretch")

# ---------- Results ----------
if analyze_btn:
    errors = []
    if not drug_name.strip():
        errors.append(t("app.err_drug"))
    if not disease_name.strip():
        errors.append(t("app.err_disease"))

    if errors:
        st.warning(t("app.warn_missing", fields=t("app.and_join").join(errors)))
    else:
        with st.spinner(t("app.spinner_query", drug=drug_name, disease=disease_name)):
            # 药物靶点来自 ChEMBL + DGIdb,疾病靶点来自 Open Targets + DisGeNET,均已去重
            drug_df = Utils.get_drug_targets_merged(drug_name)
            drug_targets = drug_df["approvedSymbol"].tolist() if not drug_df.empty else []
            disease_df = Utils.get_disease_targets_merged(disease_name)
            # 从表格取出基因符号列表,供韦恩图和统计使用
            disease_targets = disease_df["approvedSymbol"].tolist() if not disease_df.empty else []

        st.markdown('<div class="result-card">', unsafe_allow_html=True)

        drug_set = set(drug_targets)
        disease_set = set(disease_targets)

        col_a, col_b, col_c = st.columns([1, 1, 1])
        col_a.metric(t("app.metric_drug"), drug_name)
        col_b.metric(t("app.metric_targets"), len(drug_set))
        col_c.metric(t("app.metric_disease_targets"), len(disease_set))

        common = drug_set & disease_set

        st.metric(t("app.metric_common"), len(common))

        # ----- Venn diagram -----
        fig, ax = plt.subplots(figsize=(6, 4))# 创建 6x4 英寸画布
        # 图内文字保持英文:matplotlib 默认字体不支持中文,中文会显示成方框。
        v = venn2(
            subsets=(len(drug_set - disease_set),
                     len(disease_set - drug_set),
                     len(common)),
            set_labels=(f"{drug_name} targets", f"{disease_name} targets"),
            ax=ax,
        )# 绘制韦恩图
        if v and v.get_label_by_id("11"):
            v.get_label_by_id("11").set_fontsize(12)# 设置交集区域字号
        ax.set_title(f"{drug_name} vs {disease_name} — Venn Diagram", fontsize=14)# 设置图标题
        st.pyplot(fig)# 显示韦恩图

        with st.expander(t("app.expander_targets"), expanded=False):
            # ----- 药物靶点表 -----
            st.markdown(t("app.drug_targets_md", name=drug_name, count=len(drug_set)))
            if not drug_df.empty:
                st.dataframe(
                    drug_df.sort_values(["source_count", "score_pchembl"],
                                        ascending=False, na_position="last"),
                    width="stretch",
                )
            else:
                st.write(t("app.no_drug_data"))

            st.divider()

            # ----- 疾病靶点表 -----
            st.markdown(t("app.disease_targets_md", name=disease_name, count=len(disease_set)))
            if not disease_df.empty:
                st.dataframe(
                    disease_df.sort_values(["source_count", "score_ot"],
                                           ascending=False, na_position="last"),
                    width="stretch",
                )
            else:
                st.write(t("app.no_disease_data"))

            st.divider()

            # ----- 交集靶点表 -----
            st.markdown(t("app.common_targets_md", drug=drug_name, disease=disease_name,
                          count=len(common)))
            if common:
                # 合并两侧信息：每个交集基因保留药物侧和疾病侧的数据库来源与分数
                common_df = pd.merge(
                    drug_df.drop(columns=["source_count"]),
                    disease_df.drop(columns=["source_count"]),
                    on="approvedSymbol", how="inner",
                )
                # 两侧共 4 个数据库，支持的数据库越多排序越靠前
                common_df["source_total"] = common_df[
                    ["source_chmbl", "source_dgidb", "source_ot", "source_disgenet"]
                ].sum(axis=1)
                st.dataframe(
                    common_df.sort_values(["source_total", "score_pchembl", "score_ot"],
                                          ascending=False, na_position="last"),
                    width="stretch",
                )
            else:
                st.write(t("app.no_common_data"))

        st.markdown('</div>', unsafe_allow_html=True)

        st.success(t("app.success_done"))
else:
    st.info(t("app.info_hint"))
