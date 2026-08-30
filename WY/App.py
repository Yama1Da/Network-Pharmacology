import streamlit as st
import Utils
from i18n import t, set_lang, get_lang
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# 初始化日志系统（幂等：每天一个文件、保留两年、四级全录）
Utils.setup_logging()

#页面布局
st.set_page_config(page_title=t("app.page_title"), layout="centered")

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
        gap: 0.5rem;
        pointer-events: none;
    }
    /* 容器不拦截点击，按钮自身恢复可点击 */
    div.st-key-top_right_logs .stButton {
        pointer-events: auto;
    }
    /* Logs button: #909399 background, white text */
    div.st-key-top_right_logs .stButton button {
        background-color: #909399;
        color: #FFFFFF;
        border: 1px solid #909399;
        box-shadow: none;
    }
    div.st-key-top_right_logs .stButton button:hover,
    div.st-key-top_right_logs .stButton button:active {
        background-color: #909399;
        color: #FFFFFF;
        border: 1px solid #909399;
        box-shadow: none;
    }
    /* Language sliding switch: knob left = 中文, knob right = English */
    .st-key-lang_switch [data-baseweb="checkbox"] {
        position: relative;
    }
    .st-key-lang_switch [data-baseweb="checkbox"]::before {
        content: "中文";
        position: absolute;
        right: calc(100% + 8px);
        top: 50%;
        transform: translateY(-50%);
        color: #909399;
        font-size: 0.9rem;
        white-space: nowrap;
    }
    .st-key-lang_switch [data-baseweb="checkbox"]::after {
        content: "English";
        position: absolute;
        left: calc(100% + 8px);
        top: 50%;
        transform: translateY(-50%);
        color: #909399;
        font-size: 0.9rem;
        white-space: nowrap;
    }
    /* 当前语言一侧文字高亮（深色加粗），另一侧灰色 */
    .st-key-lang_switch [data-baseweb="checkbox"]:not(:has(input:checked))::before {
        color: #1F2937;
        font-weight: 600;
    }
    .st-key-lang_switch [data-baseweb="checkbox"]:has(input:checked)::after {
        color: #1F2937;
        font-weight: 600;
    }
    /* 轨道颜色随状态：中文 = 灰，English = 蓝 */
    .st-key-lang_switch [data-baseweb="checkbox"] > div {
        background-color: #909399;
    }
    .st-key-lang_switch [data-baseweb="checkbox"]:has(input:checked) > div {
        background-color: #4A90D9;
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
    # 滑动开关：滑块居左 = 中文，居右 = English；点击一次滑动切换
    is_en = get_lang() == "en"
    switched = st.toggle("Language", value=is_en, label_visibility="collapsed",
                         key="lang_switch")
    if switched != is_en:
        set_lang("en" if switched else "zh")
        st.rerun()
    if st.button(t("app.logs_btn"), key="logs_btn"):
        st.switch_page("pages/Logs.py")

# ---------- Management buttons (fixed top-left, open full-screen pages) ----------
with st.container(key="top_left_buttons"):
    # 下拉式菜单：四个信息管理页面的前往按钮（不含 Data Management）
    with st.popover(t("app.popover_source"), key="source_mgmt_popover"):
        if st.button(t("app.btn_drug_mgmt"), key="drug_mgmt_btn", width="stretch"):
            st.switch_page("pages/Drug_Target_Information_Management.py")
        if st.button(t("app.btn_disease_mgmt"), key="disease_mgmt_btn", width="stretch"):
            st.switch_page("pages/Disease_Target_Information_Management.py")
        if st.button(t("app.btn_ppi_info"), key="ppi_info_btn", width="stretch"):
            st.switch_page("pages/PPI_Information.py")
        if st.button(t("app.btn_enrichment_info"), key="enrichment_info_btn", width="stretch"):
            st.switch_page("pages/Enrichment_Information.py")
    # 下拉式菜单：PPI/富集主分析 + 分子对接（分子对接为富集下游页面）
    with st.popover(t("app.popover_docking"), key="docking_mgmt_popover"):
        if st.button(t("app.btn_data_network"), key="data_network_btn", width="stretch"):
            st.switch_page("pages/Data_Network.py")
        if st.button(t("app.btn_molecular_docking"), key="molecular_docking_btn", width="stretch"):
            st.session_state["_docking_entry"] = True  # 通知对接页：从入口进入，重置状态
            st.switch_page("pages/Molecular_Docking.py")

st.divider()# 水平分隔线，用于分隔管理按钮和搜索栏

# ---------- Search bar ----------
col_drug, col_disease = st.columns(2)# 两个列，用于输入药物名称和疾病名称
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
            # 药物靶点：ChEMBL + DGIdb 去重整合；疾病靶点：Open Targets + DisGeNET 去重整合
            drug_df = Utils.get_drug_targets_merged(drug_name)
            drug_targets = drug_df["approvedSymbol"].tolist() if not drug_df.empty else []
            disease_df = Utils.get_disease_targets_merged(disease_name)
            # 韦恩图/统计需要基因符号列表，从表格中取出
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
        fig, ax = plt.subplots(figsize=(6, 4))# 创建一个6x4英寸的子图
        v = venn2(
            subsets=(len(drug_set - disease_set),
                     len(disease_set - drug_set),
                     len(common)),
            set_labels=(t("app.targets_of", name=drug_name), t("app.targets_of", name=disease_name)),
            ax=ax,
        )# 绘制Venn图
        if v and v.get_label_by_id("11"):
            v.get_label_by_id("11").set_fontsize(12)# 设置Venn图中“交集”的字体大小
        ax.set_title(t("app.venn_title", drug=drug_name, disease=disease_name), fontsize=14)# 设置Venn图的标题
        st.pyplot(fig)# 显示Venn图

        with st.expander(t("app.expander_targets"), expanded=False):
            # ----- 第一段：药物靶点（ChEMBL + DGIdb 去重整合结果） -----
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

            # ----- 第二段：疾病靶点（Open Targets + DisGeNET 去重整合结果） -----
            st.markdown(t("app.disease_targets_md", name=disease_name, count=len(disease_set)))
            if not disease_df.empty:
                st.dataframe(
                    disease_df.sort_values(["source_count", "score_ot"],
                                           ascending=False, na_position="last"),
                    width="stretch",
                )
            else:
                st.write(t("app.no_disease_data"))

        st.markdown('</div>', unsafe_allow_html=True)

        st.success(t("app.success_done"))
else:
    st.info(t("app.info_hint"))
