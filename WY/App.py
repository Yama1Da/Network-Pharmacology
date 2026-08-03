import streamlit as st
import Utils
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

st.set_page_config(page_title="Network Pharmacology Platform", layout="centered")

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
    /* Top-left corner management buttons (fixed overlay, no card background) */
    div.st-key-top_left_buttons {
        position: fixed;
        top: 3.5rem;
        left: 1rem;
        z-index: 9999;
        display: flex;
        flex-direction: row;
        gap: 0.6rem;
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
    @keyframes gentle-float {
        0%   { transform: translateY(0); }
        50%  { transform: translateY(-4px); }
        100% { transform: translateY(0); }
    }
</style>
""", unsafe_allow_html=True)

# ---------- Title ----------
st.markdown('<div class="title-area"><h1>🧬 Network Pharmacology Analysis Platform</h1><p>Input drug and disease names to analyze overlapping targets</p></div>',
            unsafe_allow_html=True)

# ---------- Management buttons (fixed top-left, open full-screen pages) ----------
with st.container(key="top_left_buttons"):
    if st.button("Drug Target Information Management", key="drug_mgmt_btn"):
        st.switch_page("pages/Drug_Target_Information_Management.py")
    if st.button("Disease Target Information Management", key="disease_mgmt_btn"):
        st.switch_page("pages/Disease_Target_Information_Management.py")

st.divider()

# ---------- Search bar ----------
col_drug, col_disease = st.columns(2)
with col_drug:
    drug_name = st.text_input("Drug name", value="aspirin",
                               label_visibility="collapsed",
                               placeholder="Drug name, e.g. aspirin")
with col_disease:
    disease_name = st.text_input("Disease name", value="cancer",
                                  label_visibility="collapsed",
                                  placeholder="Disease name, e.g. cancer")

analyze_btn = st.button("Analyze", type="primary", use_container_width=True)

# ---------- Results ----------
if analyze_btn:
    errors = []
    if not drug_name.strip():
        errors.append("drug name")
    if not disease_name.strip():
        errors.append("disease name")

    if errors:
        st.warning(f"Please enter {' and '.join(errors)}")
    else:
        with st.spinner(f"Querying targets for {drug_name} and {disease_name}…"):
            drug_df, _ = Utils.fetch_drug_targets_chembl(drug_name)
            drug_targets = drug_df["approvedSymbol"].tolist() if not drug_df.empty else []
            disease_df = Utils.get_disease_targets(disease_name)
            # 韦恩图/统计需要基因符号列表，从表格中取出
            disease_targets = disease_df["approvedSymbol"].tolist() if not disease_df.empty else []

        st.markdown('<div class="result-card">', unsafe_allow_html=True)

        drug_set = set(drug_targets)
        disease_set = set(disease_targets)

        col_a, col_b, col_c = st.columns([1, 1, 1])
        col_a.metric("Drug", drug_name)
        col_b.metric("Targets", len(drug_set))
        col_c.metric("Disease targets", len(disease_set))

        common = drug_set & disease_set

        st.metric("🎯 Common targets", len(common))

        # ----- Venn diagram -----
        fig, ax = plt.subplots(figsize=(6, 4))
        v = venn2(
            subsets=(len(drug_set - disease_set),
                     len(disease_set - drug_set),
                     len(common)),
            set_labels=(f"{drug_name} targets", f"{disease_name} targets"),
            ax=ax,
        )
        if v and v.get_label_by_id("11"):
            v.get_label_by_id("11").set_fontsize(12)
        ax.set_title(f"{drug_name} vs {disease_name} — Venn Diagram", fontsize=14)
        st.pyplot(fig)

        with st.expander("View target lists", expanded=False):
            # ----- 第一段：药物靶点（全量，按 pchembl 降序） -----
            st.markdown(f"**💊 {drug_name} targets** — {len(drug_set)} genes ({len(drug_df)} records, sorted by pChEMBL desc)")
            st.write(drug_df if not drug_df.empty else "（未获取到药物靶点数据）")

            st.divider()

            # ----- 第二段：疾病靶点（表格，含 Open Targets 关联指标，按分数降序） -----
            st.markdown(f"**🦠 {disease_name} targets** — {len(disease_targets)} genes")
            if not disease_df.empty:
                st.dataframe(disease_df, use_container_width=True)
            else:
                st.write("（未获取到疾病靶点数据）")

        st.markdown('</div>', unsafe_allow_html=True)

        st.success("✅ Analysis complete")
else:
    st.info("💡 Enter drug and disease names above and click Analyze")