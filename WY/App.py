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
</style>
""", unsafe_allow_html=True)

# ---------- Title ----------
st.markdown('<div class="title-area"><h1>🧬 Network Pharmacology Analysis Platform</h1><p>Input drug and disease names to analyze overlapping targets</p></div>',
            unsafe_allow_html=True)

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
            disease_targets = Utils.get_disease_targets(disease_name)

        st.markdown('<div class="result-card">', unsafe_allow_html=True)

        col_a, col_b, col_c = st.columns([1, 1, 1])
        col_a.metric("Drug", drug_name)
        col_b.metric("Targets", len(drug_targets))
        col_c.metric("Disease targets", len(disease_targets))

        drug_set = set(drug_targets)
        disease_set = set(disease_targets)
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

        with st.expander("View target list", expanded=False):
            st.write(drug_df if not drug_df.empty else [])

        st.markdown('</div>', unsafe_allow_html=True)

        st.success("✅ Analysis complete")
else:
    st.info("💡 Enter drug and disease names above and click Analyze")