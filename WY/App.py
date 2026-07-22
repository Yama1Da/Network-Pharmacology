import streamlit as st
import Utils

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
        border: 1px solid #e5e7eb;
        border-radius: 12px;
        padding: 1.5rem 2rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.04);
    }
    /* Alert styling */
    .stAlert {
        border-radius: 8px;
    }
</style>
""", unsafe_allow_html=True)

# ---------- Title ----------
st.markdown('<div class="title-area"><h1>🧬 Network Pharmacology Analysis Platform</h1><p>Input drug name to query its target information</p></div>',
            unsafe_allow_html=True)

# ---------- Search bar ----------
col1, col2 = st.columns([3, 1])
with col1:
    drug_name = st.text_input("Drug name", value="aspirin",
                               label_visibility="collapsed",
                               placeholder="e.g. aspirin")
with col2:
    analyze_btn = st.button("Analyze", type="primary", use_container_width=True)

# ---------- Results ----------
if analyze_btn:
    if not drug_name.strip():
        st.warning("Please enter a drug name")
    else:
        with st.spinner(f"Querying targets for {drug_name}…"):
            drug_targets = Utils.get_drug_targets(drug_name)

        st.markdown('<div class="result-card">', unsafe_allow_html=True)

        col_a, col_b = st.columns([1, 1])
        col_a.metric("Drug", drug_name)
        col_b.metric("Targets", len(drug_targets))

        with st.expander("View target list", expanded=False):
            st.write(drug_targets)

        st.markdown('</div>', unsafe_allow_html=True)

        st.success("✅ Analysis complete")
else:
    st.info("💡 Enter a drug name above and click Analyze")
