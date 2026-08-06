import os
import streamlit as st
import pandas as pd
import Utils

st.set_page_config(page_title="PPI Information", layout="wide")

CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "ppi",
)


# ---------- Hide auto navigation sidebar for full-screen display ----------
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        display: none;
    }
    /* 行内 Download 按钮（#207D62 底、#C8F0DF 字） */
    .st-key-file_rows [class*="st-key-download_"] button {
        background-color: #207D62;
        color: #C8F0DF;
        border: 1px solid #207D62;
    }
</style>
""", unsafe_allow_html=True)

st.title("🕸️ PPI Information")
st.caption("Cached STRING interaction files under Static/data/cache/ppi")

col_back, _ = st.columns([1, 4])
with col_back:
    if st.button("← Back to Main", key="back_btn"):
        st.switch_page("App.py")

st.divider()

ppi_files = Utils.list_cache_files("ppi")
if ppi_files.empty:
    st.info("No cached PPI files. Run a Data Network analysis first.")
else:
    # 表头
    col_h_name, col_h_time, col_h_dl = st.columns([4, 2, 1])
    col_h_name.markdown("**File name**")
    col_h_time.markdown("**Created time**")
    col_h_dl.markdown("**Download**")

    # 每行一个 Download 按钮，下载磁盘上的缓存 CSV（文件名含来源标识，如 StringPPI_aspirin_cancer_xxx.csv）
    with st.container(key="file_rows"):
        for _, row in ppi_files.iterrows():
            col_name, col_time, col_dl = st.columns([4, 2, 1])
            col_name.write(row["File name"])
            col_time.write(row["Created time"])
            with col_dl:
                file_path = os.path.join(CACHE_DIR, row["File name"])
                if os.path.exists(file_path):
                    with open(file_path, "rb") as f:
                        csv_bytes = f.read()
                    st.download_button(
                        "Download",
                        data=csv_bytes,
                        file_name=row["File name"],
                        mime="text/csv",
                        key=f"download_{row['File name']}",
                    )
