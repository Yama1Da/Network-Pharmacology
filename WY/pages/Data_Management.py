import os
import streamlit as st
import pandas as pd

import config
import Utils

st.set_page_config(page_title="Data Management", layout="wide")

# ---------- Hide auto navigation sidebar for full-screen display ----------
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        display: none;
    }
    /* Download button colors (#207D62 bg, #C8F0DF text) */
    .st-key-dir_rows button {
        background-color: #207D62;
        color: #C8F0DF;
        border: 1px solid #207D62;
    }
</style>
""", unsafe_allow_html=True)

st.title("📁 Data Management")
st.caption("Downloaded result groups under Static/data/downloads")

col_back, _ = st.columns([1, 4])
with col_back:
    if st.button("← Back to Network", key="back_btn"):
        st.switch_page("pages/Data_Network.py")

st.divider()

dirs_df = Utils.list_download_dirs()
if dirs_df.empty:
    st.info("No downloaded result groups yet. Run a PPI analysis first.")
else:
    # 表头
    col_h_name, col_h_time, col_h_btn = st.columns([3, 2, 1])
    col_h_name.markdown("**Directory**")
    col_h_time.markdown("**Created time**")
    col_h_btn.markdown("**Download**")

    # 每行一个 Download 按钮，打包并下载该行目录的 zip
    with st.container(key="dir_rows"):
        for _, row in dirs_df.iterrows():
            dir_path = os.path.join(config.DOWNLOADS_DIR, row["Directory"])
            col_name, col_time, col_btn = st.columns([3, 2, 1])
            col_name.write(row["Directory"])
            col_time.write(row["Created time"])
            with col_btn:
                st.download_button(
                    "Download",
                    Utils.make_dir_zip(dir_path),
                    f"{row['Directory']}.zip",
                    "application/zip",
                    key=f"dl_{row['Directory']}",
                    width="stretch",
                )
