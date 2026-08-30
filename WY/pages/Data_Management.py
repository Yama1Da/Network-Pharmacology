import os
import streamlit as st
import pandas as pd

import config
import Utils
from i18n import t

st.set_page_config(page_title=t("data_mgmt.page_title"), layout="wide")

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

st.title(t("data_mgmt.title"))
st.caption(t("data_mgmt.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    if st.button(t("data_mgmt.back"), key="back_btn"):
        st.switch_page("pages/Data_Network.py")

st.divider()

dirs_df = Utils.list_download_dirs()
# 排除分子对接成果区与收藏区：它们由 Molecular_Docking / Star Management 页管理，不在此提供整包下载
excluded = {os.path.basename(config.MD_RESULT_DIR), os.path.basename(config.MD_STAR_DIR)}
dirs_df = dirs_df[~dirs_df["Directory"].isin(excluded)]
if dirs_df.empty:
    st.info(t("data_mgmt.empty"))
else:
    # 表头
    col_h_name, col_h_time, col_h_btn = st.columns([3, 2, 1])
    col_h_name.markdown(t("common.col_directory"))
    col_h_time.markdown(t("common.col_created"))
    col_h_btn.markdown(t("common.col_download"))

    # 每行一个 Download 按钮，打包并下载该行目录的 zip
    with st.container(key="dir_rows"):
        for _, row in dirs_df.iterrows():
            dir_path = os.path.join(config.DOWNLOADS_DIR, row["Directory"])
            col_name, col_time, col_btn = st.columns([3, 2, 1])
            col_name.write(row["Directory"])
            col_time.write(row["Created time"])
            with col_btn:
                st.download_button(
                    t("common.download"),
                    Utils.make_dir_zip(dir_path),
                    f"{row['Directory']}.zip",
                    "application/zip",
                    key=f"dl_{row['Directory']}",
                    width="stretch",
                )
