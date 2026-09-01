import os
import streamlit as st
import pandas as pd
import Utils
from i18n import t

st.set_page_config(page_title=t("enrich_info.page_title"), layout="wide")

CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "enrichment",
)


# 隐藏侧边栏,全屏显示页面
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

st.title(t("enrich_info.title"))
st.caption(t("enrich_info.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    if st.button(t("common.back_main"), key="back_btn"):
        st.switch_page("App.py")

st.divider()

enr_files = Utils.list_cache_files("enrichment")
if enr_files.empty:
    st.info(t("enrich_info.empty"))
else:
    # 表头
    col_h_name, col_h_time, col_h_dl = st.columns([4, 2, 1])
    col_h_name.markdown(t("common.col_file_name"))
    col_h_time.markdown(t("common.col_created"))
    col_h_dl.markdown(t("common.col_download"))

    # 每行一个 Download 按钮,下载对应的富集分析结果 CSV
    with st.container(key="file_rows"):
        for _, row in enr_files.iterrows():
            col_name, col_time, col_dl = st.columns([4, 2, 1])
            col_name.write(row["File name"])
            col_time.write(row["Created time"])
            with col_dl:
                file_path = os.path.join(CACHE_DIR, row["File name"])
                if os.path.exists(file_path):
                    with open(file_path, "rb") as f:
                        csv_bytes = f.read()
                    st.download_button(
                        t("common.download"),
                        data=csv_bytes,
                        file_name=row["File name"],
                        mime="text/csv",
                        key=f"download_{row['File name']}",
                    )
