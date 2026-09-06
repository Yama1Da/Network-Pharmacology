import streamlit as st

import config
import Utils
from i18n import t

st.set_page_config(page_title=t("enrich_info.page_title"), layout="wide")

# 隐藏侧边栏,全屏显示页面;行内 Download 按钮(#207D62 底、#C8F0DF 字)
Utils.hide_sidebar("""
    .st-key-file_rows [class*="st-key-download_"] button {
        background-color: #207D62;
        color: #C8F0DF;
        border: 1px solid #207D62;
    }
""")

st.title(t("enrich_info.title"))
st.caption(t("enrich_info.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    Utils.back_button()

st.divider()

enr_files = Utils.list_cache_files("enrichment")
if enr_files.empty:
    st.info(t("enrich_info.empty"))
else:
    # 每行一个 Download 按钮,下载对应的富集分析结果 CSV
    Utils.render_file_table(enr_files, config.ENRICHMENT_CACHE_DIR, search_key="search_enrich")
