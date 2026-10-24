import streamlit as st

from WY import config
from WY.core.files import list_cache_files
from WY.ui.widgets import (
    CSS_BTN_DOWNLOAD,
    back_button,
    hide_sidebar,
    render_file_action_bar,
    render_file_table,
    show_file_op_result,
)
from WY.i18n import t

st.set_page_config(page_title=t("ppi_info.page_title"), layout="wide")

# 隐藏侧边栏,全屏显示页面;行内 Download 按钮绿底、删除按钮红底白字
hide_sidebar(CSS_BTN_DOWNLOAD + """
    .st-key-ppi_files_del_btn button {
        background-color: #F53F3F;
        color: #FFFFFF;
        border: 1px solid #F53F3F;
    }
""")

st.title(t("ppi_info.title"))
st.caption(t("ppi_info.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    back_button()

st.divider()

show_file_op_result("ppi_files")

ppi_files = list_cache_files("ppi")
if ppi_files.empty:
    st.info(t("ppi_info.empty"))
else:
    # 多选文件表 + 批量操作:下载选中(ZIP)/删除(确认弹窗)
    render_file_table(ppi_files, config.PPI_CACHE_DIR, search_key="search_ppi",
                            ns="ppi_files")
    render_file_action_bar("ppi_files", config.PPI_CACHE_DIR, zip_name="PPI_files.zip")
