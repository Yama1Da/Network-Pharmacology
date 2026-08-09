import os

import streamlit as st

import config
import Utils

# 初始化日志系统（幂等：直接访问本页 URL 时也保证日志就绪）
Utils.setup_logging()

st.set_page_config(page_title="Logs", layout="wide")

# ---------- Hide auto navigation sidebar for full-screen display ----------
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        display: none;
    }
</style>
""", unsafe_allow_html=True)

st.title("📜 Logs")
st.caption("Local log files under WY/Static/log — click a file name to view its content")

col_back, _ = st.columns([1, 4])
with col_back:
    if st.button("← Back to Main", key="back_btn"):
        st.switch_page("App.py")

st.divider()

log_files = Utils.list_log_files()
if log_files.empty:
    st.info("No log files yet. Run the app to generate daily logs.")
else:
    # 每个日志文件一个展开区：点击文件名即查看内容
    # app.log 是当前正在写入的文件（今天的日志），其余带日期后缀为历史文件
    for _, row in log_files.iterrows():
        file_path = os.path.join(config.LOG_DIR, row["File name"])
        is_current = row["File name"] == "app.log"
        tag = " 📌 current (today)" if is_current else ""
        label = f"{row['File name']}{tag}  ·  modified {row['Modified time']}"
        with st.expander(label):
            if not os.path.exists(file_path):
                st.info("File no longer exists on disk.")
            else:
                try:
                    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
                        content = f.read()
                    st.code(content, language=None)
                except OSError as e:
                    st.error(f"Failed to read {row['File name']}: {e}")
