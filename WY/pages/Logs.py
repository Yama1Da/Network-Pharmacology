import os

import streamlit as st

import config
import Utils
from i18n import t

# 初始化日志系统
Utils.setup_logging()

st.set_page_config(page_title=t("logs.page_title"), layout="wide")

# 隐藏侧边栏,全屏显示页面
Utils.hide_sidebar()

st.title(t("logs.title"))
st.caption(t("logs.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    Utils.back_button()

st.divider()

# 运行中可临时开启 DEBUG;正常使用仍默认 INFO,避免日志过大
if "logs_debug_enabled" not in st.session_state:
    st.session_state["logs_debug_enabled"] = config.LOG_LEVEL.upper() == "DEBUG"
level_name = "DEBUG" if st.session_state["logs_debug_enabled"] else "INFO"

col_level, _ = st.columns([1, 3])
with col_level:
    st.markdown(t("logs.level_title"))
    st.toggle(
        t("logs.debug_toggle"),
        key="logs_debug_enabled",
        help=t("logs.debug_toggle_help"),
    )
    Utils.set_logging_level(level_name)
    st.caption(t("logs.level_current", level=level_name))

st.divider()

log_files = Utils.list_log_files()
if log_files.empty:
    st.info(t("logs.empty"))
else:
    # 每个日志文件一个展开区,点击查看内容
    # app.log 是今天的日志,其余为历史日志
    for _, row in log_files.iterrows():
        file_path = os.path.join(config.LOG_DIR, row["File name"])
        is_current = row["File name"] == "app.log"
        tag = t("logs.tag_current") if is_current else ""
        label = f"{row['File name']}{tag}  ·  {t('logs.modified')} {row['Modified time']}"
        with st.expander(label):
            if not os.path.exists(file_path):
                st.info(t("logs.gone"))
            else:
                try:
                    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
                        content = f.read()
                    st.code(content, language=None)
                except OSError as e:
                    st.error(t("logs.read_failed", name=row["File name"], err=e))
