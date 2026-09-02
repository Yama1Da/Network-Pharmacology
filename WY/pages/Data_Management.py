import os
import time
import streamlit as st

import config
import Utils
from i18n import t
from core.interactome import (
    build_interactome,
    interactome_age_days,
    is_built,
    is_stale,
    require_meta,
)

st.set_page_config(page_title=t("data_mgmt.page_title"), layout="wide")

# 隐藏侧边栏,全屏显示页面
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

# 背景互作网络:显示构建状态,支持一键重建
# 三种状态:未构建(错误提醒) / 已构建但数据陈旧(警告) / 已构建(正常)
st.subheader(t("data_mgmt.itx_title"))

if "itx_build_result" not in st.session_state:
    st.session_state.itx_build_result = None  # 保存本次构建结果(成功或失败信息)

if st.session_state.itx_build_result is not None:
    result = st.session_state.itx_build_result
    if result[0] == "ok":
        st.success(t("data_mgmt.itx_build_done", nodes=result[1], edges=result[2], secs=result[3]))
    else:
        st.error(t("data_mgmt.itx_build_fail", err=result[1]))
    st.session_state.itx_build_result = None

itx_built = is_built("Human")
if itx_built:
    meta = require_meta("Human")
    age = interactome_age_days("Human")
    st.info(t("data_mgmt.itx_built",
              nodes=meta["n_nodes"], edges=meta["n_edges"],
              ver=meta["string_version"], conf=meta["confidence"],
              date=str(meta.get("built_at", ""))[:10], age=age if age is not None else "-"))
    if is_stale("Human"):
        st.warning(t("data_mgmt.itx_stale", age=age))
else:
    st.error(t("data_mgmt.itx_missing"))

if st.button(t("data_mgmt.itx_build_btn"), key="itx_build_btn", type="primary", width="stretch"):
    progress = st.progress(0.0, text=t("data_mgmt.itx_stage_processing"))
    t0 = time.time()

    def _on_progress(stage: str, fraction, detail: str) -> None:
        # 把构建阶段名翻译成当前语言,显示在进度条上
        base = t(f"data_mgmt.itx_stage_{stage}")
        if detail == "cached":
            detail_txt = t("data_mgmt.itx_cached")
        else:
            detail_txt = detail
        text = f"{base} {detail_txt}".strip()
        progress.progress(min(max(float(fraction or 0.0), 0.0), 1.0), text=text)

    try:
        meta_new = build_interactome(progress_cb=_on_progress)
        st.session_state.itx_build_result = (
            "ok", meta_new["n_nodes"], meta_new["n_edges"], time.time() - t0)
    except Exception as exc:
        st.session_state.itx_build_result = ("fail", str(exc), None, None)
    st.rerun()

st.divider()

dirs_df = Utils.list_download_dirs()
# 分子对接成果区和收藏区由各自页面管理,这里不提供下载
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
