import os
import time
import streamlit as st

import config
from core.files import list_download_dirs, make_dir_zip
from ui.widgets import (
    back_button,
    collect_selected,
    hide_sidebar,
    render_file_action_bar,
    render_row_checkbox,
    render_select_all_checkbox,
    search_hit,
    search_input,
    show_file_op_result,
)
from i18n import t
from core.interactome import (
    build_interactome,
    clean_old_version_caches,
    interactome_age_days,
    is_built,
    is_stale,
    require_meta,
)

st.set_page_config(page_title=t("data_mgmt.page_title"), layout="wide")

# 隐藏侧边栏,全屏显示页面;Download 按钮(#207D62 底、#C8F0DF 字),删除按钮红底白字
hide_sidebar("""
    .st-key-dir_rows button {
        background-color: #207D62;
        color: #C8F0DF;
        border: 1px solid #207D62;
    }
    .st-key-data_dirs_del_btn button {
        background-color: #F53F3F;
        color: #FFFFFF;
        border: 1px solid #F53F3F;
    }
""")

st.title(t("data_mgmt.title"))
st.caption(t("data_mgmt.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    back_button("data_mgmt.back", target="pages/Data_Network.py")

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
        if result[4]:
            st.info(t("data_mgmt.itx_old_cleaned", n=result[4]))
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
    st.error(t("data_mgmt.itx_missing", ver=config.STRING_VERSION))

if st.button(t("data_mgmt.itx_build_btn"), key="itx_build_btn", type="primary", width="stretch"):
    progress = st.progress(0.0, text=t("data_mgmt.itx_stage_processing"))
    t0 = time.time()

    def _on_progress(stage: str, fraction, detail: str) -> None:
        # 并发分段不可用回退单流:进度条文字会被后续进度覆盖,另发一条常驻警告提醒慢
        if stage == "fallback":
            st.warning(t("data_mgmt.itx_fallback_notice", detail=detail))
            return
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
        # 新版本构建成功后清理其它版本/置信度组合的三件套(齿轮确认弹窗已告知该清理;
        # 只在成功后删,保证任何时刻至少有一个可用背景网络)
        cleaned = clean_old_version_caches("Human")
        st.session_state.itx_build_result = (
            "ok", meta_new["n_nodes"], meta_new["n_edges"], time.time() - t0, cleaned)
    except Exception as exc:
        st.session_state.itx_build_result = ("fail", str(exc), None, None)
    st.rerun()

st.divider()

dirs_df_all = list_download_dirs()
# 分子对接成果区和收藏区由各自页面管理,这里不提供下载
excluded = {os.path.basename(config.MD_RESULT_DIR), os.path.basename(config.MD_STAR_DIR)}
dirs_df_all = dirs_df_all[~dirs_df_all["Directory"].isin(excluded)]
# 批量删除结果提示(删除在确认弹窗内完成,结果跨 rerun 传递)
show_file_op_result("data_dirs")

query = search_input("search_dirs")
dirs_df = (dirs_df_all[dirs_df_all["Directory"].map(
    lambda d: search_hit(d, query))] if query else dirs_df_all)
if dirs_df.empty:
    st.info(t("common.search_no_match", query=query) if query else t("data_mgmt.empty"))
else:
    all_dirs = dirs_df_all["Directory"].tolist()
    # 表头(行首多选,支持批量下载/删除)
    col_h_sel, col_h_name, col_h_time, col_h_btn = st.columns([0.5, 2.5, 2, 1])
    with col_h_sel:
        select_all = render_select_all_checkbox("data_dirs", all_dirs)
    col_h_name.markdown(t("common.col_directory"))
    col_h_time.markdown(t("common.col_modified"))
    col_h_btn.markdown(t("common.col_download"))

    # 每行一个 Download 按钮，打包并下载该行目录的 zip
    with st.container(key="dir_rows"):
        for _, row in dirs_df.iterrows():
            dir_path = os.path.join(config.DOWNLOADS_DIR, row["Directory"])
            col_sel, col_name, col_time, col_btn = st.columns([0.5, 2.5, 2, 1])
            with col_sel:
                render_row_checkbox("data_dirs", row["Directory"], select_all)
            col_name.write(row["Directory"])
            col_time.write(row["Modified time"])
            with col_btn:
                st.download_button(
                    t("common.download"),
                    make_dir_zip(dir_path),
                    f"{row['Directory']}.zip",
                    "application/zip",
                    key=f"dl_{row['Directory']}",
                    width="stretch",
                )

    # 批量操作:下载选中(各目录打成一个 ZIP)/删除选中(确认弹窗)
    collect_selected("data_dirs", all_dirs)
    render_file_action_bar("data_dirs", config.DOWNLOADS_DIR, zip_name="Data_files.zip",
                                 dirs=True)
