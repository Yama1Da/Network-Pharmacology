import os
import shutil
from datetime import datetime

import streamlit as st

import config
import Utils
from i18n import t

st.set_page_config(page_title=t("star.page_title"), layout="wide")

# 隐藏侧边栏,全屏显示页面;下载按钮(#00B42A 底白字)、删除按钮(#F53F3F 底白字)
Utils.hide_sidebar("""
    [class*="st-key-star_dl_"] button {
        background-color: #00B42A;
        color: #FFFFFF;
        border: 1px solid #00B42A;
    }
    [class*="st-key-star_del_"] button {
        background-color: #F53F3F;
        color: #FFFFFF;
        border: 1px solid #F53F3F;
    }
""")

st.title(t("star.title"))
st.caption(t("star.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    Utils.back_button("star.back", target="pages/Molecular_Docking.py")

st.divider()


@st.dialog(t("star.dialog_title"))
def confirm_delete() -> None:
    """删除确认弹窗:点"是"才删除,点"否"只关闭窗口。"""
    name = st.session_state.get("_star_pending", "")
    path = st.session_state.get("_star_pending_path", "")
    st.write(t("star.dialog_confirm", name=name))
    col_yes, col_no = st.columns(2)
    with col_yes:
        if st.button(t("common.yes"), key="dlg_yes", width="stretch"):
            shutil.rmtree(path, ignore_errors=True)
            if os.path.isdir(path):
                # rmtree 静默失败(如 Windows 目录被占用):按失败提示,不虚报成功
                st.session_state["_star_delete_failed"] = name
            else:
                st.session_state["_star_deleted"] = name
            st.rerun()  # 关闭弹窗并刷新列表
    with col_no:
        if st.button(t("common.no"), key="dlg_no", width="stretch"):
            st.rerun()  # 仅关闭弹窗


if st.session_state.pop("_star_deleted", None):
    st.success(t("star.deleted"))
failed_name = st.session_state.pop("_star_delete_failed", None)
if failed_name:
    st.error(t("common.delete_failed", names=failed_name))

star_dir = config.MD_STAR_DIR
all_names = []
if os.path.isdir(star_dir):
    all_names = sorted(n for n in os.listdir(star_dir)
                       if os.path.isdir(os.path.join(star_dir, n)))
# 批量删除结果提示(多选删除在确认弹窗内完成,结果跨 rerun 传递)
Utils.show_file_op_result("star_dirs")

query = Utils.search_input("search_star")
names = [n for n in all_names if Utils.search_hit(n, query)] if query else all_names
if not names:
    st.info(t("common.search_no_match", query=query) if query else t("star.empty"))
else:
    st.caption(t("star.total", count=len(names)))

    # 表头(行首多选,支持批量下载/删除)
    col_h_sel, col_h_name, col_h_time, col_h_dl, col_h_del = st.columns([0.5, 3.5, 2, 1, 1])
    with col_h_sel:
        select_all = Utils.render_select_all_checkbox("star_dirs", all_names)
    col_h_name.markdown(t("star.col_name"))
    col_h_time.markdown(t("star.col_time"))
    col_h_dl.markdown(t("common.col_download"))
    col_h_del.markdown(t("star.col_delete"))

    for name in names:
        path = os.path.join(star_dir, name)
        # 用 mtime 而非 getctime:后者在 Linux 上是 inode 变更时间,跨平台语义不一致;
        # 收藏目录创建后内容不再变化,mtime 即收藏时间
        created = datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")
        col_sel, col_a, col_b, col_c, col_d = st.columns([0.5, 3.5, 2, 1, 1])
        with col_sel:
            Utils.render_row_checkbox("star_dirs", name, select_all)
        col_a.write(name)
        col_b.write(created)
        with col_c:
            st.download_button("⬇️", data=Utils.make_dir_zip(path),
                               file_name=name + ".zip", mime="application/zip",
                               key=f"star_dl_{name}")
        with col_d:
            if st.button("🗑️", key=f"star_del_{name}"):
                st.session_state["_star_pending"] = name
                st.session_state["_star_pending_path"] = path
                confirm_delete()

    # 批量操作:下载选中(各目录打成一个 ZIP)/删除选中(确认弹窗)
    Utils.collect_selected("star_dirs", all_names)
    Utils.render_file_action_bar("star_dirs", star_dir, zip_name="Starred_results.zip",
                                 dirs=True, confirm_key="star.dialog_confirm_batch")
