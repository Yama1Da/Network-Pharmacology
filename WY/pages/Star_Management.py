import os
import shutil
from datetime import datetime

import streamlit as st

import config
import Utils

st.set_page_config(page_title="Star Management", layout="wide")

# ---------- Hide auto navigation sidebar for full-screen display ----------
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        display: none;
    }
    /* 下载按钮：#00B42A 底、白字 */
    [class*="st-key-star_dl_"] button {
        background-color: #00B42A;
        color: #FFFFFF;
        border: 1px solid #00B42A;
    }
    /* 删除按钮：#F53F3F 底、白字 */
    [class*="st-key-star_del_"] button {
        background-color: #F53F3F;
        color: #FFFFFF;
        border: 1px solid #F53F3F;
    }
</style>
""", unsafe_allow_html=True)

st.title("⭐ Star Management")
st.caption("收藏的对接成果（MDresult/star/ 区）— 不受自动清理影响，仅可手动删除")

col_back, _ = st.columns([1, 4])
with col_back:
    if st.button("← Back to Molecular Docking", key="back_btn"):
        st.switch_page("pages/Molecular_Docking.py")

st.divider()


@st.dialog("删除收藏")
def confirm_delete() -> None:
    """删除确认悬浮窗：选择"是"才删除，"否"关闭窗口。待删项经 session_state 传递。"""
    name = st.session_state.get("_star_pending", "")
    path = st.session_state.get("_star_pending_path", "")
    st.write(f"确定要删除收藏「{name}」吗？删除后不可恢复。")
    col_yes, col_no = st.columns(2)
    with col_yes:
        if st.button("是", key="dlg_yes", width="stretch"):
            shutil.rmtree(path, ignore_errors=True)
            st.session_state["_star_deleted"] = name
            st.rerun()  # rerun 关闭对话框并刷新列表
    with col_no:
        if st.button("否", key="dlg_no", width="stretch"):
            st.rerun()  # 仅关闭对话框


if st.session_state.pop("_star_deleted", None):
    st.success("已删除该收藏")

star_dir = config.MD_STAR_DIR
if not os.path.isdir(star_dir) or not any(
        os.path.isdir(os.path.join(star_dir, n)) for n in os.listdir(star_dir)):
    st.info("暂无收藏。在分子对接页保存成果后点击 ☆ 即可收藏。")
else:
    names = sorted(n for n in os.listdir(star_dir)
                   if os.path.isdir(os.path.join(star_dir, n)))
    st.caption(f"共 {len(names)} 份收藏。下载到本地保存，删除后不可恢复。")

    # 表头
    col_h_name, col_h_time, col_h_dl, col_h_del = st.columns([4, 2, 1, 1])
    col_h_name.markdown("**收藏名称**")
    col_h_time.markdown("**收藏时间**")
    col_h_dl.markdown("**下载**")
    col_h_del.markdown("**删除**")

    for name in names:
        path = os.path.join(star_dir, name)
        created = datetime.fromtimestamp(os.path.getctime(path)).strftime("%Y-%m-%d %H:%M:%S")
        col_a, col_b, col_c, col_d = st.columns([4, 2, 1, 1])
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
