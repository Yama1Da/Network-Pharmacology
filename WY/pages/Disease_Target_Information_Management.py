import os
import streamlit as st
import pandas as pd
import Utils

st.set_page_config(page_title="Disease Target Information Management", layout="wide")

CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "disease",
)
PREFIX = "OpenTargets_"


def refresh_file(filename: str) -> dict:
    """删除旧缓存并直接调用 Open Targets 获取方法，重新生成该文件（更新内容与创建时间）。"""
    disease_name = filename[len(PREFIX):-len(".csv")]
    # 先删除旧文件，使重新生成的文件创建时间更新
    os.remove(os.path.join(CACHE_DIR, filename))
    # 绕过 st.cache_data，强制重新请求 API
    Utils.get_disease_targets.clear()
    df = Utils.get_disease_targets(disease_name)
    if not df.empty:
        return {"Status": "OK", "Message": f"{len(df)} targets"}
    return {"Status": "FAILED", "Message": "fetch failed"}


# ---------- Hide auto navigation sidebar for full-screen display ----------
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        display: none;
    }
    /* Update buttons colors (#6347A1 bg, #E0D8F2 text) */
    .st-key-update_btn button,
    .st-key-file_rows button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }
</style>
""", unsafe_allow_html=True)

st.title("🦠 Disease Target Information Management")
st.caption("Cached files under Static/data/cache/disease")

# 单文件更新结果（st.rerun 后仍展示）
last_update = st.session_state.pop("last_update", None)
if last_update:
    if last_update["Status"] == "OK":
        st.success(f"✅ {last_update['File']} updated ({last_update['Message']})")
    else:
        st.error(f"❌ {last_update['File']} update failed: {last_update['Message']}")

col_back, col_update = st.columns([1, 1])
with col_back:
    if st.button("← Back to Main", key="back_btn"):
        st.switch_page("App.py")
with col_update:
    all_update_clicked = st.button("All Updates", key="update_btn")

if all_update_clicked:
    with st.spinner("Updating all cached disease targets…"):
        rows = []
        for filename in sorted(os.listdir(CACHE_DIR)):
            if not filename.startswith(PREFIX) or not filename.endswith(".csv"):
                continue
            try:
                rows.append({"File": filename, **refresh_file(filename)})
            except Exception:
                rows.append({"File": filename, "Status": "FAILED", "Message": "update error"})
        results = pd.DataFrame(rows, columns=["File", "Status", "Message"])
    if results.empty:
        st.info("No cached disease target files to update.")
    else:
        st.dataframe(results, use_container_width=True, hide_index=True)
        if results["Status"].eq("OK").all():
            st.success("✅ All updates complete")
        else:
            st.warning("Some files failed to update")

st.divider()

disease_files = Utils.list_cache_files("disease")
if disease_files.empty:
    st.info("No cached disease target files.")
else:
    # 表头
    col_h_name, col_h_time, col_h_btn = st.columns([3, 2, 1])
    col_h_name.markdown("**File name**")
    col_h_time.markdown("**Created time**")
    col_h_btn.markdown("**Update**")

    # 每行一个 Update 按钮，只更新该行对应的文件
    with st.container(key="file_rows"):
        for _, row in disease_files.iterrows():
            col_name, col_time, col_btn = st.columns([3, 2, 1])
            col_name.write(row["File name"])
            col_time.write(row["Created time"])
            with col_btn:
                if st.button("Update", key=f"update_{row['File name']}"):
                    with st.spinner(f"Updating {row['File name']}…"):
                        result = refresh_file(row["File name"])
                    st.session_state["last_update"] = {"File": row["File name"], **result}
                    st.rerun()
