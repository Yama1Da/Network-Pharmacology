import os
import streamlit as st
import pandas as pd
import Utils

st.set_page_config(page_title="Drug Target Information Management", layout="wide")

CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "drug",
)
PREFIX = "ChEMBL_"


def refresh_file(filename: str, cache_hours: int = 1) -> dict:
    """删除旧缓存并直接调用 ChEMBL 获取方法，重新生成该文件（更新内容与创建时间）。

    cache_hours : chembl 客户端 HTTP 缓存有效期（小时），单行更新 1h / 全量更新 12h。
    """
    drug_name = filename[len(PREFIX):-len(".csv")]
    # 按调用方式调整 chembl HTTP 缓存 TTL，让更新尽量拿到较新数据
    Utils.set_chembl_cache_ttl(cache_hours)
    # 先删除旧文件，使重新生成的文件创建时间更新
    os.remove(os.path.join(CACHE_DIR, filename))
    df, csv_path = Utils.fetch_drug_targets_chembl(drug_name)
    if not df.empty and csv_path:
        return {"Status": "OK", "Message": f"{len(df)} targets"}
    return {"Status": "FAILED", "Message": "fetch failed"}


# ---------- Hide auto navigation sidebar for full-screen display ----------
#对于更新说明一下，由于ChEMBL对速率限制+防止被封禁，所以每次更新所等待的时间约莫在5分钟以上，且单行数据更新限制为1小时，全部数据更新限制为12小时，短时间内多次更新仅改变时间戳
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        display: none;
    }
    /* 行内 Update 按钮（#6347A1 底、#E0D8F2 字） */
    .st-key-update_btn button,
    .st-key-file_rows [class*="st-key-update_"] button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }
    /* 行内 Download 按钮（#207D62 底、#C8F0DF 字） */
    .st-key-file_rows [class*="st-key-download_"] button {
        background-color: #207D62;
        color: #C8F0DF;
        border: 1px solid #207D62;
    }
</style>
""", unsafe_allow_html=True)

st.title("💊 Drug Target Information Management")
st.caption("Cached files under Static/data/cache/drug")

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
    with st.spinner("Updating all cached drug targets…"):
        rows = []
        for filename in sorted(os.listdir(CACHE_DIR)):
            if not filename.startswith(PREFIX) or not filename.endswith(".csv"):
                continue
            try:
                rows.append({"File": filename, **refresh_file(filename, cache_hours=12)})
            except Exception:
                rows.append({"File": filename, "Status": "FAILED", "Message": "update error"})
        results = pd.DataFrame(rows, columns=["File", "Status", "Message"])
    if results.empty:
        st.info("No cached drug target files to update.")
    else:
        st.dataframe(results, width="stretch", hide_index=True)
        if results["Status"].eq("OK").all():
            st.success("✅ All updates complete")
        else:
            st.warning("Some files failed to update")

st.divider()

drug_files = Utils.list_cache_files("drug")
clicked_file = None

if drug_files.empty:
    st.info("No cached drug target files.")
else:
    # 表头
    col_h_name, col_h_time, col_h_btn, col_h_dl = st.columns([3, 2, 1, 1])
    col_h_name.markdown("**File name**")
    col_h_time.markdown("**Created time**")
    col_h_btn.markdown("**Update**")
    col_h_dl.markdown("**Download**")

    # 每行：Update（只更新该文件）+ Download（下载磁盘上的 CSV，保证持久化）
    with st.container(key="file_rows"):
        for _, row in drug_files.iterrows():
            col_name, col_time, col_btn, col_dl = st.columns([3, 2, 1, 1])
            col_name.write(row["File name"])
            col_time.write(row["Created time"])
            with col_btn:
                if st.button("Update", key=f"update_{row['File name']}"):
                    clicked_file = row["File name"]
            with col_dl:
                file_path = os.path.join(CACHE_DIR, row["File name"])
                if os.path.exists(file_path):
                    with open(file_path, "rb") as f:
                        csv_bytes = f.read()
                    st.download_button(
                        "Download",
                        data=csv_bytes,
                        file_name=row["File name"],
                        mime="text/csv",
                        key=f"download_{row['File name']}",
                    )

# 更新在容器外执行后再 rerun，避免旧行以残影形式残留在页面上
if clicked_file:
    with st.spinner(f"Updating {clicked_file}…"):
        result = refresh_file(clicked_file, cache_hours=1)
    st.session_state["last_update"] = {"File": clicked_file, **result}
    st.rerun()
