import os
import streamlit as st
import pandas as pd
import Utils

st.set_page_config(page_title="Drug Target Information Management", layout="wide")

CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "drug",
)
MERGED_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "merged",
)
# 展示/更新支持的文件名前缀：两个来源缓存 + 去重整合结果
PREFIXES = ("ChEMBL_", "DGIdb_", "merged_drug_")


def _file_path(filename: str) -> str:
    """merged_drug_ 前缀的文件位于 merged 目录，其余位于 drug 目录。"""
    if filename.startswith("merged_drug_"):
        return os.path.join(MERGED_DIR, filename)
    return os.path.join(CACHE_DIR, filename)


def _invalidate_merged_drug(drug_name: str) -> None:
    """更新来源文件后删除对应 merged 缓存，避免整合结果与来源不一致。"""
    path = os.path.join(MERGED_DIR, f"merged_drug_{drug_name}.csv")
    if os.path.exists(path):
        os.remove(path)


def refresh_file(filename: str, cache_hours: int = 1) -> dict:
    """按文件名前缀路由到对应抓取/整合函数，重新生成该文件。

    cache_hours : 仅 ChEMBL 来源生效（chembl 客户端 HTTP 缓存有效期）。
    """
    os.remove(_file_path(filename))

    if filename.startswith("ChEMBL_"):
        drug_name = filename[len("ChEMBL_"):-len(".csv")]
        Utils.set_chembl_cache_ttl(cache_hours)
        df, csv_path = Utils.fetch_drug_targets_chembl(drug_name)
        _invalidate_merged_drug(drug_name)
        if not df.empty and csv_path:
            return {"Status": "OK", "Message": f"{len(df)} targets"}
    elif filename.startswith("DGIdb_"):
        drug_name = filename[len("DGIdb_"):-len(".csv")]
        df = Utils.fetch_drug_targets_dgidb(drug_name)
        _invalidate_merged_drug(drug_name)
        if not df.empty:
            return {"Status": "OK", "Message": f"{len(df)} targets"}
    elif filename.startswith("merged_drug_"):
        drug_name = filename[len("merged_drug_"):-len(".csv")]
        df = Utils.get_drug_targets_merged(drug_name)
        if not df.empty:
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
            if not filename.startswith(PREFIXES) or not filename.endswith(".csv"):
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
merged_files = Utils.list_cache_files("merged")
merged_files = merged_files[merged_files["File name"].str.startswith("merged_drug_")]
drug_files = pd.concat([drug_files, merged_files], ignore_index=True)
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
                file_path = _file_path(row["File name"])
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
