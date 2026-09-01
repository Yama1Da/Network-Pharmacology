import os
import logging
import streamlit as st
import pandas as pd
import Utils
from i18n import t

st.set_page_config(page_title=t("drug_mgmt.page_title"), layout="wide")

CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "drug",
)
MERGED_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "merged",
)


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

    - 来源文件（ChEMBL_ / DGIdb_）：先备份旧文件，抓取成功后才替换；
      失败时恢复旧文件并返回 FAILED，避免数据丢失。
    - merged 文件（merged_drug_）：先更新两个来源文件，再重建整合结果；
      任一来源更新失败则保留旧 merged 并返回 FAILED。

    cache_hours : 仅 ChEMBL 来源生效（chembl 客户端 HTTP 缓存有效期）。
    """
    # ===== merged 文件：先更新来源，再重建 =====
    if filename.startswith("merged_drug_"):
        drug_name = filename[len("merged_drug_"):-len(".csv")]
        for source in (f"ChEMBL_{drug_name}.csv", f"DGIdb_{drug_name}.csv"):
            result = refresh_file(source, cache_hours)
            if result["Status"] != "OK":
                return {"Status": "FAILED",
                        "Message": f"source '{source}' update failed: {result['Message']}"}
        merged_path = _file_path(filename)
        if os.path.exists(merged_path):
            os.remove(merged_path)  # 来源更新时可能已失效删除，这里兜底
        df = Utils.get_drug_targets_merged(drug_name)
        if not df.empty:
            return {"Status": "OK", "Message": f"{len(df)} targets"}
        return {"Status": "FAILED", "Message": "fetch failed"}

    # ===== 来源文件：备份 → 抓取 → 成功替换 / 失败恢复 =====
    path = _file_path(filename)
    backup = path + ".bak"
    if os.path.exists(path):
        os.replace(path, backup)

    try:
        if filename.startswith("ChEMBL_"):
            drug_name = filename[len("ChEMBL_"):-len(".csv")]
            Utils.set_chembl_cache_ttl(cache_hours)
            df, csv_path = Utils.fetch_drug_targets_chembl(drug_name)
            ok = not df.empty and csv_path and os.path.exists(csv_path)
        elif filename.startswith("DGIdb_"):
            drug_name = filename[len("DGIdb_"):-len(".csv")]
            df = Utils.fetch_drug_targets_dgidb(drug_name)
            # 新缓存一律以净化命名落盘，存在性判断须与 Utils.safe_name 一致
            ok = not df.empty and os.path.exists(
                _file_path(f"DGIdb_{Utils.safe_name(drug_name)}.csv"))
        else:
            raise ValueError(f"unknown file prefix: {filename}")
    except Exception:
        logging.exception("refresh_file failed for '%s'", filename)
        ok = False

    if ok:
        if os.path.exists(backup):
            os.remove(backup)
        _invalidate_merged_drug(drug_name)
        return {"Status": "OK", "Message": f"{len(df)} targets"}
    # 失败：恢复旧文件，保持旧整合结果有效
    if os.path.exists(backup):
        os.replace(backup, path)
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

st.title(t("drug_mgmt.title"))
st.caption(t("drug_mgmt.caption"))

# 单文件更新结果（st.rerun 后仍展示）
last_update = st.session_state.pop("last_update", None)
if last_update:
    if last_update["Status"] == "OK":
        st.success(t("common.file_updated", file=last_update["File"], msg=last_update["Message"]))
    else:
        st.error(t("common.file_update_failed", file=last_update["File"], msg=last_update["Message"]))

col_back, col_update = st.columns([1, 1])
with col_back:
    if st.button(t("common.back_main"), key="back_btn"):
        st.switch_page("App.py")
with col_update:
    all_update_clicked = st.button(t("common.all_updates"), key="update_btn")

if all_update_clicked:
    with st.spinner(t("drug_mgmt.updating_all")):
        rows = []
        # 第一阶段：先更新所有来源文件（ChEMBL_ / DGIdb_）
        drug_names = set()
        for filename in sorted(os.listdir(CACHE_DIR)):
            if not filename.endswith(".csv"):
                continue
            for prefix in ("ChEMBL_", "DGIdb_"):
                if filename.startswith(prefix):
                    drug_names.add(filename[len(prefix):-len(".csv")])
                    try:
                        rows.append({"File": filename, **refresh_file(filename, cache_hours=12)})
                    except Exception:
                        rows.append({"File": filename, "Status": "FAILED", "Message": "update error"})
                    break
        # 第二阶段：再重建所有 merged 文件（来源已在第一阶段刷新，此处直接重建）
        for drug_name in sorted(drug_names):
            merged_file = f"merged_drug_{drug_name}.csv"
            merged_path = os.path.join(MERGED_DIR, merged_file)
            if os.path.exists(merged_path):
                os.remove(merged_path)
            try:
                df = Utils.get_drug_targets_merged(drug_name)
                if not df.empty:
                    rows.append({"File": merged_file, "Status": "OK", "Message": f"{len(df)} targets"})
                else:
                    rows.append({"File": merged_file, "Status": "FAILED", "Message": "fetch failed"})
            except Exception:
                rows.append({"File": merged_file, "Status": "FAILED", "Message": "update error"})
        results = pd.DataFrame(rows, columns=["File", "Status", "Message"])
    if results.empty:
        st.info(t("drug_mgmt.no_files_to_update"))
    else:
        st.dataframe(results, width="stretch", hide_index=True)
        if results["Status"].eq("OK").all():
            st.success(t("common.all_complete"))
        else:
            st.warning(t("common.some_failed"))

st.divider()

drug_files = Utils.list_cache_files("drug")
merged_files = Utils.list_cache_files("merged")
merged_files = merged_files[merged_files["File name"].str.startswith("merged_drug_")]
drug_files = pd.concat([drug_files, merged_files], ignore_index=True)
clicked_file = None

if drug_files.empty:
    st.info(t("drug_mgmt.no_files"))
else:
    # 表头
    col_h_name, col_h_time, col_h_btn, col_h_dl = st.columns([3, 2, 1, 1])
    col_h_name.markdown(t("common.col_file_name"))
    col_h_time.markdown(t("common.col_created"))
    col_h_btn.markdown(t("common.col_update"))
    col_h_dl.markdown(t("common.col_download"))

    # 每行：Update（只更新该文件）+ Download（下载磁盘上的 CSV，保证持久化）
    with st.container(key="file_rows"):
        for _, row in drug_files.iterrows():
            col_name, col_time, col_btn, col_dl = st.columns([3, 2, 1, 1])
            col_name.write(row["File name"])
            col_time.write(row["Created time"])
            with col_btn:
                if st.button(t("common.update"), key=f"update_{row['File name']}"):
                    clicked_file = row["File name"]
            with col_dl:
                file_path = _file_path(row["File name"])
                if os.path.exists(file_path):
                    with open(file_path, "rb") as f:
                        csv_bytes = f.read()
                    st.download_button(
                        t("common.download"),
                        data=csv_bytes,
                        file_name=row["File name"],
                        mime="text/csv",
                        key=f"download_{row['File name']}",
                    )

# 更新在容器外执行后再 rerun，避免旧行以残影形式残留在页面上
# try/except 兜底：refresh_file 内部异常（如整合函数抛错）也会记录结果并 rerun，
# 避免页面停留在"运行中"的半渲染状态
if clicked_file:
    with st.spinner(t("common.updating_file", file=clicked_file)):
        try:
            result = refresh_file(clicked_file, cache_hours=1)
        except Exception:
            logging.exception("update failed for '%s'", clicked_file)
            result = {"Status": "FAILED", "Message": "update error"}
    st.session_state["last_update"] = {"File": clicked_file, **result}
    st.rerun()
