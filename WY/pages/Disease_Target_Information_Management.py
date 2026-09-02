import os
import logging
import streamlit as st
import pandas as pd
import Utils
from i18n import t

st.set_page_config(page_title=t("disease_mgmt.page_title"), layout="wide")

CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "disease",
)
MERGED_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "merged",
)


def _file_path(filename: str) -> str:
    """merged_disease_ 前缀的文件位于 merged 目录，其余位于 disease 目录。"""
    if filename.startswith("merged_disease_"):
        return os.path.join(MERGED_DIR, filename)
    return os.path.join(CACHE_DIR, filename)


def _invalidate_merged_disease(disease_name: str) -> None:
    """更新来源文件后删除对应 merged 缓存，避免整合结果与来源不一致。"""
    path = os.path.join(MERGED_DIR, f"merged_disease_{disease_name}.csv")
    if os.path.exists(path):
        os.remove(path)


def refresh_file(filename: str) -> dict:
    """按文件名重新抓取或整合数据,更新指定文件。

    来源文件(OpenTargets_/DisGeNET_)先备份旧文件,抓取成功才替换,失败则恢复旧文件,避免数据丢失;
    merged 整合文件先更新两个来源文件再重建,任一来源失败则保留旧的整合结果。
    """
    # merged 整合文件:先更新两个来源文件,再重建
    if filename.startswith("merged_disease_"):
        disease_name = filename[len("merged_disease_"):-len(".csv")]
        for source in (f"OpenTargets_{disease_name}.csv", f"DisGeNET_{disease_name}.csv"):
            result = refresh_file(source)
            if result["Status"] != "OK":
                return {"Status": "FAILED",
                        "Message": f"source '{source}' update failed: {result['Message']}"}
        merged_path = _file_path(filename)
        if os.path.exists(merged_path):
            os.remove(merged_path)  # 可能已在上一步删除,这里兜底再删一次
        df = Utils.get_disease_targets_merged(disease_name)
        if not df.empty:
            return {"Status": "OK", "Message": f"{len(df)} targets"}
        return {"Status": "FAILED", "Message": "fetch failed"}

    # 来源文件:先备份,抓取成功才替换,失败则恢复旧文件
    path = _file_path(filename)
    backup = path + ".bak"
    if os.path.exists(path):
        os.replace(path, backup)

    try:
        if filename.startswith("OpenTargets_"):
            disease_name = filename[len("OpenTargets_"):-len(".csv")]
            # 不使用缓存,强制重新请求 API
            Utils.get_disease_targets.clear()
            df = Utils.get_disease_targets(disease_name)
            # 文件名经过规范化处理,按同样规则检查新文件是否生成
            ok = not df.empty and os.path.exists(
                _file_path(f"OpenTargets_{Utils.safe_name(disease_name)}.csv"))
        elif filename.startswith("DisGeNET_"):
            disease_name = filename[len("DisGeNET_"):-len(".csv")]
            df = Utils.fetch_disease_targets_disgenet(disease_name)
            ok = not df.empty and os.path.exists(
                _file_path(f"DisGeNET_{Utils.safe_name(disease_name)}.csv"))
        else:
            raise ValueError(f"unknown file prefix: {filename}")
    except Exception:
        logging.exception("refresh_file failed for '%s'", filename)
        ok = False

    if ok:
        if os.path.exists(backup):
            os.remove(backup)
        _invalidate_merged_disease(disease_name)
        return {"Status": "OK", "Message": f"{len(df)} targets"}
    # 失败:恢复旧文件
    if os.path.exists(backup):
        os.replace(backup, path)
    return {"Status": "FAILED", "Message": "fetch failed"}


# 隐藏侧边栏,全屏显示页面
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

st.title(t("disease_mgmt.title"))
st.caption(t("disease_mgmt.caption"))

# 展示上一次单文件更新的结果
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
    with st.spinner(t("disease_mgmt.updating_all")):
        rows = []
        # 第一步:先更新所有来源文件
        disease_names = set()
        for filename in sorted(os.listdir(CACHE_DIR)):
            if not filename.endswith(".csv"):
                continue
            for prefix in ("OpenTargets_", "DisGeNET_"):
                if filename.startswith(prefix):
                    disease_names.add(filename[len(prefix):-len(".csv")])
                    try:
                        rows.append({"File": filename, **refresh_file(filename)})
                    except Exception:
                        rows.append({"File": filename, "Status": "FAILED", "Message": "update error"})
                    break
        # 第二步:再重建所有 merged 整合文件
        for disease_name in sorted(disease_names):
            merged_file = f"merged_disease_{disease_name}.csv"
            merged_path = os.path.join(MERGED_DIR, merged_file)
            if os.path.exists(merged_path):
                os.remove(merged_path)
            try:
                df = Utils.get_disease_targets_merged(disease_name)
                if not df.empty:
                    rows.append({"File": merged_file, "Status": "OK", "Message": f"{len(df)} targets"})
                else:
                    rows.append({"File": merged_file, "Status": "FAILED", "Message": "fetch failed"})
            except Exception:
                rows.append({"File": merged_file, "Status": "FAILED", "Message": "update error"})
        results = pd.DataFrame(rows, columns=["File", "Status", "Message"])
    if results.empty:
        st.info(t("disease_mgmt.no_files_to_update"))
    else:
        st.dataframe(results, width="stretch", hide_index=True)
        if results["Status"].eq("OK").all():
            st.success(t("common.all_complete"))
        else:
            st.warning(t("common.some_failed"))

st.divider()

disease_files = Utils.list_cache_files("disease")
merged_files = Utils.list_cache_files("merged")
merged_files = merged_files[merged_files["File name"].str.startswith("merged_disease_")]
disease_files = pd.concat([disease_files, merged_files], ignore_index=True)
clicked_file = None

if disease_files.empty:
    st.info(t("disease_mgmt.no_files"))
else:
    # 表头
    col_h_name, col_h_time, col_h_btn, col_h_dl = st.columns([3, 2, 1, 1])
    col_h_name.markdown(t("common.col_file_name"))
    col_h_time.markdown(t("common.col_created"))
    col_h_btn.markdown(t("common.col_update"))
    col_h_dl.markdown(t("common.col_download"))

    # 每行:Update 只更新该文件,Download 下载该 CSV
    with st.container(key="file_rows"):
        for _, row in disease_files.iterrows():
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

# 单文件更新在页面渲染完成后执行,完成或出错都会记录结果并刷新页面
if clicked_file:
    with st.spinner(t("common.updating_file", file=clicked_file)):
        try:
            result = refresh_file(clicked_file)
        except Exception:
            logging.exception("update failed for '%s'", clicked_file)
            result = {"Status": "FAILED", "Message": "update error"}
    st.session_state["last_update"] = {"File": clicked_file, **result}
    st.rerun()
