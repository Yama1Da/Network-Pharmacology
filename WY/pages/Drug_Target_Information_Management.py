import logging
import os

import pandas as pd
import streamlit as st

import config
import Utils
from i18n import t

st.set_page_config(page_title=t("drug_mgmt.page_title"), layout="wide")

# 隐藏侧边栏,全屏显示页面;行内 Update 按钮(#6347A1 底、#E0D8F2 字)、Download 按钮(#207D62 底、#C8F0DF 字)
# 说明:ChEMBL 请求已改为批量分页,单个文件更新通常几十秒内完成;短时间内重复更新只会刷新时间戳
Utils.hide_sidebar("""
    .st-key-update_btn button,
    .st-key-file_rows [class*="st-key-update_"] button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }
    .st-key-file_rows [class*="st-key-download_"] button {
        background-color: #207D62;
        color: #C8F0DF;
        border: 1px solid #207D62;
    }
""")

# 各来源文件的抓取逻辑:返回 (表, 是否成功),供 Utils.refresh_cache_file 调用
def _fetch_chembl(name: str, cache_hours):
    Utils.set_chembl_cache_ttl(cache_hours)
    df, csv_path = Utils.fetch_drug_targets_chembl(name)
    ok = not df.empty and bool(csv_path) and os.path.exists(csv_path)
    return df, ok


def _fetch_dgidb(name: str, cache_hours):
    df = Utils.fetch_drug_targets_dgidb(name)
    # 文件名经过规范化处理,按同样规则检查新文件是否生成
    ok = not df.empty and os.path.exists(
        os.path.join(config.DRUG_CACHE_DIR, f"DGIdb_{Utils.safe_name(name)}.csv"))
    return df, ok


def _invalidate_merged_drug(drug_name: str) -> None:
    """更新来源文件后删除对应 merged 缓存，避免整合结果与来源不一致。"""
    # 传入的是原始药物名,文件名需走同样的 safe_name 规范化才能对上
    path = os.path.join(config.MERGED_CACHE_DIR,
                        f"{config.MERGED_DRUG_PREFIX}{Utils.safe_name(drug_name)}.csv")
    if os.path.exists(path):
        os.remove(path)


_SPEC = {
    "cache_dir": config.DRUG_CACHE_DIR,
    "merged_dir": config.MERGED_CACHE_DIR,
    "merged_prefix": config.MERGED_DRUG_PREFIX,
    "source_prefixes": ("ChEMBL_", "DGIdb_"),
    "fetch_source": {"ChEMBL_": _fetch_chembl, "DGIdb_": _fetch_dgidb},
    "fetch_merged": Utils.get_drug_targets_merged,
    "invalidate_merged": _invalidate_merged_drug,
}

st.title(t("drug_mgmt.title"))
st.caption(t("drug_mgmt.caption"))

# 展示上一次单文件更新的结果
last_update = st.session_state.pop("last_update", None)
if last_update:
    if last_update["Status"] == "OK":
        st.success(t("common.file_updated", file=last_update["File"], msg=last_update["Message"]))
    else:
        st.error(t("common.file_update_failed", file=last_update["File"], msg=last_update["Message"]))

col_back, col_update = st.columns([1, 1])
with col_back:
    Utils.back_button()
with col_update:
    all_update_clicked = st.button(t("common.all_updates"), key="update_btn")

if all_update_clicked:
    with st.spinner(t("drug_mgmt.updating_all")):
        results = Utils.update_all_cache_files(_SPEC, cache_hours=12)
    if results.empty:
        st.info(t("drug_mgmt.no_files_to_update"))
    else:
        st.dataframe(results, width="stretch", hide_index=True)
        if results["Status"].eq("OK").all():
            st.success(t("common.all_complete"))
        else:
            st.warning(t("common.some_failed"))

st.divider()

source_files = Utils.list_cache_files("drug")
merged_files = Utils.list_cache_files("merged")
merged_files = merged_files[merged_files["File name"].str.startswith("merged_drug_")]

# 数据来源筛选(默认全开):ChEMBL/DGIdb 作用于来源文件;merged 为两来源的整合缓存,单独可筛
st.caption(t("common.source_filter"))
c_ch, c_dg, c_m = st.columns(3)
with c_ch:
    show_chembl = st.checkbox("ChEMBL", value=True, key="flt_chembl")
with c_dg:
    show_dgidb = st.checkbox("DGIdb", value=True, key="flt_dgidb")
with c_m:
    show_merged = st.checkbox(t("common.filter_merged"), value=True, key="flt_merged")

def _match_source(n: str) -> bool:
    if n.startswith("ChEMBL_"):
        return show_chembl
    if n.startswith("DGIdb_"):
        return show_dgidb
    return True

n_source = len(source_files)
n_merged = len(merged_files)
if n_source:
    source_files = source_files[source_files["File name"].map(_match_source)]
if not show_merged:
    merged_files = merged_files.iloc[0:0]
drug_files = pd.concat([source_files, merged_files], ignore_index=True)
clicked_file = None

if drug_files.empty:
    # 有文件但全被筛掉时,按"筛选无结果"提示,与"目录本来就空"区分开
    st.info(t("common.filter_no_match") if (n_source or n_merged) else t("drug_mgmt.no_files"))
else:
    clicked_file = Utils.render_updateable_file_table(drug_files, _SPEC, search_key="search_drug")

# 单文件更新在页面渲染完成后执行,完成或出错都会记录结果并刷新页面
if clicked_file:
    with st.spinner(t("common.updating_file", file=clicked_file)):
        try:
            result = Utils.refresh_cache_file(clicked_file, _SPEC, cache_hours=1)
        except Exception:
            logging.exception("update failed for '%s'", clicked_file)
            result = {"Status": "FAILED", "Message": "update error"}
    st.session_state["last_update"] = {"File": clicked_file, **result}
    st.rerun()
