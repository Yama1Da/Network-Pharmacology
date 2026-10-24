"""ML 结果管理页(方案 §10 增补:结果保留语义与治理)。

ML 的每次分析都会生成一个 <tag>_<时间戳> 结果目录,只增不减 —— 本页提供
治理能力:双根扫描列出全部结果(CLI 的 results/mlscreen + 页面的 GEOresult)、
关键词筛选、勾选删除(二次确认)、打包下载(整个结果目录 zip)。

数据层在 core/ml_results(列表/删除/打包全部带校验与留痕),本页只做渲染与
调用,不含文件逻辑 —— 前端后续改写 TS 时以 core 为 API 边界。
"""
import os
import time

import pandas as pd
import streamlit as st

from WY.core.ml_results import delete_ml_runs, list_ml_runs, zip_ml_runs
from WY.i18n import t
from WY.ui.widgets import back_button, hide_sidebar

st.set_page_config(page_title=t("mri.title"), layout="wide")

hide_sidebar()

st.title(t("mri.title"))
st.caption(t("mri.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    back_button()

st.divider()

# ---------- 全量列表(数据层已按时间新→旧排序) ----------
runs = list_ml_runs()
total_mb = sum(r["size_bytes"] for r in runs) / 1048576
st.caption(t("mri.total", n=len(runs), size=f"{total_mb:.0f}"))

if not runs:
    st.info(t("mri.empty"))
    st.stop()

# 关键词筛选(运行名/来源,大小写不敏感)
query = st.text_input(t("mri.search"), key="mri_search").strip().lower()
if query:
    runs = [r for r in runs
            if query in r["name"].lower()
            or query in str(r.get("gse") or "").lower()
            or query in str(r.get("source") or "").lower()]

# 组装展示表(列头复用 eml 词条;路径映射表与展示行严格同序)
rows, run_paths = [], []
for r in runs:
    _s = r["summary"] or {}
    _m = _s.get("ml_meta") or {}
    rows.append({
        t("eml.metric_run"): r["gse"] or os.path.basename(str(_s.get("expr") or r["name"])),
        t("eml.metric_samples"): _s.get("n_samples", ""),
        t("eml.metric_auc"): _m.get("cv_auc") if _m.get("cv_auc") is not None else "N/A",
        t("eml.metric_shortlist"): _s.get("shortlist_count", ""),
        t("eml.metric_version"): f"v{_s.get('ml_version', '?')}",
        t("mri.col_stale"): t("mri.stale_yes", n=len(r["stale"])) if r["stale"] else "—",
        t("mri.col_time"): time.strftime("%Y-%m-%d %H:%M", time.localtime(r["mtime"])),
        t("mri.col_size"): f"{r['size_bytes'] / 1048576:.1f} MB",
    })
    run_paths.append(r["path"])

df = pd.DataFrame(rows)
event = st.dataframe(df, selection_mode="multi-row", on_select="rerun",
                     height=420, hide_index=True, use_container_width=True)
sel_paths = [run_paths[i] for i in event.selection.rows]

# ---------- 操作栏:勾选 → 确认 → 删除;打包下载(先打包后删的治理路径) ----------
col_del, col_zip, _ = st.columns([1, 1, 2])
with col_del:
    confirm = st.checkbox(t("mri.confirm_cb", n=len(sel_paths)), key="mri_confirm",
                          disabled=not sel_paths)
    if st.button(t("mri.delete_btn"), key="mri_delete", type="primary",
                 disabled=not sel_paths or not confirm):
        deleted = delete_ml_runs(sel_paths)
        st.success(t("mri.deleted_ok", n=len(deleted)))
        st.rerun()
with col_zip:
    if st.button(t("mri.zip_btn"), key="mri_zip", disabled=not sel_paths):
        st.download_button(
            t("mri.zip_btn"),
            data=zip_ml_runs(sel_paths),
            file_name=f"ml_results_{time.strftime('%Y%m%d_%H%M%S')}.zip",
            mime="application/zip", key="mri_zip_dl")
if not sel_paths:
    st.caption(t("mri.none_selected"))
