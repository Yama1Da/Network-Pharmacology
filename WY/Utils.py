"""Utils 工具入口。

具体功能(靶点获取、PPI、富集分析、分子对接等)都在 core/ 包中实现,
这里统一转出,页面代码仍按 import Utils + Utils.xxx 的方式调用。
界面公共小组件(隐藏侧边栏、返回按钮、下载按钮等)也在这里实现。
"""

import os

import pandas as pd
import streamlit as st

from i18n import t
from core.chembl import fetch_drug_targets_chembl, set_chembl_cache_ttl
from core.dgidb import fetch_drug_targets_dgidb
from core.disgenet import fetch_disease_targets_disgenet
from core.docking import (
    clean_docking_cache,
    fetch_ligand_sdf_3d,
    fetch_receptor_candidates_batch,
    fetch_receptor_pdb_candidates,
    infer_box_from_cofactor,
    parse_vina_scores,
    prepare_ligand_pdbqt,
    prepare_receptor_pdbqt,
    prune_docking_results,
    run_vina,
    split_vina,
    star_docking_result,
    write_docking_metadata,
)
from core.enrichment import render_enrichment_plots, run_enrichment
from core.http_client import clear_http_cache
from core.files import (
    list_cache_files,
    list_download_dirs,
    list_log_files,
    make_dir_zip,
    make_ppi_zip,
)
from core.logging_setup import setup_logging
from core.opentargets import get_disease_targets
from core.ppi import (
    build_ppi_graph,
    compute_ppi_metrics,
    fetch_string_ppi,
    mark_hub_genes,
    render_pyvis_html,
    render_static_plot,
    save_ppi_results,
)
from core.targets import (
    get_common_targets,
    get_disease_targets_merged,
    get_drug_targets_merged,
    get_targets_concurrent,
)
from core.text_utils import resolve_original_name, safe_name


# ---------- 界面公共组件 ----------


def hide_sidebar(extra_css: str = "") -> None:
    """隐藏 Streamlit 侧边栏,全屏显示页面;extra_css 为该页面额外的样式规则。"""
    st.markdown(f"""
<style>
    [data-testid="stSidebar"] {{
        display: none;
    }}
{extra_css}</style>
""", unsafe_allow_html=True)


def back_button(label_key: str = "common.back_main", *, key: str = "back_btn",
                target: str = "App.py") -> None:
    """"返回"按钮:点击后切到 target 页面;页面自行把它放进列布局。"""
    if st.button(t(label_key), key=key):
        st.switch_page(target)


def csv_download_button(file_path: str, file_name: str, *, key: str) -> None:
    """行内 CSV 下载按钮:文件存在时读取字节并渲染按钮,不存在时静默跳过。"""
    if not os.path.exists(file_path):
        return
    with open(file_path, "rb") as f:
        st.download_button(t("common.download"), data=f.read(), file_name=file_name,
                           mime="text/csv", key=key)


def download_button(label_key: str, data, file_name: str, mime: str) -> None:
    """通宽下载按钮:label 经 i18n 翻译,data 为内存字节/DataFrame 内容等。"""
    st.download_button(t(label_key), data, file_name, mime, width="stretch")


def download_file_button(label_key: str, path: str, file_name: str, mime: str) -> None:
    """通宽下载按钮(磁盘文件版):文件存在时才渲染。"""
    if not os.path.exists(path):
        return
    with open(path, "rb") as f:
        download_button(label_key, f.read(), file_name, mime)


def df_download_button(df, file_name: str, *, key: str) -> None:
    """行内 CSV 下载按钮:直接下载内存中的 DataFrame(utf-8-sig 编码,Excel 可直接打开)。"""
    st.download_button(t("common.download"), data=df.to_csv(index=False).encode("utf-8-sig"),
                       file_name=file_name, mime="text/csv", key=key)


def render_file_table(files_df, cache_dir: str) -> None:
    """"文件名 | 创建时间 | 下载"三列文件列表:表头 + 每行一个下载按钮。"""
    col_h_name, col_h_time, col_h_dl = st.columns([4, 2, 1])
    col_h_name.markdown(t("common.col_file_name"))
    col_h_time.markdown(t("common.col_created"))
    col_h_dl.markdown(t("common.col_download"))

    with st.container(key="file_rows"):
        for _, row in files_df.iterrows():
            col_name, col_time, col_dl = st.columns([4, 2, 1])
            col_name.write(row["File name"])
            col_time.write(row["Created time"])
            with col_dl:
                csv_download_button(os.path.join(cache_dir, row["File name"]),
                                    row["File name"], key=f"download_{row['File name']}")


# ---------- 数据管理页公共逻辑(药物/疾病两个管理页共用) ----------
# 页面通过 spec 字典描述各自的数据源差异:
#   cache_dir        来源文件所在目录
#   merged_dir       merged 整合文件所在目录
#   merged_prefix    整合文件前缀("merged_drug_" / "merged_disease_")
#   source_prefixes  来源文件前缀元组,如 ("ChEMBL_", "DGIdb_")
#   fetch_source     {前缀: fn(name, cache_hours) -> (表, 是否成功)}:重新抓取一个来源文件
#   fetch_merged     fn(name) -> 表:重建整合文件
#   invalidate_merged fn(name):删除对应整合缓存


def _mgmt_file_path(filename: str, spec: dict) -> str:
    """整合前缀的文件位于 merged 目录,其余位于来源缓存目录。"""
    if filename.startswith(spec["merged_prefix"]):
        return os.path.join(spec["merged_dir"], filename)
    return os.path.join(spec["cache_dir"], filename)


def refresh_cache_file(filename: str, spec: dict, cache_hours: int | None = None) -> dict:
    """按文件名重新抓取或整合数据,更新指定文件(药物/疾病管理页共用)。

    来源文件先备份旧文件,抓取成功才替换,失败则恢复旧文件,避免数据丢失;
    merged 整合文件先更新两个来源文件再重建,任一来源失败则保留旧的整合结果。
    cache_hours 透传给 fetch_source(ChEMBL 更新需区分单文件/全部更新的 TTL)。

    文件名里的名字是 safe_name 有损转换的结果(如 Alzheimer_s_Disease),
    先按映射还原成原始查询名(Alzheimer's Disease)再抓取,否则查不到数据。
    """
    import logging

    # merged 整合文件:先更新两个来源文件,再重建
    if filename.startswith(spec["merged_prefix"]):
        stem = filename[len(spec["merged_prefix"]):-len(".csv")]
        name = resolve_original_name(stem) or stem
        for prefix in spec["source_prefixes"]:
            # 来源文件名用 safe 名(name 已还原为原始名,含空格/撇号等)
            source = f"{prefix}{safe_name(name)}.csv"
            result = refresh_cache_file(source, spec, cache_hours)
            if result["Status"] != "OK":
                return {"Status": "FAILED",
                        "Message": f"source '{source}' update failed: {result['Message']}"}
        merged_path = _mgmt_file_path(filename, spec)
        if os.path.exists(merged_path):
            os.remove(merged_path)  # 可能已在上一步删除,这里兜底再删一次
        df = spec["fetch_merged"](name)
        if not df.empty:
            return {"Status": "OK", "Message": f"{len(df)} targets"}
        return {"Status": "FAILED", "Message": "fetch failed"}

    # 来源文件:先备份,抓取成功才替换,失败则恢复旧文件
    path = _mgmt_file_path(filename, spec)
    backup = path + ".bak"
    if os.path.exists(path):
        os.replace(path, backup)

    # 清空 HTTP 响应缓存,避免强制更新时命中旧响应导致更新无效
    clear_http_cache()

    ok = False
    df = pd.DataFrame()
    try:
        for prefix in spec["source_prefixes"]:
            if filename.startswith(prefix):
                stem = filename[len(prefix):-len(".csv")]
                name = resolve_original_name(stem) or stem
                df, ok = spec["fetch_source"][prefix](name, cache_hours)
                break
        else:
            raise ValueError(f"unknown file prefix: {filename}")
    except Exception:
        logging.exception("refresh_file failed for '%s'", filename)
        ok = False

    if ok:
        if os.path.exists(backup):
            os.remove(backup)
        spec["invalidate_merged"](name)
        return {"Status": "OK", "Message": f"{len(df)} targets"}
    # 失败:恢复旧文件
    if os.path.exists(backup):
        os.replace(backup, path)
    return {"Status": "FAILED", "Message": "fetch failed"}


def update_all_cache_files(spec: dict, cache_hours: int | None = None):
    """更新某类别全部来源文件并重建全部整合文件,返回结果表(File/Status/Message)。"""
    import logging

    rows: list[dict] = []
    # 第一步:先更新所有来源文件
    names = set()
    for filename in sorted(os.listdir(spec["cache_dir"])):
        if not filename.endswith(".csv"):
            continue
        for prefix in spec["source_prefixes"]:
            if filename.startswith(prefix):
                stem = filename[len(prefix):-len(".csv")]
                names.add(resolve_original_name(stem) or stem)
                try:
                    rows.append({"File": filename,
                                 **refresh_cache_file(filename, spec, cache_hours)})
                except Exception:
                    rows.append({"File": filename, "Status": "FAILED", "Message": "update error"})
                break
    # 第二步:再重建所有 merged 整合文件(name 已还原为原始名,文件名需再走 safe_name)
    for name in sorted(names):
        merged_file = f"{spec['merged_prefix']}{safe_name(name)}.csv"
        merged_path = os.path.join(spec["merged_dir"], merged_file)
        if os.path.exists(merged_path):
            os.remove(merged_path)
        try:
            df = spec["fetch_merged"](name)
            if not df.empty:
                rows.append({"File": merged_file, "Status": "OK", "Message": f"{len(df)} targets"})
            else:
                rows.append({"File": merged_file, "Status": "FAILED", "Message": "fetch failed"})
        except Exception:
            rows.append({"File": merged_file, "Status": "FAILED", "Message": "update error"})
    return pd.DataFrame(rows, columns=["File", "Status", "Message"])


def render_updateable_file_table(files_df, spec: dict) -> str | None:
    """"文件名 | 创建时间 | 更新 | 下载"四列文件列表;返回被点击 Update 的文件名。"""
    col_h_name, col_h_time, col_h_btn, col_h_dl = st.columns([3, 2, 1, 1])
    col_h_name.markdown(t("common.col_file_name"))
    col_h_time.markdown(t("common.col_created"))
    col_h_btn.markdown(t("common.col_update"))
    col_h_dl.markdown(t("common.col_download"))

    # 每行:Update 只更新该文件,Download 下载该 CSV
    clicked = None
    with st.container(key="file_rows"):
        for _, row in files_df.iterrows():
            col_name, col_time, col_btn, col_dl = st.columns([3, 2, 1, 1])
            col_name.write(row["File name"])
            col_time.write(row["Created time"])
            with col_btn:
                if st.button(t("common.update"), key=f"update_{row['File name']}"):
                    clicked = row["File name"]
            with col_dl:
                csv_download_button(_mgmt_file_path(row["File name"], spec),
                                    row["File name"], key=f"download_{row['File name']}")
    return clicked
