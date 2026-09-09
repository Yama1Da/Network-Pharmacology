"""Utils 工具入口。

具体功能(靶点获取、PPI、富集分析、分子对接等)都在 core/ 包中实现,
这里统一转出,页面代码仍按 import Utils + Utils.xxx 的方式调用。
界面公共小组件(隐藏侧边栏、返回按钮、下载按钮等)也在这里实现。
"""

import logging
import os
import shutil

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
from core.enrichment import (
    migrate_legacy_cache_names,
    parse_enrichment_filename,
    render_enrichment_plots,
    run_enrichment,
)
from core.http_client import clear_http_cache
from core.files import (
    list_cache_files,
    list_download_dirs,
    list_log_files,
    make_dir_zip,
    make_dirs_zip,
    make_ppi_zip,
)
from core.logging_setup import set_logging_level, setup_logging
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

# 常用按钮配色(与各页 hide_sidebar 的 extra_css 拼接;红/删除按钮因键名各页不同,留在页内)
CSS_BTN_DOWNLOAD = """    .st-key-file_rows [class*="st-key-download_"] button {
        background-color: #207D62;
        color: #C8F0DF;
        border: 1px solid #207D62;
    }"""
CSS_BTN_UPDATE = """    [class*="update_btn"] button,
    .st-key-file_rows [class*="st-key-update_"] button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }"""


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


# ---------- 多选小件(文件表与目录表共用) ----------
# 选中状态按 "{ns}_sel::{条目名}" 存在会话里,全量选中项写入 "{ns}_selected",
# 供 render_file_action_bar 批量操作读取。ns 同一页面内必须一致。


def render_select_all_checkbox(ns: str, names: list[str]) -> bool:
    """表头"全选"复选框;切换时批量改写各行复选框状态(须在行复选框渲染前调用)。"""
    select_all = st.checkbox("", key=f"{ns}_all", help=t("common.select_all"))
    if st.session_state.get(f"{ns}_all_prev") != select_all:
        st.session_state[f"{ns}_all_prev"] = select_all
        for n in names:
            st.session_state[f"{ns}_sel::{n}"] = select_all
    return select_all


def render_row_checkbox(ns: str, name: str, select_all: bool) -> None:
    """行首复选框;value 仅在新键首次创建时生效,新出现的条目跟随当前全选状态。"""
    st.checkbox("", key=f"{ns}_sel::{name}", value=select_all)


def collect_selected(ns: str, names: list[str]) -> list[str]:
    """按全量条目列表汇总选中项(被搜索过滤掉的行不影响已勾选状态)并写入会话。"""
    selected = [n for n in names if st.session_state.get(f"{ns}_sel::{n}")]
    st.session_state[f"{ns}_selected"] = selected
    return selected


def render_file_table(files_df, cache_dir: str, *, search_key: str,
                      display_map: dict[str, str] | None = None,
                      ns: str | None = None) -> list[str] | None:
    """"文件名 | 创建时间 | 下载"三列文件列表:顶部搜索框过滤文件名,表头 + 每行一个下载按钮。

    display_map 提供 {真实文件名: 展示名}:列表展示展示名(友好名称),搜索同时匹配
    展示名与真实文件名,下载按钮仍下载真实文件。

    ns 不为 None 时进入多选模式:行首加复选框、表头加全选;选中项写入
    st.session_state[f"{ns}_selected"] 并作为返回值,供 render_file_action_bar 批量操作。
    """
    display_map = display_map or {}

    def _disp(name: str) -> str:
        return display_map.get(name, name)

    all_names = files_df["File name"].tolist()
    query = st.text_input(t("common.search_placeholder"), key=search_key).strip().lower()
    if query:
        # 展示名/真实文件名任一命中即保留;子串大小写不敏感匹配
        files_df = files_df[files_df["File name"].map(
            lambda n: query in n.lower() or query in _disp(n).lower())]
    if files_df.empty:
        st.info(t("common.search_no_match", query=query))
        return None

    select_all = False
    if ns is None:
        col_h_name, col_h_time, col_h_dl = st.columns([4, 2, 1])
    else:
        col_h_sel, col_h_name, col_h_time, col_h_dl = st.columns([0.5, 3.5, 2, 1])
        with col_h_sel:
            select_all = render_select_all_checkbox(ns, all_names)
    col_h_name.markdown(t("common.col_file_name"))
    col_h_time.markdown(t("common.col_created"))
    col_h_dl.markdown(t("common.col_download"))

    with st.container(key="file_rows"):
        for _, row in files_df.iterrows():
            name = row["File name"]
            if ns is None:
                col_name, col_time, col_dl = st.columns([4, 2, 1])
            else:
                col_sel, col_name, col_time, col_dl = st.columns([0.5, 3.5, 2, 1])
                with col_sel:
                    render_row_checkbox(ns, name, select_all)
            col_name.write(_disp(name))
            col_time.write(row["Created time"])
            with col_dl:
                csv_download_button(os.path.join(cache_dir, name),
                                    name, key=f"download_{name}")

    if ns is None:
        return None
    return collect_selected(ns, all_names)


# ---------- 缓存/结果条目的批量操作(多选删除/打包下载) ----------
# 删除交互与收藏管理页一致:点删除先弹确认框,点"是"才真正删除。


@st.dialog(t("common.delete_dialog_title"))
def _confirm_batch_delete(ns: str, base_dir: str, *, dirs: bool = False,
                          confirm_key: str = "common.delete_dialog_confirm") -> None:
    """批量删除确认弹窗:点"是"才删除并刷新列表,点"否"只关闭窗口。

    dirs=False 删除文件(os.remove);True 删除目录及内容(shutil.rmtree)。
    confirm_key 指定确认文案的 i18n 键,目录类页面可用更贴切的提示。
    """
    names = st.session_state.get(f"{ns}_pending_delete", [])
    st.write(t(confirm_key, n=len(names)))
    if names:
        st.caption(", ".join(names[:10]) + (" …" if len(names) > 10 else ""))
    col_yes, col_no = st.columns(2)
    with col_yes:
        if st.button(t("common.yes"), key=f"{ns}_dlg_yes", width="stretch"):
            failed: list[str] = []
            for n in names:
                path = os.path.join(base_dir, n)
                try:
                    if dirs:
                        shutil.rmtree(path, ignore_errors=True)
                        if os.path.isdir(path):  # rmtree 静默失败时目录仍在,按失败处理
                            raise OSError(path)
                    else:
                        os.remove(path)
                except OSError:
                    logging.exception("Failed to delete %s", path)
                    failed.append(n)
                else:
                    st.session_state.pop(f"{ns}_sel::{n}", None)
            st.session_state.pop(f"{ns}_pending_delete", None)
            st.session_state[f"{ns}_deleted"] = (len(names) - len(failed), failed)
            st.rerun()  # 关闭弹窗并刷新列表
    with col_no:
        if st.button(t("common.no"), key=f"{ns}_dlg_no", width="stretch"):
            st.session_state.pop(f"{ns}_pending_delete", None)
            st.rerun()  # 仅关闭弹窗


def show_file_op_result(ns: str) -> None:
    """显示最近一次批量删除的结果(删除在确认弹窗内完成,结果跨 rerun 传递到页面顶部)。"""
    res = st.session_state.pop(f"{ns}_deleted", None)
    if not res:
        return
    n, failed = res
    if n:
        st.success(t("common.delete_done", n=n))
    if failed:
        st.error(t("common.delete_failed", names=", ".join(failed)))


def render_file_action_bar(ns: str, base_dir: str, *, zip_name: str,
                           dirs: bool = False,
                           confirm_key: str = "common.delete_dialog_confirm") -> None:
    """选中条目的批量操作栏:下载 ZIP / 删除(确认弹窗)。

    选中项来自多选模式写入的 {ns}_selected(render_file_table 或
    render_select_all_checkbox + render_row_checkbox + collect_selected 组合)。
    dirs=False 条目为 base_dir 下的文件;True 为子目录(删除整目录、zip 内保留目录结构)。
    """
    selected = [n for n in st.session_state.get(f"{ns}_selected", [])
                if os.path.exists(os.path.join(base_dir, n))]
    col_cnt, col_dl, col_del = st.columns([2, 2, 2])
    col_cnt.caption(t("common.selected_count", n=len(selected)))
    with col_dl:
        if dirs:
            data = make_dirs_zip({n: os.path.join(base_dir, n) for n in selected}) \
                if selected else b""
        else:
            data = make_ppi_zip({n: os.path.join(base_dir, n) for n in selected}) \
                if selected else b""
        st.download_button(t("common.download_selected"), data=data,
                           file_name=zip_name, mime="application/zip",
                           key=f"{ns}_dl_btn", disabled=not selected, width="stretch")
    with col_del:
        if st.button(t("common.delete_selected"), key=f"{ns}_del_btn",
                     disabled=not selected, width="stretch"):
            st.session_state[f"{ns}_pending_delete"] = selected
            _confirm_batch_delete(ns, base_dir, dirs=dirs, confirm_key=confirm_key)


# ---------- 富集结果文件的友好命名展示 ----------
# 缓存文件名(数据源_标识_短哈希.csv)面向程序;下面两个函数把它翻译成用户可读的名称,
# 富集信息页与分子对接页的文件选择共用,保证两处看到的名字一致。


def humanize_analysis_label(label: str) -> str:
    """"药物_疾病"分析标识转成人读文本:首个下划线换成"→",其余下划线换回空格。"""
    label = (label or "").strip()
    if not label:
        return ""
    if "_" in label:
        drug, disease = label.split("_", 1)
        return f"{drug} → {disease.replace('_', ' ')}"
    return label.replace("_", " ")


def enrichment_display_name(file_name: str) -> str:
    """富集缓存文件名 → 用户友好名称,如 gProfiler_aspirin_cancer_<哈希>.csv
    → "aspirin → cancer · gProfiler";模块缓存的 `_moduleN` 后缀单独译成"模块 N",
    避免被埋进疾病名文本里难以辨认。无法识别的文件名原样返回。"""
    info = parse_enrichment_filename(file_name)
    if info is None:
        return file_name
    label = info["label"]
    module_desc = ""
    if "_module" in label:
        head, _, num = label.rpartition("_module")
        if num.isdigit():
            module_desc = t("enrich_file.module_n", n=int(num))
            label = head
    desc = humanize_analysis_label(label)
    if info["kind"] == "gsea":
        side = ""
        if info["side"] == "drug":
            side = t("enrich_file.side_drug")
        elif info["side"] == "disease":
            side = t("enrich_file.side_disease")
        return " · ".join(x for x in (desc, "GSEA", side) if x)
    return " · ".join(x for x in (desc, info["source"], module_desc) if x)


# ---------- 数据管理页公共逻辑(药物/疾病两个管理页共用) ----------
# 页面通过 spec 字典描述各自的数据源差异:
#   cache_dir        来源文件所在目录
#   merged_dir       merged 整合文件所在目录
#   merged_prefix    整合文件前缀("merged_drug_v1_" 等,带过滤规则版本号,取自 config)
#   source_prefixes  来源文件前缀元组,如 ("ChEMBL_", "DGIdb_")
#   fetch_source     {前缀: fn(name, cache_hours) -> (表, 是否成功)}:重新抓取一个来源文件
#   fetch_merged     fn(name) -> 表:重建整合文件
#   invalidate_merged fn(name):删除对应整合缓存


def _mgmt_file_path(filename: str, spec: dict) -> str:
    """整合前缀的文件位于 merged 目录,其余位于来源缓存目录。"""
    if filename.startswith(spec["merged_prefix"]):
        return os.path.join(spec["merged_dir"], filename)
    return os.path.join(spec["cache_dir"], filename)


def refresh_cache_file(filename: str, spec: dict, cache_hours: int | None = None,
                       clear_http: bool = True) -> dict:
    """按文件名重新抓取或整合数据,更新指定文件(药物/疾病管理页共用)。

    来源文件先备份旧文件,抓取成功才替换,失败则恢复旧文件,避免数据丢失;
    merged 整合文件先更新两个来源文件再重建,任一来源失败则保留旧的整合结果。
    cache_hours 透传给 fetch_source(ChEMBL 更新需区分单文件/全部更新的 TTL)。
    clear_http=False 供并发批处理使用:由 update_all_cache_files 在入口统一清一次,
    避免批量中途清缓存干扰其他线程的在途请求。

    文件名里的名字是 safe_name 有损转换的结果(如 Alzheimer_s_Disease),
    先按映射还原成原始查询名(Alzheimer's Disease)再抓取,否则查不到数据。
    """
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

    # 清空 HTTP 响应缓存,避免强制更新时命中旧响应导致更新无效;
    # 并发批处理改为入口统一清一次,见 clear_http 参数说明
    if clear_http:
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


def update_all_cache_files(spec: dict, cache_hours: int | None = None, *,
                           progress: bool = True, max_workers: int = 4):
    """更新某类别全部来源文件并重建全部整合文件,返回结果表(File/Status/Message)。

    来源文件互相独立,线程池并发抓取提速;抓取链路自带按源限流(throttle 与
    DisGeNET 配额均有锁保护),并发安全。HTTP 响应缓存在入口只清一次:
    每文件中途清理会干扰其他线程的在途请求,且批量开始后本批抓取是唯一的
    缓存写入方,清一次即可保证全部抓到新数据。
    merged 整合依赖全部来源的最新结果,保持串行收尾。
    progress=True 时渲染进度条(仅限在 Streamlit 页面线程调用);max_workers 为并发上限。
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    bar = st.progress(0.0) if progress else None

    def _tick(frac: float, key: str, **fmt) -> None:
        if bar is not None:
            bar.progress(frac, text=t(key, **fmt))

    # 第一步:先并发更新所有来源文件
    jobs: list[str] = []
    names = set()
    for filename in sorted(os.listdir(spec["cache_dir"])):
        if not filename.endswith(".csv"):
            continue
        for prefix in spec["source_prefixes"]:
            if filename.startswith(prefix):
                stem = filename[len(prefix):-len(".csv")]
                names.add(resolve_original_name(stem) or stem)
                jobs.append(filename)
                break

    rows: list[dict] = []
    if jobs:
        clear_http_cache()

        def _run_one(filename: str) -> dict:
            return refresh_cache_file(filename, spec, cache_hours, clear_http=False)

        results: dict[str, dict] = {}
        workers = max(1, min(max_workers, len(jobs)))
        if workers == 1:
            for filename in jobs:
                try:
                    results[filename] = _run_one(filename)
                except Exception:
                    results[filename] = {"Status": "FAILED", "Message": "update error"}
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(_run_one, f): f for f in jobs}
                for done, fut in enumerate(as_completed(futures), 1):
                    filename = futures[fut]
                    try:
                        results[filename] = fut.result()
                    except Exception:
                        results[filename] = {"Status": "FAILED", "Message": "update error"}
                    _tick(0.8 * done / len(jobs), "mgmt.prog.sources",
                          done=done, total=len(jobs))
        # 按文件名排序输出,与并发完成顺序无关
        rows = [{"File": f, **results[f]} for f in sorted(results)]

    # 第二步:再重建所有 merged 整合文件(name 已还原为原始名,文件名需再走 safe_name)
    merged_names = sorted(names)
    for idx, name in enumerate(merged_names, 1):
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
        _tick(0.8 + 0.2 * idx / len(merged_names), "mgmt.prog.merged",
              done=idx, total=len(merged_names))
    if bar is not None:
        bar.progress(1.0, text=t("mgmt.prog.done"))
    return pd.DataFrame(rows, columns=["File", "Status", "Message"])


def render_updateable_file_table(files_df, spec: dict, *, search_key: str) -> str | None:
    """"文件名 | 创建时间 | 更新 | 下载"四列文件列表;返回被点击 Update 的文件名。顶部搜索框过滤文件名。"""
    query = st.text_input(t("common.search_placeholder"), key=search_key).strip().lower()
    if query:
        files_df = files_df[files_df["File name"].str.lower().str.contains(query, regex=False)]
    if files_df.empty:
        st.info(t("common.search_no_match", query=query))
        return None
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


def render_source_management_page(
    ns: str,
    spec: dict,
    *,
    source_filters: list[tuple[str, str]],
    cache_hours_all: int | None = None,
    cache_hours_single: int | None = None,
) -> None:
    """药物/疾病靶点信息管理页的完整渲染流程(两个页面文件只保留各自的
    _fetch/_invalidate/_SPEC 差异,页面调用本函数前不得执行任何 Streamlit 命令)。

    ns     会话键与文案键前缀("drug"/"disease"):会话键写成 "{ns}_update_btn" 等,
           避免两页共用同名键导致更新提示跨页面串台;文案键拼成 "{ns}_mgmt.title" 等。
    source_filters  来源筛选勾选框 [(显示名, widget key)],与 merged 筛选开关同排一行;
           顺序必须与 spec["source_prefixes"] 一一对应。
    cache_hours_all/cache_hours_single  "全部更新"与"单文件更新"的 TTL:
           药物页固定 12/1(ChEMBL 批量分页较慢,短时重复更新只刷时间戳),疾病页维持默认。
    """
    st.set_page_config(page_title=t(f"{ns}_mgmt.page_title"), layout="wide")

    # 隐藏侧边栏,全屏显示页面;行内 Update 按钮紫底、Download 按钮绿底
    hide_sidebar(CSS_BTN_UPDATE + "\n" + CSS_BTN_DOWNLOAD)

    st.title(t(f"{ns}_mgmt.title"))
    st.caption(t(f"{ns}_mgmt.caption"))

    # 展示上一次单文件更新的结果
    last_update = st.session_state.pop(f"{ns}_last_update", None)
    if last_update:
        if last_update["Status"] == "OK":
            st.success(t("common.file_updated", file=last_update["File"], msg=last_update["Message"]))
        else:
            st.error(t("common.file_update_failed", file=last_update["File"], msg=last_update["Message"]))

    col_back, col_update = st.columns([1, 1])
    with col_back:
        back_button()
    with col_update:
        all_update_clicked = st.button(t("common.all_updates"), key=f"{ns}_update_btn")

    if all_update_clicked:
        with st.spinner(t(f"{ns}_mgmt.updating_all")):
            results = update_all_cache_files(spec, cache_hours=cache_hours_all)
        if results.empty:
            st.info(t(f"{ns}_mgmt.no_files_to_update"))
        else:
            st.dataframe(results, width="stretch", hide_index=True)
            if results["Status"].eq("OK").all():
                st.success(t("common.all_complete"))
            else:
                st.warning(t("common.some_failed"))

    st.divider()

    source_files = list_cache_files(ns)
    merged_files = list_cache_files("merged")
    # 按类别前缀(不带版本号)过滤,旧版本号的整合文件也一并列出
    merged_files = merged_files[merged_files["File name"].str.startswith(f"merged_{ns}_")]

    # 数据来源筛选(默认全开):来源筛选作用于来源文件;merged 为两来源的整合缓存,单独可筛
    st.caption(t("common.source_filter"))
    cols = st.columns(len(source_filters) + 1)
    flags: list[bool] = []
    for col, (label, key) in zip(cols, source_filters):
        with col:
            flags.append(st.checkbox(label, value=True, key=key))
    with cols[-1]:
        show_merged = st.checkbox(t("common.filter_merged"), value=True,
                                  key=f"{ns}_flt_merged")

    # 来源前缀与勾选框按传入顺序一一对应
    show_by_prefix = dict(zip(spec["source_prefixes"], flags))

    def _match_source(n: str) -> bool:
        for prefix, show in show_by_prefix.items():
            if n.startswith(prefix):
                return show
        return True

    n_source = len(source_files)
    n_merged = len(merged_files)
    if n_source:
        source_files = source_files[source_files["File name"].map(_match_source)]
    if not show_merged:
        merged_files = merged_files.iloc[0:0]
    files = pd.concat([source_files, merged_files], ignore_index=True)
    clicked_file = None

    if files.empty:
        # 有文件但全被筛掉时,按"筛选无结果"提示,与"目录本来就空"区分开
        st.info(t("common.filter_no_match") if (n_source or n_merged) else t(f"{ns}_mgmt.no_files"))
    else:
        clicked_file = render_updateable_file_table(files, spec, search_key=f"search_{ns}")

    # 单文件更新在页面渲染完成后执行,完成或出错都会记录结果并刷新页面
    if clicked_file:
        with st.spinner(t("common.updating_file", file=clicked_file)):
            try:
                result = refresh_cache_file(clicked_file, spec, cache_hours=cache_hours_single)
            except Exception:
                logging.exception("update failed for '%s'", clicked_file)
                result = {"Status": "FAILED", "Message": "update error"}
        st.session_state[f"{ns}_last_update"] = {"File": clicked_file, **result}
        st.rerun()


# ---------- 设置面板(主页面板齿轮)公共逻辑 ----------

def save_ui_settings(overrides: dict) -> None:
    """把齿轮面板/对接页修改的分析参数持久化到 ui_settings.json。

    config.py 启动时会读取该文件套用覆盖值,因此设置重启后仍生效;
    只写白名单内的键(数值参数与 ChEMBL 黑名单),其余忽略。
    与文件中已有内容按键合并(不同入口各自改各自的参数,互不覆盖)。
    """
    import json

    import config

    clean: dict = {}
    for k, v in overrides.items():
        if k == "CHEMBL_TARGET_BLACKLIST":
            if isinstance(v, (list, tuple)):
                clean[k] = sorted({str(x).strip().upper() for x in v if str(x).strip()})
        elif k in config.UI_SETTINGS_KEYS and isinstance(v, bool):
            clean[k] = v
        elif k in config.UI_SETTINGS_KEYS and isinstance(v, int):
            clean[k] = v  # 步进器参数按整数保存,避免 +- 调节后出现小数
        elif k in config.UI_SETTINGS_KEYS and isinstance(v, (int, float)) \
                and not isinstance(v, bool):
            clean[k] = float(v)
    if not clean:
        return
    merged: dict = {}
    if os.path.exists(config.UI_SETTINGS_FILE):
        try:
            with open(config.UI_SETTINGS_FILE, "r", encoding="utf-8") as f:
                merged = json.load(f) or {}
        except Exception:
            merged = {}  # 旧文件损坏时按空文件处理,以本次写入为准
    merged.update(clean)
    os.makedirs(os.path.dirname(config.UI_SETTINGS_FILE), exist_ok=True)
    with open(config.UI_SETTINGS_FILE, "w", encoding="utf-8") as f:
        json.dump(merged, f, ensure_ascii=False, indent=2)


def purge_merged_cache() -> int:
    """清空 merged 整合缓存目录,返回删除的文件数。

    降噪阈值等参数变更后,旧参数算出的整合结果不再可信,统一清除,
    下次分析按新参数从原始缓存重建(原始缓存不受影响,无需重新联网)。
    """
    import glob

    import config

    n = 0
    for path in glob.glob(os.path.join(config.MERGED_CACHE_DIR, "*.csv")):
        try:
            os.remove(path)
            n += 1
        except OSError:
            pass
    if n:
        logging.info("Purged %d merged cache files (settings changed)", n)
    return n
