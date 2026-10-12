"""药物/疾病靶点信息管理页的整页渲染流程。

页面文件只保留各自的 _fetch/_invalidate/_SPEC 差异,渲染统一走
render_source_management_page;数据刷新编排本体在 core.management,
本模块负责 Streamlit 呈现(进度条、结果提示、两段式按钮状态机)。"""

import logging

import pandas as pd
import streamlit as st

from i18n import t
from core.files import list_cache_files
from core.management import (_mgmt_file_path, refresh_cache_file,
                             update_all_cache_files)
from ui.widgets import (CSS_BTN_DOWNLOAD, CSS_BTN_UPDATE, CSS_HIDE_IDLE_PROGRESS,
                        back_button, busy_pending, enter_busy, finish_busy,
                        hide_sidebar, is_busy, render_file_action_bar,
                        render_file_table, show_busy_notice, show_file_op_result)

# refresh_cache_file 结果 Message 里的"未配置 API 密钥"标记:
# fetch_source 返回 "failed:no_api_key",core 原样并入 Message,见 core.management
_NO_KEY_MARK = "no_api_key"


def _has_no_key_mark(messages: pd.Series) -> bool:
    return messages.fillna("").astype(str).str.contains(_NO_KEY_MARK).any()


def show_update_summary(results: pd.DataFrame) -> None:
    """更新结果表的状态汇总:全 OK 成功;有 FAILED 警告"部分失败",其中因未配置
    API 密钥导致的失败给出指明原因的专用提示;仅 0 靶点警告检查查询物;
    结果为空(全新安装、缓存目录尚未生成)提示先去首页完成一次分析。"""
    if results.empty:
        st.info(t("common.update_no_files"))
        return
    if results["Status"].eq("OK").all():
        st.success(t("common.all_complete"))
    elif results["Status"].eq("FAILED").any():
        if _has_no_key_mark(results["Message"]):
            st.warning(t("common.some_failed_no_key"))
        else:
            st.warning(t("common.some_failed"))
    else:
        st.warning(t("common.some_empty"))


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

    # 批量任务运行期间 busy=True:所有按钮禁用,防止点击排队到任务结束后再次触发
    # (两段式 busy 状态机的完整说明见 ui.widgets)
    busy = is_busy(f"{ns}_busy")

    # 隐藏侧边栏,全屏显示页面;行内 Update 按钮紫底、Download 按钮绿底、批量删除按钮红底
    hide_sidebar(CSS_BTN_UPDATE + "\n" + CSS_BTN_DOWNLOAD + "\n" + CSS_HIDE_IDLE_PROGRESS + f"""
    .st-key-{ns}_files_del_btn button {{
        background-color: #F53F3F;
        color: #FFFFFF;
        border: 1px solid #F53F3F;
    }}""")

    st.title(t(f"{ns}_mgmt.title"))
    st.caption(t(f"{ns}_mgmt.caption"))

    # 展示上一次单文件更新的结果
    last_update = st.session_state.pop(f"{ns}_last_update", None)
    if last_update:
        if last_update["Status"] == "OK":
            st.success(t("common.file_updated", file=last_update["File"], msg=last_update["Message"]))
        elif last_update["Status"] == "EMPTY":
            # 抓取成功但 0 靶点:警告而非错误,提示检查查询物
            st.warning(t("common.file_update_empty", file=last_update["File"]))
        elif _NO_KEY_MARK in str(last_update.get("Message", "")):
            # 未配置 API 密钥的失败:指明原因的专用提示,避免误判为网络故障
            st.warning(t("common.file_update_no_key", file=last_update["File"]))
        else:
            st.error(t("common.file_update_failed", file=last_update["File"], msg=last_update["Message"]))

    # 批量删除结果(删除在确认弹窗内完成,跨 rerun 传回页面顶部)
    show_file_op_result(f"{ns}_files")

    # 上一次批量任务被在途点击打断时给出提示(见 ui.widgets 两段式 busy 状态机)
    show_busy_notice(f"{ns}_busy")

    col_back, col_update = st.columns([1, 1])
    with col_back:
        back_button(disabled=busy)
    with col_update:
        all_update_clicked = st.button(t("common.all_updates"), key=f"{ns}_update_btn",
                                       disabled=busy)

    if all_update_clicked:
        # 两段式:点击轮只置位并立刻重渲染(本轮按钮的 disabled 与上一轮相同),
        # 下一轮 busy=True,全部按钮以禁用态绘制后再执行任务。
        enter_busy(f"{ns}_busy", pending="all")

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

    # 批量更新的进度条:render_file_table 经插槽把进度条画到搜索框正下方。
    # 进度条元素必须每轮常驻、只做原地更新(改值/文案):若该位置只在更新轮才插入
    # 元素(或往空容器里填内容),Streamlit 前端对页面中部的增删帧会错位协调,
    # 把上一帧的列表行残留成可交互的"残影"。空闲轮文案为空串,由页面 CSS
    # (CSS_HIDE_IDLE_PROGRESS)把无文案的进度条隐藏。文件列表为空(没走插槽)时
    # 执行块兜底就地创建,故容器须在分支外先初始化。
    bars: list = []
    if files.empty:
        # 有文件但全被筛掉时,按"筛选无结果"提示,与"目录本来就空"区分开
        st.info(t("common.filter_no_match") if (n_source or n_merged) else t(f"{ns}_mgmt.no_files"))
    else:
        # 多选模式:行首复选框 + 表头全选,path_of 让来源文件与 merged 整合文件各按其目录解析
        def _progress_slot() -> None:
            bars.append(st.progress(0.0, text=""))

        render_file_table(files, spec["cache_dir"], search_key=f"search_{ns}",
                          ns=f"{ns}_files", path_of=lambda n: _mgmt_file_path(n, spec),
                          disabled=busy, pre_rows=_progress_slot)
    clicked_file = st.session_state.pop(f"{ns}_files_update_clicked", None)

    # 选中条目的批量操作:更新选中 / 下载选中(ZIP) / 删除选中(确认弹窗);
    # "更新选中"复用 update_all_cache_files 的并发抓取(filenames 限定选中清单)
    batch_update = render_file_action_bar(
        f"{ns}_files", None, zip_name=f"{ns}_files_selected.zip",
        path_of=lambda n: _mgmt_file_path(n, spec),
        update_label=t("common.update_selected"), disabled=busy)
    if batch_update:
        # 两段式点击轮:置位 busy 与待办清单后立刻重渲染,禁用轮再执行(见 ui.widgets)
        enter_busy(f"{ns}_busy", pending=batch_update)

    # 批量执行轮(busy=True):本轮所有按钮已按禁用态渲染,跑完复位并重渲染展示结果。
    # 待办为 "all"(全部更新)或选中文件名清单(更新选中);
    # 选中更新复用 update_all_cache_files 的并发抓取(filenames 限定选中清单)。
    pending_update = busy_pending(f"{ns}_busy")
    if busy and pending_update is not None:
        is_all = pending_update == "all"
        completed = False
        try:
            # 进度条本体已由插槽每轮常驻在搜索框下方,这里只做原地更新;
            # 文件列表为空没走插槽时,兜底就地创建(此时必然在页面底部)
            bar = bars[0] if bars else st.progress(0.0)
            results = update_all_cache_files(
                spec,
                cache_hours=cache_hours_all if is_all else cache_hours_single,
                filenames=None if is_all else pending_update,
                on_progress=lambda frac, text: bar.progress(frac, text=text))
            completed = True
        finally:
            finish_busy(f"{ns}_busy", completed=completed)
        st.session_state[f"{ns}_batch_results"] = results
        st.rerun()  # 复位 busy 后重渲染:按钮恢复可用,结果由下方 pop 展示
    batch_results = st.session_state.pop(f"{ns}_batch_results", None)
    if batch_results is not None:
        if batch_results.empty:
            st.info(t(f"{ns}_mgmt.no_files_to_update"))
        else:
            st.dataframe(batch_results, width="stretch", hide_index=True)
            show_update_summary(batch_results)

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
