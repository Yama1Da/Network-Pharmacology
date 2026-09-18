"""跨页面复用的界面小组件:布局辅助、下载按钮、多选小件、搜索框与历史、
文件表、列表分页、批量操作栏(删除确认/打包下载/更新选中)。只做展示与交互,
数据抓取与刷新编排在 core/(见 core.management)。"""

import logging
import os
import shutil
import threading
from collections.abc import Callable

import streamlit as st

from i18n import t
from core.files import make_dirs_zip, make_ppi_zip
from core.text_utils import load_name_map

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

# 空闲进度条隐藏:批量更新/网络构建的进度条元素每轮常驻、只做原地更新
# (元素结构跨帧稳定,否则 Streamlit 前端对页面中部的增删帧会错位残留"残影",
# 详见 ui.source_management)。进度条空闲时恒为 0%,该规则把 0% 的进度条整体
# 隐藏;更新开始后数值大于 0 自然恢复显示,完成后重渲染回到 0% 重新隐藏。
# aria-label 是 react-aria 固定值"progress",与文案无关,不能用作判据。
CSS_HIDE_IDLE_PROGRESS = """    .stApp div.stProgress:has([role="progressbar"][aria-valuenow="0"]) {
        display: none;
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
                target: str = "App.py", disabled: bool = False) -> None:
    """"返回"按钮:点击后切到 target 页面;页面自行把它放进列布局。
    disabled 供批量任务运行时禁用,避免点击排队到任务结束后触发跳转、丢掉结果展示。"""
    if st.button(t(label_key), key=key, disabled=disabled):
        st.switch_page(target)


def csv_download_button(file_path: str, file_name: str, *, key: str,
                        disabled: bool = False) -> None:
    """行内 CSV 下载按钮:文件存在时读取字节并渲染按钮,不存在时静默跳过。"""
    if not os.path.exists(file_path):
        return
    with open(file_path, "rb") as f:
        st.download_button(t("common.download"), data=f.read(), file_name=file_name,
                           mime="text/csv", key=key, disabled=disabled)


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
    select_all = st.checkbox(t("common.select_all"), key=f"{ns}_all",
                             label_visibility="collapsed")
    if st.session_state.get(f"{ns}_all_prev") != select_all:
        st.session_state[f"{ns}_all_prev"] = select_all
        for n in names:
            st.session_state[f"{ns}_sel::{n}"] = select_all
    return select_all


def render_row_checkbox(ns: str, name: str, select_all: bool,
                        label: str | None = None) -> None:
    """行首复选框;value 仅在新键首次创建时生效,新出现的条目跟随当前全选状态。

    label 折叠隐藏、仅供读屏(建议传展示名),缺省用真实名。
    """
    st.checkbox(label or name, key=f"{ns}_sel::{name}", value=select_all,
                label_visibility="collapsed")


def collect_selected(ns: str, names: list[str]) -> list[str]:
    """按全量条目列表汇总选中项(被搜索过滤掉的行不影响已勾选状态)并写入会话。"""
    selected = [n for n in names if st.session_state.get(f"{ns}_sel::{n}")]
    st.session_state[f"{ns}_selected"] = selected
    return selected


# ---------- 搜索框 + 搜索历史 ----------

# search_history.json 的进程内读写锁:读-改-写是事务,多标签页/多会话并发时互斥,
# 否则并发写入互相整文件覆盖、丢历史
_search_history_lock = threading.Lock()


def load_search_history(search_key: str) -> list[str]:
    """读取指定搜索框的历史记录(最新的在前,最多 SEARCH_HISTORY_MAX 条)。

    全部历史按搜索框 key 分组存于 Static/data/cache/search_history.json,
    重启 Streamlit 后仍在;文件缺失或损坏时按空历史处理。
    """
    import json

    import config

    if not os.path.exists(config.SEARCH_HISTORY_FILE):
        return []
    try:
        with open(config.SEARCH_HISTORY_FILE, "r", encoding="utf-8") as f:
            data = json.load(f) or {}
    except Exception:
        return []  # 文件损坏时按空历史处理,下次写入会重建
    hist = data.get(search_key, [])
    return [h for h in hist if isinstance(h, str) and h][:config.SEARCH_HISTORY_MAX]


def add_search_history(search_key: str, query: str) -> None:
    """记录一次搜索:去重后提到最前,超出 SEARCH_HISTORY_MAX 的旧记录淘汰。

    与最近一次记录相同时不重写文件,避免页面内其他操作触发的 rerun 反复落盘。
    """
    import json

    import config

    query = (query or "").strip().lower()
    if not query:
        return
    with _search_history_lock:  # 读-改-写整段互斥,并发时后写者不能覆盖先写者的条目
        hist = load_search_history(search_key)
        if hist and hist[0] == query:
            return
        data: dict = {}
        if os.path.exists(config.SEARCH_HISTORY_FILE):
            try:
                with open(config.SEARCH_HISTORY_FILE, "r", encoding="utf-8") as f:
                    data = json.load(f) or {}
            except Exception:
                data = {}  # 旧文件损坏时按空文件处理,以本次写入为准
        data[search_key] = [query] + [h for h in hist if h != query]
        os.makedirs(os.path.dirname(config.SEARCH_HISTORY_FILE), exist_ok=True)
        # 原子写(临时文件+替换),中断不会留下半截 JSON
        from core.caching import save_json_file
        save_json_file(config.SEARCH_HISTORY_FILE, data)


def search_input(search_key: str) -> str:
    """带历史记录的搜索框:每个搜索框按自己的 key 维护独立历史。

    输入框下方以 pills 展示历史(最多 6 条),点击历史项回填搜索框并重渲染;
    每次提交的非空搜索写入 JSON 历史,重启 Streamlit 后仍有效。返回小写查询串。
    """
    hist = load_search_history(search_key)
    # 点击历史项后先回填搜索框的会话值,必须在 widget 创建前设置才会生效
    picked = st.session_state.pop(f"{search_key}_hist_pick", None)
    if picked:
        st.session_state[search_key] = picked
    # pills 的选中态会跨 rerun 保留:点击历史项后残留的选中值与新查询不一致,
    # 会被下面误判成又一次点击,把输入框弹回历史词、丢掉用户刚输入的查询。
    # 故在 pills 实例化前清掉选中(创建后改控件键会抛异常);本次点击的浏览器
    # 事件在控件创建时写入,不受影响。
    st.session_state.pop(f"{search_key}_hist_pills", None)
    query = st.text_input(t("common.search_placeholder"), key=search_key).strip().lower()
    if query:
        # 先记录本次搜索再渲染历史,保证提交的这一轮就能看到自己的历史项
        add_search_history(search_key, query)
        hist = load_search_history(search_key)
    if hist:
        choice = st.pills(t("common.search_history"), hist,
                          key=f"{search_key}_hist_pills", label_visibility="collapsed")
        if choice and choice != (st.session_state.get(search_key) or ""):
            st.session_state[f"{search_key}_hist_pick"] = choice
            st.rerun()
    return query


def loose_name(text: str) -> str:
    """搜索用归一化:转小写并去掉空白、下划线与其他标点。

    让 "Type 2 diabetes mellitus" 能匹配 Type_2_diabetes_mellitus、
    "Alzheimer's" 能匹配 Alzheimer_s 这类安全文件名与原始输入的差异。
    """
    import re

    return re.sub(r"[\W_]+", "", str(text), flags=re.UNICODE).lower()


def search_hit(name: str, query: str, name_map: dict[str, str] | None = None) -> bool:
    """判断文件/目录名是否命中搜索词(归一化子串匹配,分隔符与标点差异不敏感)。

    除直接匹配外,还借 name_map.json 把文件名内嵌的安全名(如 Type_2_diabetes_mellitus)
    还原成原始名(Type 2 diabetes mellitus)后再匹配一次,覆盖原始名含被
    safe_name 替换掉的字符(如 "My+Drug@2024" ↔ My_Drug_2024)的场景。
    name_map 传入时直接使用(避免文件列表逐个重读),缺省现读一次。
    """
    q = loose_name(query)
    if not q:
        return True
    if q in loose_name(name):
        return True
    if name_map is None:
        from core.text_utils import load_name_map

        name_map = load_name_map()
    expanded = str(name)
    for key, original in name_map.items():
        if key and key in expanded:
            expanded = expanded.replace(key, str(original))
    return q in loose_name(expanded)


# ---------- 列表分页 ----------
# 文件/目录列表统一按 PAGE_SIZE 行分页;页码与最近一次搜索词按 search_key 存会话,
# 与搜索历史等键一样跨 rerun 保留。翻页条画在行列表之后,点击按钮改写页码后整页重绘。
PAGE_SIZE = 10


def page_slice(search_key: str, n_total: int, *, query: str = "") -> tuple[int, int, int, int]:
    """计算过滤后列表当前页的 [start, end) 切片,返回 (start, end, 页码, 总页数)。

    只做状态与边界处理,不渲染:搜索词变化时页码归 1(用户新搜索应从头浏览);
    删除/筛选导致当前页越界时收拢到最后一页,避免停在空页。
    """
    page_key = f"{search_key}_page"
    last_q_key = f"{search_key}_last_query"
    if st.session_state.get(last_q_key) != query:
        st.session_state[last_q_key] = query
        st.session_state[page_key] = 0
    n_pages = max(1, -(-n_total // PAGE_SIZE))  # n_total=0 时也给出 1 页,调用方负责空态
    page = min(max(int(st.session_state.get(page_key, 0)), 0), n_pages - 1)
    st.session_state[page_key] = page
    return page * PAGE_SIZE, (page + 1) * PAGE_SIZE, page, n_pages


def render_pager(search_key: str, page: int, n_pages: int, *, disabled: bool = False) -> None:
    """行列表下方的翻页条:◀ 上一页 | 第 x / n 页 | 下一页 ▶;不足两页时不渲染。

    按钮只改写会话页码并 st.rerun,由下一轮 page_slice 按新页码切片;
    disabled 供批量任务运行轮禁用,与行内按钮一致,防止点击排队到任务后触发。
    """
    if n_pages <= 1:
        return
    col_prev, col_status, col_next = st.columns([1, 2, 1])
    with col_prev:
        if st.button(t("common.page_prev"), key=f"{search_key}_page_prev",
                     disabled=disabled or page <= 0, width="stretch"):
            st.session_state[f"{search_key}_page"] = page - 1
            st.rerun()
    col_status.caption(t("common.page_status", cur=page + 1, total=n_pages))
    with col_next:
        if st.button(t("common.page_next"), key=f"{search_key}_page_next",
                     disabled=disabled or page >= n_pages - 1, width="stretch"):
            st.session_state[f"{search_key}_page"] = page + 1
            st.rerun()


def render_file_table(files_df, cache_dir: str, *, search_key: str,
                      display_map: dict[str, str] | None = None,
                      ns: str | None = None,
                      path_of: Callable[[str], str] | None = None,
                      disabled: bool = False,
                      pre_rows: Callable[[], None] | None = None) -> list[str] | None:
    """文件名 | 最后修改时间 | 下载"三列文件列表:顶部搜索框过滤文件名,表头 + 每行一个下载按钮。

    display_map 提供 {真实文件名: 展示名}:列表展示展示名(友好名称),搜索同时匹配
    展示名与真实文件名,下载按钮仍下载真实文件。

    ns 不为 None 时进入多选模式:行首加复选框、表头加全选;选中项写入
    st.session_state[f"{ns}_selected"] 并作为返回值,供 render_file_action_bar 批量操作。
    "全选"作用于全量条目(含被搜索滤掉、不在当前页的行),与分页无关,由批量操作栏
    的已选计数向用户明示。

    path_of 提供时(ns 多选模式下)每行多一个 Update 按钮(管理页用),行下载/更新路径
    都按它解析(来源文件与 merged 整合文件分属不同目录);单文件更新点击记录在
    st.session_state[f"{ns}_update_clicked"],由页面取出处理。

    disabled 供批量任务运行时禁用行内按钮与翻页条,避免点击排队到任务结束后再次触发。

    pre_rows 可选插槽:在搜索框之后、表头与行之前渲染(如批量更新的进度条,使其
    出现在搜索框正下方而非整页底部)。插槽内创建的 widget 由调用方持有引用并后续更新。
    """
    display_map = display_map or {}
    resolve = path_of or (lambda n: os.path.join(cache_dir, n))

    def _disp(name: str) -> str:
        return display_map.get(name, name)

    all_names = files_df["File name"].tolist()
    query = search_input(search_key)
    if query:
        # 展示名/真实文件名/还原后的原始名任一命中即保留;分隔符与标点差异不敏感
        name_map = load_name_map()
        files_df = files_df[files_df["File name"].map(
            lambda n: search_hit(n, query, name_map) or search_hit(_disp(n), query, name_map))]
    if files_df.empty:
        st.info(t("common.search_no_match", query=query))
        return None

    if pre_rows is not None:
        pre_rows()
    start, end, page, n_pages = page_slice(search_key, len(files_df), query=query)

    with_update = bool(ns and path_of)
    select_all = False
    if ns is None:
        col_h_name, col_h_time, col_h_dl = st.columns([4, 2, 1])
    elif with_update:
        col_h_sel, col_h_name, col_h_time, col_h_upd, col_h_dl = st.columns([0.5, 3, 2, 1, 1])
        with col_h_sel:
            select_all = render_select_all_checkbox(ns, all_names)
        col_h_upd.markdown(t("common.col_update"))
    else:
        col_h_sel, col_h_name, col_h_time, col_h_dl = st.columns([0.5, 3.5, 2, 1])
        with col_h_sel:
            select_all = render_select_all_checkbox(ns, all_names)
    col_h_name.markdown(t("common.col_file_name"))
    col_h_time.markdown(t("common.col_modified"))
    col_h_dl.markdown(t("common.col_download"))

    with st.container(key="file_rows"):
        for _, row in files_df.iloc[start:end].iterrows():
            name = row["File name"]
            if ns is None:
                col_name, col_time, col_dl = st.columns([4, 2, 1])
            elif with_update:
                col_sel, col_name, col_time, col_upd, col_dl = st.columns([0.5, 3, 2, 1, 1])
                with col_sel:
                    render_row_checkbox(ns, name, select_all, label=_disp(name))
                with col_upd:
                    if st.button(t("common.update"), key=f"update_{name}",
                                 disabled=disabled):
                        st.session_state[f"{ns}_update_clicked"] = name
            else:
                col_sel, col_name, col_time, col_dl = st.columns([0.5, 3.5, 2, 1])
                with col_sel:
                    render_row_checkbox(ns, name, select_all, label=_disp(name))
            col_name.write(_disp(name))
            col_time.write(row["Modified time"])
            with col_dl:
                csv_download_button(resolve(name), name, key=f"download_{name}",
                                    disabled=disabled)

    render_pager(search_key, page, n_pages, disabled=disabled)

    if ns is None:
        return None
    return collect_selected(ns, all_names)


# ---------- 缓存/结果条目的批量操作(多选删除/打包下载) ----------
# 删除交互与收藏管理页一致:点删除先弹确认框,点"是"才真正删除。


def _confirm_batch_delete(ns: str, base_dir: str | None, *, dirs: bool = False,
                          confirm_key: str = "common.delete_dialog_confirm",
                          path_of=None) -> None:
    """批量删除确认弹窗:点"是"才删除并刷新列表,点"否"只关闭窗口。

    dirs=False 删除文件(os.remove);True 删除目录及内容(shutil.rmtree)。
    confirm_key 指定确认文案的 i18n 键,目录类页面可用更贴切的提示。
    path_of 提供时按它把条目名解析成完整路径(条目分属多个目录的页面用,如管理页的
    来源文件与 merged 整合文件);缺省用 base_dir 拼接。
    st.dialog 装饰器放在函数内而不是模块级:@st.dialog 的参数只在装饰时求值,
    本模块是进程级单次导入,放外面会把弹窗标题冻结成首个会话的语言。
    """
    @st.dialog(t("common.delete_dialog_title"))
    def _dialog() -> None:
        names = st.session_state.get(f"{ns}_pending_delete", [])
        st.write(t(confirm_key, n=len(names)))
        if names:
            st.caption(", ".join(names[:10]) + (" …" if len(names) > 10 else ""))
        col_yes, col_no = st.columns(2)
        with col_yes:
            if st.button(t("common.yes"), key=f"{ns}_dlg_yes", width="stretch"):
                failed: list[str] = []
                resolve = path_of or (lambda n: os.path.join(base_dir, n))
                for n in names:
                    path = resolve(n)
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

    _dialog()


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


def render_file_action_bar(ns: str, base_dir: str | None, *, zip_name: str,
                           dirs: bool = False,
                           confirm_key: str = "common.delete_dialog_confirm",
                           path_of: Callable[[str], str] | None = None,
                           update_label: str | None = None,
                           disabled: bool = False) -> list[str] | None:
    """选中条目的批量操作栏:下载 ZIP / 删除(确认弹窗);update_label 提供时多一列
    "更新选中"按钮(管理页用),被点击时返回选中的文件名列表交由页面执行,否则恒返回 None。

    选中项来自多选模式写入的 {ns}_selected(render_file_table 或
    render_select_all_checkbox + render_row_checkbox + collect_selected 组合)。
    dirs=False 条目为 base_dir 下的文件;True 为子目录(删除整目录、zip 内保留目录结构)。
    path_of 提供时按它把条目名解析成完整路径(条目分属多目录的页面用,如管理页的
    来源文件与 merged 整合文件);缺省用 base_dir 拼接。
    disabled 供批量任务运行时禁用整栏按钮。注意"更新选中"按钮自身在点击当轮不得
    变为 disabled:disabled 是 widget 身份的一部分,参数变化会让点击值丢失,
    页面须走"置位 busy → st.rerun → 禁用轮执行"的两段式(见 ui.source_management)。
    """
    resolve = path_of or (lambda n: os.path.join(base_dir, n))
    selected = [n for n in st.session_state.get(f"{ns}_selected", [])
                if os.path.exists(resolve(n))]
    if update_label:
        col_cnt, col_upd, col_dl, col_del = st.columns([2, 2, 2, 2])
    else:
        col_cnt, col_dl, col_del = st.columns([2, 2, 2])
    col_cnt.caption(t("common.selected_count", n=len(selected)))
    update_clicked = None
    if update_label:
        with col_upd:
            if st.button(update_label, key=f"{ns}_batch_update_btn",
                         disabled=disabled or not selected, width="stretch"):
                update_clicked = selected
    with col_dl:
        if dirs:
            data = make_dirs_zip({n: resolve(n) for n in selected}) if selected else b""
        else:
            data = make_ppi_zip({n: resolve(n) for n in selected}) if selected else b""
        st.download_button(t("common.download_selected"), data=data,
                           file_name=zip_name, mime="application/zip",
                           key=f"{ns}_dl_btn",
                           disabled=disabled or not selected, width="stretch")
    with col_del:
        if st.button(t("common.delete_selected"), key=f"{ns}_del_btn",
                     disabled=disabled or not selected, width="stretch"):
            st.session_state[f"{ns}_pending_delete"] = selected
            _confirm_batch_delete(ns, base_dir, dirs=dirs, confirm_key=confirm_key,
                                  path_of=path_of)
    return update_clicked
