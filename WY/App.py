import pandas as pd
import streamlit as st
import config
import Utils
import i18n
from i18n import DEFAULT_LANG, t, get_lang
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# 语言读取注入:i18n 包本身不依赖 streamlit(core 也能安全引用),
# 由宿主在此注册"当前语言"的存取;语言切换直接写 st.session_state["lang"]。
# 必须在本文件任何 t()/get_lang() 调用之前注册。
i18n.set_lang_provider(lambda: st.session_state.get("lang", DEFAULT_LANG))

# 初始化日志:每天一个日志文件,保留两年。
Utils.setup_logging()

st.set_page_config(page_title=t("app.page_title"), layout="centered")

is_en = get_lang() == "en"

# ---------- Custom Minimal CSS ----------
st.markdown("""
<style>
    /* Center container */
    .main > div {
        max-width: 720px;
        margin: 0 auto;
    }
    /* Title area */
    .title-area {
        text-align: center;
        padding: 2.5rem 0 1rem 0;
    }
    .title-area h1 {
        font-size: 2rem;
        font-weight: 600;
        letter-spacing: 1px;
    }
    .title-area p {
        color: #6b7280;
        font-size: 0.95rem;
        margin-top: 0.3rem;
    }
    /* Enlarge input */
    .stTextInput input {
        font-size: 1.2rem;
        padding: 0.7rem 1rem;
    }
    /* Enlarge button */
    .stButton button {
        font-size: 1.1rem;
        padding: 0.5rem 1rem;
    }
    /* Result card */
    .result-card {
        background: #ffffff;
        border-radius: 12px;
        padding: 0.3rem 0.3rem;
    }
    /* Alert styling */
    .stAlert {
        border-radius: 8px;
    }
    /* Hide auto multipage navigation sidebar */
    [data-testid="stSidebar"] {
        display: none;
    }
    /* Top-right corner language + Logs buttons (fixed overlay) */
    div.st-key-top_right_logs {
        position: fixed;
        top: 3.5rem;
        right: 1rem;
        left: 0;
        z-index: 9999;
        display: flex;
        flex-direction: row;
        justify-content: flex-end;
        align-items: center;
        gap: 0.5rem;
        pointer-events: none;
    }
    /* 容器不拦截点击，按钮自身恢复可点击 */
    div.st-key-top_right_logs .stButton {
        pointer-events: auto;
    }
    /* Logs button: #909399 background, white text；高度与蓝色按钮一致（40px、单行） */
    div.st-key-top_right_logs .stButton button {
        background-color: #909399;
        color: #FFFFFF;
        border: 1px solid #909399;
        box-shadow: none;
        height: 40px;
        padding: 4px 12px;
        font-size: 16px;
        line-height: 1;
        white-space: nowrap;
    }
    div.st-key-top_right_logs .stButton button:hover,
    div.st-key-top_right_logs .stButton button:active {
        background-color: #909399;
        color: #FFFFFF;
        border: 1px solid #909399;
        box-shadow: none;
    }
    /* 左下角齿轮 FAB:无底色,只显示一个白色描边齿轮(带中孔);点击弹出设置悬浮窗 */
    div.st-key-gear_fab {
        position: fixed;
        bottom: 1.5rem;
        left: 1.5rem;
        z-index: 9999;
        width: 18px;
        height: 18px;
        pointer-events: auto;
    }
    div.st-key-gear_fab .stButton button {
        width: 18px;
        height: 18px;
        min-height: 18px;
        padding: 0;
        border-radius: 50%;
        background-color: rgba(255, 255, 255, 0.08);
        color: transparent;
        border: none;
        box-shadow: none;
        line-height: 1;
        display: flex;
        align-items: center;
        justify-content: center;
        /* 齿轮用内联 SVG 绘制(白色描边、空心中心,40px):Streamlit 打包的
           图标字体是静态实例,无可变轴,CSS 调不出空心,故不再依赖字体 */
        background-image: url("data:image/svg+xml,%3Csvg xmlns=%22http://www.w3.org/2000/svg%22 viewBox=%220 0 24 24%22 fill=%22none%22 stroke=%22%23ffffff%22 stroke-width=%221.6%22 stroke-linecap=%22round%22 stroke-linejoin=%22round%22%3E%3Cpath d=%22M10.325 4.317c.426-1.756 2.924-1.756 3.35 0a1.724 1.724 0 0 0 2.573 1.066c1.543-.94 3.31.826 2.37 2.37a1.724 1.724 0 0 0 1.065 2.572c1.756.426 1.756 2.924 0 3.35a1.724 1.724 0 0 0-1.066 2.573c.94 1.543-.826 3.31-2.37 2.37a1.724 1.724 0 0 0-2.572 1.065c-.426 1.756-2.924 1.756-3.35 0a1.724 1.724 0 0 0-2.573-1.066c-1.543.94-3.31-.826-2.37-2.37a1.724 1.724 0 0 0-1.065-2.572c-1.756-.426-1.756-2.924 0-3.35a1.724 1.724 0 0 0 1.066-2.573c-.94-1.543.826-3.31 2.37-2.37c1 .608 2.296.07 2.572-1.065z%22/%3E%3Ccircle cx=%2212%22 cy=%2212%22 r=%223%22/%3E%3C/svg%3E");
        background-repeat: no-repeat;
        background-position: center;
        background-size: 18px 18px;
        filter: drop-shadow(0 1px 3px rgba(0, 0, 0, 0.3));
    }
    /* Material 字体齿轮文字隐藏(齿轮由按钮背景 SVG 绘制) */
    div.st-key-gear_fab .stButton button span[role="img"] {
        color: transparent !important;
    }
    div.st-key-gear_fab .stButton button:hover,
    div.st-key-gear_fab .stButton button:active {
        background-color: rgba(255, 255, 255, 0.18);
        color: transparent;
        box-shadow: none;
        filter: drop-shadow(0 2px 6px rgba(0, 0, 0, 0.5));
    }
    /* 语言切换：iOS 风格滑动开关，文本直接嵌在开关轨道内
       （左=中文红 / 右=English蓝，文字始终位于滑块对侧） */
    .st-key-lang_switch_wrap {
        display: flex;
        flex-direction: row;
        align-items: center;
        justify-content: flex-end;
        gap: 0.5rem !important;
        /* 在上层纵向布局中靠右对齐，且宽度收缩为内容宽 */
        align-self: flex-end;
    }
    /* fixed 容器默认不拦截点击，语言开关整体恢复可点击 */
    div.st-key-top_right_logs .st-key-lang_switch_wrap {
        pointer-events: auto;
    }
    /* 轨道：显式放大为带文字的宽开关。Streamlit 1.62 开关结构：
       轨道是 label 内第一个 div，选中态在 label[data-selected] 上 */
    .st-key-lang_switch_wrap label > div:first-of-type {
        position: relative;
        width: 80px !important;
        height: 26px !important;
        padding: 0 !important;
    }
    /* 轨道配色：未选中（中文模式）=红；选中（English模式）=蓝 */
    .st-key-lang_switch_wrap label:not([data-selected="true"]) > div:first-of-type {
        background-color: #E53E3E !important;
        border-color: #E53E3E !important;
    }
    .st-key-lang_switch_wrap label[data-selected="true"] > div:first-of-type {
        background-color: #3B82F6 !important;
        border-color: #3B82F6 !important;
    }
    /* 滑块：绝对定位（20x20，3px 内边距），选中时滑到右端 */
    .st-key-lang_switch_wrap label > div:first-of-type > div {
        position: absolute !important;
        top: 3px !important;
        left: 3px !important;
        width: 20px !important;
        height: 20px !important;
        transform: translateX(0) !important;
    }
    .st-key-lang_switch_wrap label[data-selected="true"] > div:first-of-type > div {
        transform: translateX(54px) !important;
    }
    /* 轨道内文字：显示当前语言，位于滑块对侧 */
    .st-key-lang_switch_wrap label > div:first-of-type::after {
        position: absolute;
        top: 50%;
        transform: translateY(-50%);
        font-size: 12px;
        font-weight: 700;
        line-height: 1;
        letter-spacing: 0.3px;
        color: #FFFFFF;
        pointer-events: none;
    }
    .st-key-lang_switch_wrap label:not([data-selected="true"]) > div:first-of-type::after {
        content: "中文";
        right: 8px;
    }
    .st-key-lang_switch_wrap label[data-selected="true"] > div:first-of-type::after {
        content: "English";
        left: 8px;
    }
    /* Top-left corner management buttons (fixed overlay, no card background) */
    div.st-key-top_left_buttons {
        position: fixed;
        top: 3.5rem;
        left: 1rem;
        z-index: 9999;
        display: flex;
        flex-direction: row;
        gap: 0.6rem;
        width: fit-content;
    }
    /* Management buttons: #167A9C background, #C6E7F5 text, gently float up/down */
    div.st-key-top_left_buttons .stButton button {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
        box-shadow: none;
        animation: gentle-float 2s ease-in-out infinite;
    }
    /* On hover: stop floating, keep colors unchanged */
    div.st-key-top_left_buttons .stButton button:hover,
    div.st-key-top_left_buttons .stButton button:active {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
        box-shadow: none;
        animation: none;
    }
    /* Popover 触发按钮（Source Information Management / PPI-Enrichment & Molecular Docking）同色 */
    .st-key-source_mgmt_popover button,
    .st-key-docking_mgmt_popover button {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
        box-shadow: none;
    }
    .st-key-source_mgmt_popover button:hover,
    .st-key-source_mgmt_popover button:active,
    .st-key-docking_mgmt_popover button:hover,
    .st-key-docking_mgmt_popover button:active {
        background-color: #167A9C;
        color: #C6E7F5;
        border: 1px solid #167A9C;
        box-shadow: none;
    }
    /* 下拉面板选项：浅灰底、深灰 hover 变色，醒目可辨 */
    .st-key-drug_mgmt_btn button,
    .st-key-disease_mgmt_btn button,
    .st-key-ppi_info_btn button,
    .st-key-enrichment_info_btn button,
    .st-key-data_network_btn button,
    .st-key-molecular_docking_btn button {
        background-color: #F0F0F0;
        color: #333333;
        border: 1px solid #D0D0D0;
        width: 100%;
    }
    .st-key-drug_mgmt_btn button:hover,
    .st-key-disease_mgmt_btn button:hover,
    .st-key-ppi_info_btn button:hover,
    .st-key-enrichment_info_btn button:hover,
    .st-key-data_network_btn button:hover,
    .st-key-molecular_docking_btn button:hover {
        background-color: #B8B8B8;
        color: #000000;
        border-color: #888888;
    }
    @keyframes gentle-float {
        0%   { transform: translateY(0); }
        50%  { transform: translateY(-4px); }
        100% { transform: translateY(0); }
    }
</style>
""", unsafe_allow_html=True)

# ---------- Settings floating window (modal dialog, opened by the bottom-left gear) ----------
@st.dialog(t("app.settings.title"))
def _open_settings_dialog() -> None:
    """设置悬浮窗:靶点降噪阈值与黑名单、富集/模块参数,改动即时生效并持久化。"""
    st.caption(t("app.settings.caption"))
    # ChEMBL 黑名单:可在常用选项里勾选,也支持手输任意基因符号
    new_bl = st.multiselect(
        t("app.settings.blacklist"),
        options=list(config.CHEMBL_BLACKLIST_OPTIONS),
        default=list(config.CHEMBL_TARGET_BLACKLIST),
        accept_new_options=True,
        key="set_blacklist",
    )
    def _stepper(label: str, cfg_key: str, lo: int, hi: int, step: int, skey: str) -> int:
        """+- 步进器:只允许按钮加减(不提供键入入口),到达上下限后对应按钮立即置灰。
        用 on_click 回调改值(回调先于本次渲染执行),保证同一次点击里数值与按钮置灰状态同步刷新。"""
        if skey not in st.session_state:
            st.session_state[skey] = int(getattr(config, cfg_key))

        def _dec() -> None:
            st.session_state[skey] = max(lo, int(st.session_state[skey]) - step)

        def _inc() -> None:
            st.session_state[skey] = min(hi, int(st.session_state[skey]) + step)

        val = int(st.session_state[skey])
        col_dec, col_val, col_inc = st.columns([1, 2, 1])
        with col_dec:
            st.button("−", key=f"{skey}_dec", on_click=_dec,
                      disabled=val <= lo, use_container_width=True)
        with col_val:
            st.markdown(f"<p style='text-align:center;margin:0.35rem 0 0'>"
                        f"<b>{val}</b></p>",
                        unsafe_allow_html=True)
        with col_inc:
            st.button("+", key=f"{skey}_inc", on_click=_inc,
                      disabled=val >= hi, use_container_width=True)
        st.caption(label)
        return val

    ot_topn = _stepper(t("app.settings.ot_topn"), "OT_TOP_N",
                       lo=200, hi=1000, step=50, skey="set_ot_topn")
    dg_pages = _stepper(t("app.settings.dg_pages"), "DISGENET_PAGES",
                        lo=1, hi=10, step=1, skey="set_dg_pages")
    new_vals = {
        "CHEMBL_PCHEMBL_THRESHOLD": st.number_input(
            t("app.settings.pchembl"), min_value=0.0, max_value=12.0, step=0.1,
            value=float(config.CHEMBL_PCHEMBL_THRESHOLD), key="set_pchembl"),
        "DGIDB_MIN_SCORE": st.number_input(
            t("app.settings.dgidb"), min_value=0.0, max_value=1.0, step=0.01,
            value=float(config.DGIDB_MIN_SCORE), key="set_dgidb"),
        "CHEMBL_CORROBORATION_REQUIRED": st.checkbox(
            t("app.settings.corroboration"),
            value=bool(config.CHEMBL_CORROBORATION_REQUIRED), key="set_corroboration"),
        "OT_MIN_SCORE": st.number_input(
            t("app.settings.ot"), min_value=0.0, max_value=1.0, step=0.05,
            value=float(config.OT_MIN_SCORE), key="set_ot"),
        "OT_TOP_N": ot_topn,
        "DISGENET_PAGES": dg_pages,
        "ENRICHMENT_PVALUE": st.number_input(
            t("app.settings.enr_p"), min_value=0.001, max_value=0.5, step=0.005,
            value=float(config.ENRICHMENT_PVALUE), key="set_enr_p"),
        "LOUVAIN_RESOLUTION": st.number_input(
            t("app.settings.louvain"), min_value=0.3, max_value=3.0, step=0.1,
            value=float(config.LOUVAIN_RESOLUTION), key="set_louvain"),
        "POCKET_CUTOFF": st.number_input(
            t("app.settings.pocket"), min_value=3.0, max_value=6.0, step=0.5,
            value=float(config.POCKET_CUTOFF), key="set_pocket"),
        "PPI_CONFIDENCE": st.slider(
            t("app.settings.ppi_conf"), min_value=0.1, max_value=0.9, step=0.05,
            value=float(config.PPI_CONFIDENCE), key="set_ppi_conf",
            help=t("app.settings.ppi_conf_help")),
    }
    # 有变化才动手:改 config 运行时值 + 持久化 + 清 merged 整合缓存
    changed = {k: v for k, v in new_vals.items() if v != getattr(config, k)}
    if sorted(new_bl) != sorted(config.CHEMBL_TARGET_BLACKLIST):
        changed["CHEMBL_TARGET_BLACKLIST"] = new_bl
    if changed:
        for _k, _v in changed.items():
            setattr(config, _k, _v)
        # PPI 置信度变更时同步派生阈值:背景互作组(RWR/邻近度)的筛选与缓存文件名都依赖它
        config.INTERACTOME_CONFIDENCE = int(config.PPI_CONFIDENCE * 1000)
        Utils.save_ui_settings(changed)
        Utils.purge_merged_cache()
        st.toast(t("app.settings.toast"))

    # ---------- STRING 数据版本切换 ----------
    # 查询 API 与背景互作组按版本同源派生;切换需二次确认,保存前经下载站探活,
    # 新版本发布后在此输入版本号即可,无需改代码。齿轮弹窗内不能嵌套弹窗,
    # 确认步骤用两段式按钮实现(先「切换到」再「确认切换」)。
    st.divider()
    st.caption(t("app.settings.string_ver_caption"))
    _known = list(config.STRING_KNOWN_VERSIONS)
    _custom_opt = t("app.settings.string_ver_custom_opt")
    _choice = st.selectbox(
        t("app.settings.string_ver"), _known + [_custom_opt],
        index=(_known.index(config.STRING_VERSION)
               if config.STRING_VERSION in _known else len(_known)),
        key="set_string_ver_choice", label_visibility="collapsed")
    if _choice == _custom_opt:
        _raw = st.text_input(
            t("app.settings.string_ver_custom"),
            value=config.STRING_VERSION if config.STRING_VERSION not in _known else "",
            placeholder="v13.0", key="set_string_ver_input")
        _target = config.normalized_string_version(_raw)
        if _raw.strip() and _target is None:
            # 非法格式给内联提示,不再静默禁用按钮让用户摸不着头脑
            st.error(t("app.settings.string_ver_bad_format", raw=_raw.strip()))
    else:
        _target = _choice
    if st.button(
            t("app.settings.string_ver_apply", ver=_target or "?"),
            key="set_string_ver_apply",
            disabled=not _target or _target == config.STRING_VERSION,
            use_container_width=True):
        from core.interactome import string_version_exists

        if string_version_exists(_target):
            st.session_state["set_string_ver_pending"] = _target
        else:
            st.session_state["set_string_ver_pending"] = None
            st.error(t("app.settings.string_ver_missing", ver=_target))
    _pending = st.session_state.get("set_string_ver_pending")
    if _pending and _pending != _target:
        # 待确认期间又改了所选版本:旧确认作废
        st.session_state["set_string_ver_pending"] = None
        _pending = None
    if _pending:
        st.warning(t("app.settings.string_ver_confirm", ver=_pending))
        if st.button(t("app.settings.string_ver_confirm_btn", ver=_pending),
                     key="set_string_ver_confirm", type="primary",
                     use_container_width=True):
            config.STRING_VERSION = _pending
            # 与 config.py 末尾的派生保持同一变换:版本化子域名把 v12.0 写作 version-12-0
            config.STRING_BASE_URL = (
                f"https://version-{_pending[1:].replace('.', '-')}.string-db.org/api/json")
            from core.interactome import clear_runtime_caches

            clear_runtime_caches()  # 图/度表/LCC 上下文按版本隔离,清一次防旧版本残留
            Utils.save_ui_settings({"STRING_VERSION": _pending})
            # 清掉版本控件自身的会话状态:重开齿轮时按新当前版本重新定位选中项,
            # 避免停留在切换前的旧选项上被误点切回
            for _wk in ("set_string_ver_choice", "set_string_ver_input"):
                st.session_state.pop(_wk, None)
            st.session_state["set_string_ver_pending"] = None
            st.toast(t("app.settings.string_ver_done", ver=_pending))
            st.rerun()


# ---------- Title ----------
def _string_ver_status_line() -> str:
    """首页常驻的 STRING 版本行(版本 + 背景网络构建状态):切换版本后不至于遗忘当前用的是哪版数据。"""
    from core.interactome import is_built

    status = t("app.string_ver_built") if is_built("Human") else t("app.string_ver_not_built")
    return t("app.string_ver_line", ver=config.STRING_VERSION, status=status)


st.markdown(f'<div class="title-area"><h1>{t("app.title")}</h1><p>{t("app.subtitle")}</p>'
            f'<p>{_string_ver_status_line()}</p></div>',
            unsafe_allow_html=True)

# ---------- Bottom-left gear: opens the settings modal ----------
with st.container(key="gear_fab"):
    if st.button(":material/settings:", key="gear_fab_btn"):
        _open_settings_dialog()

# ---------- Language sliding switch + Logs button (fixed top-right) ----------
with st.container(key="top_right_logs"):
    # 语言切换开关:左为中文(红色),右为英文(蓝色),文字显示在开关轨道上。
    with st.container(key="lang_switch_wrap"):
        # 若语言在其它页面被切换过,让开关跟随当前语言。
        if st.session_state.get("_lang_switch_synced") != get_lang():
            st.session_state["_lang_switch"] = is_en
            st.session_state["_lang_switch_synced"] = get_lang()
        lang_is_en = st.toggle("Language", key="_lang_switch",
                               label_visibility="collapsed")
        if lang_is_en != is_en:
            st.session_state["lang"] = "en" if lang_is_en else "zh"
            st.session_state["_lang_switch_synced"] = get_lang()
            st.rerun()
    if st.button(t("app.logs_btn"), key="logs_btn"):
        st.switch_page("pages/Logs.py")

# ---------- Management buttons (fixed top-left, open full-screen pages) ----------
with st.container(key="top_left_buttons"):
    # 下拉菜单:跳转到四个信息管理页面。
    with st.popover(t("app.popover_source"), key="source_mgmt_popover"):
        if st.button(t("app.btn_drug_mgmt"), key="drug_mgmt_btn", width="stretch"):
            st.switch_page("pages/Drug_Target_Information_Management.py")
        if st.button(t("app.btn_disease_mgmt"), key="disease_mgmt_btn", width="stretch"):
            st.switch_page("pages/Disease_Target_Information_Management.py")
        if st.button(t("app.btn_ppi_info"), key="ppi_info_btn", width="stretch"):
            st.switch_page("pages/PPI_Information.py")
        if st.button(t("app.btn_enrichment_info"), key="enrichment_info_btn", width="stretch"):
            st.switch_page("pages/Enrichment_Information.py")
    # 下拉菜单:主分析(PPI/富集)与分子对接入口。
    with st.popover(t("app.popover_docking"), key="docking_mgmt_popover"):
        if st.button(t("app.btn_data_network"), key="data_network_btn", width="stretch"):
            st.switch_page("pages/Data_Network.py")
        if st.button(t("app.btn_molecular_docking"), key="molecular_docking_btn", width="stretch"):
            st.session_state["_docking_entry"] = True  # 通知对接页:从入口进入并重置状态
            st.switch_page("pages/Molecular_Docking.py")

st.divider()

# ---------- Search bar ----------
col_drug, col_disease = st.columns(2)
with col_drug:
    drug_name = st.text_input(t("app.drug_label"), value="metformin",
                               label_visibility="collapsed",
                               placeholder=t("app.drug_placeholder"))
with col_disease:
    disease_name = st.text_input(t("app.disease_label"), value="Type 2 diabetes mellitus",
                                  label_visibility="collapsed",
                                  placeholder=t("app.disease_placeholder"))

analyze_btn = st.button(t("app.analyze_btn"), type="primary", width="stretch")

# ---------- Results ----------
if analyze_btn:
    errors = []
    if not drug_name.strip():
        errors.append(t("app.err_drug"))
    if not disease_name.strip():
        errors.append(t("app.err_disease"))

    if errors:
        st.warning(t("app.warn_missing", fields=t("app.and_join").join(errors)))
    else:
        with st.spinner(t("app.spinner_query", drug=drug_name, disease=disease_name)):
            # 四路数据库并发请求(ChEMBL/DGIdb/OpenTargets/DisGeNET);每路写各自的来源 CSV,
            # merged 表由主线程最后写出,并发不会写错文件;进度条随各路完成逐格推进
            progress = st.progress(0.0, text=t("app.progress_prepare"))

            def _on_fetch_progress(done: list[str], total: int) -> None:
                waiting = "、".join(
                    t(f"source.{key}")
                    for key in ("chembl", "dgidb", "opentargets", "disgenet")
                    if key not in done)
                progress.progress(
                    (len(done) / total) if total else 1.0,
                    text=t("app.progress_fetch", done=len(done), total=total,
                           pending=waiting or t("app.progress_done_word")))

            drug_df, disease_df, failed_sources = Utils.get_targets_concurrent(
                drug_name, disease_name, progress_cb=_on_fetch_progress)
            progress.empty()
            # 从表格取出基因符号列表,供韦恩图和统计使用
            drug_targets = drug_df["approvedSymbol"].tolist() if not drug_df.empty else []
            disease_targets = disease_df["approvedSymbol"].tolist() if not disease_df.empty else []

        if failed_sources:
            # 单路失败时结果可能不完整:此时该侧 merged 缓存不会被写入,重试只补拉失败源
            st.warning(t("app.warn_source_failed",
                         sources=t("app.and_join").join(
                             t(f"source.{key}") for key in failed_sources)))

        st.markdown('<div class="result-card">', unsafe_allow_html=True)

        drug_set = set(drug_targets)
        disease_set = set(disease_targets)

        col_a, col_b, col_c = st.columns([1, 1, 1])
        col_a.metric(t("app.metric_drug"), drug_name)
        col_b.metric(t("app.metric_targets"), len(drug_set))
        col_c.metric(t("app.metric_disease_targets"), len(disease_set))

        common = drug_set & disease_set

        st.metric(t("app.metric_common"), len(common))

        # ----- Venn diagram -----
        fig, ax = plt.subplots(figsize=(6, 4))
        # 图内文字保持英文:matplotlib 默认字体不支持中文,中文会显示成方框。
        v = venn2(
            subsets=(len(drug_set - disease_set),
                     len(disease_set - drug_set),
                     len(common)),
            set_labels=(f"{drug_name} targets", f"{disease_name} targets"),
            ax=ax,
        )
        if v and v.get_label_by_id("11"):
            v.get_label_by_id("11").set_fontsize(12)# 设置交集区域字号
        ax.set_title(f"{drug_name} vs {disease_name} — Venn Diagram", fontsize=14)
        st.pyplot(fig)

        with st.expander(t("app.expander_targets"), expanded=False):
            # ----- 药物靶点表 -----
            st.markdown(t("app.drug_targets_md", name=drug_name, count=len(drug_set)))
            if not drug_df.empty:
                st.dataframe(
                    drug_df.sort_values(["source_count", "score_pchembl"],
                                        ascending=False, na_position="last"),
                    width="stretch",
                )
            else:
                st.write(t("app.no_drug_data"))

            st.divider()

            # ----- 疾病靶点表 -----
            st.markdown(t("app.disease_targets_md", name=disease_name, count=len(disease_set)))
            if not disease_df.empty:
                st.dataframe(
                    disease_df.sort_values(["source_count", "score_ot"],
                                           ascending=False, na_position="last"),
                    width="stretch",
                )
            else:
                st.write(t("app.no_disease_data"))

            st.divider()

            # ----- 交集靶点表 -----
            st.markdown(t("app.common_targets_md", drug=drug_name, disease=disease_name,
                          count=len(common)))
            if common:
                # 合并两侧信息：每个交集基因保留药物侧和疾病侧的数据库来源与分数
                common_df = pd.merge(
                    drug_df.drop(columns=["source_count"]),
                    disease_df.drop(columns=["source_count"]),
                    on="approvedSymbol", how="inner",
                )
                # 两侧共 4 个数据库，支持的数据库越多排序越靠前
                common_df["source_total"] = common_df[
                    ["source_chmbl", "source_dgidb", "source_ot", "source_disgenet"]
                ].sum(axis=1)
                st.dataframe(
                    common_df.sort_values(["source_total", "score_pchembl", "score_ot"],
                                          ascending=False, na_position="last"),
                    width="stretch",
                )
            else:
                st.write(t("app.no_common_data"))

        st.markdown('</div>', unsafe_allow_html=True)

        st.success(t("app.success_done"))
else:
    st.info(t("app.info_hint"))
