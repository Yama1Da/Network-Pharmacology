"""轻量 i18n：zh/en 文案字典 + t() 取词。

用法：
    from i18n import t
    st.title(t("app.title"))
    t("star.dialog_confirm", name=name)   # 带参数模板
"""
import streamlit as st

# 默认语言：英文（会话未切换过语言时使用）
DEFAULT_LANG = "en"
# 支持的语言（供后续语言切换按钮渲染）
SUPPORTED_LANGS = {"en": "English", "zh": "中文"}

TEXT: dict[str, dict[str, str]] = {
    # ==================== English (default) ====================
    "en": {
        # ---- 公共 ----
        "common.yes": "Yes",
        "common.no": "No",
        "common.download": "Download",
        "common.back_main": "← Back to Main",
        "common.col_directory": "**Directory**",
        "common.col_created": "**Created time**",
        "common.col_download": "**Download**",
        # ---- App.py（主页） ----
        "app.page_title": "Network Pharmacology Platform",
        "app.title": "🧬 Network Pharmacology Analysis Platform",
        "app.subtitle": "Input drug and disease names to analyze overlapping targets",
        "app.logs_btn": "📜 Logs",
        "app.popover_source": "Source Information Management",
        "app.btn_drug_mgmt": "Drug Target Information Management",
        "app.btn_disease_mgmt": "Disease Target Information Management",
        "app.btn_ppi_info": "PPI Information",
        "app.btn_enrichment_info": "Enrichment Information",
        "app.popover_docking": "PPI-Enrichment & Molecular Docking",
        "app.btn_data_network": "Data Network Pharmacology Analysis",
        "app.btn_molecular_docking": "Molecular Docking",
        "app.drug_label": "Drug name",
        "app.drug_placeholder": "Drug name, e.g. aspirin",
        "app.disease_label": "Disease name",
        "app.disease_placeholder": "Disease name, e.g. cancer",
        "app.analyze_btn": "Analyze",
        "app.err_drug": "drug name",
        "app.err_disease": "disease name",
        "app.and_join": " and ",
        "app.warn_missing": "Please enter {fields}",
        "app.spinner_query": "Querying targets for {drug} and {disease}…",
        "app.metric_drug": "Drug",
        "app.metric_targets": "Targets",
        "app.metric_disease_targets": "Disease targets",
        "app.metric_common": "🎯 Common targets",
        "app.targets_of": "{name} targets",
        "app.venn_title": "{drug} vs {disease} — Venn Diagram",
        "app.expander_targets": "View target lists",
        "app.drug_targets_md": "**💊 {name} targets** — {count} genes (ChEMBL + DGIdb, deduplicated)",
        "app.no_drug_data": "(No drug target data retrieved)",
        "app.disease_targets_md": "**🦠 {name} targets** — {count} genes (Open Targets + DisGeNET, deduplicated)",
        "app.no_disease_data": "(No disease target data retrieved)",
        "app.success_done": "✅ Analysis complete",
        "app.info_hint": "💡 Enter drug and disease names above and click Analyze",
        # ---- Logs.py ----
        "logs.page_title": "Logs",
        "logs.title": "📜 Logs",
        "logs.caption": "Local log files under WY/Static/log — click a file name to view its content",
        "logs.empty": "No log files yet. Run the app to generate daily logs.",
        "logs.tag_current": " 📌 current (today)",
        "logs.modified": "modified",
        "logs.gone": "File no longer exists on disk.",
        "logs.read_failed": "Failed to read {name}: {err}",
        # ---- Data_Management.py ----
        "data_mgmt.page_title": "Data Management",
        "data_mgmt.title": "📁 Data Management",
        "data_mgmt.caption": "Downloaded result groups under Static/data/downloads",
        "data_mgmt.back": "← Back to Network",
        "data_mgmt.empty": "No downloaded result groups yet. Run a PPI analysis first.",
        # ---- Star_Management.py ----
        "star.page_title": "Star Management",
        "star.title": "⭐ Star Management",
        "star.caption": "Starred docking results (MDresult/star/) — exempt from auto-cleanup, manual deletion only",
        "star.back": "← Back to Molecular Docking",
        "star.dialog_title": "Delete Starred Result",
        "star.dialog_confirm": 'Are you sure you want to delete "{name}"? This cannot be undone.',
        "star.deleted": "Starred result deleted",
        "star.empty": "No starred results yet. Save a result on the Molecular Docking page and click ☆ to star it.",
        "star.total": "{count} starred result(s) in total. Download to keep local copies; deletion is irreversible.",
        "star.col_name": "**Name**",
        "star.col_time": "**Starred time**",
        "star.col_delete": "**Delete**",
    },
    # ==================== 中文 ====================
    "zh": {
        # ---- 公共 ----
        "common.yes": "是",
        "common.no": "否",
        "common.download": "下载",
        "common.back_main": "← 返回主页",
        "common.col_directory": "**目录**",
        "common.col_created": "**创建时间**",
        "common.col_download": "**下载**",
        # ---- App.py（主页） ----
        "app.page_title": "网络药理学平台",
        "app.title": "🧬 网络药理学分析平台",
        "app.subtitle": "输入药物与疾病名称，分析交集靶点",
        "app.logs_btn": "📜 日志",
        "app.popover_source": "信息来源管理",
        "app.btn_drug_mgmt": "药物靶点信息管理",
        "app.btn_disease_mgmt": "疾病靶点信息管理",
        "app.btn_ppi_info": "PPI 信息",
        "app.btn_enrichment_info": "富集信息",
        "app.popover_docking": "PPI-富集与分子对接",
        "app.btn_data_network": "数据网络药理学分析",
        "app.btn_molecular_docking": "分子对接",
        "app.drug_label": "药物名称",
        "app.drug_placeholder": "药物名称，如 aspirin",
        "app.disease_label": "疾病名称",
        "app.disease_placeholder": "疾病名称，如 cancer",
        "app.analyze_btn": "分析",
        "app.err_drug": "药物名称",
        "app.err_disease": "疾病名称",
        "app.and_join": " 和 ",
        "app.warn_missing": "请输入{fields}",
        "app.spinner_query": "正在查询 {drug} 与 {disease} 的靶点…",
        "app.metric_drug": "药物",
        "app.metric_targets": "药物靶点",
        "app.metric_disease_targets": "疾病靶点",
        "app.metric_common": "🎯 交集靶点",
        "app.targets_of": "{name} 的靶点",
        "app.venn_title": "{drug} vs {disease} — 韦恩图",
        "app.expander_targets": "查看靶点列表",
        "app.drug_targets_md": "**💊 {name} 靶点** — {count} 个基因（ChEMBL + DGIdb，去重）",
        "app.no_drug_data": "（未获取到药物靶点数据）",
        "app.disease_targets_md": "**🦠 {name} 靶点** — {count} 个基因（Open Targets + DisGeNET，去重）",
        "app.no_disease_data": "（未获取到疾病靶点数据）",
        "app.success_done": "✅ 分析完成",
        "app.info_hint": "💡 在上方输入药物与疾病名称，然后点击「分析」",
        # ---- Logs.py ----
        "logs.page_title": "日志",
        "logs.title": "📜 日志",
        "logs.caption": "WY/Static/log 下的本地日志文件 — 点击文件名查看内容",
        "logs.empty": "暂无日志文件。运行应用后将生成按天的日志。",
        "logs.tag_current": " 📌 当前（今天）",
        "logs.modified": "修改于",
        "logs.gone": "文件已不在磁盘上。",
        "logs.read_failed": "读取 {name} 失败：{err}",
        # ---- Data_Management.py ----
        "data_mgmt.page_title": "数据管理",
        "data_mgmt.title": "📁 数据管理",
        "data_mgmt.caption": "Static/data/downloads 下已下载的结果分组",
        "data_mgmt.back": "← 返回网络分析",
        "data_mgmt.empty": "暂无已下载的结果分组。请先运行 PPI 分析。",
        # ---- Star_Management.py ----
        "star.page_title": "收藏管理",
        "star.title": "⭐ 收藏管理",
        "star.caption": "收藏的对接成果（MDresult/star/ 区）— 不受自动清理影响，仅可手动删除",
        "star.back": "← 返回分子对接",
        "star.dialog_title": "删除收藏",
        "star.dialog_confirm": "确定要删除收藏「{name}」吗？删除后不可恢复。",
        "star.deleted": "已删除该收藏",
        "star.empty": "暂无收藏。在分子对接页保存成果后点击 ☆ 即可收藏。",
        "star.total": "共 {count} 份收藏。下载到本地保存，删除后不可恢复。",
        "star.col_name": "**收藏名称**",
        "star.col_time": "**收藏时间**",
        "star.col_delete": "**删除**",
    },
}


def get_lang() -> str:
    """当前语言（session 级，跨页面持久）；默认英文。仅读 session_state，可在 set_page_config 前调用。"""
    return st.session_state.get("lang", DEFAULT_LANG)


def set_lang(lang: str) -> None:
    """切换语言（仅接受 TEXT 中已定义的语言）。"""
    if lang in TEXT:
        st.session_state["lang"] = lang


def t(key: str, **kwargs) -> str:
    """按键取当前语言文案；缺失时回退默认语言，再缺失返回 key 本身。"""
    text = TEXT.get(get_lang(), {}).get(key)
    if text is None:
        text = TEXT[DEFAULT_LANG].get(key, key)
    return text.format(**kwargs) if kwargs else text
