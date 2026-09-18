import streamlit as st

import config
from core.enrichment import migrate_legacy_cache_names, parse_enrichment_filename
from core.files import list_cache_files
from core.management import enrichment_display_name
from ui.widgets import (
    CSS_BTN_DOWNLOAD,
    back_button,
    hide_sidebar,
    render_file_action_bar,
    render_file_table,
    show_file_op_result,
)
from i18n import t

st.set_page_config(page_title=t("enrich_info.page_title"), layout="wide")

# 隐藏侧边栏,全屏显示页面;行内 Download 按钮绿底、删除按钮红底白字
hide_sidebar(CSS_BTN_DOWNLOAD + """
    .st-key-enr_files_del_btn button {
        background-color: #F53F3F;
        color: #FFFFFF;
        border: 1px solid #F53F3F;
    }
""")

st.title(t("enrich_info.title"))
st.caption(t("enrich_info.caption"))

col_back, _ = st.columns([1, 4])
with col_back:
    back_button()

st.divider()

# 旧版完整哈希命名的历史文件先改成短哈希命名(幂等),再列出,保证列表命名统一可读
migrate_legacy_cache_names()

show_file_op_result("enr_files")

enr_files = list_cache_files("enrichment", exts=(".csv",))
if enr_files.empty:
    st.info(t("enrich_info.empty"))
else:
    # 数据来源筛选(默认全开):GSEA 按侧别拆成两项;无法识别命名的文件不受筛选影响
    st.caption(t("common.source_filter"))
    c_gd, c_gr, c_gp, c_er = st.columns(4)
    with c_gd:
        show_gsea_disease = st.checkbox(f"GSEA · {t('enrich_file.side_disease')}",
                                        value=True, key="flt_gsea_disease")
    with c_gr:
        show_gsea_drug = st.checkbox(f"GSEA · {t('enrich_file.side_drug')}",
                                     value=True, key="flt_gsea_drug")
    with c_gp:
        show_gprofiler = st.checkbox("gProfiler", value=True, key="flt_gprofiler")
    with c_er:
        show_enrichr = st.checkbox("Enrichr", value=True, key="flt_enrichr")

    def _match_source(n: str) -> bool:
        info = parse_enrichment_filename(n)
        if info is None:
            return True
        if info["kind"] == "gsea":
            if info["side"] == "drug":
                return show_gsea_drug
            if info["side"] == "disease":
                return show_gsea_disease
            return show_gsea_drug or show_gsea_disease  # 无侧别标记的旧文件:任一 GSEA 开着即显示
        if info["source"] == "gProfiler":
            return show_gprofiler
        if info["source"] == "Enrichr":
            return show_enrichr
        return True

    enr_files = enr_files[enr_files["File name"].map(_match_source)]
    if enr_files.empty:
        st.info(t("common.filter_no_match"))
    else:
        # 最新分析排在最前;文件名展示为"药物 → 疾病 · 数据源"的友好名称,搜索仍可匹配原始文件名
        enr_files = enr_files.sort_values("Modified time", ascending=False)
        display_map = {n: enrichment_display_name(n) for n in enr_files["File name"]}
        # 多选文件表 + 批量操作:下载选中(ZIP)/删除(确认弹窗)
        render_file_table(enr_files, config.ENRICHMENT_CACHE_DIR,
                                search_key="search_enrich", display_map=display_map,
                                ns="enr_files")
        render_file_action_bar("enr_files", config.ENRICHMENT_CACHE_DIR,
                                     zip_name="Enrichment_files.zip")
