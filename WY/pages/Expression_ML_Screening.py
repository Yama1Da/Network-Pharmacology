"""GEO 表达谱与 ML 共识筛选页面(M4,模块收官 GUI)。

四区布局(方案 §M4):
① 输入区:疾病名检索 GEO 数据集 / 直接输入 GSE 编号 / 本地矩阵直跑(方案 §6)
② 配置区:分组正则 + 平台排序 CSV 多选(按 analysis_meta.json 元数据展示与过滤,
   只列排序表 *_query_genes.csv)+ top-K/预过滤/bootstrap 参数
③ 运行区:一键执行(DEG → ML 共识 → shortlist),产物经共用管线落盘
   Static/data/downloads/GEOresult/<GSE>_<时间戳>/,与 PPI(_data)/对接(MDresult)
   同属 Web 端产物区;CLI 默认仍落 results/mlscreen(方案 §6),历史列表双根同扫;
   GEO 下载失败降级用旧缓存时,显示黄色警告条 + 强制更新按钮(失效场景②)
④ 结果区:共识表/ROC/shortlist 下载/历史结果(近 3 次,读 summary.json;
   陈旧检测=失效场景③;治理与批量操作在结果管理页)

本页面不新增能力、只换交互形态:全链流程在 core/ml_pipeline.run_screen_chain,
CLI 与页面共用同一套编排,防止两处漂移。
"""
import os
import warnings

import pandas as pd
import streamlit as st

from WY import config
from WY.core.expression import assess_expression_dataset, normalize_gpl
from WY.core.geo import fetch_geo_dataset
from WY.core.ml_pipeline import STAGES, run_screen_chain
from WY.core.ppi import (
    base_disease_of,
    list_analyzed_diseases,
    list_analyzed_drugs,
    list_ranked_csv_candidates,
)
from WY.i18n import t
from WY.ui.widgets import back_button, hide_sidebar, search_selectbox

st.set_page_config(page_title=t("eml.page_title"), layout="wide")

# 隐藏侧边栏;「ML 结果管理」入口仿数据网络页「数据管理」:同位(返回键旁)、同配色(#6347A1)
hide_sidebar("""
    .st-key-ml_mgmt_btn button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }
""")

st.title(t("eml.title"))
st.caption(t("eml.caption"))

col_back, col_mgmt, _ = st.columns([1, 1, 4])
with col_back:
    back_button()
with col_mgmt:
    if st.button(t("mri.title"), key="ml_mgmt_btn"):
        st.switch_page("pages/ML_Results_Information.py")

st.divider()

# ---------- 会话状态 ----------
ss = st.session_state
ss.setdefault("eml", {})
ss.setdefault("eml_busy", False)


def _fmt_size(n: float) -> str:
    """字节数 → KB/MB 显示(下载耗时预估用)。"""
    return f"{n / 1024:.0f} KB" if n < 1048576 else f"{n / 1048576:.1f} MB"


def _fmt_eta(seconds: float) -> str:
    """秒数 → 秒 / 分钟 / 小时+分钟 显示。"""
    if seconds < 90:
        return t("eml.eta_sec", n=max(1, round(seconds)))
    if seconds < 3600:
        return t("eml.eta_min", n=max(1, round(seconds / 60)))
    return t("eml.eta_hr_min", h=int(seconds // 3600), m=round(seconds % 3600 / 60))

# ---------- ① 输入区 ----------
st.subheader(t("eml.step1"))

tab_search, tab_direct, tab_local = st.tabs(
    [t("eml.tab_search"), t("eml.tab_direct"), t("eml.tab_local")])

with tab_search:
    # 疾病名优先从本机 App 分析记录取(用户跑过靶点的疾病):下拉列出干净的
    # 基准疾病名(剥掉 rwrNN 变体后缀),accept_new_options 仍允许手输新名称;
    # 无记录时退回文本框直接输入
    disease_options = list_analyzed_diseases()
    col_d, col_l = st.columns([3, 1])
    with col_d:
        if disease_options:
            disease = search_selectbox(t("eml.disease_label"), disease_options,
                                       key="eml_disease_sel")
            st.caption(t("eml.analyzed_hint", n=len(disease_options)))
        else:
            disease = st.text_input(t("eml.disease_label"), key="eml_disease_txt",
                                    placeholder=t("eml.disease_ph"))
            st.caption(t("eml.analyzed_empty"))
    with col_l:
        limit = st.number_input(t("eml.limit_label"), min_value=5, max_value=50,
                                value=15, key="eml_limit")
    if st.button(t("eml.search_btn"), key="eml_search", type="primary"):
        if not disease.strip():
            st.warning(t("eml.disease_empty"))
        else:
            with st.spinner(t("eml.searching")):
                try:
                    from WY.core.geo import search_geo_datasets
                    results = search_geo_datasets(disease.strip(), limit=int(limit))
                    ss.eml["search_results"] = results
                    ss.eml["search_disease"] = disease.strip()  # 结果带检索上下文
                    # 新检索 = 新候选列表:上一轮的 GSE 选择与全部衍生探测/预检失效,
                    # 不清的话选择框会残留旧 GSE(甚至已判定不可用的那个)
                    ss.eml.pop("selected_gse", None)
                    ss.eml.pop("probe", None)
                    ss.eml.pop("dataset", None)
                    for _wk in ("eml_selected_gse_sel", "eml_selected_gse"):
                        ss.pop(_wk, None)
                except Exception as exc:
                    st.error(t("eml.search_failed", exc=exc))
    if "search_results" in ss.eml and not ss.eml["search_results"].empty:
        # 检索结果跨页持久是 session 设计;但当前疾病与检索时不一致(改选/清空,
        # 或切页回来)要标明来历,避免对着一张对不上号的旧结果表
        _sd = ss.eml.get("search_disease", "")
        if disease.strip().lower() != _sd.lower():
            st.warning(t("eml.search_stale", disease=_sd))
        st.dataframe(ss.eml["search_results"], use_container_width=True,
                     hide_index=True)
        st.caption(t("eml.search_hint"))
    elif "search_results" in ss.eml:
        st.info(t("eml.no_results"))

with tab_direct:
    st.caption(t("eml.direct_hint"))

with tab_local:
    st.caption(t("eml.local_hint"))
    col_le, col_ll = st.columns(2)
    with col_le:
        expr_path = st.text_input(t("eml.local_expr"), key="eml_local_expr",
                                  placeholder=r"D:\data\GSE12345_expr.csv")
    with col_ll:
        labels_path = st.text_input(t("eml.local_labels"), key="eml_local_labels",
                                    placeholder=r"D:\data\labels.csv")
    case_label = st.text_input(t("eml.local_case"), key="eml_local_case",
                               placeholder=t("eml.local_case_ph"))

# 数据集选择:检索后有结果 → 下拉点选(accept_new_options 仍可直填其他编号);
# 未检索 → 文本框直填(tab_direct 路径)
_srch = ss.eml.get("search_results")
_gse_opts = []
if _srch is not None and not _srch.empty and "gse" in _srch.columns:
    _gse_opts = [str(g).strip() for g in _srch["gse"].tolist() if str(g).strip()]
if _gse_opts:
    selected_gse = search_selectbox(t("eml.select_gse"), _gse_opts,
                                    key="eml_selected_gse_sel",
                                    placeholder="GSE63067").upper()
else:
    selected_gse = st.text_input(t("eml.select_gse"), key="eml_selected_gse",
                                 placeholder="GSE63067").strip().upper()
if selected_gse:
    ss.eml["selected_gse"] = selected_gse
    st.info(t("eml.gse_confirmed", gse=selected_gse))
    # 下载耗时预估(方案 §4.1:跨境带宽波动大,选型阶段先知道要等多久)。
    # 每个 GSE 每会话只探测一次(存 session),网络变化后换 GSE 或刷新页面重测;
    # 本地矩阵模式不探测。单个小 Range 请求,不落盘、不写缓存、不影响正式下载。
    if not expr_path.strip() and ss.eml.get("probe", {}).get("gse") != selected_gse:
        with st.spinner(t("eml.est_probing")):
            try:
                from WY.core.geo import probe_series_matrix
                ss.eml["probe"] = {"gse": selected_gse,
                                   "r": probe_series_matrix(selected_gse), "err": ""}
            except Exception as _exc:
                ss.eml["probe"] = {"gse": selected_gse, "r": None, "err": str(_exc)}
        # 平台注释预检(方案 §11):检索结果带平台且为芯片型 → HEAD 一次注释文件,
        # 404 在下载前就亮"不可用"灯。结果随 probe 缓存,每 GSE 每会话只探一次。
        _sr = ss.eml.get("search_results")
        if (_sr is not None and not _sr.empty
                and {"gse", "platform", "type"} <= set(_sr.columns)):
            _row = _sr[_sr["gse"].astype(str).str.upper() == selected_gse]
            if len(_row) and str(_row["type"].iloc[0]) == "芯片":
                _pgpl = str(_row["platform"].iloc[0]).strip()
                if _pgpl:
                    from WY.core.expression import probe_platform_annot
                    ss.eml["probe"]["annot"] = probe_platform_annot(_pgpl)
                    ss.eml["probe"]["gpl"] = normalize_gpl(_pgpl)
        # 矩阵已在本地缓存 → 直接解析做权威判定(基因级矩阵不需要注释文件,
        # 上面基于检索元数据的注释警告对它就是误报 —— GSE281291 即此形态);
        # 未缓存则不做(预检不该替正式下载),保持条件措辞
        from WY.core.geo import _geo_cache_path
        if os.path.isfile(_geo_cache_path(selected_gse)):
            try:
                ss.eml["probe"]["assess"] = assess_expression_dataset(
                    fetch_geo_dataset(selected_gse))
            except Exception:
                pass  # 解析失败不阻断预估展示,警告回退到元数据条件版
        elif (ss.eml["probe"].get("r")
                and ss.eml["probe"]["r"].get("size_bytes", 10 ** 9) < 20 * 1024):
            # 小体积 series matrix(≈3KB)几乎必然是元数据型空表:预检 counts
            # 可用性 —— 有则放行(运行时自动回退 DESeq2),没有才判死
            from WY.core.geo import has_counts_supplement
            ss.eml["probe"]["assess"] = {
                "matrix_empty": True,
                "counts_available": has_counts_supplement(selected_gse),
                "symbol_level": None, "annot": None, "gpl": "", "group_fields": [],
                "n_samples": None, "usable": False}
    _probe = ss.eml.get("probe", {})
    if not expr_path.strip() and _probe.get("gse") == selected_gse:
        _cache_dir = os.path.join(config.CACHE_DIR, config.GEO_SERIES_CACHE_SUBDIR)
        if _probe.get("r"):
            _pr = _probe["r"]
            if _pr["speed_bps"] < 1024:
                # < 1 KB/s:字节仍在流但已无可用带宽,按"近乎断流"处理,
                # 不渲染体积/速度换算出的夸张 ETA
                st.warning(t("eml.est_dead",
                             speed=f"{_pr['speed_bps']:.0f} B/s",
                             path=_cache_dir))
            else:
                _eta_raw = _pr["size_bytes"] / _pr["speed_bps"]
                _eta = min(_eta_raw, 5400)  # ETA 封顶展示,超时部分交给慢速提示
                st.caption(t("eml.est_line",
                             size=_fmt_size(_pr["size_bytes"]),
                             speed=f"{_pr['speed_bps'] / 1024:.0f} KB/s",
                             eta=_fmt_eta(_eta)))
                if not _pr["reliable"] or _eta_raw > 300:  # 预计超 5 分钟或测速被时限截断
                    st.caption(t("eml.est_slow", path=_cache_dir))
        elif _probe.get("err"):
            st.caption(t("eml.est_failed", path=_cache_dir))
        # 警告两级:矩阵已缓存 → 权威判定(探针级且注释缺失才算不可用;
        # 基因级矩阵不需要注释,基于元数据的警告对它就是误报);未缓存 → 元数据条件版
        _as = _probe.get("assess")
        if _as is not None:
            if _as.get("matrix_empty"):
                # 空矩阵 ≠ 判死:有 counts 补充文件时运行时会自动回退 DESeq2
                if _as.get("counts_available"):
                    st.info(t("eml.usa_matrix_empty_counts"))
                else:
                    st.error(t("eml.usa_matrix_empty"))
            elif not _as["symbol_level"] and _as["annot"] == "missing":
                st.error(t("eml.usa_annot_missing", gpl=_as["gpl"]))
        elif (_probe.get("annot") == "missing"
              and ss.eml.get("overview_done_gse") != selected_gse):
            # 总览已跑过 → 权威判定已存,此处的条件版警告不再重复渲染
            # (否则同一页会出现"⚠️大概率不可用 + ❌不可用"两套话)
            st.warning(t("eml.annot_missing_warn", gpl=_probe.get("gpl", "")))

# 数据源判定:本地矩阵(两路都填)与 GSE 二选一(方案 §6,与管线防呆一致)
use_local = bool(expr_path.strip()) and bool(labels_path.strip())
local_partial = bool(expr_path.strip()) != bool(labels_path.strip())
has_gse = bool(selected_gse.strip())
source_conflict = use_local and has_gse

# 可用性锁定(方案 §11:确定性不可用 → 不允许进入配置与运行)。
# 只认权威判定(空矩阵 / 探针级且注释缺失);分组字段启发式不算 —— 分组写在
# 样本标题里的数据集(如 GSE218929)预检看不到,锁掉就误杀了。
_probe_d = ss.eml.get("probe") or {}
_as = _probe_d.get("assess")
unusable = bool(has_gse and _as and _probe_d.get("gse") == selected_gse
                and ((_as.get("matrix_empty") and not _as.get("counts_available", False))
                     or (not _as.get("symbol_level")
                         and _as.get("annot") == "missing")))

st.divider()

# ---------- ② 配置区(数据集被判定不可用时整体锁定) ----------
st.subheader(t("eml.step2"))
if unusable:
    st.error(t("eml.step2_blocked"))

col_case, col_ctrl = st.columns(2)
with col_case:
    case_pat = st.text_input(t("eml.case_pattern"), key="eml_case",
                             placeholder="steatosis|non-alcoholic",
                             disabled=unusable)
with col_ctrl:
    ctrl_pat = st.text_input(t("eml.ctrl_pattern"), key="eml_ctrl",
                             placeholder="healthy", disabled=unusable)

# 分组总览按钮(运行中禁用:重复触发会对同一数据集并发下载,文件锁只能保正确、不能省时间)
if selected_gse.strip() and st.button(t("eml.overview_btn"), key="eml_overview",
                                      disabled=ss.eml_busy):
    ss.eml_busy = True
    try:
        with st.spinner(t("eml.loading")):
            try:
                from WY.core.geo import fetch_geo_dataset, group_fields_overview
                ds = fetch_geo_dataset(selected_gse.strip())
                ss.eml["dataset"] = ds
                if ds.issues and not (ds.table.empty
                                      and ((ss.eml.get("probe") or {}).get("assess")
                                           or {}).get("counts_available")):
                    # 「数据表为空」对 counts 类数据集是 RNA-seq 常态而非问题:
                    # 运行时会自动回退 DESeq2,ℹ️ 提示与下方 counts 列名已经说明
                    # 情况 —— 此时不再渲染黄色警告,避免与 ℹ️ 自相矛盾
                    st.warning(t("eml.dataset_issues",
                                 issues="; ".join(ds.issues)))

                # 分组字段取值分布(counts 数据集同样展示:元数据字段仍是重要参考;
                # 但"正则匹配对象=样本注释全文"的说法对 counts 数据集不成立 ——
                # 运行时匹配的是 counts 列名,提示由下方 counts_overview_hint 给出)
                overview = group_fields_overview(ds)
                st.markdown(t("eml.overview_title"))
                st.dataframe(overview, use_container_width=True, hide_index=True)
                if not ds.table.empty:
                    st.caption(t("eml.overview_hint"))

                # counts 类数据集(空 series matrix):追加样本列清单 —— 正则的
                # 匹配对象;并按样本标题把元数据分组信息(如 treatment)合并进来,
                # 用户能直接看出每个列属于哪个组,正则照着写即可
                if ds.issues and ds.table.empty and selected_gse.strip():
                    try:
                        from WY.core.geo import fetch_gse_count_matrix
                        cnt = fetch_gse_count_matrix(selected_gse.strip())
                        # 候选分组列:取值 2-8 类的 characteristics 列(剔除属性列)
                        _grp_cols = [c for c in ds.samples.columns
                                     if c.startswith("characteristics")
                                     and 2 <= ds.samples[c].astype(str).nunique() <= 8
                                     and not any(k in c.lower() for k in
                                                 ("age", "sex", "gender", "date",
                                                  "batch"))]
                        _title = (ds.samples["title"].astype(str)
                                  if "title" in ds.samples.columns else None)
                        _rows = []
                        for _c in [str(c) for c in cnt["counts"].columns]:
                            _row = {t("eml.col_sample"): _c}
                            if _title is not None:
                                _m = ds.samples[_title == _c]
                                if not _m.empty:
                                    for _gc in _grp_cols:
                                        _row[_gc] = str(_m.iloc[0][_gc])[:64]
                            _rows.append(_row)
                        st.markdown(t("eml.counts_overview_title",
                                      n=len(_rows)))
                        st.dataframe(pd.DataFrame(_rows), use_container_width=True,
                                     hide_index=True)
                        st.caption(t("eml.counts_overview_hint"))
                    except Exception as exc:
                        st.error(t("eml.overview_failed", exc=exc))

                # 可用性预检(方案 §11:空矩阵/平台注释缺失/无病例对照设计,
                # 在选型阶段排除而不是点「开始分析」后才发现)。
                # 结果写回 probe 会话键:Step 2/3 的锁定与警告同一来源
                a = assess_expression_dataset(ds)
                ss.eml.setdefault("probe", {})["assess"] = a
                ss.eml["overview_done_gse"] = selected_gse
                if a["matrix_empty"]:
                    if a.get("counts_available"):
                        st.info(t("eml.usa_matrix_empty_counts"))
                    else:
                        st.error(t("eml.usa_matrix_empty"))
                else:
                    st.caption(t("eml.usa_matrix",
                                 shape=f"{a['shape'][0]}×{a['shape'][1]}",
                                 idkind=t("eml.usa_id_symbol" if a["symbol_level"]
                                          else "eml.usa_id_probe")))
                    if not a["symbol_level"]:
                        if a["annot"] in ("cached", "ok"):
                            st.caption(t("eml.usa_annot_ok", gpl=a["gpl"]))
                        elif a["annot"] == "missing":
                            st.error(t("eml.usa_annot_missing", gpl=a["gpl"]))
                        else:
                            st.caption(t("eml.usa_annot_unknown", gpl=a["gpl"]))
                if a["group_fields"]:
                    st.caption(t("eml.usa_groups",
                                 fields="; ".join(a["group_fields"][:3])))
                elif not (ds.issues and ds.table.empty):
                    # counts 数据集的分组信息在样本列清单里,不渲染"未发现分组字段"
                    st.warning(t("eml.usa_groups_missing"))
                if a["n_samples"] < 2 * config.ML_MIN_SAMPLES_PER_GROUP:
                    st.warning(t("eml.usa_samples_low",
                                 n=a["n_samples"],
                                 min=2 * config.ML_MIN_SAMPLES_PER_GROUP))
            except Exception as exc:
                st.error(t("eml.overview_failed", exc=exc))
    finally:
        ss.eml_busy = False

# 平台排序 CSV 自动发现
col_rank, col_param = st.columns([2, 1])
with col_rank:
    st.markdown(t("eml.ranked_csv_title"))
    # 自动发现:与疾病下拉同源(list_analyzed_pairs),只收排序表 *_query_genes.csv
    # (nodes/edges 是网络结构表,混入会被当排序用,语义错误)。列表跟随上方选定
    # 的疾病自动限定 —— rwr 等参数变体的排序以文件维度并列展示,跨疾病的排序
    # 不再出现,从源头排除"把别的疾病的排序勾进交集"的可能(方案 §M4 按元数据
    # 匹配,不用文件名猜)
    # 排序文件发现与药物过滤下拉的数据源均在 core(list_ranked_csv_candidates /
    # list_analyzed_drugs),页面只做筛选与勾选
    ranked_files = list_ranked_csv_candidates()
    # 当前疾病上下文(下拉或手输,同样剥掉可能手输的 rwrNN 后缀);为空则展示全部。
    # 匹配只允许"用户输的更泛 → 命中更具体的记录"(如输 cancer 带出 gastric cancer),
    # 不做反向:选 gastric cancer 不能把泛癌 cancer 的排序带进交集(疾病混淆)
    _ctx = base_disease_of(disease).lower() if disease.strip() else ""
    if ranked_files:
        # 药物过滤下拉:选项来自本机已分析药物(同疾病下拉的数据源模式),
        # accept_new_options 保留旧的自由输入过滤能力;留空 = 不过滤
        _drugs = list_analyzed_drugs()
        f_drug = ""
        if _drugs:
            f_drug = search_selectbox(t("eml.filter_drug"), _drugs,
                                      key="eml_ranked_drug",
                                      placeholder=t("eml.filter_all")).lower()
        shown = [rf for rf in ranked_files
                 if (not _ctx
                     or _ctx == rf["base_disease"].lower()
                     or _ctx in rf["base_disease"].lower())
                 and (not f_drug or f_drug in rf["drug"].lower())]
        if shown:
            selected_files = st.multiselect(
                t("eml.ranked_csv_select"),
                options=[rf["label"] for rf in shown],
                key="eml_ranked", disabled=unusable)
            ss.eml["ranked_paths"] = [
                rf["path"] for rf in shown if rf["label"] in selected_files]
        else:
            st.caption(t("eml.no_ranked_for_disease"))
    else:
        st.caption(t("eml.no_ranked_csv"))

with col_param:
    top_k = st.number_input(t("eml.topk_label"), min_value=5, max_value=100,
                            value=config.ML_SHORTLIST_TOP_K, key="eml_topk",
                            disabled=unusable)
    # 默认随 config(500 基因 × 30 轮,页面/CLI 统一口径):本机实测 500×30 ML 段约 2 分钟
    # (LOOCV 小样本数据集约 2.5 分钟);SVM-RFE 对特征数近平方级,调大候选池会显著变慢
    prefilter_n = st.number_input(
        t("eml.prefilter_label"), min_value=50, max_value=5000,
        value=config.ML_DEG_PREFILTER_TOP_N, key="eml_prefilter",
        disabled=unusable)
    bootstrap_n = st.number_input(
        t("eml.bootstrap_label"), min_value=5, max_value=200,
        value=config.ML_BOOTSTRAP_ROUNDS, key="eml_bootstrap",
        disabled=unusable)
    st.caption(t("eml.perf_hint"))

st.divider()

# ---------- 运行进度面板 ----------
class _MLProgress:
    """Step3 运行进度面板:总进度条 + 各阶段状态行(⏳ 待执行 / 🔄 进行中 / ✅ 完成)。

    与数据网络页 Build Report 的 _ReportProgress 同款;区别在 ML 管线于脚本
    主线程同步执行,回调里可直接更新组件,无需线程轮询。管线只在阶段开始时
    发键,前序阶段随之置为完成;补充信息行(stage=None)不改变面板状态。
    """

    _ICON = {"pending": "⏳", "running": "🔄", "done": "✅"}

    def __init__(self, steps: list[str]):
        self._steps = list(steps)
        self._states = {k: "pending" for k in self._steps}
        with st.container(border=True):
            st.caption(t("eml.prog.caption"))
            self._bar = st.progress(0.0)
            self._rows = {k: st.empty() for k in self._steps}
        for k in self._steps:
            self._redraw(k)

    def _redraw(self, key: str) -> None:
        self._rows[key].markdown(
            f"{self._ICON[self._states[key]]} {t(f'eml.prog.{key}')}")

    def _set(self, key: str, state: str) -> None:
        self._states[key] = state
        self._redraw(key)

    def start(self, key: str) -> None:
        """某阶段开始:前序全部置为完成,本阶段置为进行中,总进度推进到阶段起点;
        键不在面板里时忽略(管线发出的补充信息行 stage=None 走这里)。"""
        if key not in self._states:
            return
        idx = self._steps.index(key)
        for prev in self._steps[:idx]:
            self._set(prev, "done")
        self._set(key, "running")
        self._bar.progress(idx / len(self._steps))

    def finish(self) -> None:
        """全链完成:全部阶段置为完成,总进度拉满。"""
        for k in self._steps:
            self._set(k, "done")
        self._bar.progress(1.0)


# ---------- ③ 运行区 ----------
st.subheader(t("eml.step3"))
can_run = ((has_gse and case_pat.strip() and ctrl_pat.strip()) or use_local) \
    and not source_conflict and not unusable
if unusable:
    pass  # 锁定原因已在 Step 2 顶部渲染,此处不重复
elif source_conflict:
    st.error(t("eml.source_conflict"))
elif local_partial:
    st.warning(t("eml.local_missing"))
elif not can_run:
    st.info(t("eml.fill_required"))
if ss.eml_busy:
    st.caption(t("eml.busy_hint"))

if st.button(t("eml.run_btn"), key="eml_run", type="primary",
             disabled=(not can_run) or ss.eml_busy):
    ss.eml_busy = True
    # 进度面板:管线发 STAGES 七阶段键,阶段开始即推进入面板;补充信息行
    # (message 为中文)不进面板,避免与界面语言冲突
    prog = _MLProgress(STAGES)

    def _on_progress(stage, _message):
        prog.start(stage)

    try:
        with st.spinner(t("eml.running")):
            try:
                ranked_csvs = {os.path.splitext(os.path.basename(p))[0]: p
                               for p in ss.eml.get("ranked_paths", [])}
                common = dict(prefilter_top=int(prefilter_n),
                              bootstrap_rounds=int(bootstrap_n),
                              ranked_csvs=ranked_csvs, top_k=int(top_k),
                              out_base=config.GEO_RESULT_DIR,
                              out_subdir=None,
                              progress=_on_progress)
                if use_local:
                    res = run_screen_chain(
                        expr=expr_path.strip(), labels=labels_path.strip(),
                        case_label=case_label.strip() or None, **common)
                else:
                    res = run_screen_chain(
                        gse=selected_gse.strip().upper(),
                        case_pattern=case_pat.strip(),
                        control_pattern=ctrl_pat.strip(), **common)
                prog.finish()
                ss.eml["result"] = {
                    "consensus": res["consensus"],
                    "shortlist": res["shortlist"],
                    "funnel": res["funnel"],
                    "deg": res["deg"],
                    "meta": res["meta"],
                    "roc": res["roc"],
                    "roc_png": res.get("roc_png"),
                    "gse": res["gse"] or res["tag"],
                    "out_dir": res["out_dir"],
                    "elapsed": res["summary"]["elapsed_total_s"],
                }
                st.success(t("eml.run_success",
                             gse=res["gse"] or res["tag"],
                             elapsed=res["summary"]["elapsed_total_s"]))
                st.caption(t("eml.saved_to", dir=res["out_dir"]))
                # 降级标记入会话(方案失效场景②):强制更新按钮必须渲染在运行按钮
                # 块之外 —— 块内按钮在本次 run 结束后不再重建,点击会被 Streamlit
                # 静默丢弃(审查 BUG-12)
                if res["summary"].get("geo_cache_fallback"):
                    ss.eml["geo_fallback"] = {
                        "gse": selected_gse.strip().upper(),
                        "date": res["summary"].get("geo_cache_date", "")}
                else:
                    ss.eml.pop("geo_fallback", None)
            except Exception as exc:
                st.error(t("eml.run_failed", exc=exc))
    finally:
        ss.eml_busy = False

# ---------- GEO 下载降级提示(方案失效场景②) ----------
# 常驻渲染:强制更新按钮任何时候可点(点击触发重跑,本块会重新执行,点击不会丢)
if ss.eml.get("geo_fallback"):
    _fb = ss.eml["geo_fallback"]
    st.warning(t("eml.cache_fallback", date=_fb["date"]))
    if st.button(t("eml.force_redownload"), key="eml_force_dl",
                 disabled=ss.eml_busy):
        with st.spinner(t("eml.loading")):
            try:
                fetch_geo_dataset(_fb["gse"], force_refresh=True)
                ss.eml.pop("geo_fallback", None)
                ss.eml.pop("probe", None)  # 清掉探测/预检缓存:注释文件已更新,下一步全部重判
                st.success(t("eml.force_dl_ok"))
            except Exception as exc:
                st.error(t("eml.run_failed", exc=exc))

st.divider()

# ---------- ④ 结果区 ----------
st.subheader(t("eml.step4"))

if "result" not in ss.eml:
    st.info(t("eml.no_result"))
else:
    r = ss.eml["result"]
    m = r["meta"]

    # 指标卡
    col1, col2, col3, col4 = st.columns(4)
    col1.metric(t("eml.metric_samples"), m["n_samples"])
    col2.metric(t("eml.metric_genes"), m["n_genes"])
    col3.metric(t("eml.metric_auc"), m["cv_auc"] if m["cv_auc"] is not None else "N/A")
    col4.metric(t("eml.metric_spearman"), m["algo_spearman_mean"])

    # ROC 曲线
    st.markdown(t("eml.roc_title"))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.metrics import roc_curve, roc_auc_score
        fpr, tpr, _ = roc_curve(r["roc"]["y_true"], r["roc"]["scores"])
        auc_val = roc_auc_score(r["roc"]["y_true"], r["roc"]["scores"])
        fig, ax = plt.subplots(figsize=(4.2, 4.0), dpi=120)
        ax.plot(fpr, tpr, lw=1.6, color="#1f6fb4")
        ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color="#999")
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_title(f"CV-ROC (AUC={auc_val:.3f})", fontsize=10)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)
    except Exception:
        st.warning(t("eml.roc_failed"))
    # ROC PNG 下载(产物目录里的同名文件;页面图表是矢量重绘,论文用 PNG 从这里拿)
    if r.get("roc_png") and os.path.isfile(r["roc_png"]):
        with open(r["roc_png"], "rb") as _f:
            st.download_button(t("eml.dl_roc"), data=_f.read(),
                               file_name=os.path.basename(r["roc_png"]),
                               mime="image/png", key="eml_dl_roc")

    # 共识表(top 50,可排序)
    st.markdown(t("eml.consensus_title"))
    top50 = r["consensus"].head(50).copy()
    top50.insert(0, "gene", top50.index)
    st.dataframe(top50, use_container_width=True, hide_index=True)

    # 漏斗
    st.markdown(t("eml.funnel_title"))
    st.dataframe(r["funnel"], use_container_width=True, hide_index=True)

    # shortlist 下载
    st.markdown(t("eml.shortlist_title",
                  n=int(r["shortlist"]["in_shortlist"].sum())))
    st.dataframe(r["shortlist"], use_container_width=True, hide_index=True)

    # CSV 下载按钮
    col_dl1, col_dl2, col_dl3 = st.columns(3)
    with col_dl1:
        st.download_button(
            t("eml.dl_consensus"),
            data=r["consensus"].to_csv(index_label="gene").encode("utf-8-sig"),
            file_name=f"ml_consensus_{r['gse']}.csv",
            mime="text/csv", key="eml_dl_consensus")
    with col_dl2:
        st.download_button(
            t("eml.dl_shortlist"),
            data=r["shortlist"].to_csv(index=False).encode("utf-8-sig"),
            file_name=f"shortlist_{r['gse']}.csv",
            mime="text/csv", key="eml_dl_shortlist")
    with col_dl3:
        st.download_button(
            t("eml.dl_deg"),
            data=r["deg"].to_csv(index=False).encode("utf-8-sig"),
            file_name=f"deg_{r['gse']}.csv",
            mime="text/csv", key="eml_dl_deg")

st.divider()

# ---------- 历史结果(近 3 次;翻页/搜索/勾选删除/打包下载在结果管理页) ----------
from WY.core.ml_results import list_ml_runs

st.subheader(t("eml.hist_recent"))
_runs = list_ml_runs()[:3]
if not _runs:
    st.caption(t("eml.hist_empty"))
else:
    _pick = st.selectbox(t("eml.hist_recent"),
                         options=[r["name"] for r in _runs],
                         index=0, key="eml_hist3",
                         label_visibility="collapsed")
    _r = next(r for r in _runs if r["name"] == _pick)
    _s = _r["summary"] or {}
    _m = _s.get("ml_meta") or {}
    st.caption(t("eml.hist_line",
                 auc=_m.get("cv_auc") if _m.get("cv_auc") is not None else "N/A",
                 n=_s.get("n_samples", "?"),
                 sl=_s.get("shortlist_count", "?"),
                 ver=_s.get("ml_version", "?")))
    if _r["stale"]:
        _msgs = []
        for _rs in _r["stale"]:
            if _rs["type"] == "version":
                _msgs.append(t("eml.stale_version", old=_rs.get("old"), new=_rs.get("new")))
            elif _rs["type"] == "ranked_changed":
                _msgs.append(t("eml.stale_ranked", name=_rs["name"]))
            else:
                _msgs.append(t("eml.stale_ranked_missing", name=_rs["name"]))
        st.warning(f"{t('eml.stale_title')}: " + "; ".join(_msgs))
    st.caption(f"{t('eml.history_files')}: {_r['path']}")

st.divider()
st.caption(t("eml.footer"))
