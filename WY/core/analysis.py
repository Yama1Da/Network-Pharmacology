"""Build Report 用例编排：一次"PPI 网络药理学分析"的完整流程、并发任务与结果落盘。

模块与依赖方向（上层引下层，禁止反向）：
    analysis → targets / ppi / network_algorithms / enrichment
             → config（路径与参数常量；按 core 惯例函数内延迟导入）

界面层（Streamlit 页面 / 未来 API 端点）只做三件事：收集参数 → 调 run_ppi_analysis →
展示/下发结果。进度经 AnalysisProgress 协议回调，编排层不依赖任何 UI 框架；
matplotlib 的 pyplot 全局状态非线程安全，所有绘图统一留在调用线程执行。

与 cli.py 的关系：CLI 的 cmd_analyze 等命令目前各自编排 core 调用，后续可逐步
收敛到本模块（网页端 Build Report 与 CLI 共享同一份分析语义）。
"""

import logging
import re
from concurrent.futures import ThreadPoolExecutor, wait
from typing import Protocol, TypedDict

import pandas as pd

from .enrichment import (
    build_ranked_list,
    fetch_enrichment_gsea_prerank,
    gsea_gene_sets,
    render_drug_target_pathway_plot,
    render_enrichment_plots,
    run_enrichment,
)
from .network_algorithms import (
    ModuleAnalysis,
    compute_module_analysis,
    run_proximity_cached,
    run_rwr_cached,
)
from .ppi import (
    build_main_figures,
    build_ppi_graph,
    compute_mcc_scores,
    compute_ppi_metrics,
    fetch_string_ppi,
    mark_hub_genes,
    save_gsea_results,
    save_module_computation,
    save_ppi_results,
    save_proximity_results,
    save_query_genes,
    save_rwr_results,
)
from .targets import (
    get_common_targets,
    get_disease_targets_merged,
    get_drug_targets_merged,
)


class AnalysisProgress(Protocol):
    """进度回调协议：编排层经它汇报阶段状态与总进度，不依赖任何 UI 框架。
    Streamlit 页面的 _ReportProgress 与未来 API 的 job 进度记录都实现这两个方法。"""

    def set(self, key: str, state: str) -> None:
        """更新某阶段状态（key 为 targets/string/graph/... 阶段名，state 为 pending/running/done）。"""
        ...

    def advance(self, fraction: float) -> None:
        """推进总进度条，fraction ∈ [0, 1]。"""
        ...


class GseaSideInfo(TypedDict):
    """单侧（药物/疾病）GSEA prerank 的结果条目。"""
    df: pd.DataFrame   # 该侧 GSEA 结果表（失败/无结果时为空表）
    stats: dict        # build_ranked_list 的排秩统计
    reason: str        # 失败/缓存原因，成功为空串


class PpiAnalysisResult(TypedDict):
    """run_ppi_analysis 的返回契约：一次完整分析的内存结果与落盘路径。
    展示层只读本结构，不重算；改字段名须同步展示层所有读取点。"""
    genes: list[str]
    edges_df: pd.DataFrame
    graph: "nx.Graph"
    metrics_df: pd.DataFrame
    hub_df: pd.DataFrame
    hub_rule: str
    hub_method_eff: str
    expansion_df: pd.DataFrame | None
    html_content: str
    static_fig: "matplotlib.figure.Figure | None"
    string_fig: "matplotlib.figure.Figure | None"
    cyto_fig: "matplotlib.figure.Figure | None"
    mcc_map: dict[str, float]
    core_fig: "matplotlib.figure.Figure | None"
    core_genes: list[str]
    enrichment_df: pd.DataFrame
    enrichment_source: str
    enrichment_reason: str
    enrichment_bar_fig: "matplotlib.figure.Figure | None"
    enrichment_dot_fig: "matplotlib.figure.Figure | None"
    dtp_fig: "matplotlib.figure.Figure | None"
    gsea: dict[str, GseaSideInfo]        # 仅 names 模式非空
    proximity: dict | None               # 仅 names 模式有值
    rwr: dict                            # 仅 names 模式非空
    modules: ModuleAnalysis | None       # 模块分析失败时为 None
    saved: dict                          # save_analysis_outputs 的落盘路径字典
    prefix: str
    params: tuple                        # (mode, drug, disease, genes_text, score,
                                         #  network_type, rwr_expand, hub_top_n, hub_method)


def pair_prefix(drug: str, disease: str) -> str:
    """药病对的缓存 label / 结果目录前缀(唯一定义点,网页 Build Report 与 CLI run_pair 共用):
    非单词字符统一压成 _。PPI 边表/富集/邻近度/RWR 缓存文件名都带它,
    网页跑过的药病对 CLI 直接命中缓存。"""
    return re.sub(r"[^\w\-]+", "_", f"{drug}_{disease}")


class RwrExpansion(TypedDict):
    """expand_targets_with_rwr 的信息字典契约。"""
    attempted: bool    # 是否真正尝试了游走(未启用或基因集为空时为 False)
    ok: bool           # attempted 且成功
    added: list[dict]  # 成功时的扩展条目 [{"gene": str, "score": float}],已按 (-score, 基因名) 排序
    cached: bool       # 本次游走是否命中缓存
    error: str | None  # attempted 但失败的原因


def expand_targets_with_rwr(genes: list[str], rwr_expand: int,
                            prefix: str) -> tuple[list[str], str, RwrExpansion]:
    """RWR 查询集扩展(Step 0.5 共享语义,网页 Build Report 与 CLI run_pair 共用):
    以当前基因集为种子在背景互作组上游走,取得分最高的前 rwr_expand 个非种子并入查询集。

    游走结果只依赖种子与 α,与 N 无关,游走缓存 label 用 f"expand_{prefix}" 不带 _rwr{N}
    (改 N 复用同一份游走);只有返回的 prefix 带 _rwr{N},供下游(STRING/富集)命名。
    未启用(<=0)或基因集为空时不尝试,原样返回;失败时按原基因集返回并把原因写入
    info["error"],调用方自行降级,不让分析挂掉。prefix 后缀标记"请求了扩展":
    成功时查询集确已扩展;失败时查询集仍是原交集,但下游缓存/目录沿用同一命名规则
    便于对照(旧 CLI 失败路径不带后缀,收敛时已统一到网页规则)。"""
    info: RwrExpansion = {"attempted": False, "ok": False, "added": [],
                          "cached": False, "error": None}
    if int(rwr_expand) <= 0 or not genes:
        return list(genes), prefix, info
    rwr_res = run_rwr_cached(genes, label=f"expand_{prefix}")
    info["attempted"] = True
    if rwr_res.get("ok"):
        # 平局按基因名稳定排序,保证扩展清单可复现
        top = sorted(rwr_res.get("top_nonseed", []),
                     key=lambda r: (-float(r["score"]), str(r["gene"])))[:int(rwr_expand)]
        info["ok"] = True
        info["added"] = top
        info["cached"] = bool(rwr_res.get("cached"))
        genes = sorted(set(genes) | {str(r["gene"]) for r in top})
    else:
        info["error"] = rwr_res.get("error") or "unknown"
    return genes, f"{prefix}_rwr{int(rwr_expand)}", info


def compute_ppi_modules(graph, organism, directory, prefix) -> dict:
    """模块检测 + 富集 + 结果落盘（交互区手动重跑入口）。目录须已存在。"""
    return save_module_computation(directory, prefix,
                                   compute_module_analysis(graph, organism, label_prefix=prefix))


def gsea_prerank_side(side: str, drug_name: str, disease_name: str, prefix: str):
    """对药物/疾病一侧构建排秩表并运行 GSEA prerank（Build Report 自动流程与手动重跑共用）；
    返回 (结果表, 排秩统计, 状态原因)。"""
    if side == "drug":
        src_df = get_drug_targets_merged(drug_name.strip())
        score_col = "score_pchembl"
    else:
        src_df = get_disease_targets_merged(disease_name.strip())
        score_col = "score_ot"
    rnk, rstats = build_ranked_list(src_df, "approvedSymbol", score_col)
    gsea_status: dict = {}
    # 缓存文件名里侧别(drug/disease)放在最前,避免与疾病名中可能出现的"disease"字样混淆
    gdf = fetch_enrichment_gsea_prerank(
        rnk, libs=gsea_gene_sets(), label=f"{side}_{prefix}", status=gsea_status)
    return gdf, rstats, gsea_status.get("reason", "")


def run_enrichment_task(genes: list[str], organism: str, prefix: str):
    """功能富集(g:Profiler 主 + Enrichr 备,输入为交集靶点),返回 (结果表, 数据源, 状态原因)。"""
    from WY import config

    enr_status: dict = {}
    df, source = run_enrichment(genes, organism=organism,
                                      pvalue=config.ENRICHMENT_PVALUE, label=prefix,
                                      status=enr_status)
    return df, source, enr_status.get("reason", "")


def run_gsea_task(drug_name: str, disease_name: str, prefix: str):
    """GSEA 排序富集(仅 names 模式):药物侧按 pChEMBL 活性分排序、疾病侧按 Open Targets
    关联分排序,检验通路基因是否整体偏向高排名。失败只返回空结果,不影响主报告。"""
    gsea_results: dict = {}
    gsea_drug_df = gsea_disease_df = None
    try:
        for side_label in ("drug", "disease"):
            gdf, rstats, reason = gsea_prerank_side(side_label, drug_name, disease_name, prefix)
            gsea_results[side_label] = {"df": gdf, "stats": rstats, "reason": reason}
            if side_label == "drug":
                gsea_drug_df = gdf if not gdf.empty else None
            else:
                gsea_disease_df = gdf if not gdf.empty else None
    except Exception:
        logging.exception("GSEA step failed during report build — skipped gracefully")
        gsea_results, gsea_drug_df, gsea_disease_df = {}, None, None
    return gsea_results, gsea_drug_df, gsea_disease_df


def run_network_tests_task(drug_name: str, disease_name: str, prefix: str):
    """网络邻近度 + RWR(仅 names 模式):检验药物靶点与疾病基因在互作网络上的距离,
    并以 RWR 给非种子蛋白排序。两者共用背景网络上下文,放在同一任务里顺序执行。
    各步失败只跳过自身(邻近度返回哨兵,RWR 置空),不影响主报告。"""
    proximity_result: dict | None = None
    rwr_results: dict = {}
    drug_df_full = pd.DataFrame()
    disease_df_full = pd.DataFrame()
    try:
        drug_df_full = get_drug_targets_merged(drug_name.strip())
        disease_df_full = get_disease_targets_merged(disease_name.strip())
        proximity_result = run_proximity_cached(
            drug_df_full["approvedSymbol"].tolist(),
            disease_df_full["approvedSymbol"].tolist(),
            label=prefix,
        )
    except Exception:
        logging.exception("Proximity step failed during report build — skipped gracefully")
        proximity_result = {"ok": False, "error": "exception", "coverage": {}, "cached": False}
    if not (drug_df_full.empty or disease_df_full.empty):
        try:
            drug_syms = {str(g).strip().upper() for g in drug_df_full["approvedSymbol"]}
            disease_syms = {str(g).strip().upper() for g in disease_df_full["approvedSymbol"]}
            rwr_results["disease_seeds"] = run_rwr_cached(
                list(disease_syms), label=f"{prefix}_disease_seeds", annotate=drug_syms)
            rwr_results["drug_seeds"] = run_rwr_cached(
                list(drug_syms), label=f"{prefix}_drug_seeds", annotate=disease_syms)
        except Exception:
            logging.exception("RWR step failed during report build — skipped gracefully")
            rwr_results = {}
    return proximity_result, rwr_results


def run_modules_task(graph, organism: str, prefix: str):
    """模块检测+富集的纯计算(不落盘):结果目录由 save_ppi_results 重建(会先删除旧目录),
    模块文件必须在主线程保存完主结果后再写,否则会被清掉。失败返回 None,只跳过模块分析。"""
    try:
        return compute_module_analysis(graph, organism, label_prefix=prefix)
    except Exception:
        logging.exception("Module detection step failed during report build — skipped gracefully")
        return None


def save_analysis_outputs(*, graph, static_fig, metrics_df, edges_df, html_content, prefix,
                          enrichment_df, enrichment_bar_fig, enrichment_dot_fig,
                          core_fig, dtp_fig, gsea_drug_df, gsea_disease_df,
                          string_fig, cyto_fig,
                          expansion_df, rwr_results, proximity_result, modules_data,
                          drug_name: str = "", disease_name: str = "",
                          progress=None) -> tuple[dict, dict | None]:
    """Step 5/5.5:全部结果保存到下载目录(同名分析旧目录会被替换,不会堆积)。

    模块文件必须在主目录重建后落盘,否则会被 save_ppi_results 的旧目录清理删掉。
    drug_name/disease_name 写入 analysis_meta.json,供 ML 筛选页按元数据发现排序 CSV。
    返回 (saved 路径字典, modules_result)。"""
    def _prog(step: str | None = None, state: str | None = None,
              fraction: float | None = None) -> None:
        if progress is not None:
            if step and state:
                progress.set(step, state)
            if fraction is not None:
                progress.advance(fraction)

    _prog("save", "running", 0.80)
    saved = save_ppi_results(
        graph, static_fig, metrics_df, edges_df, html_content, prefix=prefix,
        enrichment_df=enrichment_df if not enrichment_df.empty else None,
        enrichment_bar_fig=enrichment_bar_fig,
        enrichment_dot_fig=enrichment_dot_fig,
        core_fig=core_fig,
        target_pathway_fig=dtp_fig,
        string_fig=string_fig,
        cyto_fig=cyto_fig,
        drug=drug_name, disease=disease_name,
    )
    # GSEA 结果 CSV 保存进分析目录（空结果不保存），打包下载时自动包含
    gsea_saved = save_gsea_results(saved["directory"], prefix, gsea_drug_df, gsea_disease_df)
    saved.update(gsea_saved)
    # 查询集清单：RWR 扩展启用时记录种子/扩展来源与得分,未启用不落盘
    saved.update(save_query_genes(saved["directory"], prefix, expansion_df))
    # RWR 结果保存：种子基因（全量）与非种子候选蛋白（前若干名）合并为一张表，is_seed 列区分；
    # 空结果或失败时不保存。
    if rwr_results:

        def _rwr_df(r: dict):
            if not r.get("ok"):
                return None
            frames = []
            if r.get("seeds"):
                frames.append(pd.DataFrame(r["seeds"]).assign(is_seed=True))
            if r.get("top_nonseed"):
                frames.append(pd.DataFrame(r["top_nonseed"]).assign(is_seed=False))
            if not frames:
                return None
            return pd.concat(frames, ignore_index=True)

        rwr_saved = save_rwr_results(
            saved["directory"], prefix,
            disease_df=_rwr_df(rwr_results.get("disease_seeds", {})),
            drug_df=_rwr_df(rwr_results.get("drug_seeds", {})),
        )
        saved.update(rwr_saved)

    # 邻近度结果保存:报告页只作展示,这份 CSV 是它的唯一持久记录;失败/未计算不保存
    if proximity_result and proximity_result.get("ok"):
        saved.update(save_proximity_results(saved["directory"], prefix, proximity_result))

    # Step 5.5: 模块检测结果落盘（主目录已重建完毕,此时写入不会被清掉）
    modules_result: dict | None = None
    if modules_data is not None:
        try:
            modules_result = save_module_computation(saved["directory"], prefix, modules_data)
            saved.update(modules_result["saved_paths"])
        except Exception:
            logging.exception("Module detection step failed during report build — skipped gracefully")
            modules_result = None
    _prog("save", "done", 1.0)
    return saved, modules_result


def run_ppi_analysis(mode, drug_name, disease_name, gene_text, score,
                     network_type: str = "full",
                     rwr_expand: int = 0, hub_top_n: int = 10, hub_method: str = "degree",
                     progress: AnalysisProgress | None = None) -> PpiAnalysisResult:
    """执行完整分析（交集→[RWR扩展]→STRING→建图→指标→富集→可视化→保存），返回结果字典。
    物种固定为人类（疾病靶点数据源仅有人类数据）；PPI 建网只隐藏孤立节点、不做最低连接度剪枝
    （会抵消 RWR 扩展拉入的外周靶点）；RWR 扩展只作用于建网/拓扑/hub/模块/富集,邻近度与 GSEA 仍基于原始靶点集。"""
    from WY import config

    species_id = config.SPECIES_MAP["Human"]["ncbi_taxon"]
    organism = config.SPECIES_MAP["Human"]["gprofiler"]
    def _prog(step: str | None = None, state: str | None = None,
              fraction: float | None = None) -> None:
        if progress is not None:
            if step and state:
                progress.set(step, state)
            if fraction is not None:
                progress.advance(fraction)

    # 疾病靶点目前只有 Open Targets 人类数据，因此整个分析按人类处理，
    # 界面也只开放 Human 选项；多物种支持待后续接入相应疾病数据库后再开放。
    # Step 0: 确定输入基因列表
    _prog("targets", "running", 0.02)
    if mode == "names":
        genes = get_common_targets(drug_name, disease_name, top_n=config.STRING_QUERY_LIMIT)
        prefix = pair_prefix(drug_name, disease_name)
    else:
        tokens = re.split(r"[\s,;]+", gene_text)
        genes = sorted({t.strip().upper() for t in tokens if t.strip()})
        prefix = "manual"

    # Step 0.5: RWR 扩展(可选),共享语义见 expand_targets_with_rwr;扩展失败按原始交集
    # 继续并保留纯种子清单,不让整个分析挂掉
    expansion_df: pd.DataFrame | None = None
    if int(rwr_expand) > 0 and genes:
        seed_rows = [{"gene": g, "source": "seed", "rwr_score": ""} for g in genes]
        genes, prefix, rwr_exp = expand_targets_with_rwr(genes, int(rwr_expand), prefix)
        if rwr_exp["ok"]:
            expansion_df = pd.DataFrame(
                seed_rows + [{"gene": str(r["gene"]), "source": "expanded",
                              "rwr_score": float(r["score"])} for r in rwr_exp["added"]])
        else:
            logging.warning("RWR expansion skipped: %s", rwr_exp["error"])
            expansion_df = pd.DataFrame(seed_rows)
    _prog("targets", "done", 0.08)

    # Step 1: 从 STRING 获取蛋白相互作用（按置信度分数过滤,孤立节点在下游展示时隐藏）
    _prog("string", "running", 0.08)
    edges_df = fetch_string_ppi(genes, species=species_id, score=score,
                                      network_type=network_type, label=prefix)
    _prog("string", "done", 0.30)

    # Step 2-4: 建图 + 拓扑指标 + Hub 筛选(排名制 top-N + 小网络保护)。
    # hub_method="mcc" 时 MCC 须在标记前算好,后面核心靶点图复用同一份结果。
    _prog("graph", "running", 0.30)
    graph = build_ppi_graph(edges_df, genes)
    metrics_df = compute_ppi_metrics(graph)
    mcc_map: dict[str, float] = {}
    if hub_method == "mcc":
        try:
            mcc_map = compute_mcc_scores(graph)
        except Exception:
            logging.exception("MCC computation failed — hub ranking falls back to degree")
    metrics_df = mark_hub_genes(metrics_df,
                                      top_n=hub_top_n, method=hub_method,
                                      mcc_map=mcc_map,
                                      small_net_factor=config.HUB_SMALL_NET_FACTOR)
    hub_df = metrics_df[metrics_df["is_hub"]]
    hub_rule = str(metrics_df["hub_rule"].iloc[0]) if len(metrics_df) else ""
    # MCC 排名被 mark_hub_genes 回退为 degree 时,metrics_df 不会出现 mcc 列
    hub_method_eff = ("degree"
                      if hub_rule.startswith("top") and hub_method == "mcc"
                      and "mcc" not in metrics_df.columns else hub_method)
    _prog("graph", "done", 0.36)

    # ---- Step 4.5-5.5(数据部分):富集/GSEA/邻近度+RWR/模块检测四个互相独立的任务并发执行 ----
    # 只把 API 请求与本地计算放进线程池,墙钟时间从"各步之和"变成"最慢一步"。
    # matplotlib 的 pyplot 是全局状态、非线程安全,所有绘图统一留在调用线程做。
    with ThreadPoolExecutor(max_workers=4, thread_name_prefix="report") as pool:
        fut_enr = pool.submit(run_enrichment_task, genes, organism, prefix)
        fut_gsea = (pool.submit(run_gsea_task, drug_name, disease_name, prefix)
                    if mode == "names" else None)
        fut_net = (pool.submit(run_network_tests_task, drug_name, disease_name, prefix)
                   if mode == "names" else None)
        fut_mod = pool.submit(run_modules_task, graph, organism, prefix)

        fut_key = {f: k for f, k in ((fut_enr, "enrichment"), (fut_gsea, "gsea"),
                                     (fut_net, "network"), (fut_mod, "modules")) if f is not None}
        for k in fut_key.values():
            _prog(k, "running")
        _prog(fraction=0.40)

        # 调用线程在等待期间做 pyvis 交互图、静态图、STRING/Cytoscape 风格图与核心靶点图渲染(纯本地计算);
        # MCC(CytoHubba 默认算法)与核心靶点图失败只跳过自身,不影响主报告。
        # hub_method="mcc" 时已在 Hub 筛选阶段算过 MCC,此处直接复用。
        (html_content, static_fig, mcc_map, hub_df, core_fig, core_genes,
         string_fig, cyto_fig) = build_main_figures(
            graph, metrics_df, edges_df, hub_method, mcc_map)
        _prog(fraction=0.48)

        # 轮询任务完成情况并刷新进度回调(如 Streamlit 页面则只能在脚本主线程更新组件)
        pending = set(fut_key)
        per_task = 0.28 / len(pending)  # 任务全完成时进度从 0.48 推进到 0.76
        while pending:
            done, pending = wait(pending, timeout=0.25)
            for f in done:
                _prog(fut_key[f], "done")
            if done:
                _prog(fraction=0.48 + per_task * (len(fut_key) - len(pending)))

        enrichment_df, enrichment_source, enrichment_reason = fut_enr.result()
        enrichment_bar_fig = enrichment_dot_fig = None
        if not enrichment_df.empty:
            enrichment_bar_fig, enrichment_dot_fig = render_enrichment_plots(enrichment_df)
        # 药物-靶点-通路网络图(纯本地绘图,基于富集结果):无 KEGG 通路时为 None,
        # 失败只跳过自身
        dtp_fig = None
        if not enrichment_df.empty:
            try:
                dtp_fig = render_drug_target_pathway_plot(
                    drug_name.strip() if mode == "names" else "", genes, enrichment_df)
            except Exception:
                logging.exception("Drug-target-pathway network figure failed — skipped")

        gsea_results: dict = {}
        gsea_drug_df = gsea_disease_df = None
        if fut_gsea is not None:
            gsea_results, gsea_drug_df, gsea_disease_df = fut_gsea.result()

        proximity_result: dict | None = None
        rwr_results: dict = {}
        if fut_net is not None:
            proximity_result, rwr_results = fut_net.result()

        modules_data = fut_mod.result()

    # Step 5/5.5: 全部结果保存到下载目录
    saved, modules_result = save_analysis_outputs(
        graph=graph, static_fig=static_fig, metrics_df=metrics_df, edges_df=edges_df,
        html_content=html_content, prefix=prefix,
        enrichment_df=enrichment_df, enrichment_bar_fig=enrichment_bar_fig,
        enrichment_dot_fig=enrichment_dot_fig, core_fig=core_fig, dtp_fig=dtp_fig,
        gsea_drug_df=gsea_drug_df, gsea_disease_df=gsea_disease_df,
        string_fig=string_fig, cyto_fig=cyto_fig,
        expansion_df=expansion_df, rwr_results=rwr_results,
        proximity_result=proximity_result, modules_data=modules_data,
        drug_name=drug_name, disease_name=disease_name,
        progress=progress)

    return {
        "genes": genes,
        "edges_df": edges_df,
        "graph": graph,
        "metrics_df": metrics_df,
        "hub_df": hub_df,
        "hub_rule": hub_rule,
        "hub_method_eff": hub_method_eff,
        "expansion_df": expansion_df,
        "html_content": html_content,
        "static_fig": static_fig,
        "string_fig": string_fig,
        "cyto_fig": cyto_fig,
        "mcc_map": mcc_map,
        "core_fig": core_fig,
        "core_genes": core_genes,
        "enrichment_df": enrichment_df,
        "enrichment_source": enrichment_source,
        "enrichment_reason": enrichment_reason,
        "enrichment_bar_fig": enrichment_bar_fig,
        "enrichment_dot_fig": enrichment_dot_fig,
        "dtp_fig": dtp_fig,
        "gsea": gsea_results,
        "proximity": proximity_result,
        "rwr": rwr_results,
        "modules": modules_result,
        "saved": saved,
        "prefix": prefix,
        "params": (mode, drug_name, disease_name, gene_text, score,
                   network_type, int(rwr_expand), int(hub_top_n), hub_method),
    }
