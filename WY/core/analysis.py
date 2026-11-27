"""Build Report 用例编排：一次"PPI 网络药理学分析"的完整流程、并发任务与结果落盘。

模块与依赖方向（上层引下层，禁止反向）：
    analysis → targets / ppi / network_algorithms / enrichment
             → config（路径与参数常量；按 core 惯例函数内延迟导入）

界面层（Vue SPA 经 WY/api/、CLI 子命令、tools 脚本）只做三件事：收集参数 → 调 run_ppi_analysis →
展示/下发结果。进度经 AnalysisProgress 协议回调，编排层不依赖任何 UI 框架；
matplotlib 的 pyplot 全局状态非线程安全，所有绘图统一留在调用线程执行。

与 cli.py 的关系：网页端 Build Report 与 CLI 共享同一份分析语义,共享点集中在
本模块——pair_prefix（缓存命名）、expand_targets_with_rwr（RWR 扩展）、
build_ppi_stack（抓边→建图→指标→MCC→hub→核心靶点）、run_enrichment_step（富集）;
CLI 各步骤只做产物落盘与状态字典,网页侧只做并发调度与绘图。GSEA/模块分析/图形
渲染是网页专属产物形态,不进 CLI;邻近度/RWR 双方直接调 core.network_algorithms
的同名函数,天然共享。改任何计算规则只改本模块,两侧自动一致。
"""

import hashlib
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
    pick_top_mcc_genes,
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
    现由 api/ppi.py 的适配器实现(jobs 进度回调)。"""

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


def build_ppi_stack(genes: list[str], *, score: float, label: str | None,
                    hub_top_n: int, hub_method: str) -> dict:
    """PPI 数据栈（纯计算，不落盘不画图；CLI _ppi_step 与网页 run_ppi_analysis 共用）：
    STRING_QUERY_LIMIT 截断 → STRING 抓边 → 建图 → 拓扑指标 → MCC → hub 标记 → 核心靶点。

    截断规则与 fetch_string_ppi 内部一致，此处先行截断并如实返回 truncated 标记与
    截断后的 genes（网页手动模式原先静默依赖 fetch 内部截断，收敛后两端标记一致）。
    MCC 无条件尝试（核心靶点无论 hub_method 取哪种都按 MCC 挑选；失败降级为空映射继续）；
    MCC 成功且 metrics_df 尚无 mcc 列时补一列（fillna 0.0，保证指标表 CSV 完整——
    网页侧原先在 build_main_figures 里补列且零值落成 NaN，收敛后统一为 0.0）。

    返回 dict：genes（截断后）/ truncated / edges_df（fetch 失败归一化为空表，可能为空）/
    graph / metrics_df / mcc_map / hub_rule / hub_method_eff / core_genes。
    hub_method_eff：请求 "mcc" 但 mark_hub_genes 实际退回 degree 时为 "degree"（与网页
    原判定一致），否则等于 hub_method。边表为空时 graph/metrics 仍会构建（全孤立节点图），
    由调用方按需分支（CLI 判失败，网页沿用全孤立节点流程）。
    """
    from WY import config

    truncated = len(genes) > config.STRING_QUERY_LIMIT
    if truncated:
        genes = genes[:config.STRING_QUERY_LIMIT]
    edges_df = fetch_string_ppi(
        genes, species=config.SPECIES_MAP["Human"]["ncbi_taxon"],
        score=score, network_type="full", label=label)
    if edges_df is None:
        # fetch_string_ppi 失败可返回 None:归一化为空表,调用方统一按 .empty 分支
        edges_df = pd.DataFrame()
    graph = build_ppi_graph(edges_df, genes)
    metrics_df = compute_ppi_metrics(graph)
    try:
        mcc_map: dict[str, float] = compute_mcc_scores(graph)
    except Exception:
        logging.exception("MCC computation failed; core targets skipped")
        mcc_map = {}
    if mcc_map and "mcc" not in metrics_df.columns:
        # hub_method 为 degree 时 mark_hub_genes 不带 mcc 列,提前补上保证指标表完整
        metrics_df["mcc"] = metrics_df["gene"].map(mcc_map).fillna(0.0)
    metrics_df = mark_hub_genes(
        metrics_df, top_n=hub_top_n, method=hub_method,
        mcc_map=mcc_map, small_net_factor=config.HUB_SMALL_NET_FACTOR)
    hub_rule = str(metrics_df["hub_rule"].iloc[0]) if len(metrics_df) else ""
    # MCC 排名被 mark_hub_genes 回退为 degree 时,metrics_df 不会出现 mcc 列
    hub_method_eff = ("degree"
                      if hub_rule.startswith("top") and hub_method == "mcc"
                      and "mcc" not in metrics_df.columns else hub_method)
    core_genes = pick_top_mcc_genes(metrics_df, mcc_map, config.MCC_TOP_N)
    return {"genes": genes, "truncated": truncated, "edges_df": edges_df,
            "graph": graph, "metrics_df": metrics_df, "mcc_map": mcc_map,
            "hub_rule": hub_rule, "hub_method_eff": hub_method_eff,
            "core_genes": core_genes}


def run_enrichment_step(genes: list[str], *, organism: str | None = None,
                        pvalue: float | None = None, sources: tuple | None = None,
                        label: str | None = None) -> tuple:
    """功能富集计算（纯计算，不落盘；CLI _enrichment_step 与网页 run_enrichment_task 共用）：
    g:Profiler 主 + Enrichr 自动降级。organism/pvalue 缺省取 config 默认值，
    sources 缺省 None（走 run_enrichment 内部默认）。返回 (结果表, 数据源, 状态原因)。"""
    from WY import config

    enr_status: dict = {}
    df, source = run_enrichment(
        genes, organism=organism or config.SPECIES_MAP["Human"]["gprofiler"],
        pvalue=pvalue if pvalue is not None else config.ENRICHMENT_PVALUE,
        sources=sources, label=label, status=enr_status)
    return df, source, enr_status.get("reason", "")


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
    """功能富集任务(网页 Build Report 的线程池并发入口),返回 (结果表, 数据源, 状态原因)。
    计算统一走 run_enrichment_step(与 CLI _enrichment_step 共用同一份)。"""
    return run_enrichment_step(genes, organism=organism, label=prefix)


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
                          genes: list | None = None,
                          progress=None) -> tuple[dict, dict | None]:
    """Step 5/5.5:全部结果保存到下载目录(同名分析旧目录会被替换,不会堆积)。

    模块文件必须在主目录重建后落盘,否则会被 save_ppi_results 的旧目录清理删掉。
    drug_name/disease_name 写入 analysis_meta.json,供 ML 筛选页按元数据发现排序 CSV;
    genes 同样入 meta(手动基因模式的哈希桶名靠它提供人类可读的清单)。
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
        genes=genes,
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
        # 手动模式按基因列表指纹分桶(与 fetch_string_ppi 的基因指纹同款 md5):
        # 不同基因列表各自成目录互不覆盖,同一列表幂等复用同名目录与富集缓存;
        # 边表取数仍传 "manual"(见下方 fetch_string_ppi),与 CLI 共享同一特判与缓存。
        digest = hashlib.md5("\r".join(genes).encode("utf-8")).hexdigest()[:8]
        prefix = f"manual_{digest}"
        # 表单里残留的药病名与本分析无关,置空:不进 meta(ML 页按药病名发现排序
        # CSV,残留值会让手动桶冒充某个药病对),params 快照也置空(前端 genes 模式
        # 的 snapshot 同步传空,两侧一致才不会误报参数过期)
        drug_name = disease_name = ""

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
    # Step 1-4: STRING 抓边 → 建图 → 拓扑指标 → MCC → hub 标记 → 核心靶点,
    # 统一走 build_ppi_stack(与 CLI _ppi_step 共用同一份计算)。
    # 手动模式固定传 "manual":命中 fetch_string_ppi 的 label=="manual" 特判
    # (文件名附加基因指纹),与 CLI run_ppi --genes 共享同一份边表缓存。
    stack = build_ppi_stack(
        genes, score=score, label="manual" if mode != "names" else prefix,
        hub_top_n=hub_top_n, hub_method=hub_method)
    genes = stack["genes"]  # STRING_QUERY_LIMIT 截断后(names 模式不会截,手动模式可能截)
    edges_df = stack["edges_df"]
    graph = stack["graph"]
    metrics_df = stack["metrics_df"]
    mcc_map = stack["mcc_map"]
    hub_rule = stack["hub_rule"]
    hub_method_eff = stack["hub_method_eff"]
    _prog("string", "done", 0.30)
    _prog("graph", "running", 0.30)
    hub_df = metrics_df[metrics_df["is_hub"]]
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

        # 轮询任务完成情况并刷新进度回调(回调在调用线程内触发,实现需自行保证线程安全)
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
        genes=genes,
        progress=progress)

    # ---- Step 5.6(可选,W7 招牌图,2026-11-27):暴露约束前后排序对比 ----
    # 仅 names 模式;数据源 = RWR 药侧种子得分(一期等权)× w(R) 档位权重
    # (w_unknown=0.5 为 1.7 校准前临时口径);渲染门槛 max(5, ceil(5%×N)) 见
    # 方案 §3.2。失败/低于门槛只跳过自身或回结论行摘要,不推翻主报告。
    reach_compare: dict | None = None
    if mode == "names" and drug_name.strip():
        try:
            from .reach_compare import save_reach_compare
            from .reachability import compute_reachability
            drug_rwr = (rwr_results or {}).get("drug_seeds") or {}
            if drug_rwr.get("ok") and drug_rwr.get("seeds"):
                seed_scores = {str(r["gene"]): float(r["score"])
                               for r in drug_rwr["seeds"]}
                # 逐靶点档位直连计算层(core→core,分层红线:core 不 import api);
                # 归一化缓存未就绪 → RuntimeError → 本步整体跳过(结论行也不出)
                res = compute_reachability(
                    get_drug_targets_merged(drug_name.strip()), drug_name.strip())
                tiers = {str(g): str(t) for g, t in
                         zip(res["targets"]["gene"], res["targets"]["tier"])}
                reach_compare = save_reach_compare(
                    saved["directory"], drug_name.strip(), seed_scores, tiers)
                png = reach_compare.pop("reach_compare_png", None)
                if png:
                    saved["reach_compare_png"] = png  # 进 artifacts → 前端/产物页
        except Exception:
            logging.exception("reach compare step failed — skipped gracefully")
            reach_compare = None

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
        "reach_compare": reach_compare,
        "prefix": prefix,
        "params": (mode, drug_name, disease_name, gene_text, score,
                   network_type, int(rwr_expand), int(hub_top_n), hub_method),
    }
