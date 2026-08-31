"""网络算法模块（纯计算）：Louvain 模块检测 + 模块级批量富集。

职责边界：
- 只接收 networkx.Graph / list[dict] 等数据入参，不 import streamlit、
  不自行发起网络请求（富集委托 core.enrichment，由其走既有缓存与熔断逻辑）；
- 可视化渲染仍归 core.ppi，本模块只输出数据。

失败语义（调用方需遵守）：
- 单模块富集失败时返回空表并记录失败原因；批量流程完成后会主动删除
  失败模块可能产生的缓存文件，保证"失败/异常产生的空表不会被后续命中"，
  用户看到的空表永远只可能是"真无显著项"。
"""

import logging
import os

import pandas as pd

from .enrichment import (
    ENRICHR_LIB_MAP,
    _ENRICHR_ORGANISM_MAP,
    _enrichment_cache_key,
    fetch_enrichment_enrichr,
    fetch_enrichment_gprofiler,
)

# fetch_enrichment_* 的 status.reason 中代表"调用成功"的取值集合。
# 注意 "circuit_open"（被熔断跳过，未真正请求）与 "request_failed" 均视为不可信。
_OK_REASONS = ("success", "no_significant", "cache")


def detect_louvain_modules(graph: "nx.Graph", resolution: float = 1.0, seed: int = 0) -> list[dict]:
    """
    用 Louvain 算法把 PPI 网络划分为功能模块。

    Parameters
    ----------
    graph : networkx.Graph
        PPI 网络（可含孤立节点，孤立节点各自成单点模块）。
    resolution : float
        Louvain 分辨率，越大模块越细碎；退化（巨模块）时建议调高。
    seed : int
        随机种子。固定种子保证同一网络两次检测结果完全一致。

    Returns
    -------
    list[dict]
        每个模块一个字典：module_id（按规模降序重排为 1..n）、genes（排序后的
        基因列表）、size、avg_degree。排序规则：size 降序，同规模按基因名字典序，
        因此输出完全确定。
    """
    from networkx.algorithms.community import louvain_communities

    if graph is None or graph.number_of_nodes() == 0:
        return []
    communities = louvain_communities(graph, resolution=resolution, seed=seed)
    # 确定性排序：先按规模降序，同规模按（排序后）基因名字典序，
    # 再重排 module_id，保证同参数两次运行输出逐字段一致。
    comm_lists = sorted((sorted(c) for c in communities), key=lambda nodes: (-len(nodes), nodes))
    modules = []
    for i, nodes in enumerate(comm_lists, start=1):
        avg_degree = sum(graph.degree(n) for n in nodes) / len(nodes)
        modules.append({
            "module_id": i,
            "genes": nodes,
            "size": len(nodes),
            "avg_degree": round(float(avg_degree), 2),
        })
    logging.info("Louvain detected %d modules (resolution=%s, seed=%d, nodes=%d)",
                 len(modules), resolution, seed, graph.number_of_nodes())
    return modules


def _module_cache_files(genes: list[str], organism: str, sources: tuple,
                        pvalue: float, label: str) -> list[str]:
    """推算某模块富集可能产生的缓存文件路径（g:Profiler 与 Enrichr 两份）。

    命名规则与 core.enrichment 内部保持一致：label 只进文件名前缀、不进 md5 键，
    md5 由（排序去重基因 + 物种 + 数据源 + 阈值）决定。
    """
    import config

    genes_clean = sorted({g.strip() for g in genes if g and g.strip()})
    files: list[str] = []
    if not genes_clean:
        return files
    md5_gp = _enrichment_cache_key(genes_clean, organism, sources, pvalue)
    files.append(os.path.join(config.ENRICHMENT_CACHE_DIR, f"gProfiler_{label}_{md5_gp}.csv"))
    enrichr_org = _ENRICHR_ORGANISM_MAP.get(organism)
    if enrichr_org:
        libs = [ENRICHR_LIB_MAP[s] for s in sources if s in ENRICHR_LIB_MAP]
        if libs:
            md5_enr = _enrichment_cache_key(genes_clean, enrichr_org, libs, pvalue)
            files.append(os.path.join(config.ENRICHMENT_CACHE_DIR, f"Enrichr_{label}_{md5_enr}.csv"))
    return files


def _purge_module_cache(files: list[str], module_id: int, reason: str) -> None:
    """删除某模块富集的缓存文件（失败/异常场景的防污染清理，文件不存在时为无害空操作）。"""
    for f in files:
        try:
            if os.path.exists(f):
                os.remove(f)
                logging.warning("Removed suspect enrichment cache of module %d (%s): %s",
                                module_id, reason, os.path.basename(f))
        except Exception:
            logging.exception("Failed to remove suspect cache file %s", f)


def _enrich_single_module(genes: list[str], organism: str, sources: tuple,
                          pvalue: float, label: str) -> tuple[pd.DataFrame, str | None, bool, str]:
    """
    单模块富集：g:Profiler 优先 → Enrichr 降级，并区分"真无显著项"与"服务失败"。

    Returns
    -------
    tuple
        (富集结果表, 实际数据源, 是否失败, 失败/状态说明)。
        failed=True 时结果表必为空表，调用方应清理该模块缓存并向用户提示错误；
        failed=False 且结果表为空表示两条链路都成功但确实无显著条目。
    """
    st_gp: dict = {}
    df = fetch_enrichment_gprofiler(genes, organism=organism, sources=sources,
                                    pvalue=pvalue, label=label, status=st_gp)
    if not df.empty:
        return df, "gProfiler", False, ""

    enrichr_org = _ENRICHR_ORGANISM_MAP.get(organism)
    libs = [ENRICHR_LIB_MAP[s] for s in sources if s in ENRICHR_LIB_MAP] if enrichr_org else []
    st_enr: dict = {}
    if enrichr_org and libs:
        df = fetch_enrichment_enrichr(genes, libs=libs, organism=enrichr_org,
                                      pvalue=pvalue, label=label, status=st_enr)
        if not df.empty:
            return df, "Enrichr", False, ""

    # 两条链路都拿到空表：只有每条"已启用"链路的 status 都可信（请求成功但无显著项 /
    # 正常缓存）才判定为真无显著项；g:Profiler 被熔断跳过、请求失败、Enrichr 请求失败，
    # 任何一项不可信即整体视为失败（大鼠无 Enrichr 库时视为"可信缺省"）。
    gp_ok = bool(st_gp.get("ok")) and st_gp.get("reason") in _OK_REASONS
    enrichr_applicable = bool(enrichr_org and libs)
    enrichr_ok = (not enrichr_applicable) or (
        bool(st_enr.get("ok")) and st_enr.get("reason") in _OK_REASONS
    )
    if gp_ok and enrichr_ok:
        return pd.DataFrame(), None, False, "no_significant"
    detail = f"gProfiler:{st_gp.get('reason', 'unknown')}; Enrichr:{st_enr.get('reason', 'not_applicable')}"
    logging.warning("Module enrichment failed (%s) for %d genes", detail, len(genes))
    return pd.DataFrame(), None, True, detail


def enrich_modules(
    modules: list[dict],
    organism: str,
    min_module_size: int = 5,
    top_n_modules: int = 10,
    sources: tuple | None = None,
    pvalue: float = 0.05,
) -> tuple[pd.DataFrame, list[dict], list[int]]:
    """
    批量模块富集：对满足规模条件的模块逐个做 GO/KEGG 富集并汇总为长表。

    失败处理约定（防止用户使用错误信息）：
    - 某模块失败 → 该模块计入 failures（页面提示错误），且批量结束后主动删除
      该模块可能产生的缓存文件；这样后续重试不会命中失败期间的任何空表/半截表；
    - 某模块两链路均成功但无显著条目 → 计入 no_sig（页面作一般性提示，不报错）。

    Parameters
    ----------
    modules : list[dict]
        detect_louvain_modules 的输出（需已按 size 降序）。
    organism : str
        g:Profiler 物种名（hsapiens / mmusculus / rnorvegicus）。
    min_module_size : int
        参与富集的模块规模下限（基因过少统计功效不足）。
    top_n_modules : int
        最多富集的模块数量（取规模最大的前 N 个）。
    sources : tuple | None
        g:Profiler 数据源，默认取 config.ENRICHMENT_SOURCES。
    pvalue : float
        显著性阈值。

    Returns
    -------
    tuple[pd.DataFrame, list[dict], list[int]]
        (富集长表: 统一列 + module_id + module_size, 失败列表, 无显著项模块ID列表)。
        长表按 module_id 升序、p_value 升序排列，便于页面直接 groupby.head 展示。
    """
    import config

    sources = sources or config.ENRICHMENT_SOURCES
    candidates = [m for m in modules if m["size"] >= min_module_size][:top_n_modules]

    frames: list[pd.DataFrame] = []
    failures: list[dict] = []
    no_sig: list[int] = []

    for m in candidates:
        label = f"module{m['module_id']}"
        cache_files = _module_cache_files(m["genes"], organism, sources, pvalue, label)
        try:
            df, used_source, failed, err = _enrich_single_module(
                m["genes"], organism, sources, pvalue, label)
        except Exception:
            logging.exception("Module %d enrichment crashed", m["module_id"])
            _purge_module_cache(cache_files, m["module_id"], "exception")
            failures.append({"module_id": m["module_id"], "size": m["size"], "error": "exception"})
            continue

        if failed:
            # 失败：空表 + 记录 + 清缓存（防止失败期间的空表被 TTL 内的后续运行命中）
            _purge_module_cache(cache_files, m["module_id"], err)
            failures.append({"module_id": m["module_id"], "size": m["size"], "error": err})
        elif df.empty:
            # 真无显著项：正常实现不会写空缓存，这里防御性清理一次以绝后患
            _purge_module_cache(cache_files, m["module_id"], "no_significant")
            no_sig.append(m["module_id"])
        else:
            out = df.copy()
            out.insert(0, "module_id", m["module_id"])
            out.insert(1, "module_size", m["size"])
            frames.append(out)

    if frames:
        long_df = pd.concat(frames, ignore_index=True)
        long_df = long_df.sort_values(["module_id", "p_value"]).reset_index(drop=True)
    else:
        long_df = pd.DataFrame()

    logging.info("Module enrichment done: %d modules, %d succeeded, %d failed, %d no-significant",
                 len(candidates), len(frames), len(failures), len(no_sig))
    return long_df, failures, no_sig
