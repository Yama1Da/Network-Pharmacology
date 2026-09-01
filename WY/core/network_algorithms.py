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

import functools
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


# ---------- 网络邻近度检验（M3，依赖 core.interactome 的背景互作组） ----------
# 方法学：Guney et al. 2016（Nat Commun）的 proximity 框架 / 论文 2.1.4 节。
# d_c(A,B) = mean_{b∈B} min_{a∈A} d(a,b)；Z = (d_obs - μ_rand) / σ_rand，
# Z<0 表示药物靶点与疾病蛋白的网络距离显著近于度匹配随机期望。
# 实现要点：
# - 全部计算限制在背景互作组的最大连通分量（LCC）内——跨分量距离无定义，
#   这是 Guney 框架的标准处理，同时天然规避不可达特判；
# - 置换随机化用 scipy 稀疏矩阵"级别同步 BFS"（每次毫秒级），networkx 纯 Python
#   参考实现 network_proximity() 保留给小图与自检交叉验证。


def _multi_source_bfs(Mt, source_idx) -> "np.ndarray":
    """稀疏邻接（已转置 CSR）上的多源 BFS 距离；节点全在 LCC 内，必然可达（>=0）。"""
    import numpy as np

    n = Mt.shape[0]
    dist = np.full(n, -1, dtype=np.int64)
    reached = np.zeros(n, dtype=bool)
    frontier = np.zeros(n, dtype=bool)
    source_idx = np.asarray(source_idx, dtype=int)
    frontier[source_idx] = True
    dist[source_idx] = 0  # 源节点自身距离为 0——缺失会让 A∩B 的节点以 -1 混入均值
    reached |= frontier
    depth = 0
    while True:
        depth += 1
        nxt = Mt.dot(frontier) > 0  # 布尔向量过稀疏矩阵 → 一步可达集合
        nxt &= ~reached
        if not nxt.any():
            break
        dist[nxt] = depth
        reached |= nxt
        frontier = nxt
    return dist


@functools.lru_cache(maxsize=2)
def _lcc_context(species: str):
    """背景互作组的 LCC 计算上下文（进程内缓存）。

    Returns
    -------
    tuple
        (nodes: 排序后基因元组, idx: 基因→行号, Mt: 转置 0-1 邻接 CSR,
         degree_map: 基因→LCC内度数)
    """
    import networkx as nx

    from . import interactome as itx

    graph = itx.load_interactome(species)
    lcc_nodes = max(nx.connected_components(graph), key=len)
    sub = graph.subgraph(lcc_nodes)
    nodes = tuple(sorted(sub.nodes()))
    idx = {g: i for i, g in enumerate(nodes)}
    adj = nx.to_scipy_sparse_array(sub, nodelist=nodes, weight=None, format="csr")
    Mt = adj.T.tocsr()
    degree_map = {g: int(d) for g, d in sub.degree()}
    logging.info("Proximity LCC context ready: %d/%d nodes in LCC",
                 len(nodes), graph.number_of_nodes())
    return nodes, idx, Mt, degree_map


def _build_degree_bins(degree_map: dict, bin_size: int) -> dict:
    """按度值区间 [k*bin, (k+1)*bin) 分箱（论文 Algorithm 1 的 bins）。"""
    bins: dict[int, list[str]] = {}
    for gene, deg in degree_map.items():
        bins.setdefault(int(deg) // bin_size, []).append(gene)
    for k in bins:
        bins[k].sort()
    return bins


def degree_matched_sample(degree_map: dict, seeds: list, bin_size: int, rng) -> list[str]:
    """
    度匹配随机抽样（论文 Algorithm 1）：每个种子从其度数所在箱内随机取一个替换节点。

    - 排除种子自身；箱内无候选时按 |度差| 向最近非空箱扩散；
    - 不同种子可能取到同一替换节点（与参考实现一致，随机集内部去重不影响
      距离计算——集合语义）；
    - 同一 rng 序列下结果确定。
    """
    bins = _build_degree_bins(degree_map, bin_size)
    max_k = max(bins) if bins else 0
    seeds_u = sorted({str(s) for s in seeds if s in degree_map})
    picked: list[str] = []
    for seed_gene in seeds_u:
        k = int(degree_map[seed_gene]) // bin_size
        chosen = None
        for offset in range(0, max_k + 2):
            for k_try in ((k - offset, k + offset) if offset else (k,)):
                pool = [g for g in bins.get(k_try, ()) if g != seed_gene]
                if pool:
                    chosen = pool[int(rng.integers(len(pool)))]
                    break
            if chosen is not None:
                break
        if chosen is not None:
            picked.append(chosen)
    return picked


def network_proximity(graph: "nx.Graph", set_a: list, set_b: list) -> float:
    """
    参考实现（networkx 级别同步多源 BFS）：d_c(A,B) = mean_{b∈B} min_{a∈A} d(a,b)。

    仅供小图与自检交叉验证；页面链路走 proximity_zscore 的稀疏快路径。
    图中不存在的基因被忽略；B 中存在从 A 不可达的节点时抛 ValueError。
    """
    a = {str(g) for g in set_a if g in graph}
    b = {str(g) for g in set_b if g in graph}
    if not a or not b:
        raise ValueError("network_proximity: empty set after filtering to graph")
    dist: dict[str, int] = {g: 0 for g in a}
    frontier = set(a)
    depth = 0
    while frontier:
        depth += 1
        nxt: set[str] = set()
        for u in frontier:
            for v in graph[u]:
                if v not in dist:
                    dist[v] = depth
                    nxt.add(v)
        frontier = nxt
    unreachable = b - dist.keys()
    if unreachable:
        raise ValueError(f"network_proximity: {len(unreachable)} nodes in B unreachable from A")
    return sum(dist[g] for g in b) / len(b)


def _proximity_zscore_from_context(ctx: tuple, drug_genes: list, disease_genes: list,
                                   n_perm: int, seed: int, bin_size: int,
                                   drug_ratio: float, disease_ratio: float) -> dict:
    """给定 LCC 上下文与已映射的基因集，计算 d_c / Z / 经验 p（纯计算，便于自检注入）。"""
    import numpy as np

    nodes, idx, Mt, degree_map = ctx
    result: dict = {
        "ok": False, "error": None, "d_c": None, "z": None, "p_empirical": None,
        "n_perm": int(n_perm), "seed": int(seed),
        "mean_random": None, "std_random": None,
        "lcc_nodes": len(nodes), "cached": False,
    }
    d_lcc = sorted({g for g in drug_genes if g in idx})
    s_lcc = sorted({g for g in disease_genes if g in idx})
    result["coverage"] = {
        "drug_in_lcc": len(d_lcc), "drug_total": len(drug_genes),
        "disease_in_lcc": len(s_lcc), "disease_total": len(disease_genes),
    }
    if not d_lcc or not s_lcc:
        result["error"] = "empty_after_filter"
        return result
    if len(d_lcc) / max(len(drug_genes), 1) < drug_ratio:
        result["error"] = "low_coverage_drug"
        return result
    if len(s_lcc) / max(len(disease_genes), 1) < disease_ratio:
        result["error"] = "low_coverage_disease"
        return result

    rng = np.random.default_rng(seed)
    a_idx = np.array([idx[g] for g in d_lcc])
    b_idx = np.array([idx[g] for g in s_lcc])

    dist_obs = _multi_source_bfs(Mt, a_idx)
    d_obs = float(dist_obs[b_idx].mean())

    rand = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        ra = degree_matched_sample(degree_map, d_lcc, bin_size, rng)
        rb = degree_matched_sample(degree_map, s_lcc, bin_size, rng)
        ra_idx = np.array([idx[g] for g in ra])
        rb_idx = np.array([idx[g] for g in rb])
        rand[i] = float(_multi_source_bfs(Mt, ra_idx)[rb_idx].mean())

    mu = float(rand.mean())
    sd = float(rand.std())
    z = (d_obs - mu) / sd if sd > 0 else 0.0
    p_emp = (1 + int((rand <= d_obs).sum())) / (n_perm + 1)
    result.update(ok=True, d_c=d_obs, z=float(z), p_empirical=float(p_emp),
                  mean_random=mu, std_random=sd)
    return result


def proximity_zscore(drug_genes: list, disease_genes: list, species: str = "Human",
                     n_perm: int | None = None, seed: int | None = None) -> dict:
    """
    网络邻近度 Z-score 主入口（背景互作组 LCC 内计算，度匹配随机化置换检验）。

    覆盖率守卫：基因先经全网络口径归一化，再限制到 LCC；任一集合的入网率低于
    config.PROXIMITY_MIN_COVERAGE 或过滤后为空时，返回 ok=False + error，
    绝不输出"带病结果"。

    Raises
    ------
    core.interactome.InteractomeNotBuiltError
        背景互作组缓存未构建时上抛（调用方翻译为用户提示）。
    """
    import config
    from . import interactome as itx

    n_perm = int(n_perm if n_perm is not None else config.PROXIMITY_N_PERM)
    seed = int(seed if seed is not None else config.PROXIMITY_SEED)
    ctx = _lcc_context(species)
    d_in, _ = itx.get_network_coverage(drug_genes, species)
    s_in, _ = itx.get_network_coverage(disease_genes, species)
    return _proximity_zscore_from_context(
        ctx, d_in, s_in, n_perm=n_perm, seed=seed,
        bin_size=config.PROXIMITY_DEGREE_BIN_SIZE,
        drug_ratio=config.PROXIMITY_MIN_COVERAGE,
        disease_ratio=config.PROXIMITY_MIN_COVERAGE,
    )


def run_proximity_cached(drug_genes: list, disease_genes: list,
                         species: str = "Human", label: str | None = None) -> dict:
    """
    邻近度分析入口（带文件缓存）：成功才落盘 JSON，失败/未构建即时返回不缓存。

    缓存键 = 排序去重后的两个基因集 + 物种 + 置换参数 + 背景网络版本/阈值，
    任一项变化自然换文件；命中 TTL 内缓存时返回结果并置 cached=True。
    """
    import hashlib
    import json
    import os

    import config

    from .caching import _is_cache_fresh
    from . import interactome as itx

    def _norm(genes: list) -> list[str]:
        return sorted({str(g).strip().upper() for g in genes if g and str(g).strip()})

    payload = ("\r".join(_norm(drug_genes)) + "||" + "\r".join(_norm(disease_genes))
               + f"|{species}|{config.PROXIMITY_N_PERM}|{config.PROXIMITY_SEED}"
               + f"|{config.PROXIMITY_DEGREE_BIN_SIZE}"
               + f"|{config.STRING_VERSION}|{config.INTERACTOME_CONFIDENCE}")
    md5 = hashlib.md5(payload.encode("utf-8")).hexdigest()
    name = f"proximity_{label}_{md5}.json" if label else f"proximity_{md5}.json"
    os.makedirs(config.PROXIMITY_CACHE_DIR, exist_ok=True)
    cache_path = os.path.join(config.PROXIMITY_CACHE_DIR, name)

    if _is_cache_fresh(cache_path):
        try:
            with open(cache_path, "r", encoding="utf-8") as f:
                cached = json.load(f)
            cached["cached"] = True
            logging.info("Proximity cache hit: %s", name)
            return cached
        except Exception:
            # 损坏缓存按"不存在"处理并清除（同 M1 的空缓存清理语义）
            logging.exception("Proximity cache unreadable, removing: %s", cache_path)
            try:
                os.remove(cache_path)
            except OSError:
                pass

    try:
        result = proximity_zscore(drug_genes, disease_genes, species=species)
    except itx.InteractomeNotBuiltError:
        logging.warning("Proximity skipped: interactome not built")
        return {"ok": False, "error": "not_built", "d_c": None, "z": None,
                "p_empirical": None, "coverage": {}, "cached": False}
    if result.get("ok"):
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False)
    return result
