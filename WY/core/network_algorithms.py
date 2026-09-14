"""网络算法计算模块：Louvain 模块检测、模块级批量富集、网络邻近度 Z-score、RWR 靶点排序。

单个模块富集失败时返回空表并记录原因，且会清除失败可能产生的缓存文件；
因此调用方看到的空表只可能意味着"真无显著项"。
"""

import functools
import logging
import os

import pandas as pd

from .enrichment import (
    _ENRICHR_ORGANISM_MAP,
    _enrichment_cache_key,
    _enrichr_libs,
    enrichment_cache_filename,
    fetch_enrichment_enrichr,
    fetch_enrichment_gprofiler,
    TRUSTED_REASONS,
)

# fetch_enrichment_* 返回状态中代表"调用成功"的取值复用 core.enrichment.TRUSTED_REASONS；
# "circuit_open"（被熔断跳过）、"probe_failed"（冷却到期的探测失败）与 "request_failed" 均视为不可信。
_UNTRUSTED_REASONS = ("request_failed", "probe_failed", "circuit_open")


def detect_louvain_modules(graph: "nx.Graph", resolution: float = 1.0, seed: int = 0) -> list[dict]:
    """用 Louvain 算法按 PPI 连接紧密程度把网络自动分成若干群，每群作为一个功能模块。
    resolution 越大模块越细碎;固定 seed 保证同一网络两次结果完全一致;孤立节点各自成单点模块,module_id 按基因数从多到少编号 1..n。
    """
    from networkx.algorithms.community import louvain_communities

    if graph is None or graph.number_of_nodes() == 0:
        return []
    communities = louvain_communities(graph, resolution=resolution, seed=seed)
    # 模块按基因数从多到少排序编号（同规模按基因名字典序）；同一网络每次运行结果完全一致。
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
    """推算某模块富集可能产生的缓存文件路径（g:Profiler 与 Enrichr 两份），命名规则与 core.enrichment 内部一致。"""
    import config

    genes_clean = sorted({g.strip() for g in genes if g and g.strip()})
    files: list[str] = []
    if not genes_clean:
        return files
    md5_gp = _enrichment_cache_key(genes_clean, organism, sources, pvalue)
    files.append(os.path.join(config.ENRICHMENT_CACHE_DIR,
                              enrichment_cache_filename("gProfiler", label, md5_gp)))
    enrichr_org = _ENRICHR_ORGANISM_MAP.get(organism)
    if enrichr_org:
        libs = _enrichr_libs(sources)
        if libs:
            md5_enr = _enrichment_cache_key(genes_clean, enrichr_org, libs, pvalue)
            files.append(os.path.join(config.ENRICHMENT_CACHE_DIR,
                                      enrichment_cache_filename("Enrichr", label, md5_enr)))
    return files


def _purge_module_cache(files: list[str], module_id: int, reason: str) -> None:
    """删除某模块富集的缓存文件，防止失败产生的空表被后续运行误用；文件不存在时为无害空操作。"""
    for f in files:
        try:
            if os.path.exists(f):
                os.remove(f)
                logging.warning("Removed suspect enrichment cache of module %d (%s): %s",
                                module_id, reason, os.path.basename(f))
        except Exception:
            logging.exception("Failed to remove suspect cache file %s", f)


def _purge_orphan_module_caches(module_ids: set[int], label_prefix: str) -> None:
    """删除同分析标识下编号已不属于本次模块集的旧模块缓存(上次模块数更多时留下的孤儿文件)。
    缓存命中本就有基因集 md5 键校验,孤儿文件不会被误用,删除只为目录整洁;只匹配
    gProfiler_/Enrichr_ 前缀下 `<标识>_module<纯数字>.csv` 命名,不碰其他文件。"""
    import config

    cache_dir = config.ENRICHMENT_CACHE_DIR
    if not os.path.isdir(cache_dir):
        return
    stem = f"{label_prefix}_module" if label_prefix else "module"
    for source in ("gProfiler", "Enrichr"):
        head = f"{source}_{stem}"
        for name in os.listdir(cache_dir):
            if not name.startswith(head) or not name.endswith(".csv"):
                continue
            num = name[len(head):-len(".csv")]
            if not num.isdigit() or int(num) in module_ids:
                continue
            try:
                os.remove(os.path.join(cache_dir, name))
                logging.info("Orphan module cache removed: %s", name)
            except OSError:
                logging.exception("Failed to remove orphan module cache %s", name)


def _enrich_single_module(genes: list[str], organism: str, sources: tuple,
                          pvalue: float, label: str) -> tuple[pd.DataFrame, bool, str]:
    """单模块富集：先 g:Profiler，不行再降级用 Enrichr。返回 (结果表, 是否失败, 失败/状态说明)。
    failed=True 时结果表必为空表;failed=False 的空表表示两条链路都成功但确实无显著条目。
    """
    st_gp: dict = {}
    df = fetch_enrichment_gprofiler(genes, organism=organism, sources=sources,
                                    pvalue=pvalue, label=label, status=st_gp)
    if not df.empty:
        return df, False, ""

    enrichr_org = _ENRICHR_ORGANISM_MAP.get(organism)
    libs = _enrichr_libs(sources) if enrichr_org else []
    st_enr: dict = {}
    if enrichr_org and libs:
        df = fetch_enrichment_enrichr(genes, libs=libs, organism=enrichr_org,
                                      pvalue=pvalue, label=label, status=st_enr)
        if not df.empty:
            return df, False, ""

    # 两条链路都拿到空表时：只有每条已启用链路都真实成功（无显著项/正常缓存）才算真无显著项；
    # 被熔断跳过或请求失败一律按失败处理（大鼠无 Enrichr 库时视为可信缺省）。
    gp_ok = bool(st_gp.get("ok")) and st_gp.get("reason") in TRUSTED_REASONS
    enrichr_applicable = bool(enrichr_org and libs)
    enrichr_ok = (not enrichr_applicable) or (
        bool(st_enr.get("ok")) and st_enr.get("reason") in TRUSTED_REASONS
    )
    if gp_ok and enrichr_ok:
        return pd.DataFrame(), False, "no_significant"
    if (st_gp.get("reason") in _UNTRUSTED_REASONS
            and enrichr_applicable and st_enr.get("reason") in _UNTRUSTED_REASONS):
        detail = "both_sources_down"  # 双源都不可用：页面据此显著告警,提示稍后重试
    else:
        detail = f"gProfiler:{st_gp.get('reason', 'unknown')}; Enrichr:{st_enr.get('reason', 'not_applicable')}"
    logging.warning("Module enrichment failed (%s) for %d genes", detail, len(genes))
    return pd.DataFrame(), True, detail


def enrich_modules(
    modules: list[dict],
    organism: str,
    min_module_size: int = 5,
    top_n_modules: int = 10,
    sources: tuple | None = None,
    pvalue: float = 0.05,
    label_prefix: str = "",
) -> tuple[pd.DataFrame, list[dict], list[int]]:
    """批量模块富集：对满足规模条件的模块（取规模最大的前 top_n_modules 个）逐个做 GO/KEGG 富集并汇总为长表。
    label_prefix 为所属分析的"药物_疾病"标识,拼进缓存文件名让模块结果能认出归属;某模块失败时计入 failures 并清除
    该模块可能产生的缓存文件（保证空表不是失败残留）,两链路成功但无显著条目的模块计入 no_sig。"""
    import config

    sources = sources or config.ENRICHMENT_SOURCES
    candidates = [m for m in modules if m["size"] >= min_module_size][:top_n_modules]
    # 先清掉上次分析遗留的孤儿模块缓存(编号超出本次模块集),目录不随多次分析无限增长
    _purge_orphan_module_caches({m["module_id"] for m in candidates}, label_prefix)

    def _enrich_one(m: dict) -> tuple[str, dict]:
        """单个模块的富集 + 失败处理,返回 (类别, 数据):ok/failed/no_sig。
        线程安全:各模块的缓存文件名互不相同,状态都记在局部变量里。
        """
        label = f"{label_prefix}_module{m['module_id']}" if label_prefix else f"module{m['module_id']}"
        cache_files = _module_cache_files(m["genes"], organism, sources, pvalue, label)
        try:
            df, failed, err = _enrich_single_module(
                m["genes"], organism, sources, pvalue, label)
        except Exception:
            logging.exception("Module %d enrichment crashed", m["module_id"])
            _purge_module_cache(cache_files, m["module_id"], "exception")
            return "failed", {"module_id": m["module_id"], "size": m["size"], "error": "exception"}

        if failed:
            # 失败：记录并清缓存，防止失败期间的空表被后续运行命中
            _purge_module_cache(cache_files, m["module_id"], err)
            return "failed", {"module_id": m["module_id"], "size": m["size"], "error": err}
        if df.empty:
            # 真无显著项：正常不会写空缓存，这里防御性清理一次
            _purge_module_cache(cache_files, m["module_id"], "no_significant")
            return "no_sig", {"module_id": m["module_id"]}
        out = df.copy()
        out.insert(0, "module_id", m["module_id"])
        out.insert(1, "module_size", m["size"])
        return "ok", {"df": out}

    # 并发富集:模块间是独立的 API 请求,串行时 10 个模块可能要等几分钟;
    # 限制并发数,避免同时打到 g:Profiler/Enrichr 的请求过多。
    if candidates:
        from concurrent.futures import ThreadPoolExecutor

        max_workers = min(4, len(candidates))
        if max_workers > 1:
            with ThreadPoolExecutor(max_workers=max_workers,
                                    thread_name_prefix="modenrich") as pool:
                outcomes = list(pool.map(_enrich_one, candidates))
        else:
            outcomes = [_enrich_one(m) for m in candidates]
    else:
        outcomes = []

    frames: list[pd.DataFrame] = []
    failures: list[dict] = []
    no_sig: list[int] = []
    for (bucket, data) in outcomes:  # pool.map 保持原模块顺序
        if bucket == "ok":
            frames.append(data["df"])
        elif bucket == "failed":
            failures.append(data)
        else:
            no_sig.append(data["module_id"])

    if frames:
        long_df = pd.concat(frames, ignore_index=True)
        long_df = long_df.sort_values(["module_id", "p_value"]).reset_index(drop=True)
    else:
        long_df = pd.DataFrame()

    logging.info("Module enrichment done: %d modules, %d succeeded, %d failed, %d no-significant",
                 len(candidates), len(frames), len(failures), len(no_sig))
    return long_df, failures, no_sig


# ---------- 网络邻近度检验（M3，依赖 core.interactome 的背景互作组） ----------
# 方法学：Guney et al. 2016（Nat Commun）的 proximity 框架。
# d_c(A,B) = mean_{b∈B} min_{a∈A} d(a,b)；Z = (d_obs - μ_rand) / σ_rand，
# Z<0 表示药物靶点与疾病蛋白的网络距离显著近于度匹配随机期望。
# 全部计算限制在背景互作组最大连通分量（LCC）内——跨分量距离无定义，标准处理。


def _multi_source_bfs(Mt, source_idx) -> "np.ndarray":
    """稀疏邻接矩阵上的多源 BFS，返回各节点距离；节点全在 LCC 内，必然可达（>=0）。"""
    import numpy as np

    n = Mt.shape[0]
    dist = np.full(n, -1, dtype=np.int64)
    reached = np.zeros(n, dtype=bool)
    frontier = np.zeros(n, dtype=bool)
    source_idx = np.asarray(source_idx, dtype=int)
    frontier[source_idx] = True
    dist[source_idx] = 0  # 源节点自身距离置 0，否则 A∩B 的节点会以 -1 混入均值
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
    """背景互作组 LCC 的计算上下文（进程内缓存）：节点表、基因→行号、邻接矩阵、度数表。
    缓存键含磁盘构建记录的置信度与 STRING 版本:重建为不同参数后不会命中旧图;
    同参数重建则由 interactome.clear_runtime_caches() 一并清除本缓存。"""
    from . import interactome as itx

    meta = itx.require_meta(species)  # 未构建时由此抛 InteractomeNotBuiltError
    return _lcc_context_cached(species, int(meta["confidence"]), str(meta["string_version"]))


@functools.lru_cache(maxsize=2)
def _lcc_context_cached(species: str, confidence: int, string_version: str):
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
    """按度值区间 [k*bin, (k+1)*bin) 给基因分箱（度匹配抽样用）。"""
    bins: dict[int, list[str]] = {}
    for gene, deg in degree_map.items():
        bins.setdefault(int(deg) // bin_size, []).append(gene)
    for k in bins:
        bins[k].sort()
    return bins


def degree_matched_sample(degree_map: dict, seeds: list, bin_size: int, rng,
                          bins: dict | None = None) -> list[str]:
    """度匹配随机抽样（置换检验的零分布）：为每个种子从其度数箱内随机取一个替换节点,使随机集与真实集
    度数分布相当;排除种子自身,箱内无候选时向最近非空箱扩散,同一随机序列下结果确定。bins 可传预构建分箱（循环内复用）。
    性能:低度箱可达上万成员,按箱建一次"基因→位置"映射后用索引跳过种子,等价于过滤再抽样但避免每种子复制整箱（随机数消耗一致）。"""
    if bins is None:
        bins = _build_degree_bins(degree_map, bin_size)
    max_k = max(bins) if bins else 0
    seeds_u = sorted({str(s) for s in seeds if s in degree_map})
    pos_maps: dict[int, dict[str, int]] = {}  # 箱号 → {基因: 在箱内的位置},按需构建一次
    picked: list[str] = []
    for seed_gene in seeds_u:
        k = int(degree_map[seed_gene]) // bin_size
        chosen = None
        for offset in range(0, max_k + 2):
            for k_try in ((k - offset, k + offset) if offset else (k,)):
                pool = bins.get(k_try)
                if not pool:
                    continue
                pos = pos_maps.get(k_try)
                if pos is None:
                    pos = pos_maps[k_try] = {g: i for i, g in enumerate(pool)}
                if seed_gene in pos:
                    if len(pool) == 1:
                        continue  # 箱内只有种子自己,视为无候选,转向邻箱
                    # 从"去掉种子的箱"里均匀取一个:索引 j 跳过种子的位置
                    j = int(rng.integers(len(pool) - 1))
                    chosen = pool[j + 1] if j >= pos[seed_gene] else pool[j]
                else:
                    chosen = pool[int(rng.integers(len(pool)))]
                break
            if chosen is not None:
                break
        if chosen is not None:
            picked.append(chosen)
    return picked


def _proximity_zscore_from_context(ctx: tuple, drug_genes: list, disease_genes: list,
                                   n_perm: int, seed: int, bin_size: int,
                                   drug_ratio: float, disease_ratio: float) -> dict:
    """给定 LCC 上下文与已映射的基因集，计算 d_c / Z / 经验 p（纯计算，便于自检注入）。"""
    import numpy as np

    nodes, idx, Mt, degree_map = ctx
    result: dict = {
        "ok": False, "error": None, "d_c": None, "z": None, "p_empirical": None,
        "d_c_excl": None, "z_excl": None, "p_excl": None, "n_excl": None,
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

    # ---- 去交集版本：把同时出现在药物靶点里的疾病基因从 B 中剔除后再算 d_c ----
    # 共享靶点距离为 0 会拉低 d_c；剔除后若 Z 仍 < 0，说明邻近性并非只由
    # 药物直接作用于疾病基因贡献。
    a_set = set(d_lcc)
    b_only = [g for g in s_lcc if g not in a_set]
    result["n_excl"] = len(b_only)
    d_obs_excl = float(dist_obs[[idx[g] for g in b_only]].mean()) if b_only else None
    rand_excl = np.full(n_perm, np.nan)

    rand = np.empty(n_perm, dtype=float)
    bins = _build_degree_bins(degree_map, bin_size)  # 置换循环内复用,避免每轮重建分箱
    for i in range(n_perm):
        ra = degree_matched_sample(degree_map, d_lcc, bin_size, rng, bins=bins)
        rb = degree_matched_sample(degree_map, s_lcc, bin_size, rng, bins=bins)
        ra_idx = np.array([idx[g] for g in ra])
        rb_idx = np.array([idx[g] for g in rb])
        dist_r = _multi_source_bfs(Mt, ra_idx)
        rand[i] = float(dist_r[rb_idx].mean())
        if b_only:
            ra_set = set(ra)
            rb_only_idx = [idx[g] for g in rb if g not in ra_set]
            if rb_only_idx:
                rand_excl[i] = float(dist_r[np.array(rb_only_idx)].mean())

    mu = float(rand.mean())
    sd = float(rand.std())
    p_emp = (1 + int((rand <= d_obs).sum())) / (n_perm + 1)
    if sd > 0:
        result.update(ok=True, d_c=d_obs, z=float((d_obs - mu) / sd),
                      p_empirical=float(p_emp), mean_random=mu, std_random=sd)
    else:
        # 零分布退化(所有随机距离相同):z 无意义,输出 None(与去交集分支一致),
        # 不把可能的显著结果伪装成 z=0 的"不显著"
        result.update(ok=True, d_c=d_obs, z=None, p_empirical=float(p_emp),
                      mean_random=mu, std_random=sd)

    if d_obs_excl is not None:
        valid = rand_excl[~np.isnan(rand_excl)]
        if len(valid) > 0 and float(valid.std()) > 0:
            mu_e, sd_e = float(valid.mean()), float(valid.std())
            z_e = (d_obs_excl - mu_e) / sd_e
            p_e = (1 + int((valid <= d_obs_excl).sum())) / (len(valid) + 1)
            result.update(d_c_excl=d_obs_excl, z_excl=float(z_e), p_excl=float(p_e),
                          mean_random_excl=mu_e)
        else:
            # 零分布退化（如去交集后随机对全为同距离）：不输出误导性的 z/p
            result.update(d_c_excl=d_obs_excl, z_excl=None, p_excl=None)
    return result


def proximity_zscore(drug_genes: list, disease_genes: list, species: str = "Human",
                     n_perm: int | None = None, seed: int | None = None) -> dict:
    """网络邻近度 Z-score 主入口（背景互作组 LCC 内计算，度匹配随机化置换检验）。
    任一基因集入网率低于 config.PROXIMITY_MIN_COVERAGE 或过滤后为空时返回 ok=False,绝不输出不可靠结果;
    背景互作组未构建时抛 InteractomeNotBuiltError,由调用方转为用户提示。"""
    import config

    n_perm = int(n_perm if n_perm is not None else config.PROXIMITY_N_PERM)
    seed = int(seed if seed is not None else config.PROXIMITY_SEED)
    ctx = _lcc_context(species)

    def _norm(genes: list) -> list[str]:
        """规范化去重(大写、去空白),不做网络过滤:入网率保护与 coverage 统计
        必须以原始输入数为分母。旧实现先用 get_network_coverage 过滤成在网基因
        再传入,导致低入网率保护失效、coverage 分母虚高(只统计 LCC 占比)。"""
        seen: set[str] = set()
        out: list[str] = []
        for raw in genes:
            g = str(raw).strip().upper()
            if g and g not in seen:
                seen.add(g)
                out.append(g)
        return out

    return _proximity_zscore_from_context(
        ctx, _norm(drug_genes), _norm(disease_genes), n_perm=n_perm, seed=seed,
        bin_size=config.PROXIMITY_DEGREE_BIN_SIZE,
        drug_ratio=config.PROXIMITY_MIN_COVERAGE,
        disease_ratio=config.PROXIMITY_MIN_COVERAGE,
    )


def _load_json_cache(cache_path: str, label: str, validate=None) -> dict | None:
    """读取 JSON 结果缓存：未命中/过期返回 None；损坏或未通过校验则清除后返回 None。"""
    import json

    from .caching import _is_cache_fresh

    if not _is_cache_fresh(cache_path):
        return None
    try:
        with open(cache_path, "r", encoding="utf-8") as f:
            cached = json.load(f)
        if validate is not None:
            validate(cached)
        return cached
    except Exception:
        # 损坏/过期缓存按"不存在"处理并清除
        logging.exception("%s cache unreadable, removing: %s", label, cache_path)
        try:
            os.remove(cache_path)
        except OSError:
            pass
        return None


def _save_json_cache(cache_path: str, result: dict) -> None:
    """计算成功后把结果写入 JSON 缓存(临时文件+原子替换,并发写不会留下半截文件)。"""
    from .caching import save_json_file

    save_json_file(cache_path, result)


def _clean_expired_caches(directory: str) -> None:
    """删除目录里超过缓存有效期(CACHE_TTL_SECONDS)的 JSON 缓存;调用方进入缓存读写前调用一次即可。
    这些结果已永久保存在分析结果目录(<prefix>_data)中,过期计算缓存删除即可,下次分析会重算并覆盖。
    """
    import time

    import config

    if not os.path.isdir(directory):
        return
    cutoff = time.time() - config.CACHE_TTL_SECONDS
    for name in os.listdir(directory):
        if not name.endswith(".json"):
            continue
        path = os.path.join(directory, name)
        try:
            if os.path.getmtime(path) < cutoff:
                os.remove(path)
                logging.info("Expired cache removed: %s", name)
        except OSError:
            logging.exception("Failed to remove expired cache %s", path)


def run_proximity_cached(drug_genes: list, disease_genes: list,
                         species: str = "Human", label: str | None = None) -> dict:
    """邻近度分析入口（带文件缓存）：成功才写缓存,失败/背景网络未构建即时返回不缓存。
    同名覆盖式缓存:文件名只由标识(label)决定,命中由文件内缓存键校验,基因集/参数变了重算覆盖;
    过期缓存文件(含其它分析的)自动删除,结果已在分析结果目录永久保存。"""
    import hashlib

    import config

    from . import interactome as itx

    def _norm(genes: list) -> list[str]:
        return sorted({str(g).strip().upper() for g in genes if g and str(g).strip()})

    payload = ("\r".join(_norm(drug_genes)) + "||" + "\r".join(_norm(disease_genes))
               + f"|{species}|{config.PROXIMITY_N_PERM}|{config.PROXIMITY_SEED}"
               + f"|{config.PROXIMITY_DEGREE_BIN_SIZE}"
               + f"|{config.STRING_VERSION}|{config.INTERACTOME_CONFIDENCE}"
               + "|v2-excl")  # v2：结果新增去交集字段，旧缓存（v1）自动失效重算
    md5 = hashlib.md5(payload.encode("utf-8")).hexdigest()
    name = f"proximity_{label}.json" if label else f"proximity_{md5[:12]}.json"
    os.makedirs(config.PROXIMITY_CACHE_DIR, exist_ok=True)
    _clean_expired_caches(config.PROXIMITY_CACHE_DIR)
    cache_path = os.path.join(config.PROXIMITY_CACHE_DIR, name)

    cached = _load_json_cache(cache_path, "Proximity")
    if cached is not None and cached.get("cache_key") == md5:
        cached["cached"] = True
        logging.info("Proximity cache hit: %s", name)
        return cached
    if cached is not None:
        # 同名但输入(基因集/参数)变了:按未命中处理,重算后覆盖
        logging.info("Proximity cache key changed, recompute: %s", name)

    try:
        result = proximity_zscore(drug_genes, disease_genes, species=species)
    except itx.InteractomeNotBuiltError:
        logging.warning("Proximity skipped: interactome not built")
        return {"ok": False, "error": "not_built", "d_c": None, "z": None,
                "p_empirical": None, "coverage": {}, "cached": False}
    if result.get("ok"):
        result["cache_key"] = md5
        _save_json_cache(cache_path, result)
    return result


# ---------- RWR 靶点优先级排序（M4，依赖 core.interactome 的背景互作组） ----------
# 方法学：Random Walk with Restart——从种子基因出发在网络上随机游走，每步以概率 α
# 继续游走、以 (1-α) 概率回到种子；稳态得分高的非种子基因即种子网络邻域内的候选靶点。


def _rwr_power_iteration(W, seed_idx: list, alpha: float, tol: float, max_iter: int) -> tuple:
    """幂迭代求解 r = α·W·r + (1-α)·p。返回 (得分向量, 迭代轮数, 是否收敛)。"""
    import numpy as np

    n = W.shape[0]
    p = np.zeros(n)
    p[np.asarray(seed_idx, dtype=int)] = 1.0 / len(seed_idx)
    r = (1.0 - alpha) * p.copy()
    converged = False
    iters = 0
    for iters in range(1, max_iter + 1):
        r_new = alpha * (W @ r) + (1.0 - alpha) * p
        if float(np.abs(r_new - r).sum()) < tol:
            r = r_new
            converged = True
            break
        r = r_new
    return r, iters, converged


def _rwr_from_context(ctx: tuple, seed_genes: list, alpha: float,
                      tol: float, max_iter: int, top_nonseed: int | None = None,
                      annotate: set | None = None, null_perms: int | None = None,
                      null_seed: int = 0, bin_size: int = 100) -> dict:
    """给定 LCC 上下文的 RWR 快路径（复用背景网络的邻接与度数表）。
    annotate 非空时附加非种子行的 "annotated" 标注与度匹配校正的富集检验——比较标注集与度匹配随机基因集的
    平均 RWR 得分,输出 z 与单侧经验 p;z 显著 > 0 说明标注集高分超越其度数所能解释的水平,不是枢纽度数假象。"""
    import numpy as np
    import scipy.sparse as sp

    nodes, idx, Mt, degree_map = ctx
    result: dict = {
        "ok": False, "error": None, "top": [], "alpha": float(alpha),
        "n_iter": None, "converged": None, "seed_mass": None,
        "n_seeds": None, "lcc_nodes": len(nodes), "cached": False,
        "coverage": {}, "enrichment": None,
    }
    # 列和 = 度 → 按列归一化得转移矩阵
    col_deg = np.array([degree_map[g] for g in nodes], dtype=float)
    W = Mt @ sp.diags(1.0 / col_deg)

    seeds_in = sorted({g for g in seed_genes if g in idx})
    result["coverage"] = {"seed_in_lcc": len(seeds_in), "seed_total": len(seed_genes)}
    if not seeds_in:
        result["error"] = "empty_after_filter"
        return result

    seed_idx = [idx[g] for g in seeds_in]
    scores, iters, converged = _rwr_power_iteration(W, seed_idx, alpha, tol, max_iter)
    seed_set = set(seeds_in)
    order = np.argsort(-scores)
    # 种子全量保留、非种子按分数截断 top_nonseed：分开返回，避免大种子集
    # 把非种子全部挤出列表，导致页面"非种子 Top-N"表变空。
    seeds_rows = [{"gene": nodes[i], "score": float(scores[i])}
                  for i in order if nodes[i] in seed_set]
    nonseed_rows = [{"gene": nodes[i], "score": float(scores[i])}
                    for i in order if nodes[i] not in seed_set]
    if top_nonseed is not None:
        nonseed_rows = nonseed_rows[:int(top_nonseed)]
    result.update(
        ok=True,
        n_seeds=len(seeds_in),
        n_iter=int(iters),
        converged=bool(converged),
        seed_mass=float(scores[seed_idx].sum()),   # 质量留在种子上的比例（信息量指标）
        total_mass=float(scores.sum()),            # 质量守恒校验（连通图应≈1.0）
        seeds=seeds_rows,
        top_nonseed=nonseed_rows,
    )

    # ---- 度匹配校正：标注集富集检验（零分布 = 度匹配随机基因集的平均得分） ----
    if annotate:
        ann = sorted({str(g).strip().upper() for g in annotate
                      if g and str(g).strip().upper() in idx})
        if ann:
            ann_idx = np.array([idx[g] for g in ann])
            obs_mean = float(scores[ann_idx].mean())
            rng = np.random.default_rng(null_seed)
            perms = int(null_perms) if null_perms else 0
            null_means = np.full(max(perms, 0), np.nan)
            bins = _build_degree_bins(degree_map, bin_size)  # 置换循环内复用,避免每轮重建分箱
            for i in range(perms):
                sample = degree_matched_sample(degree_map, ann, bin_size, rng, bins=bins)
                s_idx = np.array([idx[g] for g in sample])
                null_means[i] = float(scores[s_idx].mean())
            valid = null_means[~np.isnan(null_means)]
            enrichment = {
                "n_annotated": len(ann),
                "obs_mean": obs_mean,
                "null_mean": float(valid.mean()) if len(valid) else None,
                "n_perms": int(len(valid)),
                "z": None, "p": None,
            }
            if len(valid) > 0 and float(valid.std()) > 0:
                mu_n, sd_n = float(valid.mean()), float(valid.std())
                enrichment["z"] = float((obs_mean - mu_n) / sd_n)
                # 单侧经验 p：极端方向是 obs 偏大（富集），故数 null >= obs（与邻近度的 d_c 偏小相反）。
                enrichment["p"] = float((1 + int((valid >= obs_mean).sum())) / (len(valid) + 1))
            result["enrichment"] = enrichment
            # 行级标注（种子行不含——种子本身即查询集合，标注无信息量）
            ann_set = set(ann)
            for row in nonseed_rows:
                row["annotated"] = row["gene"] in ann_set
    return result


def run_rwr_cached(seed_genes: list, species: str = "Human", label: str | None = None,
                   annotate: set | None = None) -> dict:
    """RWR 入口（带文件缓存）：成功才写缓存,失败/背景网络未构建即时返回不缓存,非种子行截断为 config.RWR_CACHE_TOP_N。
    标注集必须参与缓存键:度匹配校正依赖全量得分向量,只能在计算期完成,不能命中旧缓存后补算。
    同名覆盖式缓存:文件名只由标识(label,含药物/疾病侧)决定,命中由文件内缓存键校验,种子/参数变了重算覆盖;过期缓存自动删除,结果已永久保存。"""
    import hashlib

    import config

    from . import interactome as itx

    seeds_norm = sorted({str(g).strip().upper() for g in seed_genes if g and str(g).strip()})
    ann_norm = sorted({str(g).strip().upper() for g in annotate
                       if g and str(g).strip()}) if annotate else []
    payload = ("\r".join(seeds_norm) + "||" + "\r".join(ann_norm)
               + f"|{species}|{config.RWR_ALPHA}|{config.RWR_TOL}|{config.RWR_MAX_ITER}"
               + f"|{config.RWR_NULL_PERMS}|{config.RWR_NULL_SEED}"
               + f"|{config.PROXIMITY_DEGREE_BIN_SIZE}"
               + f"|{config.STRING_VERSION}|{config.INTERACTOME_CONFIDENCE}")
    md5 = hashlib.md5(payload.encode("utf-8")).hexdigest()
    name = f"rwr_{label}.json" if label else f"rwr_{md5[:12]}.json"
    os.makedirs(config.RWR_CACHE_DIR, exist_ok=True)
    _clean_expired_caches(config.RWR_CACHE_DIR)
    cache_path = os.path.join(config.RWR_CACHE_DIR, name)

    def _valid_schema(d: dict) -> None:
        if "top_nonseed" not in d or "seeds" not in d:
            # 旧版格式的残留缓存：按损坏处理，清除后重算
            raise ValueError("outdated RWR cache schema")

    result = _load_json_cache(cache_path, "RWR", validate=_valid_schema)
    if result is not None and result.get("cache_key") == md5:
        result["cached"] = True
        logging.info("RWR cache hit: %s", name)
        return result
    if result is not None:
        # 同名但输入(种子/参数)变了:按未命中处理,重算后覆盖
        logging.info("RWR cache key changed, recompute: %s", name)

    try:
        ctx = _lcc_context(species)
        result = _rwr_from_context(
            ctx, seeds_norm,
            alpha=config.RWR_ALPHA, tol=config.RWR_TOL, max_iter=config.RWR_MAX_ITER,
            top_nonseed=config.RWR_CACHE_TOP_N,
            annotate=annotate,
            null_perms=config.RWR_NULL_PERMS, null_seed=config.RWR_NULL_SEED,
            bin_size=config.PROXIMITY_DEGREE_BIN_SIZE)
    except itx.InteractomeNotBuiltError:
        logging.warning("RWR skipped: interactome not built")
        return {"ok": False, "error": "not_built", "seeds": [], "top_nonseed": [],
                "coverage": {}, "cached": False}
    if result.get("ok"):
        result["cache_key"] = md5
        _save_json_cache(cache_path, result)
        result["cached"] = False
    return result
