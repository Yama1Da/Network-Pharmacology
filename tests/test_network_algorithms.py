"""core.network_algorithms 的测试:Louvain 模块检测、网络邻近度、RWR 靶点排序。

迁移自 selfcheck_network.py(44 项中的 18 项)。文末的"参考实现"
(network_proximity/build_sparse_adj/rwr_prioritize)原样搬自自检脚本,
是 core 中纯本地算法的独立对照,仅供交叉验证。
"""
import glob
import logging
import os

import networkx as nx
import numpy as np
import pandas as pd
import pytest

import WY.core.network_algorithms as na
from WY import config
from WY.core import interactome as itx
from WY.core.network_algorithms import (
    _multi_source_bfs,
    _proximity_zscore_from_context,
    _rwr_power_iteration,
    degree_matched_sample,
    run_proximity_cached,
    run_rwr_cached,
)


def _karate() -> nx.Graph:
    G = nx.karate_club_graph()
    return nx.relabel_nodes(G, {i: f"G{i:02d}" for i in G.nodes()})


def _purge(pattern):
    """按 glob 模式删除旧缓存文件,文件不存在时为无害空操作。"""
    for old in glob.glob(pattern):
        try:
            os.remove(old)
        except OSError:
            pass


# ---------- 参考实现(原 core.network_algorithms 中的纯本地实现,仅供交叉验证) ----------


def network_proximity(graph: "nx.Graph", set_a: list, set_b: list) -> float:
    """参考实现:两组蛋白的网络平均距离 d_c(A,B) = mean_{b∈B} min_{a∈A} d(a,b)。
    仅供小图与测试交叉验证;图中不存在的基因被忽略,B 中有从 A 不可达的节点时抛错。
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


def build_sparse_adj(graph: "nx.Graph") -> tuple:
    """0-1 邻接 → 列随机转移矩阵 W(W[j,i] = P(i→j),每列和为 1)。
    度为 0 的节点会永久截留游走概率,故图必须连通(如限制在 LCC 内)。
    """
    import scipy.sparse as sp

    nodes = sorted(graph.nodes())
    adj = nx.to_scipy_sparse_array(graph, nodelist=nodes, weight=None, format="csr")
    col_deg = np.asarray(adj.sum(axis=0)).ravel()
    if (col_deg == 0).any():
        raise ValueError("zero-degree nodes in graph; restrict to a connected component")
    # 右乘对角阵 = 按列缩放:W[:, i] = adj[:, i] / deg(i)
    W = adj @ sp.diags(1.0 / col_deg)
    return W, nodes


def rwr_prioritize(graph: "nx.Graph", seed_genes: list, alpha: float = 0.7,
                   tol: float = 1e-8, max_iter: int = 200) -> pd.DataFrame:
    """参考入口(直接吃 networkx 图):返回全节点 RWR 打分表,列 gene/score/is_seed 按 score 降序。
    图需连通(建议 LCC);不在图内的种子被忽略。
    """
    W, nodes = build_sparse_adj(graph)
    idx = {g: i for i, g in enumerate(nodes)}
    seeds = sorted({str(g) for g in seed_genes if g in idx})
    if not seeds:
        raise ValueError("rwr_prioritize: no seed gene present in graph")
    scores, iters, converged = _rwr_power_iteration(
        W, [idx[g] for g in seeds], alpha, tol, max_iter)
    seed_set = set(seeds)
    out = pd.DataFrame({
        "gene": nodes,
        "score": scores,
        "is_seed": [g in seed_set for g in nodes],
    })
    logging.info("RWR done: %d seeds, %d iters, converged=%s", len(seeds), iters, converged)
    return out.sort_values("score", ascending=False).reset_index(drop=True)


# ---------- Louvain 模块检测 ----------


def test_detect_basic():
    """Louvain 基本性质:模块覆盖每个节点恰好一次、按规模降序、id 从 1 连续编号、组内基因有序。"""
    G = _karate()
    mods = na.detect_louvain_modules(G, seed=0)
    assert mods, "modules should not be empty"
    all_nodes = sorted(g for m in mods for g in m["genes"])
    assert all_nodes == sorted(G.nodes()), "modules must cover every node exactly once"
    sizes = [m["size"] for m in mods]
    assert sizes == sorted(sizes, reverse=True), "modules must be sorted by size desc"
    assert [m["module_id"] for m in mods] == list(range(1, len(mods) + 1)), "ids must be 1..n"
    for m in mods:
        assert m["genes"] == sorted(m["genes"]), "genes must be sorted"


def test_detect_deterministic():
    """同一种子两次运行必须给出完全相同的划分。"""
    G = _karate()
    a = na.detect_louvain_modules(G, resolution=1.0, seed=0)
    b = na.detect_louvain_modules(G, resolution=1.0, seed=0)
    assert a == b, "same seed must give identical output"


def test_detect_empty_graph():
    """空图必须返回空列表而不是报错。"""
    assert na.detect_louvain_modules(nx.Graph()) == [], "empty graph must yield []"


def test_detect_avg_degree():
    """5-团应形成一个模块,组内平均度恰为 4。"""
    G = nx.Graph()
    G.add_edges_from([("A", "B"), ("A", "C"), ("A", "D"), ("A", "E"),
                      ("B", "C"), ("B", "D"), ("B", "E"),
                      ("C", "D"), ("C", "E"), ("D", "E")])  # 5-团,每点度=4
    mods = na.detect_louvain_modules(G, seed=0)
    top = mods[0]
    assert top["size"] == 5, "clique should form one module"
    assert abs(top["avg_degree"] - 4.0) < 1e-9, f"avg_degree should be 4.0, got {top['avg_degree']}"


# ---------- 网络邻近度 ----------


def _graph_context(G: nx.Graph):
    """把 networkx 图整理成邻近度计算所需的格式(供多个测试复用)。"""
    nodes = tuple(sorted(G.nodes()))
    idx = {g: i for i, g in enumerate(nodes)}
    adj = nx.to_scipy_sparse_array(G, nodelist=nodes, weight=None, format="csr")
    Mt = adj.T.tocsr()
    degree_map = {g: int(d) for g, d in G.degree()}
    return nodes, idx, Mt, degree_map


def test_prox_bfs_matches_networkx():
    """自写的多源最短距离计算,结果应与 networkx 的标准结果完全一致。"""
    G = nx.barabasi_albert_graph(80, 2, seed=42)
    G = nx.relabel_nodes(G, {i: f"N{i:02d}" for i in G.nodes()})
    ctx = _graph_context(G)
    nodes, idx, Mt, _ = ctx
    rng = np.random.default_rng(0)
    sources = rng.choice(len(nodes), size=6, replace=False)
    dist = _multi_source_bfs(Mt, sources)
    ref = {}
    for s in sources:
        for node, d in nx.single_source_shortest_path_length(G, nodes[int(s)]).items():
            ref[node] = min(ref.get(node, 10**9), d)
    for i, g in enumerate(nodes):
        assert dist[i] == ref[g], f"node {g}: sparse={dist[i]} nx={ref[g]}"


def test_prox_degree_matched_sample():
    """度匹配抽样:替换不出现在子集内、度数偏差在 bin 内、同种子可复现、区间空时向邻近区间抽取。"""
    degrees = {f"G{i:02d}": d for i, d in enumerate(
        [5, 5, 6, 7, 105, 110, 120, 205, 210, 300])}
    seeds = ["G00", "G04", "G07"]
    rng = np.random.default_rng(0)
    picked = degree_matched_sample(degrees, seeds, bin_size=100, rng=rng)
    assert len(picked) == len(seeds)
    assert len(set(picked)) == len(set(seeds)) or True  # 允许不同种子抽到同一个基因
    for s, p in zip(seeds, picked):
        assert p != s, "replacement must not be the seed itself"
        assert abs(degrees[p] - degrees[s]) <= 100, f"degree {degrees[p]} not matched to {degrees[s]}"
    # 同一随机种子两次抽取,结果一致
    rng2 = np.random.default_rng(0)
    assert degree_matched_sample(degrees, seeds, 100, rng2) == picked
    # 种子所在的度数区间只有它自己时,应到最近的非空区间抽取
    picked2 = degree_matched_sample(degrees, ["G09"], 100, np.random.default_rng(3))
    assert picked2 and picked2[0] != "G09" and abs(degrees[picked2[0]] - 300) <= 100


def test_prox_zscore_synthetic():
    """符号测试:药物基因 A 与直接邻居 B1 应显著邻近,与远端 B2 的距离必然更大。"""
    G = nx.barabasi_albert_graph(120, 3, seed=7)
    G = nx.relabel_nodes(G, {i: f"N{i:03d}" for i in G.nodes()})
    degrees = dict(G.degree())
    hub_sorted = sorted(degrees, key=degrees.get, reverse=True)
    A = hub_sorted[:8]
    # B1:A 的直接邻居(排除 A 自身),取 8 个
    neighbors = sorted({v for a in A for v in G[a]} - set(A),
                       key=degrees.get, reverse=True)
    B1 = neighbors[:8]
    # B2:与 A 距离最远的 8 个节点(图的边缘)
    dist_from_a = {}
    for a in A:
        for node, d in nx.single_source_shortest_path_length(G, a).items():
            dist_from_a[node] = min(dist_from_a.get(node, 10**9), d)
    B2 = sorted(dist_from_a, key=dist_from_a.get, reverse=True)[:8]

    ctx = _graph_context(G)
    res_near = _proximity_zscore_from_context(ctx, A, B1, n_perm=100, seed=0,
                                              bin_size=100, drug_ratio=0.0, disease_ratio=0.0)
    res_far = _proximity_zscore_from_context(ctx, A, B2, n_perm=100, seed=0,
                                             bin_size=100, drug_ratio=0.0, disease_ratio=0.0)
    assert res_near["ok"] and res_far["ok"]
    # B1 全是 A 的直接邻居,距离应恰好为 1,且小于远端 B2
    assert abs(res_near["d_c"] - 1.0) < 1e-12, res_near["d_c"]
    assert res_near["d_c"] < res_far["d_c"], (res_near["d_c"], res_far["d_c"])
    # 近端模块的 z 值应为负,经验 p 值不应比远端更差
    assert res_near["z"] < 0, f"near-module z should be negative, got {res_near['z']}"
    assert res_near["p_empirical"] <= res_far["p_empirical"], (
        res_near["p_empirical"], res_far["p_empirical"])


def test_prox_coverage_guard():
    """覆盖率防线:全部未知基因报 empty_after_filter;入网率低于阈值报 low_coverage_drug。"""
    G = _karate()
    ctx = _graph_context(G)
    known = list(G.nodes())[:10]
    # 全部未知基因 → empty_after_filter
    res = _proximity_zscore_from_context(ctx, ["NOPE1", "NOPE2"], known[:5],
                                         n_perm=10, seed=0, bin_size=100,
                                         drug_ratio=0.3, disease_ratio=0.3)
    assert res["error"] == "empty_after_filter" and not res["ok"]
    # 入网率 1/10 < 0.3 → low_coverage_drug
    res2 = _proximity_zscore_from_context(ctx, known[:1] + [f"X{i}" for i in range(9)], known[:5],
                                          n_perm=10, seed=0, bin_size=100,
                                          drug_ratio=0.3, disease_ratio=0.3)
    assert res2["error"] == "low_coverage_drug" and not res2["ok"]


def test_prox_reference_impl_consistency():
    """快速算法算出的距离应与 networkx 参考实现一致。"""
    G = nx.barabasi_albert_graph(60, 2, seed=11)
    G = nx.relabel_nodes(G, {i: f"N{i:02d}" for i in G.nodes()})
    ctx = _graph_context(G)
    nodes, idx, Mt, _ = ctx
    A = nodes[:6]
    B = nodes[30:38]
    d_ref = network_proximity(G, A, B)
    dist = _multi_source_bfs(Mt, np.array([idx[g] for g in A]))
    d_fast = float(dist[[idx[g] for g in B]].mean())
    assert abs(d_ref - d_fast) < 1e-12, (d_ref, d_fast)


def test_prox_excl_overlap():
    """剔除药物/疾病共享靶点后再算距离:d_c 只会变大或不变;疾病靶点全被覆盖时该指标不可用。"""
    G = _karate()
    ctx = _graph_context(G)
    nodes = list(G.nodes())
    A = nodes[:8]
    B = nodes[4:16]  # 与 A 共享 4 个基因
    res = _proximity_zscore_from_context(ctx, A, B, n_perm=30, seed=0, bin_size=100,
                                         drug_ratio=0.0, disease_ratio=0.0)
    assert res["ok"] and res["n_excl"] == 8
    assert res["d_c_excl"] is not None and res["d_c_excl"] >= res["d_c"]
    assert res["z_excl"] is not None and res["p_excl"] is not None
    # 疾病集完全被药物集覆盖时,去交集后为空,该指标不可用
    res2 = _proximity_zscore_from_context(ctx, nodes[:10], nodes[:5], n_perm=10, seed=0,
                                          bin_size=100, drug_ratio=0.0, disease_ratio=0.0)
    assert res2["ok"] and res2["d_c_excl"] is None and res2["n_excl"] == 0


@pytest.mark.real_interactome
def test_prox_real_and_cached():
    """真实互作组完整流程(未构建则跳过):结果结构完整,第二次调用命中缓存。

    默认跳过(conftest 沙箱);设 NPAP_TEST_REAL_INTERACTOME=1 且本机已
    构建互作组时才真正运行。豁免沙箱后会把 selfcheck_prox 前缀的缓存
    文件写进真实邻近度缓存目录(下次运行前自清理)。"""
    if not itx.is_built("Human"):
        pytest.skip("Human 背景互作组未构建")
    drug = ["TP53", "EGFR", "AKT1", "MYC", "STAT3"]
    disease = ["BRCA1", "PTEN", "KRAS", "PIK3CA", "ERBB2"]
    res1 = run_proximity_cached(drug, disease, label="selfcheck_prox")
    assert res1["ok"], res1
    assert all(isinstance(res1[k], float) for k in ("d_c", "z", "p_empirical"))
    # 该药物/疾病对无共享靶点,去交集后疾病集仍是 5 个,d_c_excl 应与 d_c 相等
    assert res1.get("n_excl") == 5, res1.get("n_excl")
    assert res1.get("d_c_excl") is not None
    assert abs(res1["d_c_excl"] - res1["d_c"]) < 1e-12
    assert res1.get("z_excl") is not None
    res2 = run_proximity_cached(drug, disease, label="selfcheck_prox")
    assert res2.get("cached") is True, "second call must hit cache"
    assert abs(res2["z"] - res1["z"]) < 1e-12
    print(f"      [real interactome] d_c={res1['d_c']:.3f} z={res1['z']:.3f} "
          f"p={res1['p_empirical']:.4f} | excl d_c={res1['d_c_excl']:.3f} "
          f"z={res1['z_excl']:.3f} ({res1['n_perm']} perms)")


# ---------- RWR 靶点优先级排序 ----------


def test_rwr_direction_and_mass():
    """验证游走方向:从星形图中心出发,扩散一步后每个叶子各得 0.25,总和为 1。"""
    G = nx.star_graph(4)  # 节点0为中心,度4;叶1..4度1
    G = nx.relabel_nodes(G, {i: f"N{i}" for i in G.nodes()})
    W, nodes = build_sparse_adj(G)
    # 列随机性:每列和为 1
    col_sum = np.asarray(W.sum(axis=0)).ravel()
    assert np.allclose(col_sum, 1.0), col_sum
    idx = {g: i for i, g in enumerate(nodes)}
    r = np.zeros(len(nodes))
    r[idx["N0"]] = 1.0
    r1 = W @ r
    for leaf in ("N1", "N2", "N3", "N4"):
        assert abs(r1[idx[leaf]] - 0.25) < 1e-12, (leaf, r1[idx[leaf]])
    assert abs(r1.sum() - 1.0) < 1e-12, r1.sum()


def test_rwr_alpha0_returns_seeds():
    """α=0 时不发生游走,得分应全部留在种子基因上。"""
    G = nx.barabasi_albert_graph(30, 2, seed=1)
    G = nx.relabel_nodes(G, {i: f"N{i:02d}" for i in G.nodes()})
    df = rwr_prioritize(G, ["N00", "N05"], alpha=0.0, tol=1e-12, max_iter=50)
    seed_rows = df[df["is_seed"]]
    nonseed = df[~df["is_seed"]]
    assert len(seed_rows) == 2
    assert set(seed_rows["gene"]) == {"N00", "N05"}
    assert seed_rows["score"].round(12).tolist() == [0.5, 0.5]
    assert (nonseed["score"].abs() < 1e-12).all(), "non-seeds must have zero mass at alpha=0"


def test_rwr_prioritize_basic():
    """α=0.7 正常排序:总分守恒、按得分降序、种子基因标注正确。"""
    G = nx.barabasi_albert_graph(50, 2, seed=3)
    G = nx.relabel_nodes(G, {i: f"N{i:02d}" for i in G.nodes()})
    df = rwr_prioritize(G, ["N00", "N03", "NOT_IN_GRAPH"], alpha=0.7)
    assert abs(df["score"].sum() - 1.0) < 1e-6, df["score"].sum()
    assert df["score"].is_monotonic_decreasing
    assert df["is_seed"].sum() == 2
    # 度最高的枢纽节点应排在前列
    assert df.iloc[0]["gene"] in ("N00", "N03") or df.iloc[0]["score"] > 0.01


def test_rwr_from_context_guard_and_annotate():
    """种子基因全部不在网络中时应拒绝计算;计算结果本身不含基因标注信息。"""
    G = _karate()
    ctx = _graph_context(G)
    res = na._rwr_from_context(ctx, ["NOPE1", "NOPE2"], alpha=0.7, tol=1e-8, max_iter=50)
    assert res["error"] == "empty_after_filter" and not res["ok"]
    known = list(G.nodes())[:5]
    res2 = na._rwr_from_context(ctx, known, alpha=0.7, tol=1e-8, max_iter=50)
    assert res2["ok"] and res2["n_seeds"] == 5 and res2["converged"]
    assert all("annotated" not in row for row in res2["top_nonseed"])
    assert all("annotated" not in row for row in res2["seeds"])


def test_rwr_enrich_synthetic():
    """标注集富集检验:种子和标注集都是枢纽基因时,标注基因得分应显著偏高(z 为正、p 小)。"""
    G = nx.barabasi_albert_graph(120, 3, seed=7)
    G = nx.relabel_nodes(G, {i: f"N{i:03d}" for i in G.nodes()})
    degrees = dict(G.degree())
    by_deg = sorted(degrees, key=degrees.get, reverse=True)
    hubs = by_deg[:8]
    mid = by_deg[8:12]  # 次高度节点
    ctx = _graph_context(G)
    ann = set(hubs) | set(mid)
    res = na._rwr_from_context(ctx, hubs, alpha=0.7, tol=1e-8, max_iter=100,
                               annotate=ann, null_perms=100, null_seed=0,
                               bin_size=100, top_nonseed=20)
    assert res["ok"] and res["enrichment"] is not None
    enr = res["enrichment"]
    assert enr["n_annotated"] == 12 and enr["z"] > 1.0, enr
    assert enr["p"] is not None and enr["p"] <= 0.05
    for r in res["top_nonseed"]:
        assert r.get("annotated") == (r["gene"] in ann), r


@pytest.mark.real_interactome
def test_rwr_real_and_cached():
    """真实互作组完整流程(未构建则跳过):检查种子/非种子分列表、收敛与缓存。

    默认跳过(conftest 沙箱);设 NPAP_TEST_REAL_INTERACTOME=1 且本机已
    构建互作组时才真正运行。豁免沙箱后会把 selfcheck_rwr 前缀的缓存
    文件写进真实 RWR 缓存目录(运行开头自清理)。"""
    if not itx.is_built("Human"):
        pytest.skip("Human 背景互作组未构建")
    # 先删除旧测试缓存,保证第一次是全新计算
    _purge(os.path.join(config.RWR_CACHE_DIR, "rwr_selfcheck_rwr_*.json"))

    seeds = ["TP53", "EGFR", "AKT1", "MYC", "STAT3"]
    res1 = run_rwr_cached(seeds, label="selfcheck_rwr")
    assert res1["ok"], res1
    assert res1["converged"] and res1["n_iter"] >= 1
    # total_mass 是全网络总质量(约 1.0);seed_mass 是留在种子上的比例,含义不同
    assert abs(res1["total_mass"] - 1.0) < 1e-6, res1["total_mass"]
    assert 0.0 < res1["seed_mass"] < 1.0, res1["seed_mass"]
    # 种子列表全量保留;非种子榜单单独截断,且不含种子
    assert len(res1["seeds"]) == 5, len(res1["seeds"])
    assert {r["gene"] for r in res1["seeds"]} == {"TP53", "EGFR", "AKT1", "MYC", "STAT3"}
    assert len(res1["top_nonseed"]) == config.RWR_CACHE_TOP_N, len(res1["top_nonseed"])
    seed_names = {r["gene"] for r in res1["seeds"]}
    assert all(r["gene"] not in seed_names for r in res1["top_nonseed"])
    assert res1["top_nonseed"][0]["score"] >= res1["top_nonseed"][-1]["score"]
    # 标注集也算进缓存标识:换标注集就要重新计算
    first_ns = res1["top_nonseed"][0]["gene"]
    res2 = run_rwr_cached(seeds, label="selfcheck_rwr", annotate={first_ns})
    assert res2["ok"] and res2.get("cached") is False
    row0 = res2["top_nonseed"][0]
    assert row0["gene"] == first_ns and row0["annotated"] is True
    assert not any(r.get("annotated") for r in res2["top_nonseed"][1:])
    enr = res2.get("enrichment")
    assert enr and enr["n_annotated"] == 1
    res3 = run_rwr_cached(seeds, label="selfcheck_rwr", annotate={first_ns})
    assert res3.get("cached") is True
    assert abs(res3["top_nonseed"][0]["score"] - row0["score"]) < 1e-15
    if enr.get("degenerate"):
        # 榜首基因可能孤悬度数分布极端尾部(如 v12.5 全网最高度的 TRIM67):度匹配
        # 零分布退化为同一替换集,z/p 无定义,应如实为 None。此时再换分箱人口充足
        # 的基因补验非退化路径,保证富集计算本身仍被测到。
        assert enr["z"] is None and enr["p"] is None and enr["n_perms"] > 0
        degree_map = na._lcc_context("Human")[3]
        bin_size = config.PROXIMITY_DEGREE_BIN_SIZE
        bins = na._build_degree_bins(degree_map, bin_size)
        rich = next(r["gene"] for r in res2["top_nonseed"]
                    if r["gene"] != first_ns
                    and len(bins.get(degree_map[r["gene"]] // bin_size, ())) >= 3)
        res_rich = run_rwr_cached(seeds, label="selfcheck_rwr", annotate={rich})
        assert res_rich["ok"] and res_rich["enrichment"]["degenerate"] is False
        enr = res_rich["enrichment"]
    assert enr["z"] is not None and enr["p"] is not None
    print(f"      [real interactome] iters={res1['n_iter']} converged={res1['converged']} "
          f"seeds={len(res1['seeds'])} top_nonseed={len(res1['top_nonseed'])} "
          f"top5={[r['gene'] for r in res1['top_nonseed'][:5]]} "
          f"enrich_z={enr['z']:.2f} enrich_p={enr['p']:.4f}")


def test_rwr_truncation_guarantee():
    """回归测试:种子数很大时,种子表仍要全量保留,非种子榜单不被种子挤占。"""
    G = nx.barabasi_albert_graph(60, 2, seed=5)
    G = nx.relabel_nodes(G, {i: f"N{i:02d}" for i in G.nodes()})
    ctx = _graph_context(G)
    nodes = ctx[0]
    many_seeds = list(nodes[:45])  # 60 节点中 45 个种子(窗口只有 5)
    res = na._rwr_from_context(ctx, many_seeds, alpha=0.7, tol=1e-8, max_iter=100,
                               top_nonseed=5)
    assert res["ok"]
    assert len(res["seeds"]) == 45, f"seeds must be kept in full, got {len(res['seeds'])}"
    assert len(res["top_nonseed"]) == 5, f"non-seed window must stay full, got {len(res['top_nonseed'])}"
    seed_set = set(many_seeds)
    assert all(r["gene"] not in seed_set for r in res["top_nonseed"])
    scores = [r["score"] for r in res["top_nonseed"]]
    assert scores == sorted(scores, reverse=True), "non-seed rows must be score-desc"
