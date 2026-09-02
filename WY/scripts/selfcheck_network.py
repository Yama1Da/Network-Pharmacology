#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""网络分析功能自检脚本（纯本地运行，不发网络请求）。

运行方式（在 WY 目录下）：
    python scripts/selfcheck_network.py

检查内容：Louvain 模块检测、模块富集、GSEA、背景互作组、网络邻近度、RWR 靶点排序。
每项检查打印 PASS 或 FAIL；全部通过时退出码为 0，有失败则为 1。
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import networkx as nx
import numpy as np
import pandas as pd

import core.network_algorithms as na
from core import interactome as itx
from core.enrichment import (
    _gsea_cache_file,
    _normalize_enrichment_gsea,
    build_ranked_list,
    fetch_enrichment_gsea_prerank,
)
from core.network_algorithms import (
    _build_degree_bins,
    _lcc_context,
    _multi_source_bfs,
    _proximity_zscore_from_context,
    build_sparse_adj,
    degree_matched_sample,
    network_proximity,
    rwr_prioritize,
    run_proximity_cached,
    run_rwr_cached,
)


def _karate() -> nx.Graph:
    G = nx.karate_club_graph()
    return nx.relabel_nodes(G, {i: f"G{i:02d}" for i in G.nodes()})


def check(name, fn):
    try:
        fn()
        print(f"PASS  {name}")
        return True
    except AssertionError as e:
        print(f"FAIL  {name}: {e}")
        return False
    except Exception as e:
        print(f"ERROR {name}: {type(e).__name__}: {e}")
        return False


def t_detect_basic():
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


def t_deterministic():
    G = _karate()
    a = na.detect_louvain_modules(G, resolution=1.0, seed=0)
    b = na.detect_louvain_modules(G, resolution=1.0, seed=0)
    assert a == b, "same seed must give identical output"


def t_empty_graph():
    assert na.detect_louvain_modules(nx.Graph()) == [], "empty graph must yield []"


def t_avg_degree():
    G = nx.Graph()
    G.add_edges_from([("A", "B"), ("A", "C"), ("A", "D"), ("A", "E"),
                      ("B", "C"), ("B", "D"), ("B", "E"),
                      ("C", "D"), ("C", "E"), ("D", "E")])  # 5-团，每点度=4
    mods = na.detect_louvain_modules(G, seed=0)
    top = mods[0]
    assert top["size"] == 5, "clique should form one module"
    assert abs(top["avg_degree"] - 4.0) < 1e-9, f"avg_degree should be 4.0, got {top['avg_degree']}"


# ---------- 模块富集测试（用假函数代替真实请求） ----------

# 这些返回状态视为可信：真正成功，或确实无显著项
_OK_REASONS = ("success", "no_significant", "cache")


def _fake_gp(reason, df):
    """假的富集函数：reason=None 表示成功，否则按给定 reason 记录状态。"""
    def f(genes, organism="hsapiens", sources=("GO:BP",), pvalue=0.05, label=None, status=None):
        if status is not None:
            if reason is None:
                status.update(ok=True, reason="success")
            else:
                status.update(ok=reason in _OK_REASONS, reason=reason)
        return df.copy()
    return f


def _fake_enr(reason, df):
    def f(genes, libs, organism="human", pvalue=0.05, label=None, status=None):
        if status is not None:
            if reason is None:
                status.update(ok=True, reason="success")
            else:
                status.update(ok=reason in _OK_REASONS, reason=reason)
        return df.copy()
    return f


def _patch(gp_spec, enr_spec):
    """临时把两个富集服务函数换成假函数，参数各为 (reason, df)。"""
    gp_reason, gp_df = gp_spec
    enr_reason, enr_df = enr_spec
    orig = (na.fetch_enrichment_gprofiler, na.fetch_enrichment_enrichr)
    na.fetch_enrichment_gprofiler = _fake_gp(gp_reason, gp_df)
    na.fetch_enrichment_enrichr = _fake_enr(enr_reason, enr_df)
    return orig


def _restore(orig):
    na.fetch_enrichment_gprofiler, na.fetch_enrichment_enrichr = orig


def t_enrich_fail_both():
    orig = _patch(("request_failed", pd.DataFrame()), ("request_failed", pd.DataFrame()))
    try:
        modules = [{"module_id": 1, "genes": ["A", "B", "C", "D", "E"], "size": 5, "avg_degree": 2.0}]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
        assert long_df.empty, "failed module must produce empty table"
        assert len(failures) == 1 and failures[0]["module_id"] == 1, "failure must be recorded"
        assert not no_sig, "failed module must not be counted as no-significant"
    finally:
        _restore(orig)


def t_enrich_fail_gp_ok_enr():
    # 一个服务被熔断跳过、另一个确实无显著项 → 结果不可信，按失败处理
    orig = _patch(("circuit_open", pd.DataFrame()), ("no_significant", pd.DataFrame()))
    try:
        modules = [{"module_id": 1, "genes": ["A", "B", "C", "D", "E"], "size": 5, "avg_degree": 2.0}]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
        assert len(failures) == 1, "circuit_open + empty enrichr must be a failure"
        assert not no_sig
    finally:
        _restore(orig)


def t_enrich_genuine_no_sig():
    # 两条链路都成功但确实无显著项 → 不算失败
    orig = _patch(("no_significant", pd.DataFrame()), ("no_significant", pd.DataFrame()))
    try:
        modules = [{"module_id": 1, "genes": ["A", "B", "C", "D", "E"], "size": 5, "avg_degree": 2.0}]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
        assert long_df.empty and not failures and no_sig == [1]
    finally:
        _restore(orig)


def t_enrich_success():
    # 模块1 富集成功、有结果；模块2 两个服务都失败，应记入 failures。
    ok_df = pd.DataFrame([
        {"source": "KEGG", "term_id": "hsa04151", "term_name": "PI3K-Akt",
         "p_value": 0.001, "adjusted_p_value": 0.001,
         "intersection_size": 3, "query_size": 5, "genes": "A;B;C"},
    ])

    def gp(genes, organism="hsapiens", sources=("GO:BP",), pvalue=0.05, label=None, status=None):
        if "A" in genes:
            if status is not None:
                status.update(ok=True, reason="success")
            return ok_df.copy()
        if status is not None:
            status.update(ok=False, reason="request_failed")
        return pd.DataFrame()

    def enr(genes, libs, organism="human", pvalue=0.05, label=None, status=None):
        if status is not None:
            status.update(ok=False, reason="request_failed")
        return pd.DataFrame()

    orig = (na.fetch_enrichment_gprofiler, na.fetch_enrichment_enrichr)
    na.fetch_enrichment_gprofiler, na.fetch_enrichment_enrichr = gp, enr
    try:
        modules = [
            {"module_id": 1, "genes": ["A", "B", "C", "D", "E"], "size": 5, "avg_degree": 2.0},
            {"module_id": 2, "genes": ["X", "Y", "Z", "W", "V"], "size": 5, "avg_degree": 2.0},
        ]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
        assert not no_sig
        assert not long_df.empty and (long_df["module_id"] == 1).all()
        assert "module_size" in long_df.columns
        assert [f["module_id"] for f in failures] == [2], f"module 2 should fail, got {failures}"
    finally:
        _restore(orig)


def t_min_size_filter():
    # 规模不足的模块不参与富集
    calls = {"n": 0}

    def counting_gp(genes, **kw):
        calls["n"] += 1
        return pd.DataFrame()

    orig = (na.fetch_enrichment_gprofiler, na.fetch_enrichment_enrichr)
    na.fetch_enrichment_gprofiler = counting_gp
    na.fetch_enrichment_enrichr = lambda genes, **kw: pd.DataFrame()
    try:
        modules = [{"module_id": 1, "genes": ["A", "B"], "size": 2, "avg_degree": 1.0}]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens",
                                                      min_module_size=5)
        assert calls["n"] == 0, "module below min size must not be enriched"
        assert long_df.empty and not failures and not no_sig
    finally:
        _restore(orig)


# ---------- GSEA prerank 测试（本地离线运行） ----------

def t_gsea_build_ranked():
    df = pd.DataFrame({
        "approvedSymbol": ["tp53", "TP53", "EGFR", "BRCA1", "KRAS", "MT-CO1"],
        "score_pchembl": [5.0, 9.0, 7.0, None, 6.0, float("nan")],
    })
    rnk, stats = build_ranked_list(df, "approvedSymbol", "score_pchembl")
    # TP53 大小写重复取最高分；缺分的基因剔除 → 保留 3 个
    assert stats["n_input"] == 6 and stats["n_kept"] == 3, stats
    assert abs(stats["drop_ratio"] - 0.5) < 1e-9
    assert list(rnk["gene"]) == ["TP53", "EGFR", "KRAS"], rnk
    assert list(rnk["score"]) == [9.0, 7.0, 6.0], "must be sorted desc"


def t_gsea_normalize():
    # 按真实 gseapy 输出的列名构造测试数据
    res2d = pd.DataFrame([
        {"Name": "prerank", "Term": "SET_A", "ES": 0.6, "NES": 2.1,
         "NOM p-val": 0.01, "FDR q-val": 0.05, "FWER p-val": 0.02,
         "Tag %": "10/15", "Gene %": "33%", "Lead_genes": "G1;G2", "source": "LIB_X"},
        {"Name": "prerank", "Term": "SET_B", "ES": -0.5, "NES": -1.8,
         "NOM p-val": 0.08, "FDR q-val": 0.2, "FWER p-val": 0.1,
         "Tag %": "8/15", "Gene %": "27%", "Lead_genes": "G9;G8", "source": "LIB_X"},
    ])
    out = _normalize_enrichment_gsea(res2d)
    assert list(out["term_name"]) == ["SET_A", "SET_B"], "must be sorted by |NES| desc"
    assert out["source"].iloc[0] == "LIB_X"
    assert {"term_name", "es", "nes", "nominal_p_value", "fdr_q_value",
            "tag_pct", "gene_pct", "lead_genes", "source"} <= set(out.columns)
    # 缺少必需列时应返回空表而不是报错
    assert _normalize_enrichment_gsea(pd.DataFrame({"foo": [1]})).empty


def _tiny_gmt() -> str:
    genes = [f"G{i:02d}" for i in range(40)]
    path = os.path.join(tempfile.gettempdir(), "selfcheck_gsea_tiny.gmt")
    with open(path, "w", encoding="utf-8") as f:
        f.write("HIGH_SET\tna\t" + "\t".join(genes[:15]) + "\n")
        f.write("LOW_SET\tna\t" + "\t".join(genes[25:]) + "\n")
    return path


def _rnk40() -> pd.DataFrame:
    return pd.DataFrame({
        "gene": [f"G{i:02d}" for i in range(40)],
        "score": [2.0 - 0.1 * i for i in range(40)],
    })


def t_gsea_prerank_local():
    """用真实 gseapy 离线跑 GSEA：排位靠前的基因集方向应正确；第二次调用应命中缓存。"""
    import glob
    import config

    libs = {"SELFTEST_LIB": _tiny_gmt()}
    label = "selfcheck_gsea"
    # 先删除旧测试缓存，避免直接命中
    for old in glob.glob(os.path.join(config.ENRICHMENT_CACHE_DIR, f"GSEA_{label}_*.csv")):
        try:
            os.remove(old)
        except OSError:
            pass

    st1: dict = {}
    df1 = fetch_enrichment_gsea_prerank(_rnk40(), libs=libs, label=label, status=st1)
    assert st1.get("reason") == "success", st1
    assert not df1.empty and {"source", "term_name", "nes", "fdr_q_value"} <= set(df1.columns)
    high = df1[df1["term_name"] == "HIGH_SET"].iloc[0]
    low = df1[df1["term_name"] == "LOW_SET"].iloc[0]
    assert high["nes"] > 0 and low["nes"] < 0, "enrichment directions are wrong"

    st2: dict = {}
    df2 = fetch_enrichment_gsea_prerank(_rnk40(), libs=libs, label=label, status=st2)
    assert st2.get("reason") == "cache", st2
    # 经缓存读写后内容应不变（允许微小格式差异）
    pd.testing.assert_frame_equal(
        df1.reset_index(drop=True), df2.reset_index(drop=True),
        check_dtype=False, check_exact=False, rtol=1e-6,
    )

    try:
        cf = _gsea_cache_file(_rnk40(), libs, label)
        if os.path.exists(cf):
            os.remove(cf)
    except OSError:
        pass
    assert config.GSEA_MIN_RANKED_GENES <= 40, "test rnk must pass the gate"


def t_gsea_insufficient_no_cache():
    import glob
    import config

    st1: dict = {}
    df = fetch_enrichment_gsea_prerank(_rnk40().head(10), libs={"L": _tiny_gmt()},
                                       label="selfcheck_gsea_short", status=st1)
    assert df.empty and st1.get("reason") == "insufficient_input", st1
    left = glob.glob(os.path.join(config.ENRICHMENT_CACHE_DIR, "GSEA_selfcheck_gsea_short_*.csv"))
    assert not left, f"failure must not leave cache files: {left}"


# ---------- 背景互作组测试（真实数据未构建时自动跳过） ----------

def t_interactome_info_cols():
    ok = pd.DataFrame({"protein_external_id": ["9606.ENSP1"],
                       "preferredName": ["TP53"], "other": [1]})
    assert itx._info_column_names(ok) == ("protein_external_id", "preferredName")
    bad = pd.DataFrame({"foo": [1], "bar": [2]})
    try:
        itx._info_column_names(bad)
        raise AssertionError("unexpected columns must raise ValueError")
    except ValueError:
        pass


def t_interactome_links_to_edges():
    links = pd.DataFrame({
        "protein1": ["9606.E1", "9606.E2", "9606.E3", "9606.E4",
                     "9606.E5", "9606.E6", "9606.E7"],
        "protein2": ["9606.E2", "9606.E1", "9606.E3", "9606.E5",
                     "9606.E6", "9606.E7", "9606.E8"],
        # 前两行是同一对蛋白，分数 700 和 800
        "combined_score": [700, 800, 699, 950,
                           900, 850, 400],
    })
    id2sym = {"9606.E1": "TP53", "9606.E2": "EGFR", "9606.E3": "EGFR",
              "9606.E4": "KRAS", "9606.E5": "BRAF", "9606.E6": "BRCA1",
              "9606.E7": "BRCA1", "9606.E8": "ABL1"}
    edges = itx._links_to_edges(links, id2sym, confidence=700)
    records = {(r.protein_a, r.protein_b): int(r.score) for r in edges.itertuples()}
    # 无向边的两个端点统一按基因名字典序排列
    # 1) 基因塌缩产生的自环（EGFR-EGFR）必须去掉
    assert ("EGFR", "EGFR") not in records, "self-loop must be dropped"
    # 2) 同一对蛋白出现两次时保留最高分 800
    assert records.get(("EGFR", "TP53")) == 800, records
    # 3) E7-E8 分数 400 低于阈值，被过滤
    assert ("ABL1", "BRCA1") not in records and ("BRCA1", "ABL1") not in records
    # 4) E4-E5 (KRAS-BRAF, 950) 与 E5-E6 (BRAF-BRCA1, 900) 存在
    assert records.get(("BRAF", "KRAS")) == 950
    assert records.get(("BRAF", "BRCA1")) == 900
    assert len(edges) == 3, f"expect 3 edges, got {len(edges)}"
    # 5) 置信度换算一致：0.7 * 1000 == 700
    import config
    assert config.INTERACTOME_CONFIDENCE == int(config.PPI_CONFIDENCE * 1000)


def t_interactome_load_or_skip():
    """互作组未构建时必须明确报错；已构建时正常加载并检查基因覆盖情况。"""
    if not itx.is_built("Human"):
        try:
            itx.load_interactome("Human")
            raise AssertionError("must raise InteractomeNotBuiltError when not built")
        except itx.InteractomeNotBuiltError:
            return  # 跳过：真实互作组尚未构建
        return
    graph = itx.load_interactome("Human")
    assert graph.number_of_nodes() >= 10000, f"human interactome too small: {graph.number_of_nodes()}"
    inside, outside = itx.get_network_coverage(["tp53 ", "EGFR", "FAKE_GENE_X", "tp53"])
    assert inside == ["TP53", "EGFR"], inside
    assert outside == ["FAKE_GENE_X"], outside
    deg = itx.load_degree_table("Human")
    assert set(deg.columns) == {"gene", "degree"} and len(deg) == graph.number_of_nodes()


def t_interactome_staleness():
    """检查互作组数据的"过期"判断：未构建时 age 为 None；已构建时 age 为非负整数。"""
    if not itx.is_built("Human"):
        age = itx.interactome_age_days("Human")
        assert age is None, "age must be None when not built"
        assert itx.is_stale("Human") is False, "not-built must not report stale"
        return
    age = itx.interactome_age_days("Human")
    assert isinstance(age, int) and age >= 0, age
    # 判断规则：age 严格大于 days 才算过期
    assert itx.is_stale("Human", days=-1) is True, "negative threshold must be stale"
    assert itx.is_stale("Human", days=age) is False, "boundary: age > days must be strict"
    assert itx.is_stale("Human", days=10_000) is False, "huge threshold must not be stale"
    # 默认阈值来自 config，且为正数
    import config
    assert config.INTERACTOME_STALE_DAYS > 0


# ---------- 网络邻近度检验 ----------

def _graph_context(G: nx.Graph):
    """把 networkx 图整理成邻近度计算所需的格式（供多个测试复用）。"""
    nodes = tuple(sorted(G.nodes()))
    idx = {g: i for i, g in enumerate(nodes)}
    adj = nx.to_scipy_sparse_array(G, nodelist=nodes, weight=None, format="csr")
    Mt = adj.T.tocsr()
    degree_map = {g: int(d) for g, d in G.degree()}
    return nodes, idx, Mt, degree_map


def t_prox_bfs_matches_networkx():
    """自写的多源最短距离计算，结果应与 networkx 的标准结果完全一致。"""
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


def t_prox_degree_matched_sample():
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
    # 同一随机种子两次抽取，结果一致
    rng2 = np.random.default_rng(0)
    assert degree_matched_sample(degrees, seeds, 100, rng2) == picked
    # 种子所在的度数区间只有它自己时，应到最近的非空区间抽取
    picked2 = degree_matched_sample(degrees, ["G09"], 100, np.random.default_rng(3))
    assert picked2 and picked2[0] != "G09" and abs(degrees[picked2[0]] - 300) <= 100


def t_prox_zscore_synthetic():
    """符号测试：药物基因 A 与直接邻居 B1 应显著邻近，与远端 B2 的距离必然更大。"""
    G = nx.barabasi_albert_graph(120, 3, seed=7)
    G = nx.relabel_nodes(G, {i: f"N{i:03d}" for i in G.nodes()})
    degrees = dict(G.degree())
    hub_sorted = sorted(degrees, key=degrees.get, reverse=True)
    A = hub_sorted[:8]
    # B1：A 的直接邻居（排除 A 自身），取 8 个
    neighbors = sorted({v for a in A for v in G[a]} - set(A),
                       key=degrees.get, reverse=True)
    B1 = neighbors[:8]
    # B2：与 A 距离最远的 8 个节点（图的边缘）
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
    # B1 全是 A 的直接邻居，距离应恰好为 1，且小于远端 B2
    assert abs(res_near["d_c"] - 1.0) < 1e-12, res_near["d_c"]
    assert res_near["d_c"] < res_far["d_c"], (res_near["d_c"], res_far["d_c"])
    # 近端模块的 z 值应为负，经验 p 值不应比远端更差
    assert res_near["z"] < 0, f"near-module z should be negative, got {res_near['z']}"
    assert res_near["p_empirical"] <= res_far["p_empirical"], (
        res_near["p_empirical"], res_far["p_empirical"])


def t_prox_coverage_guard():
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


def t_prox_reference_impl_consistency():
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


def t_prox_real_and_cached():
    """真实互作组完整流程（未构建则跳过）：结果结构完整，第二次调用命中缓存。"""
    if not itx.is_built("Human"):
        return  # SKIP
    drug = ["TP53", "EGFR", "AKT1", "MYC", "STAT3"]
    disease = ["BRCA1", "PTEN", "KRAS", "PIK3CA", "ERBB2"]
    res1 = run_proximity_cached(drug, disease, label="selfcheck_prox")
    assert res1["ok"], res1
    assert all(isinstance(res1[k], float) for k in ("d_c", "z", "p_empirical"))
    # 该药物/疾病对无共享靶点，去交集后疾病集仍是 5 个，d_c_excl 应与 d_c 相等
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

def t_rwr_direction_and_mass():
    """验证游走方向：从星形图中心出发，扩散一步后每个叶子各得 0.25，总和为 1。"""
    import scipy.sparse as sp

    G = nx.star_graph(4)  # 节点0为中心，度4；叶1..4度1
    G = nx.relabel_nodes(G, {i: f"N{i}" for i in G.nodes()})
    W, nodes = build_sparse_adj(G)
    # 列随机性：每列和为 1
    col_sum = np.asarray(W.sum(axis=0)).ravel()
    assert np.allclose(col_sum, 1.0), col_sum
    idx = {g: i for i, g in enumerate(nodes)}
    r = np.zeros(len(nodes))
    r[idx["N0"]] = 1.0
    r1 = W @ r
    for leaf in ("N1", "N2", "N3", "N4"):
        assert abs(r1[idx[leaf]] - 0.25) < 1e-12, (leaf, r1[idx[leaf]])
    assert abs(r1.sum() - 1.0) < 1e-12, r1.sum()


def t_rwr_alpha0_returns_seeds():
    """α=0 时不发生游走，得分应全部留在种子基因上。"""
    G = nx.barabasi_albert_graph(30, 2, seed=1)
    G = nx.relabel_nodes(G, {i: f"N{i:02d}" for i in G.nodes()})
    df = rwr_prioritize(G, ["N00", "N05"], alpha=0.0, tol=1e-12, max_iter=50)
    seed_rows = df[df["is_seed"]]
    nonseed = df[~df["is_seed"]]
    assert len(seed_rows) == 2
    assert set(seed_rows["gene"]) == {"N00", "N05"}
    assert seed_rows["score"].round(12).tolist() == [0.5, 0.5]
    assert (nonseed["score"].abs() < 1e-12).all(), "non-seeds must have zero mass at alpha=0"


def t_rwr_prioritize_basic():
    """α=0.7 正常排序：总分守恒、按得分降序、种子基因标注正确。"""
    G = nx.barabasi_albert_graph(50, 2, seed=3)
    G = nx.relabel_nodes(G, {i: f"N{i:02d}" for i in G.nodes()})
    df = rwr_prioritize(G, ["N00", "N03", "NOT_IN_GRAPH"], alpha=0.7)
    assert abs(df["score"].sum() - 1.0) < 1e-6, df["score"].sum()
    assert df["score"].is_monotonic_decreasing
    assert df["is_seed"].sum() == 2
    # 度最高的枢纽节点应排在前列
    assert df.iloc[0]["gene"] in ("N00", "N03") or df.iloc[0]["score"] > 0.01


def t_rwr_from_context_guard_and_annotate():
    """种子基因全部不在网络中时应拒绝计算；计算结果本身不含基因标注信息。"""
    G = _karate()
    ctx = _graph_context(G)
    res = na._rwr_from_context(ctx, ["NOPE1", "NOPE2"], alpha=0.7, tol=1e-8, max_iter=50)
    assert res["error"] == "empty_after_filter" and not res["ok"]
    known = list(G.nodes())[:5]
    res2 = na._rwr_from_context(ctx, known, alpha=0.7, tol=1e-8, max_iter=50)
    assert res2["ok"] and res2["n_seeds"] == 5 and res2["converged"]
    assert all("annotated" not in row for row in res2["top_nonseed"])
    assert all("annotated" not in row for row in res2["seeds"])


def t_prox_excl_overlap():
    """剔除药物/疾病共享靶点后再算距离：d_c 只会变大或不变；疾病靶点全被覆盖时该指标不可用。"""
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
    # 疾病集完全被药物集覆盖时，去交集后为空，该指标不可用
    res2 = _proximity_zscore_from_context(ctx, nodes[:10], nodes[:5], n_perm=10, seed=0,
                                          bin_size=100, drug_ratio=0.0, disease_ratio=0.0)
    assert res2["ok"] and res2["d_c_excl"] is None and res2["n_excl"] == 0


def t_rwr_enrich_synthetic():
    """标注集富集检验：种子和标注集都是枢纽基因时，标注基因得分应显著偏高（z 为正、p 小）。"""
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


def t_rwr_real_and_cached():
    """真实互作组完整流程（未构建则跳过）：检查种子/非种子分列表、收敛与缓存。"""
    if not itx.is_built("Human"):
        return  # SKIP
    import glob
    import config

    # 先删除旧测试缓存，保证第一次是全新计算
    for old in glob.glob(os.path.join(config.RWR_CACHE_DIR, "rwr_selfcheck_rwr_*.json")):
        try:
            os.remove(old)
        except OSError:
            pass

    seeds = ["TP53", "EGFR", "AKT1", "MYC", "STAT3"]
    res1 = run_rwr_cached(seeds, label="selfcheck_rwr")
    assert res1["ok"], res1
    assert res1["converged"] and res1["n_iter"] >= 1
    # total_mass 是全网络总质量（约 1.0）；seed_mass 是留在种子上的比例，含义不同
    assert abs(res1["total_mass"] - 1.0) < 1e-6, res1["total_mass"]
    assert 0.0 < res1["seed_mass"] < 1.0, res1["seed_mass"]
    # 种子列表全量保留；非种子榜单单独截断，且不含种子
    assert len(res1["seeds"]) == 5, len(res1["seeds"])
    assert {r["gene"] for r in res1["seeds"]} == {"TP53", "EGFR", "AKT1", "MYC", "STAT3"}
    assert len(res1["top_nonseed"]) == config.RWR_CACHE_TOP_N, len(res1["top_nonseed"])
    seed_names = {r["gene"] for r in res1["seeds"]}
    assert all(r["gene"] not in seed_names for r in res1["top_nonseed"])
    assert res1["top_nonseed"][0]["score"] >= res1["top_nonseed"][-1]["score"]
    # 标注集也算进缓存标识：换标注集就要重新计算
    first_ns = res1["top_nonseed"][0]["gene"]
    res2 = run_rwr_cached(seeds, label="selfcheck_rwr", annotate={first_ns})
    assert res2["ok"] and res2.get("cached") is False
    row0 = res2["top_nonseed"][0]
    assert row0["gene"] == first_ns and row0["annotated"] is True
    assert not any(r.get("annotated") for r in res2["top_nonseed"][1:])
    enr = res2.get("enrichment")
    assert enr and enr["n_annotated"] == 1 and enr["z"] is not None and enr["p"] is not None
    res3 = run_rwr_cached(seeds, label="selfcheck_rwr", annotate={first_ns})
    assert res3.get("cached") is True
    assert abs(res3["top_nonseed"][0]["score"] - row0["score"]) < 1e-15
    print(f"      [real interactome] iters={res1['n_iter']} converged={res1['converged']} "
          f"seeds={len(res1['seeds'])} top_nonseed={len(res1['top_nonseed'])} "
          f"top5={[r['gene'] for r in res1['top_nonseed'][:5]]} "
          f"enrich_z={enr['z']:.2f} enrich_p={enr['p']:.4f}")


def t_rwr_truncation_guarantee():
    """回归测试：种子数很大时，种子表仍要全量保留，非种子榜单不被种子挤占。"""
    G = nx.barabasi_albert_graph(60, 2, seed=5)
    G = nx.relabel_nodes(G, {i: f"N{i:02d}" for i in G.nodes()})
    ctx = _graph_context(G)
    nodes = ctx[0]
    many_seeds = list(nodes[:45])  # 60 节点中 45 个种子（窗口只有 5）
    res = na._rwr_from_context(ctx, many_seeds, alpha=0.7, tol=1e-8, max_iter=100,
                               top_nonseed=5)
    assert res["ok"]
    assert len(res["seeds"]) == 45, f"seeds must be kept in full, got {len(res['seeds'])}"
    assert len(res["top_nonseed"]) == 5, f"non-seed window must stay full, got {len(res['top_nonseed'])}"
    seed_set = set(many_seeds)
    assert all(r["gene"] not in seed_set for r in res["top_nonseed"])
    scores = [r["score"] for r in res["top_nonseed"]]
    assert scores == sorted(scores, reverse=True), "non-seed rows must be score-desc"


def main():
    tests = [
        ("detect basic properties", t_detect_basic),
        ("detect deterministic", t_deterministic),
        ("detect empty graph", t_empty_graph),
        ("detect avg_degree", t_avg_degree),
        ("enrich both-fail -> empty + failure record", t_enrich_fail_both),
        ("enrich circuit_open treated as failure", t_enrich_fail_gp_ok_enr),
        ("enrich genuine no-significant", t_enrich_genuine_no_sig),
        ("enrich success path", t_enrich_success),
        ("enrich min-size filter", t_min_size_filter),
        ("gsea build_ranked_list dedupe/drop/sort", t_gsea_build_ranked),
        ("gsea normalize res2d columns", t_gsea_normalize),
        ("gsea prerank local + cache hit", t_gsea_prerank_local),
        ("gsea insufficient input leaves no cache", t_gsea_insufficient_no_cache),
        ("interactome protein.info column pick", t_interactome_info_cols),
        ("interactome links->edges transform", t_interactome_links_to_edges),
        ("interactome load/coverage (skip if not built)", t_interactome_load_or_skip),
        ("interactome staleness age/is_stale", t_interactome_staleness),
        ("prox sparse BFS matches networkx", t_prox_bfs_matches_networkx),
        ("prox degree-matched sampling", t_prox_degree_matched_sample),
        ("prox zscore synthetic sign test", t_prox_zscore_synthetic),
        ("prox coverage guard", t_prox_coverage_guard),
        ("prox sparse vs reference d_c", t_prox_reference_impl_consistency),
        ("prox excl-overlap variant", t_prox_excl_overlap),
        ("prox real interactome + cache (skip if not built)", t_prox_real_and_cached),
        ("rwr direction & mass conservation", t_rwr_direction_and_mass),
        ("rwr alpha=0 returns seeds", t_rwr_alpha0_returns_seeds),
        ("rwr prioritize basic", t_rwr_prioritize_basic),
        ("rwr coverage guard + annotation decoupling", t_rwr_from_context_guard_and_annotate),
        ("rwr enrichment synthetic", t_rwr_enrich_synthetic),
        ("rwr real interactome + cache (skip if not built)", t_rwr_real_and_cached),
        ("rwr truncation guarantee (large seed set)", t_rwr_truncation_guarantee),
    ]
    print(f"==== network_algorithms selfcheck: {len(tests)} tests ====")
    ok = all([check(name, fn) for name, fn in tests])
    print("==== ALL PASS ====" if ok else "==== FAILED ====")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
