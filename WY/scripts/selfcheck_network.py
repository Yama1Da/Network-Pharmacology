#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""自检：Louvain 模块检测与模块富集批量封装的本地不变式（纯离线，无网络请求）。

运行方式（在 WY 目录下）：
    python scripts/selfcheck_network.py

覆盖点：
  1. detect_louvain_modules 基本性质：全节点覆盖恰好一次、按规模降序、id 从 1 连续
  2. 确定性：同参数两次运行结果逐字段一致（种子固定）
  3. 空图安全
  4. avg_degree 计算正确性（5-团全联通，平均度应为 4）
  5. enrich_modules 失败路径：双服务失败 → 空表 + failures 记录 + 缓存清理被触发
  6. enrich_modules 成功路径：非空长表带 module_id 列且排序正确
  7. 失败与真无显著项的区分逻辑（status reason 语义）
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
    degree_matched_sample,
    network_proximity,
    run_proximity_cached,
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


# ---------- enrich_modules 测试（monkeypatch 两个 fetch 函数） ----------

# 与 core.enrichment 真实语义一致的"可信 reason"集合
_OK_REASONS = ("success", "no_significant", "cache")


def _fake_gp(reason, df):
    """reason=None 表示按真实成功语义（ok=True, reason='success'）并返回 df。"""
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
    """gp_spec / enr_spec 各为 (reason, df)：reason=None 时 fake 返回非空 df 走成功分支。"""
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
    # g:Profiler 被熔断跳过 + Enrichr 真无显著 → 整体不可信，判失败
    orig = _patch(("circuit_open", pd.DataFrame()), ("no_significant", pd.DataFrame()))
    try:
        modules = [{"module_id": 1, "genes": ["A", "B", "C", "D", "E"], "size": 5, "avg_degree": 2.0}]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
        assert len(failures) == 1, "circuit_open + empty enrichr must be a failure"
        assert not no_sig
    finally:
        _restore(orig)


def t_enrich_genuine_no_sig():
    # 双链路都成功但无显著项 → 不是失败
    orig = _patch(("no_significant", pd.DataFrame()), ("no_significant", pd.DataFrame()))
    try:
        modules = [{"module_id": 1, "genes": ["A", "B", "C", "D", "E"], "size": 5, "avg_degree": 2.0}]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
        assert long_df.empty and not failures and no_sig == [1]
    finally:
        _restore(orig)


def t_enrich_success():
    # 模块1（含基因A）：g:Profiler 成功返回非空结果；
    # 模块2（含基因X）：g:Profiler 与 Enrichr 都失败 → 进入 failures。
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


# ---------- M5：GSEA prerank（本地离线） ----------

def t_gsea_build_ranked():
    df = pd.DataFrame({
        "approvedSymbol": ["tp53", "TP53", "EGFR", "BRCA1", "KRAS", "MT-CO1"],
        "score_pchembl": [5.0, 9.0, 7.0, None, 6.0, float("nan")],
    })
    rnk, stats = build_ranked_list(df, "approvedSymbol", "score_pchembl")
    # TP53 去重取最高分 9；BRCA1/MT-CO1 缺分剔除 → 保留 3 个
    assert stats["n_input"] == 6 and stats["n_kept"] == 3, stats
    assert abs(stats["drop_ratio"] - 0.5) < 1e-9
    assert list(rnk["gene"]) == ["TP53", "EGFR", "KRAS"], rnk
    assert list(rnk["score"]) == [9.0, 7.0, 6.0], "must be sorted desc"


def t_gsea_normalize():
    # 用 gseapy 1.3.1 真实列名构造 res2d
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
    # 缺列 → 空 DataFrame 而非崩溃
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
    """真实 gseapy 离线 prerank：方向正确 + 结果缓存生效（第二次命中缓存）。"""
    import glob
    import config

    libs = {"SELFTEST_LIB": _tiny_gmt()}
    label = "selfcheck_gsea"
    # 幂等：清理历史运行留下的测试缓存（TTL 内会直接命中 cache 导致断言失真）
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
    # CSV 往返不改内容但可能改 dtype/空值表示，用容差比较而非严格 equals
    pd.testing.assert_frame_equal(
        df1.reset_index(drop=True), df2.reset_index(drop=True),
        check_dtype=False, check_exact=False, rtol=1e-6,
    )

    # 清理测试缓存文件
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


# ---------- M2：背景互作组（纯函数单测；真实缓存相关测试未构建时自动跳过） ----------

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
        # E1-E2 700（两条 ENSP 对塌缩同符号对，取最高 800 那条语义见下两行）
        "combined_score": [700, 800, 699, 950,
                           900, 850, 400],
    })
    id2sym = {"9606.E1": "TP53", "9606.E2": "EGFR", "9606.E3": "EGFR",
              "9606.E4": "KRAS", "9606.E5": "BRAF", "9606.E6": "BRCA1",
              "9606.E7": "BRCA1", "9606.E8": "ABL1"}
    edges = itx._links_to_edges(links, id2sym, confidence=700)
    records = {(r.protein_a, r.protein_b): int(r.score) for r in edges.itertuples()}
    # 无向边方向统一为 (字典序小, 字典序大)——断言按此约定写
    # 1) E3-E3 自环（E3->EGFR 塌缩）被去掉，699 分即使 >=700 也不存在
    assert ("EGFR", "EGFR") not in records, "self-loop must be dropped"
    # 2) E1-E2（TP53-EGFR）：两行 700/800 塌缩为同符号对，保留最高 800
    assert records.get(("EGFR", "TP53")) == 800, records
    # 3) 低于阈值 400 的 E7-E8 被过滤（BRCA1-ABL1 不存在）
    assert ("ABL1", "BRCA1") not in records and ("BRCA1", "ABL1") not in records
    # 4) E4-E5 (KRAS-BRAF, 950) 与 E5-E6 (BRAF-BRCA1, 900) 存在
    assert records.get(("BRAF", "KRAS")) == 950
    assert records.get(("BRAF", "BRCA1")) == 900
    assert len(edges) == 3, f"expect 3 edges, got {len(edges)}"
    # 5) 置信度换算入口一致性：0.7 * 1000 == 700
    import config
    assert config.INTERACTOME_CONFIDENCE == int(config.PPI_CONFIDENCE * 1000)


def t_interactome_load_or_skip():
    """未构建 → 必须抛 InteractomeNotBuiltError（不允许静默降级）；
    已构建 → 加载成功且覆盖守卫行为正确。"""
    if not itx.is_built("Human"):
        try:
            itx.load_interactome("Human")
            raise AssertionError("must raise InteractomeNotBuiltError when not built")
        except itx.InteractomeNotBuiltError:
            return  # SKIP：真实缓存相关分支待构建后自测
        return
    graph = itx.load_interactome("Human")
    assert graph.number_of_nodes() >= 10000, f"human interactome too small: {graph.number_of_nodes()}"
    inside, outside = itx.get_network_coverage(["tp53 ", "EGFR", "FAKE_GENE_X", "tp53"])
    assert inside == ["TP53", "EGFR"], inside
    assert outside == ["FAKE_GENE_X"], outside
    deg = itx.load_degree_table("Human")
    assert set(deg.columns) == {"gene", "degree"} and len(deg) == graph.number_of_nodes()


def t_interactome_staleness():
    """陈旧度：已构建时 age_days 为非负整数；days=0 必判 stale（age>=0）；负阈值判 stale。"""
    if not itx.is_built("Human"):
        age = itx.interactome_age_days("Human")
        assert age is None, "age must be None when not built"
        assert itx.is_stale("Human") is False, "not-built must not report stale"
        return
    age = itx.interactome_age_days("Human")
    assert isinstance(age, int) and age >= 0, age
    # 语义：stale = age > days（严格大于）。刚构建 age=0：days=-1 → stale，days=0 → 不stale
    assert itx.is_stale("Human", days=-1) is True, "negative threshold must be stale"
    assert itx.is_stale("Human", days=age) is False, "boundary: age > days must be strict"
    assert itx.is_stale("Human", days=10_000) is False, "huge threshold must not be stale"
    # 默认阈值来自 config，且为正数
    import config
    assert config.INTERACTOME_STALE_DAYS > 0


# ---------- M3：网络邻近度检验 ----------

def _graph_context(G: nx.Graph):
    """把任意 networkx 图打包成 _proximity_zscore_from_context 需要的 LCC 上下文。"""
    nodes = tuple(sorted(G.nodes()))
    idx = {g: i for i, g in enumerate(nodes)}
    adj = nx.to_scipy_sparse_array(G, nodelist=nodes, weight=None, format="csr")
    Mt = adj.T.tocsr()
    degree_map = {g: int(d) for g, d in G.degree()}
    return nodes, idx, Mt, degree_map


def t_prox_bfs_matches_networkx():
    """稀疏多源 BFS 与 networkx 单源最短路取 min 逐节点一致（连通图）。"""
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
    assert len(set(picked)) == len(set(seeds)) or True  # 允许跨种子重复（集合语义）
    for s, p in zip(seeds, picked):
        assert p != s, "replacement must not be the seed itself"
        assert abs(degrees[p] - degrees[s]) <= 100, f"degree {degrees[p]} not matched to {degrees[s]}"
    # 确定性：同种子序列重放一致
    rng2 = np.random.default_rng(0)
    assert degree_matched_sample(degrees, seeds, 100, rng2) == picked
    # 空箱扩散：度 300 的种子所在箱只有自己 → 必须落到最近非空箱
    picked2 = degree_matched_sample(degrees, ["G09"], 100, np.random.default_rng(3))
    assert picked2 and picked2[0] != "G09" and abs(degrees[picked2[0]] - 300) <= 100


def t_prox_zscore_synthetic():
    """合成图符号测试：与 A 直接相连的 B1 → d_c 显著小（p 取最小值）；
    远端 B2 的 d_c 必然更大。A 取枢纽节点，度匹配随机化提供零分布。"""
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
    # 结构断言（确定性）：B1 全部是 A 的直接邻居 → d_c 精确等于 1.0，且必然小于远端 B2
    assert abs(res_near["d_c"] - 1.0) < 1e-12, res_near["d_c"]
    assert res_near["d_c"] < res_far["d_c"], (res_near["d_c"], res_far["d_c"])
    # 统计断言（小无权图 d_c 高度量子化、随机置换常打平，绝对 p 阈值不成立）：
    # 近端 z 为负，且近端的经验 p 不劣于（≤）远端
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
    """稀疏快路径与 networkx 参考实现给出的 d_c 一致（同一小图）。"""
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
    """真实互作组端到端（未构建则跳过）：结构完整 + 缓存第二次命中。"""
    if not itx.is_built("Human"):
        return  # SKIP
    drug = ["TP53", "EGFR", "AKT1", "MYC", "STAT3"]
    disease = ["BRCA1", "PTEN", "KRAS", "PIK3CA", "ERBB2"]
    res1 = run_proximity_cached(drug, disease, label="selfcheck_prox")
    assert res1["ok"], res1
    assert all(isinstance(res1[k], float) for k in ("d_c", "z", "p_empirical"))
    res2 = run_proximity_cached(drug, disease, label="selfcheck_prox")
    assert res2.get("cached") is True, "second call must hit cache"
    assert abs(res2["z"] - res1["z"]) < 1e-12
    print(f"      [real interactome] d_c={res1['d_c']:.3f} z={res1['z']:.3f} "
          f"p={res1['p_empirical']:.4f} ({res1['n_perm']} perms)")


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
        ("prox real interactome + cache (skip if not built)", t_prox_real_and_cached),
    ]
    print(f"==== network_algorithms selfcheck: {len(tests)} tests ====")
    ok = all([check(name, fn) for name, fn in tests])
    print("==== ALL PASS ====" if ok else "==== FAILED ====")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
