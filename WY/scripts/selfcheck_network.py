#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""网络分析功能自检脚本（纯本地运行，不发网络请求）。

运行方式（在 WY 目录下）：
    python scripts/selfcheck_network.py

检查内容：Louvain 模块检测、模块富集、GSEA、背景互作组、网络邻近度、RWR 靶点排序。
每项检查打印 PASS 或 FAIL；全部通过时退出码为 0，有失败则为 1。
"""

import glob
import logging
import os
import shutil
import sys
import tempfile
import time


from WY.scripts.selfcheck_common import patch_attr, run

import networkx as nx
import numpy as np
import pandas as pd
import requests

from WY import config
import WY.core.network_algorithms as na
from WY.core import enrichment as en
from WY.core import interactome as itx
from WY.core.caching import save_json_file
from WY.core.enrichment import (
    TRUSTED_REASONS,
    _gsea_cache_file,
    _normalize_enrichment_gsea,
    build_ranked_list,
    fetch_enrichment_gsea_prerank,
    render_drug_target_pathway_plot,
)
from WY.core.network_algorithms import (
    _multi_source_bfs,
    _proximity_zscore_from_context,
    _rwr_power_iteration,
    degree_matched_sample,
    run_proximity_cached,
    run_rwr_cached,
)
from WY.core.ppi import (
    compute_mcc_scores,
    mark_hub_genes,
    pick_top_mcc_genes,
    render_core_target_plot,
)


def _agg_backend():
    """自检为无界面运行,绘图测试前统一切到 Agg 后端(幂等)。"""
    import matplotlib

    matplotlib.use("Agg", force=True)


def _karate() -> nx.Graph:
    G = nx.karate_club_graph()
    return nx.relabel_nodes(G, {i: f"G{i:02d}" for i in G.nodes()})


# ---------- 参考实现（原 core.network_algorithms 中的纯本地实现，仅供自检交叉验证） ----------


def network_proximity(graph: "nx.Graph", set_a: list, set_b: list) -> float:
    """参考实现：两组蛋白的网络平均距离 d_c(A,B) = mean_{b∈B} min_{a∈A} d(a,b)。
    仅供小图与自检交叉验证;图中不存在的基因被忽略,B 中有从 A 不可达的节点时抛错。
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
    """0-1 邻接 → 列随机转移矩阵 W（W[j,i] = P(i→j)，每列和为 1）。
    度为 0 的节点会永久截留游走概率,故图必须连通（如限制在 LCC 内）。
    """
    import networkx as nx
    import numpy as np
    import scipy.sparse as sp

    nodes = sorted(graph.nodes())
    adj = nx.to_scipy_sparse_array(graph, nodelist=nodes, weight=None, format="csr")
    col_deg = np.asarray(adj.sum(axis=0)).ravel()
    if (col_deg == 0).any():
        raise ValueError("zero-degree nodes in graph; restrict to a connected component")
    # 右乘对角阵 = 按列缩放：W[:, i] = adj[:, i] / deg(i)
    W = adj @ sp.diags(1.0 / col_deg)
    return W, nodes


def rwr_prioritize(graph: "nx.Graph", seed_genes: list, alpha: float = 0.7,
                   tol: float = 1e-8, max_iter: int = 200) -> pd.DataFrame:
    """参考入口（直接吃 networkx 图）：返回全节点 RWR 打分表,列 gene/score/is_seed 按 score 降序。
    图需连通（建议 LCC）;不在图内的种子被忽略。
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

# 这些返回状态视为可信（即 core.enrichment.TRUSTED_REASONS）：真正成功，或确实无显著项


def _fake_fetch(reason, df):
    """假的富集服务函数工厂（g:Profiler/Enrichr 共用,调用方均按关键字传参）：reason=None 表示成功，否则按给定 reason 记录状态。"""
    def f(genes, libs=None, organism="hsapiens", sources=("GO:BP",), pvalue=0.05, label=None, status=None):
        if status is not None:
            if reason is None:
                status.update(ok=True, reason="success")
            else:
                status.update(ok=reason in TRUSTED_REASONS, reason=reason)
        return df.copy()
    return f


def _patch(gp_spec, enr_spec):
    """临时把两个富集服务函数换成假函数。参数各为 (reason, df) 元组,或直接给一个现成的假函数。"""
    def _make(spec):
        return spec if callable(spec) else _fake_fetch(*spec)

    orig = (na.fetch_enrichment_gprofiler, na.fetch_enrichment_enrichr)
    na.fetch_enrichment_gprofiler = _make(gp_spec)
    na.fetch_enrichment_enrichr = _make(enr_spec)
    return orig


def _restore(orig):
    na.fetch_enrichment_gprofiler, na.fetch_enrichment_enrichr = orig


def _module5(genes=("A", "B", "C", "D", "E"), module_id=1):
    """测试用的 5 基因模块 dict（多个富集测试共用同一形状的数据）。"""
    return {"module_id": module_id, "genes": list(genes), "size": len(genes), "avg_degree": 2.0}


def _purge(pattern):
    """按 glob 模式删除旧缓存文件,文件不存在时为无害空操作。"""
    for old in glob.glob(pattern):
        try:
            os.remove(old)
        except OSError:
            pass


def t_enrich_fail_both():
    orig = _patch(("request_failed", pd.DataFrame()), ("request_failed", pd.DataFrame()))
    try:
        modules = [_module5()]
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
        modules = [_module5()]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
        assert len(failures) == 1, "circuit_open + empty enrichr must be a failure"
        assert not no_sig
    finally:
        _restore(orig)


def t_enrich_genuine_no_sig():
    # 两条链路都成功但确实无显著项 → 不算失败
    orig = _patch(("no_significant", pd.DataFrame()), ("no_significant", pd.DataFrame()))
    try:
        modules = [_module5()]
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

    orig = _patch(gp, enr)
    try:
        modules = [
            _module5(),
            _module5(genes=("X", "Y", "Z", "W", "V"), module_id=2),
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

    orig = _patch(counting_gp, lambda genes, **kw: pd.DataFrame())
    try:
        modules = [{"module_id": 1, "genes": ["A", "B"], "size": 2, "avg_degree": 1.0}]
        long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens",
                                                      min_module_size=5)
        assert calls["n"] == 0, "module below min size must not be enriched"
        assert long_df.empty and not failures and not no_sig
    finally:
        _restore(orig)


# ---------- 富集数据源熔断测试（标记文件隔离到临时目录,不碰真实缓存） ----------


def _circuit_sandbox() -> str:
    """把富集缓存目录切到临时目录（熔断标记随目录隔离）,返回原目录供恢复。"""
    orig = config.ENRICHMENT_CACHE_DIR
    config.ENRICHMENT_CACHE_DIR = tempfile.mkdtemp(prefix="enrichcircuit_")
    return orig


def _circuit_restore(orig: str) -> None:
    shutil.rmtree(config.ENRICHMENT_CACHE_DIR, ignore_errors=True)
    config.ENRICHMENT_CACHE_DIR = orig


def _expire_circuit(circuit, stage: int) -> None:
    """写入一条"已过冷却期"的失败标记,让下一次请求作为探测放行。"""
    save_json_file(circuit._file(),
                   {"failed_at": time.time() - circuit.cooldown(stage) - 1,
                    "stage": stage})


def t_circuit_ladder_advance():
    # 阶梯语义:常规失败保持第 1 档;冷却期内的失败(并发探测撞车)不推进;
    # 到期后的探测失败推进一档;冷却时长按阶梯加长且最后一档封顶。
    orig = _circuit_sandbox()
    try:
        en._gp_circuit.trip(advance=False)
        en._gp_circuit.trip(advance=False)  # 常规失败不推进档位
        assert en._gp_circuit.read()["stage"] == 1
        en._gp_circuit.trip(advance=True)   # 冷却期内:探测撞车去抖,不推进
        assert en._gp_circuit.read()["stage"] == 1
        # 熔断冷却期内,重试提示应给出大致分钟数(向上取整,不足 1 分钟按 1 分钟)
        assert en.enrichment_retry_hint_minutes() == 1
        _expire_circuit(en._gp_circuit, 1)
        en._gp_circuit.trip(advance=True)   # 确实到期后的失败才推进
        assert en._gp_circuit.read()["stage"] == 2
        assert en._gp_circuit.state() == "open"
        cooldowns = [en._gp_circuit.cooldown(s) for s in (1, 2, 3, 4, 9)]
        assert cooldowns[0] < cooldowns[1] < cooldowns[2] <= cooldowns[3] == cooldowns[4], cooldowns
    finally:
        _circuit_restore(orig)


def t_circuit_probe_empty_result_resets():
    # 冷却到期后的探测拿到空结果（服务可用、恰好无显著项）也必须完全复位,不再保持熔断
    orig = _circuit_sandbox()
    try:
        en._gp_circuit.trip(advance=False)
        _expire_circuit(en._gp_circuit, 1)
        with patch_attr(en, "_profile_gprofiler", lambda *a, **k: pd.DataFrame()):
            st: dict = {}
            df = en.fetch_enrichment_gprofiler(["TP53", "TNF"], status=st)
            assert df.empty and st.get("reason") == "no_significant", st
            assert en._gp_circuit.state() == "closed", \
                "empty probe result must reset the breaker"
        assert en.enrichment_retry_hint_minutes() is None  # 复位后无冷却信息
    finally:
        _circuit_restore(orig)


def t_circuit_probe_fail_advances():
    # 冷却到期后的探测失败:status 记 probe_failed,冷却阶梯推进一档;
    # 冷却期内后续请求直接跳过(circuit_open),不再发起
    orig = _circuit_sandbox()
    try:
        en._gp_circuit.trip(advance=False)
        _expire_circuit(en._gp_circuit, 1)

        def _down(*a, **k):
            raise ConnectionError("g:Profiler still down")

        with patch_attr(en, "_profile_gprofiler", _down):
            st: dict = {}
            df = en.fetch_enrichment_gprofiler(["TP53", "TNF"], status=st)
            assert df.empty and st.get("reason") == "probe_failed", st
            assert en._gp_circuit.read()["stage"] == 2
            assert en._gp_circuit.state() == "open"
            st2: dict = {}
            df2 = en.fetch_enrichment_gprofiler(["TP53", "TNF"], status=st2)
            assert df2.empty and st2.get("reason") == "circuit_open", st2
    finally:
        _circuit_restore(orig)


def t_enrichr_breaker_single_attempt():
    # Enrichr 兜底源熔断:单次尝试失败即跳闸(请求内不重试);冷却期内跳过;探测成功复位
    orig = _circuit_sandbox()
    try:
        calls = {"n": 0}

        def _down(*a, **k):
            calls["n"] += 1
            raise ConnectionError("Enrichr down")

        with patch_attr(en, "_query_enrichr", _down):
            st: dict = {}
            df = en.fetch_enrichment_enrichr(["TP53", "TNF"],
                                             libs=["KEGG_2021_Human"], status=st)
            assert df.empty and st.get("reason") == "request_failed", st
            assert calls["n"] == 1, "fallback source must not retry in-request"
            assert en._enrichr_circuit.read()["stage"] == 1
            st2: dict = {}
            en.fetch_enrichment_enrichr(["TP53", "TNF"],
                                        libs=["KEGG_2021_Human"], status=st2)
            assert st2.get("reason") == "circuit_open", st2
            assert calls["n"] == 1, "breaker open must skip the request"
        _expire_circuit(en._enrichr_circuit, 1)
        with patch_attr(en, "_query_enrichr", lambda *a, **k: pd.DataFrame()):
            st3: dict = {}
            df3 = en.fetch_enrichment_enrichr(["TP53", "TNF"],
                                              libs=["KEGG_2021_Human"], status=st3)
            assert df3.empty and st3.get("reason") == "no_significant", st3
            assert en._enrichr_circuit.state() == "closed"
    finally:
        _circuit_restore(orig)


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
    libs = {"SELFTEST_LIB": _tiny_gmt()}
    label = "selfcheck_gsea"
    # 先删除旧测试缓存，避免直接命中
    _purge(os.path.join(config.ENRICHMENT_CACHE_DIR, f"GSEA_{label}_*.csv"))

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
    assert config.INTERACTOME_CONFIDENCE == int(config.PPI_CONFIDENCE * 1000)


class _FakeResp:
    """_download_gz 打桩用的最小响应对象:支持 with 上下文与流式 iter_content。"""

    def __init__(self, status, headers=None, body=b""):
        self.status_code = status
        self.headers = headers or {}
        self._body = body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def close(self):
        pass

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"http {self.status_code}")

    def iter_content(self, chunk_size=1 << 20):
        for i in range(0, len(self._body), chunk_size):
            yield self._body[i:i + chunk_size]


def _fake_string_get(payload, fail_ranges=False):
    """伪装 STRING 下载源:带 Range 头按区间切片回 206,Range 探测回总大小,
    不带 Range 回 200 全量(模拟不支持分段的源);fail_ranges=True 时区间请求一律失败。"""
    def fake_get(url, headers=None, stream=True, timeout=None, **kw):
        rng = (headers or {}).get("Range", "")
        if rng == "bytes=0-0":
            return _FakeResp(206, {"Content-Range": f"bytes 0-0/{len(payload)}"},
                             payload[:1])
        if rng.startswith("bytes="):
            if fail_ranges:
                raise IOError("simulated segment failure")
            s, e = rng[6:].split("-")
            return _FakeResp(206, {}, payload[int(s):int(e) + 1])
        return _FakeResp(200, {}, payload)

    return fake_get


def t_interactome_ranged_download():
    """Range 可用时按段并发下载并拼接出与原文件一致的产物,无 .part 残留。"""
    payload = bytes(range(256)) * 64  # 16 KiB,把分段下限调小以便离线测试
    events = []
    with patch_attr(requests, "get", _fake_string_get(payload)), \
         patch_attr(itx, "_RANGE_MIN_SIZE", 4096), \
         patch_attr(itx, "_RANGE_SEG_MIN_SIZE", 4096):
        with tempfile.TemporaryDirectory() as td:
            dest = os.path.join(td, "links.txt.gz")
            out = itx._download_gz("/x.gz", dest,
                                   progress_cb=lambda *a: events.append(a),
                                   progress_label="download_links")
            assert out == dest
            with open(dest, "rb") as f:
                assert f.read() == payload
            leftover = [n for n in os.listdir(td) if ".part" in n]
            assert not leftover, leftover
    assert any(ev[0] == "download_links" for ev in events), events


def t_interactome_fallback_notice():
    """大文件的分段在重试轮次耗尽后仍失败:回退单流,文件完整,且只发一次 fallback 提醒。"""
    payload = bytes(range(256)) * 64  # 16 KiB
    events = []
    with patch_attr(requests, "get", _fake_string_get(payload, fail_ranges=True)), \
         patch_attr(itx, "_RANGE_MIN_SIZE", 4096), \
         patch_attr(itx, "_RANGE_SEG_MIN_SIZE", 4096), \
         patch_attr(config, "API_MAX_RETRIES", 1), \
         patch_attr(config, "API_RETRY_DELAY", 0):
        with tempfile.TemporaryDirectory() as td:
            dest = os.path.join(td, "links.txt.gz")
            itx._download_gz("/x.gz", dest, progress_cb=lambda *a: events.append(a))
            with open(dest, "rb") as f:
                assert f.read() == payload
    assert sum(1 for ev in events if ev[0] == "fallback") == 1, events


def t_interactome_small_file_no_notice():
    """小文件(低于分段阈值)走单流属正常路径,不应触发 fallback 误报。"""
    payload = b"tiny" * 512  # 2 KiB

    def no_range_get(url, headers=None, stream=True, timeout=None, **kw):
        return _FakeResp(200, {}, payload)

    events = []
    with patch_attr(requests, "get", no_range_get), \
         patch_attr(itx, "_RANGE_MIN_SIZE", 4096), \
         patch_attr(itx, "_RANGE_SEG_MIN_SIZE", 4096):
        with tempfile.TemporaryDirectory() as td:
            dest = os.path.join(td, "info.txt.gz")
            itx._download_gz("/x.gz", dest, progress_cb=lambda *a: events.append(a))
            with open(dest, "rb") as f:
                assert f.read() == payload
    assert not any(ev[0] == "fallback" for ev in events), events


def t_interactome_segment_resume():
    """源每次只给出区间的一小段就提前收流:靠段内续传补完,全程不应触发 fallback。"""
    payload = bytes(range(256)) * 64  # 16 KiB
    calls = {"range": 0}

    def stingy_get(url, headers=None, stream=True, timeout=None, **kw):
        rng = (headers or {}).get("Range", "")
        if rng == "bytes=0-0":
            return _FakeResp(206, {"Content-Range": f"bytes 0-0/{len(payload)}"},
                             payload[:1])
        if rng.startswith("bytes="):
            s, e = (int(x) for x in rng[6:].split("-"))
            calls["range"] += 1
            cut = min(s + 512, e + 1)  # 每次只给 512 字节就提前收流
            return _FakeResp(206, {}, payload[s:cut])
        return _FakeResp(200, {}, payload)

    events = []
    with patch_attr(requests, "get", stingy_get), \
         patch_attr(itx, "_RANGE_MIN_SIZE", 4096), \
         patch_attr(itx, "_RANGE_SEG_MIN_SIZE", 4096), \
         patch_attr(config, "API_MAX_RETRIES", 64), \
         patch_attr(config, "API_RETRY_DELAY", 0):
        with tempfile.TemporaryDirectory() as td:
            dest = os.path.join(td, "links.txt.gz")
            itx._download_gz("/x.gz", dest, progress_cb=lambda *a: events.append(a))
            with open(dest, "rb") as f:
                assert f.read() == payload
    assert not any(ev[0] == "fallback" for ev in events), events
    assert calls["range"] > 4, calls  # 每段都要多次续传才能补完


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
    # 先删除旧测试缓存，保证第一次是全新计算
    _purge(os.path.join(config.RWR_CACHE_DIR, "rwr_selfcheck_rwr_*.json"))

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


def t_mcc_known_values():
    """MCC 在已知小图上的取值:全三角每点=4;路径端点=1、中点=2;星形中心=度;孤立点=0。"""
    tri = nx.Graph([("A", "B"), ("B", "C"), ("A", "C")])
    assert compute_mcc_scores(tri) == {"A": 4.0, "B": 4.0, "C": 4.0}

    path = nx.Graph([("A", "B"), ("B", "C")])
    mcc = compute_mcc_scores(path)
    assert mcc["A"] == 1.0 and mcc["C"] == 1.0 and mcc["B"] == 2.0, mcc

    star = nx.star_graph(4)
    star = nx.relabel_nodes(star, {0: "HUB", **{i: f"L{i}" for i in range(1, 5)}})
    mcc = compute_mcc_scores(star)
    assert mcc["HUB"] == 4.0 and all(mcc[f"L{i}"] == 1.0 for i in range(1, 5)), mcc

    iso = nx.Graph()
    iso.add_node("X")
    assert compute_mcc_scores(iso) == {"X": 0.0}


def t_mcc_pick_and_core_plot():
    """pick_top_mcc_genes 排序正确且排除 0 分节点;核心靶点子网络图可渲染。"""
    _agg_backend()
    import matplotlib.pyplot as plt

    G = _karate()
    mcc = compute_mcc_scores(G)
    df = pd.DataFrame({"gene": list(G.nodes()),
                       "degree": [d for _, d in G.degree()]})
    top = pick_top_mcc_genes(df, mcc, 10)
    assert len(top) == 10 and len(set(top)) == 10, top
    assert mcc[top[0]] >= mcc[top[-1]], "top list must be MCC-descending"
    scores = [mcc[g] for g in top]
    assert scores == sorted(scores, reverse=True), "top list must be MCC-sorted"

    # MCC=0(无互作)的节点不入选;全零图返回空列表
    G2 = G.copy()
    G2.add_node("ORPHAN")
    df2 = pd.DataFrame({"gene": list(G2.nodes()),
                        "degree": [d for _, d in G2.degree()]})
    mcc2 = compute_mcc_scores(G2)
    assert "ORPHAN" not in pick_top_mcc_genes(df2, mcc2, len(G2))
    assert pick_top_mcc_genes(df2, {n: 0.0 for n in G2.nodes()}, 10) == []

    fig = render_core_target_plot(G, top, mcc)
    assert fig is not None
    plt.close(fig)


def t_hub_rank_topn_and_fallback():
    """mark_hub_genes 排名制/小网络回退/MCC 兜底:
    大网络(>=2N 节点)按 method 取前 N;小网络回退阈值法;MCC 缺失/全零退 degree;
    并列按(分数降序,基因名升序)稳定排序;空表不崩溃。"""
    G = _karate()  # 34 节点 >= 2*10 → 排名制生效
    df = pd.DataFrame({"gene": list(G.nodes()),
                       "degree": [d for _, d in G.degree()],
                       "betweenness": 0.0, "closeness": 0.0})
    out = mark_hub_genes(df, top_n=10, method="degree")
    assert (out["hub_rule"] == "top10").all(), out["hub_rule"].unique()
    assert int(out["is_hub"].sum()) == 10
    assert list(out["hub_rank"]) == list(range(1, 35)), "rows must be rank-ordered"
    degs = list(out["degree"])
    assert degs == sorted(degs, reverse=True), "rows must be degree-descending"
    for i in range(len(out) - 1):  # 并列度数按基因名升序,保证可复现
        r1, r2 = out.iloc[i], out.iloc[i + 1]
        if r1["degree"] == r2["degree"]:
            assert str(r1["gene"]) <= str(r2["gene"]), "degree ties must sort by gene name"
    assert set(out[out["is_hub"]]["gene"]) == set(out.head(10)["gene"])

    # 小网络(<2N):回退阈值法——星形中心(度>中位数)入选,叶子(度=中位数)不入选
    star = nx.relabel_nodes(nx.star_graph(9),
                            {0: "CENTER", **{i: f"L{i}" for i in range(1, 10)}})
    sdf = pd.DataFrame({"gene": list(star.nodes()),
                        "degree": [d for _, d in star.degree()],
                        "betweenness": 0.0, "closeness": 0.0})
    s = mark_hub_genes(sdf, top_n=10, method="degree")
    assert (s["hub_rule"] == "threshold").all()
    assert list(s[s["is_hub"]]["gene"]) == ["CENTER"], s[s["is_hub"]]["gene"].tolist()

    # MCC 不可用(空 map)时回退 degree 排名,且不新增 mcc 列
    m0 = mark_hub_genes(df, top_n=10, method="mcc", mcc_map={})
    assert (m0["hub_rule"] == "top10").all()
    assert "mcc" not in m0.columns, "fallback to degree must not add mcc column"
    assert list(m0["gene"]) == list(out["gene"]), "fallback order must equal degree order"

    # MCC 可用时按 mcc 排名:把度数最低节点的 mcc 抬到最高,它必须升到第一
    mcc_map = {g: 0.0 for g in G.nodes()}
    mcc_map[out.iloc[-1]["gene"]] = 99.0
    m1 = mark_hub_genes(df, top_n=10, method="mcc", mcc_map=mcc_map)
    assert m1.iloc[0]["gene"] == out.iloc[-1]["gene"], "highest MCC must rank first"
    assert "mcc" in m1.columns and int(m1["is_hub"].sum()) == 10

    # 空表安全返回
    empty = mark_hub_genes(df.iloc[0:0])
    assert len(empty) == 0 and {"is_hub", "hub_rank", "hub_rule"} <= set(empty.columns)


def t_rwr_expansion_shared():
    """RWR 查询集扩展共享语义(core.analysis.expand_targets_with_rwr):成功路径的
    平局排序/基因并集/前缀后缀、失败路径的原集降级与 error 透传、游走缓存 label 不带
    _rwr{N}、未启用或空基因集的原样返回。网页 Build Report 与 CLI run_pair 共用此
    唯一定义;旧 CLI 失败分支曾引用仅成功分支定义的变量直接 NameError,本测试锁住该路径。"""
    import WY.core.analysis as ana

    fake = {"ok": True, "cached": True,
            "top_nonseed": [{"gene": "EGFR", "score": 0.5},
                            {"gene": "BRAF", "score": 0.5},
                            {"gene": "AKT1", "score": 0.9}]}
    with patch_attr(ana, "run_rwr_cached", lambda genes, label=None, **k: fake):
        genes, prefix, info = ana.expand_targets_with_rwr(
            ["TP53", "KRAS"], 2, "metformin_Type_2_diabetes_mellitus")
        # 平局按基因名稳定排序:AKT1(0.9) 第一,0.5 的 BRAF/EGFR 按名字升序取 BRAF
        assert [r["gene"] for r in info["added"]] == ["AKT1", "BRAF"], info["added"]
        assert genes == ["AKT1", "BRAF", "KRAS", "TP53"], genes
        assert prefix == "metformin_Type_2_diabetes_mellitus_rwr2", prefix
        assert info["ok"] and info["cached"] and info["error"] is None, info
        assert info["attempted"] is True
    # 失败:按原基因集返回、prefix 仍带后缀(标记"请求了扩展")、error 透传
    with patch_attr(ana, "run_rwr_cached",
                    lambda genes, label=None, **k: {"ok": False, "error": "boom"}):
        genes, prefix, info = ana.expand_targets_with_rwr(["TP53"], 30, "x_y")
        assert genes == ["TP53"] and info["ok"] is False, (genes, info)
        assert info["error"] == "boom" and info["attempted"] is True, info
        assert prefix == "x_y_rwr30", prefix
    # 游走缓存 label 不带 _rwr{N}(改 N 复用同一份游走)
    seen = {}
    with patch_attr(ana, "run_rwr_cached",
                    lambda genes, label=None, **k: seen.update(label=label)
                    or {"ok": False, "error": "probe"}):
        ana.expand_targets_with_rwr(["TP53"], 5, "a_b")
    assert seen["label"] == "expand_a_b", seen
    # 未启用(0)或空基因集:不尝试游走,原样返回
    for g, n in ((["TP53"], 0), ([], 30)):
        genes, prefix, info = ana.expand_targets_with_rwr(g, n, "a_b")
        assert info["attempted"] is False and prefix == "a_b" and genes == list(g), (g, info)


def t_drug_target_pathway_plot():
    """药物-靶点-通路网络图:KEGG 行参与构图、GO 行忽略;无 KEGG/无命中返回 None;
    g:Profiler 列表样式 / 分号 / 逗号三种基因格式都兼容;
    无药名与 drug_name=None 均按手动模式渲染。"""
    _agg_backend()
    import matplotlib.pyplot as plt

    targets = ["TP53", "VEGFA", "TNF", "IL6", "AKT1", "STAT3", "CASP3", "MYC",
               "EGFR", "PTGS2"]
    enr = pd.DataFrame([
        {"source": "KEGG", "term_id": "hsa04151", "term_name": "PI3K-Akt signaling pathway",
         "p_value": 1e-4, "adjusted_p_value": 1e-4, "intersection_size": 4,
         "query_size": 10, "genes": "['TP53', 'AKT1', 'MYC', 'EGFR']"},  # g:Profiler 格式
        {"source": "KEGG", "term_id": "hsa04066", "term_name": "HIF-1 signaling pathway",
         "p_value": 2e-4, "adjusted_p_value": 2e-4, "intersection_size": 3,
         "query_size": 10, "genes": "VEGFA;STAT3;TP53"},                 # Enrichr 格式
        {"source": "KEGG", "term_id": "hsa04917", "term_name": "Insulin signaling pathway",
         "p_value": 3e-4, "adjusted_p_value": 3e-4, "intersection_size": 2,
         "query_size": 10, "genes": "IL6,PTGS2"},                        # 逗号分隔
        {"source": "GO:BP", "term_id": "GO:0000001", "term_name": "apoptotic process",
         "p_value": 1e-5, "adjusted_p_value": 1e-5, "intersection_size": 2,
         "query_size": 10, "genes": "TP53;TNF"},
    ])
    fig = render_drug_target_pathway_plot("metformin", targets, enr)
    assert fig is not None
    plt.close(fig)
    assert render_drug_target_pathway_plot(
        "metformin", targets, enr[enr["source"] == "GO:BP"]) is None
    bad = enr.iloc[[0]].copy()
    bad["genes"] = "ZZZ1;ZZZ2"  # 命中基因与输入靶点无交集
    assert render_drug_target_pathway_plot("metformin", targets, bad) is None
    for no_drug in ("", None):  # 手动模式(含 drug_name 缺省容错)照常出图
        fig2 = render_drug_target_pathway_plot(no_drug, targets, enr)
        assert fig2 is not None
        plt.close(fig2)


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
        ("circuit ladder advance & debounce", t_circuit_ladder_advance),
        ("circuit probe empty result resets", t_circuit_probe_empty_result_resets),
        ("circuit probe fail advances ladder", t_circuit_probe_fail_advances),
        ("enrichr breaker single-attempt trip", t_enrichr_breaker_single_attempt),
        ("gsea build_ranked_list dedupe/drop/sort", t_gsea_build_ranked),
        ("gsea normalize res2d columns", t_gsea_normalize),
        ("gsea prerank local + cache hit", t_gsea_prerank_local),
        ("gsea insufficient input leaves no cache", t_gsea_insufficient_no_cache),
        ("interactome protein.info column pick", t_interactome_info_cols),
        ("interactome links->edges transform", t_interactome_links_to_edges),
        ("interactome ranged concurrent download", t_interactome_ranged_download),
        ("interactome fallback notice (segments fail)", t_interactome_fallback_notice),
        ("interactome small file no fallback notice", t_interactome_small_file_no_notice),
        ("interactome segment resume completes", t_interactome_segment_resume),
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
        ("rwr expansion shared semantics (web+cli)", t_rwr_expansion_shared),
        ("mcc known values", t_mcc_known_values),
        ("mcc pick + core target plot", t_mcc_pick_and_core_plot),
        ("hub rank top-N + small-net fallback", t_hub_rank_topn_and_fallback),
        ("drug-target-pathway plot", t_drug_target_pathway_plot),
    ]
    run(tests, "network_algorithms selfcheck")


if __name__ == "__main__":
    main()
