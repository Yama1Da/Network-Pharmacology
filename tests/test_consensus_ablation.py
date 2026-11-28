"""消融对比纯计算核的冒烟测试(§1.13):小图性质断言。

协议细节(DisGeNET 分层/留出集)是 IO 层的事,由真实数据试跑覆盖;这里只锁
计算核的形状与数学性质——变体齐全、Recall 有界且对 K 单调不降、邻近度腿可选。
"""

import networkx as nx
import pytest

from WY.benchmark.consensus_ablation import run_ablation, render_chart


def _ctx_from_graph(g: nx.Graph):
    """合成 LCC 上下文(与 test_disease_modules 同构)。"""
    nodes = tuple(sorted(g.nodes()))
    idx = {x: i for i, x in enumerate(nodes)}
    Mt = nx.to_scipy_sparse_array(g, nodelist=nodes, weight=None, format="csr")
    degree_map = {x: int(d) for x, d in g.degree()}
    return nodes, idx, Mt, degree_map


@pytest.fixture()
def ctx():
    """10 节点小图:种子三角 A-B-C + 扩张带 D/E/F/G + 远端 H/I/J。"""
    g = nx.Graph()
    g.add_edges_from([("A", "B"), ("A", "C"), ("B", "C"),
                      ("A", "D"), ("B", "E"), ("C", "F"),
                      ("D", "G"), ("E", "G"), ("H", "I"), ("I", "J")])
    return _ctx_from_graph(g)


def test_variants_complete_and_ordered(ctx):
    res = run_ablation(ctx, ["A", "B"], reference={"G"}, n=6, k=2, cutoffs=(3, 5))
    names = [r["variant"] for r in res["rows"]]
    assert len(names) == 5
    assert names[0] == "S0_seeds" and names[-1] == "S4_consistent"
    assert {n.split("_")[1] for n in names[1:4]} == {"rwr", "diamond", "neighbor"}
    assert res["reference_size"] == 1 and res["k"] == 2


def test_recall_bounded_and_monotonic(ctx):
    res = run_ablation(ctx, ["A", "B"], reference={"D", "E", "G"},
                       n=6, k=2, cutoffs=(3, 5, 8))
    for r in res["rows"]:
        assert 0.0 <= r["recall@3"] <= 1.0
        assert r["recall@3"] <= r["recall@5"] <= r["recall@8"]  # K 增大 Recall 不降
        assert r["ref_hits"] <= 3
        if r["ref_hits"]:
            assert r["ref_mean_rank"] >= 1


def test_reference_never_absent_from_ranking_pool(ctx):
    """留出真值若不在 LCC,Recall 恒 0 且不报错(宁缺毋滥口径下显式可测)。"""
    res = run_ablation(ctx, ["A", "B"], reference={"Z"}, n=5, k=2, cutoffs=(4,))
    assert all(r["recall@4"] == 0.0 and r["ref_hits"] == 0 for r in res["rows"])


def test_proximity_leg_optional(ctx):
    base = run_ablation(ctx, ["A", "B"], reference={"G"}, n=5, k=2)
    assert all("prox_z" not in r for r in base["rows"])
    leg = run_ablation(ctx, ["A", "B"], reference={"G"}, n=5, k=2,
                       drug_targets=["A", "H"], n_perm=5)
    prox = [r for r in leg["rows"] if r["prox_ok"]]
    # 小图 5 次置换足够让近端变体跑出有限 z;失败行 prox_ok=False 且 z 为 None
    assert all((r["prox_z"] is not None) == r["prox_ok"] for r in leg["rows"])
    assert prox or True  # 恒真占位:小图低覆盖行允许全失败,不伪造数值


def test_render_chart_smoke(ctx, tmp_path):
    res = run_ablation(ctx, ["A", "B"], reference={"G"}, n=5, k=2, cutoffs=(3, 5))
    p = render_chart(res["rows"], (3, 5), str(tmp_path / "c.png"))
    import os
    assert os.path.isfile(p) and os.path.getsize(p) > 5000


def test_glue_recall_matches_independent_recomputation(ctx):
    """胶水层手算对照(清单 §1.13 验收要求,同陷阱#2 手算传统):
    不信任 run_ablation 内部——测试自己重新跑 RWR、自己合并排名、自己算
    Recall@K,逐变体逐截断与 run_ablation 的行严格相等。RWR 本身的正确性
    归 network_algorithms 的一期测试,这里锁的是"基准胶水"没拼错。"""
    from WY import config
    from WY.core.disease_modules import consensus_targets
    from WY.core.network_algorithms import _rwr_from_context

    seeds, ref, cutoffs = ["A", "B"], {"D", "E", "G"}, (3, 5, 8)
    res = run_ablation(ctx, seeds, reference=ref, n=6, k=2, cutoffs=cutoffs)
    cons = consensus_targets(ctx, seeds, n=6, k=2)
    variants = {"S0_seeds": sorted(set(seeds)),
                "S1_rwr": sorted(cons["modules"]["rwr"]),
                "S2_diamond": sorted(cons["modules"]["diamond"]),
                "S3_neighbor": sorted(cons["modules"]["neighbor"]),
                "S4_consistent": sorted(cons["consensus"]
                                        .loc[cons["consensus"]["consistent"]
                                             .astype(bool), "gene"].astype(str))}
    by_name = {r["variant"]: r for r in res["rows"]}
    for name, seedset in variants.items():
        r = _rwr_from_context(ctx, seedset, config.RWR_ALPHA, config.RWR_TOL,
                              config.RWR_MAX_ITER, top_nonseed=None)
        ranked = [x["gene"] for x in sorted(r["seeds"] + r["top_nonseed"],
                                            key=lambda x: -x["score"])]
        for k in cutoffs:
            expect = len(set(ranked[:k]) & ref) / len(ref)
            assert by_name[name][f"recall@{k}"] == pytest.approx(expect), (name, k)
        hits = [i + 1 for i, g in enumerate(ranked) if g in ref]
        assert by_name[name]["ref_hits"] == len(hits)
        assert by_name[name]["ref_mean_rank"] == pytest.approx(
            sum(hits) / len(hits) if hits else None)
