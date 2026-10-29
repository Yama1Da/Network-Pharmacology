"""core.ppi 的测试:MCC 打分、核心靶点挑选、枢纽基因标记。

迁移自 selfcheck_network.py(44 项中的 3 项)。
"""
import networkx as nx
import pandas as pd

from WY.core.ppi import (
    compute_mcc_scores,
    mark_hub_genes,
    pick_top_mcc_genes,
    render_core_target_plot,
)


def _karate() -> nx.Graph:
    G = nx.karate_club_graph()
    return nx.relabel_nodes(G, {i: f"G{i:02d}" for i in G.nodes()})


def test_mcc_known_values():
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


def test_mcc_pick_and_core_plot():
    """pick_top_mcc_genes 排序正确且排除 0 分节点;核心靶点子网络图可渲染。"""
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


def test_hub_rank_topn_and_fallback():
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
