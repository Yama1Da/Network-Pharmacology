"""core.ppi 的图构建与拓扑指标测试(纯本地小图,不发网络)。

阶段三第 5 批新写,补自检迁移未覆盖的 build_ppi_graph / compute_ppi_metrics。
"""
import networkx as nx
import pandas as pd
import pytest

from WY.core.ppi import build_ppi_graph, compute_ppi_metrics


def _edges(*triples):
    return pd.DataFrame([{"protein_A": a, "protein_B": b, "score": s}
                         for a, b, s in triples])


def test_build_ppi_graph_keeps_isolated_nodes():
    """全部输入基因都成为节点:没有互作的基因是孤立节点,不会凭空消失。"""
    edges = _edges(("TP53", "EGFR", 900))
    g = build_ppi_graph(edges, ["TP53", "EGFR", "ORPHAN"])
    assert set(g.nodes()) == {"TP53", "EGFR", "ORPHAN"}
    assert g.number_of_edges() == 1
    assert g["TP53"]["EGFR"]["score"] == 900
    assert g["TP53"]["EGFR"]["weight"] == 0.9  # 置信度换算 0~1


def test_build_ppi_graph_dedup_keeps_highest_score():
    """同一对基因的重复边只保留置信度最高的一条(方向无关)。"""
    edges = _edges(("A", "B", 500), ("B", "A", 800), ("A", "B", 300))
    g = build_ppi_graph(edges, ["A", "B"])
    assert g.number_of_edges() == 1
    assert g["A"]["B"]["score"] == 800


def test_build_ppi_graph_drops_self_loops_and_empty():
    """基因塌缩产生的自环必须去掉;空边表 → 只有孤立节点的空图。"""
    edges = _edges(("A", "A", 900))
    g = build_ppi_graph(edges, ["A", "B"])
    assert g.number_of_edges() == 0 and set(g.nodes()) == {"A", "B"}
    empty = build_ppi_graph(pd.DataFrame(), ["X"])
    assert set(empty.nodes()) == {"X"} and empty.number_of_edges() == 0


def test_compute_ppi_metrics_sorted_and_isolated_zeros():
    """拓扑指标:按度数降序;孤立节点三项全 0;路径中点介数最高。"""
    g = nx.Graph()
    g.add_edges_from([("A", "B"), ("B", "C")])
    g.add_node("ORPHAN")
    df = compute_ppi_metrics(g)
    assert list(df["gene"]) == ["B", "A", "C", "ORPHAN"]  # 度降序(并列按节点序稳定)
    assert df.iloc[0]["degree"] == 2
    assert df.iloc[0]["betweenness"] == pytest.approx(1 / 3)  # 路径中点介数最高(归一化含 ORPHAN)
    orphan = df[df["gene"] == "ORPHAN"].iloc[0]
    assert orphan["degree"] == 0 and orphan["betweenness"] == 0.0 \
        and orphan["closeness"] == 0.0
