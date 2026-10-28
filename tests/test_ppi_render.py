"""core.ppi 的可视化渲染与结果归档测试(纯本地,Agg 后端,零网络)。

深水区冲线批 III。
"""
import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import networkx as nx
import os
import pandas as pd
import pytest

from WY import config
from WY.core.ppi import (
    _edge_evidence_map,
    render_cytoscape_style_plot,
    render_pyvis_html,
    render_static_plot,
    render_string_style_plot,
    save_ppi_results,
)


def _fixture_graph():
    """3 节点小图 + 指标表(is_hub 列齐) + 带证据通道列的边表。"""
    g = nx.Graph()
    g.add_edges_from([("TP53", "EGFR"), ("TP53", "BRCA1")])
    metrics = pd.DataFrame([
        {"gene": "TP53", "degree": 2, "betweenness": 0.6, "closeness": 0.66,
         "is_hub": True, "hub_rule": "top10"},
        {"gene": "EGFR", "degree": 1, "betweenness": 0.0, "closeness": 0.5,
         "is_hub": False, "hub_rule": "top10"},
        {"gene": "BRCA1", "degree": 1, "betweenness": 0.0, "closeness": 0.5,
         "is_hub": False, "hub_rule": "top10"},
    ])
    edges = pd.DataFrame([
        {"protein_A": "TP53", "protein_B": "EGFR", "score": 900,
         "nscore": 0.1, "fscore": 0.2, "pscore": 0.0, "ascore": 0.0,
         "escore": 0.0, "dscore": 0.0, "tscore": 0.0},
        {"protein_A": "TP53", "protein_B": "BRCA1", "score": 800,
         "nscore": 0.0, "fscore": 0.0, "pscore": 0.0, "ascore": 0.0,
         "escore": 0.0, "dscore": 0.9, "tscore": 0.0},
    ])
    return g, metrics, edges


def test_render_pyvis_html():
    """交互网页:返回含节点/连线脚本的 HTML 文本;模块着色分支不报错。"""
    g, metrics, _ = _fixture_graph()
    html = render_pyvis_html(g, metrics)
    assert isinstance(html, str) and "vis-network" in html or "network" in html.lower()
    html2 = render_pyvis_html(g, metrics, module_map={"TP53": 1, "EGFR": 1,
                                                      "BRCA1": 0})
    assert isinstance(html2, str) and html2


def test_render_static_plot():
    """静态图(论文用):返回 matplotlib Figure,颜色条/标签渲染不报错。"""
    g, metrics, _ = _fixture_graph()
    fig = render_static_plot(g, metrics)
    assert fig is not None
    plt.close(fig)


def test_render_string_style_plot():
    """STRING 风格图:证据通道映射驱动连线配色;通道列缺失退灰色不报错。"""
    g, metrics, edges = _fixture_graph()
    fig = render_string_style_plot(g, metrics, edges)
    assert fig is not None
    plt.close(fig)
    fig2 = render_string_style_plot(g, metrics, edges.drop(
        columns=["nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"]))
    assert fig2 is not None
    plt.close(fig2)


def test_render_cytoscape_style_plot():
    """Cytoscape 风格图:确定性同心圆布局,返回 Figure。"""
    g, metrics, _ = _fixture_graph()
    fig = render_cytoscape_style_plot(g, metrics)
    assert fig is not None
    plt.close(fig)


def test_edge_evidence_map_dominant_channel():
    """证据通道映射:每条边取七个通道里最高的主导通道;全零边不进映射;缺列边表空。"""
    g, metrics, edges = _fixture_graph()
    mapping = _edge_evidence_map(edges)
    assert mapping[("EGFR", "TP53")][2] == 0.2          # nscore 主导(0-1 值)
    assert mapping[("BRCA1", "TP53")][2] == 0.9
    # 通道值按 0-999 计数返回时自动换算
    edges_scaled = edges.copy()
    edges_scaled["nscore"] = 200
    mapping2 = _edge_evidence_map(edges_scaled)
    assert mapping2[("EGFR", "TP53")][2] == pytest.approx(0.2)
    # 全零通道 → 空映射;缺部分通道列 → 用剩余通道回退(0.2 从 fscore 取);全缺 → 空
    zero = edges.copy()
    for c in ("nscore", "fscore", "pscore", "ascore", "escore", "dscore", "tscore"):
        zero[c] = 0
    assert _edge_evidence_map(zero) == {}
    partial = _edge_evidence_map(edges.drop(columns=["nscore", "dscore"]))
    assert partial[("EGFR", "TP53")][1] == "Gene fusion"  # fscore 0.2 成为主导
    no_channels = edges.drop(columns=["nscore", "fscore", "pscore", "ascore",
                                      "escore", "dscore", "tscore"])
    assert _edge_evidence_map(no_channels) == {}


def test_save_ppi_results_archives_and_replaces(tmp_path, monkeypatch):
    """结果归档:建 <prefix>_data 目录、七类文件落盘、meta 带 drug/disease;
    同 prefix 旧目录被整体重建;不传富集时对应文件跳过。"""
    monkeypatch.setattr(config, "DOWNLOADS_DIR", str(tmp_path))
    g, metrics, edges = _fixture_graph()
    import matplotlib.figure

    fig = matplotlib.figure.Figure()
    paths = save_ppi_results(
        graph=g, static_fig=fig, nodes_df=metrics, edges_df=edges,
        html_content="<html/>", prefix="metformin_T2D",
        drug="metformin", disease="T2D")
    out = paths["directory"]
    for key in ("nodes", "edges", "html", "graphml", "png", "svg"):
        assert os.path.exists(paths[key]), key
    meta = pd.read_json(os.path.join(out, "analysis_meta.json"), typ="series")
    assert meta["drug"] == "metformin" and meta["string_version"] == config.STRING_VERSION
    # 旧目录整体重建:预置一个同名旧文件会被清掉
    stale = os.path.join(out, "stale_marker.txt")
    open(stale, "w").close()
    save_ppi_results(graph=g, static_fig=fig, nodes_df=metrics, edges_df=edges,
                     html_content="<html/>", prefix="metformin_T2D")
    assert not os.path.exists(stale), "旧目录应被整体重建"

