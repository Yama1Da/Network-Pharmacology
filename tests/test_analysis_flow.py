"""core.analysis 的编排层测试:run_ppi_analysis 主流程(全部重活打桩,只跑编排骨架)。

深水区补测第二批。模式 manual(手动基因)避开仅 names 模式的 GSEA/邻近度任务;
进度协议、线程池汇合、结果字典结构、RWR 扩展分支全部覆盖。
"""
import pandas as pd
import pytest

import WY.core.analysis as ana


class _FakeProgress:
    """进度协议假件:记录 set/advance 调用序列。"""

    def __init__(self):
        self.steps = []
        self.fractions = []

    def set(self, step, state):
        self.steps.append((step, state))

    def advance(self, fraction):
        self.fractions.append(fraction)


@pytest.fixture()
def stubbed_core(monkeypatch):
    """把 run_ppi_analysis 依赖的四个重活全部打桩:STRING 边、富集、模块检测、绘图与落盘。"""
    edges = pd.DataFrame([{"protein_A": "TP53", "protein_B": "EGFR", "score": 900}])
    monkeypatch.setattr(ana, "fetch_string_ppi", lambda genes, **k: edges)
    monkeypatch.setattr(ana, "run_enrichment_task",
                        lambda genes, organism, prefix: (
                            pd.DataFrame([{"source": "KEGG", "term_name": "x"}]),
                            "gProfiler", "success"))
    monkeypatch.setattr(ana, "run_modules_task",
                        lambda graph, organism, prefix: {"modules": []})
    monkeypatch.setattr(ana, "render_enrichment_plots", lambda df: (None, None))
    monkeypatch.setattr(ana, "render_drug_target_pathway_plot", lambda *a, **k: None)
    hub_df = pd.DataFrame({"gene": ["TP53", "EGFR"], "is_hub": [True, False],
                           "hub_rule": ["top10", "top10"]})
    monkeypatch.setattr(ana, "build_main_figures",
                        lambda graph, metrics_df, edges_df, hub_method, mcc_map:
                        ("<html/>", None, {}, hub_df, None, ["TP53"], None, None))
    saved_marker = {"called": 0}

    def fake_save(**kwargs):
        saved_marker["called"] += 1
        saved_marker["kwargs"] = kwargs
        return {"files": 1}, {"modules": []}

    monkeypatch.setattr(ana, "save_analysis_outputs", fake_save)
    return saved_marker


def test_run_ppi_analysis_manual_mode(stubbed_core):
    """manual 模式全链:基因解析→建图指标→富集汇合→结果字典结构完整→落盘恰好一次。"""
    res = ana.run_ppi_analysis(mode="manual", drug_name="", disease_name="",
                               gene_text="tp53, egfr, TP53", score=0.7)
    # 基因:大写去重排序
    assert res["genes"] == ["EGFR", "TP53"]
    assert res["prefix"] == "manual"
    # 建图与拓扑是真的(小图纯本地)
    assert res["graph"].number_of_nodes() == 2
    assert not res["metrics_df"].empty
    assert res["hub_rule"] == "threshold"  # 双节点小网络回退阈值法
    # 富集任务结果原样汇入
    assert res["enrichment_source"] == "gProfiler"
    assert res["enrichment_reason"] == "success"
    # 落盘恰好一次,且拿到了图/指标/边表等全部产物
    assert stubbed_core["called"] == 1
    kwargs = stubbed_core["kwargs"]
    assert kwargs["prefix"] == "manual" and not kwargs["edges_df"].empty
    # manual 模式不跑 GSEA/邻近度
    assert res["gsea"] == {} and res["proximity"] is None
    # 参数快照按位置元组记录
    assert res["params"][0] == "manual" and res["params"][4] == 0.7


def test_run_ppi_analysis_rwr_expansion_branch(stubbed_core, monkeypatch):
    """RWR 扩展:启用且成功时,查询集并入扩展基因,expansion_df 记录种子与扩展来源。"""
    monkeypatch.setattr(ana, "expand_targets_with_rwr",
                        lambda genes, n, prefix: (
                            genes + ["KRAS"], prefix + "_rwr3",
                            {"ok": True, "cached": False,
                             "added": [{"gene": "KRAS", "score": 0.5}]}))
    res = ana.run_ppi_analysis(mode="manual", drug_name="", disease_name="",
                               gene_text="TP53", score=0.7, rwr_expand=3)
    assert "KRAS" in res["genes"]
    assert res["prefix"].endswith("_rwr3")
    exp = res["expansion_df"]
    assert list(exp["source"]) == ["seed", "expanded"]
    assert exp.iloc[1]["gene"] == "KRAS"


def test_run_ppi_analysis_rwr_failure_falls_back(stubbed_core, monkeypatch):
    """扩展失败:按原始交集继续不让分析挂掉,失败原因写进 expansion 记录。"""
    monkeypatch.setattr(ana, "expand_targets_with_rwr",
                        lambda genes, n, prefix: (
                            genes, prefix, {"ok": False, "error": "no interactome"}))
    res = ana.run_ppi_analysis(mode="manual", drug_name="", disease_name="",
                               gene_text="TP53", score=0.7, rwr_expand=3)
    assert res["genes"] == ["TP53"]                       # 原始交集继续
    assert res["expansion_df"]["source"].tolist() == ["seed"]  # 纯种子清单
    assert res["params"][6] == 3


def test_run_ppi_analysis_progress_protocol(stubbed_core):
    """进度协议:阶段按序 set,分数单调推进,收在 1.0。"""
    prog = _FakeProgress()
    ana.run_ppi_analysis(mode="manual", drug_name="", disease_name="",
                         gene_text="TP53", score=0.7, progress=prog)
    step_names = [s for s, _ in prog.steps]
    for expected in ("targets", "string", "graph"):
        assert expected in step_names
    assert prog.fractions[0] == pytest.approx(0.02)
    # 编排层推进到 0.76;收尾的 1.0 由 save_analysis_outputs 负责(本桩不含)
    assert prog.fractions[-1] == pytest.approx(0.76)
    assert all(b >= a for a, b in zip(prog.fractions, prog.fractions[1:]))  # 单调
