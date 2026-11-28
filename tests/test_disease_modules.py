"""core/disease_modules 单元测试(清单 §1.4,队友名下)。

手算对照全覆盖四个陷阱:
  #2 超几何参数顺序(M=全图节点数, n=模块内节点数, N=候选度数)——
     _diamond_pvalue(2,2,2,6)=C(2,2)/C(6,2)=1/15,填错位不报错但数值不同
  小图扩张顺序手算推演(论文式示例):A,B 种子 → C(k_out=2) → F → D
  #3 交集外基因计数入报告 / #4 三方法同 N / #5 秩从 1 开始
RWR 走合成 ctx(无需真实互作组),零改动复用一期 _rwr_from_context。
"""

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from WY.core import disease_modules as dmo


def _ctx_from_graph(g: nx.Graph):
    """合成 LCC 上下文(与 network_algorithms._lcc_context 同构)。"""
    nodes = tuple(sorted(g.nodes()))
    idx = {x: i for i, x in enumerate(nodes)}
    Mt = nx.to_scipy_sparse_array(g, nodelist=nodes, weight=None, format="csr")
    degree_map = {x: int(d) for x, d in g.degree()}
    return nodes, idx, Mt, degree_map


@pytest.fixture()
def ctx():
    """6 节点小图:A-B,A-C,B-C(种子三角),A-D,D-E,B-F。
    度:A=3,B=3,C=2,D=2,E=1,F=1。扩张顺序可手算推演。"""
    g = nx.Graph()
    g.add_edges_from([("A", "B"), ("A", "C"), ("B", "C"),
                      ("A", "D"), ("D", "E"), ("B", "F")])
    return _ctx_from_graph(g)


class TestDiamond:
    def test_pvalue_param_order_hand_calc(self):
        """陷阱#2:手算对照。M=6, n_mod=2, k=2, k_out=2 →
        P(X≥2)=C(2,2)C(4,0)/C(6,2)=1/15≈0.0667;填错位会得到别的数。"""
        assert dmo._diamond_pvalue(2, 2, 2, 6) == pytest.approx(1 / 15)
        assert dmo._diamond_pvalue(1, 2, 2, 6) == pytest.approx(0.6)  # (8+1)/15
        assert dmo._diamond_pvalue(0, 2, 2, 6) == pytest.approx(1.0)  # 无连接

    def test_expansion_order_hand_calc(self, ctx):
        """小图扩张顺序手算推演(严格按原文逐步重算):
        step1 {A,B}: C: P(X≥2|6,2,2)=1/15 最小 → 入模;
        step2 {A,B,C}: F: P(X≥1|6,3,1)=0.5 < D 的 0.8 → 入模;
        step3: D: 14/15 < E 的 1 → 入模。n=5 时 join_order=[C,F,D]。"""
        res = dmo.diamond_expand(ctx, ["A", "B"], n=5)
        assert res["join_order"] == ["C", "F", "D"]
        assert res["module"] == ["A", "B", "C", "F", "D"]  # 种子在前基因序
        assert res["n_added"] == 3

    def test_seeds_out_of_lcc_counted(self, ctx):
        res = dmo.diamond_expand(ctx, ["A", "Z"], n=3)
        assert res["coverage"]["seeds_out_of_lcc"] == ["Z"]   # 不静默丢(陷阱#3)
        assert res["coverage"]["seeds_in_lcc"] == 1


class TestNeighbor:
    def test_ranking_by_seed_connection(self, ctx):
        """连接种子数降序:C(2)>D(1,度高)>F(1)>E(0);n=5 需 3 个 → C,D,F。"""
        res = dmo.neighbor_expansion(ctx, ["A", "B"], n=5)
        assert res["join_order"] == ["C", "D", "F"]
        assert res["module"] == ["A", "B", "C", "D", "F"]


class TestRwr:
    def test_reuse_and_size(self, ctx):
        res = dmo.rwr_module(ctx, ["A", "B"], n=5)
        assert res["ok"] is True
        assert len(res["module"]) == 5                       # 统一 N(陷阱#4)
        assert set(res["module"][:2]) == {"A", "B"}          # 种子在前(按 RWR 得分降序)
        assert set(res["join_order"]) <= {"C", "D", "E", "F"}

    def test_empty_seeds_reported(self, ctx):
        res = dmo.rwr_module(ctx, ["Z"], n=5)
        assert res["ok"] is False and res["error"] == "empty_after_filter"


class TestConsensus:
    def test_same_n_and_rank_semantics(self, ctx):
        out = dmo.consensus_targets(ctx, ["A", "B"], n=5, k=2)
        assert out["ok"] is True
        assert len(set(out["module_sizes"].values())) == 1   # 三方法同 N(陷阱#4)
        df = out["consensus"]
        assert df["rwr_rank"].min() == 1                     # 秩从 1(陷阱#5)
        assert df["diamond_rank"].min() == 1
        # C 在三种方法都是首个扩张成员 → n_methods=3 且 consistent;
        # 种子占秩 1-2,C 的三方法秩均=3 → 共识分 3.0,且为非种子最优
        c = df.set_index("gene").loc["C"]
        assert c["n_methods"] == 3 and bool(c["consistent"])
        assert c["consensus_score"] == pytest.approx(3.0)
        nonseed = df[~df["is_seed"]]
        assert c["consensus_score"] == nonseed["consensus_score"].min()
        assert bool(df.set_index("gene").loc["A", "is_seed"])

    def test_k_flag_and_out_of_lcc_merge(self, ctx):
        out = dmo.consensus_targets(ctx, ["A", "ZZ"], n=5, k=3)
        assert out["coverage"]["seeds_out_of_lcc"] == ["ZZ"]
        assert out["k"] == 3
        df = out["consensus"].set_index("gene")
        assert not bool(df.loc["D", "consistent"]) or df.loc["D", "n_methods"] >= 3

    def test_deterministic(self, ctx):
        a = dmo.consensus_targets(ctx, ["A", "B"], n=5)["consensus"]
        b = dmo.consensus_targets(ctx, ["A", "B"], n=5)["consensus"]
        pd.testing.assert_frame_equal(a, b)                  # 可复现红线


class TestSignificance:
    def test_smoke_small_perms(self, ctx):
        out = dmo.module_significance(ctx, ["A", "B"], n=5,
                                      n_perms=5, null_seed=0)
        assert set(out["methods"]) == {"rwr", "diamond", "neighbor"}
        for m in out["methods"].values():
            assert 0 < m["p"] <= 1.0                          # 经验单侧 p
            assert m["observed_cohesion"] >= 0
        again = dmo.module_significance(ctx, ["A", "B"], n=5,
                                        n_perms=5, null_seed=0)
        assert out["methods"] == again["methods"]             # 种子固定可复现

    def test_seed_module_cohesion_positive(self, ctx):
        out = dmo.module_significance(ctx, ["A", "B"], n=5,
                                      n_perms=3, null_seed=1)
        assert out["methods"]["diamond"]["observed_cohesion"] >= 2  # 种子三角+连接


class TestExport:
    def test_export_consensus_csv(self, ctx, tmp_path):
        out = dmo.consensus_targets(ctx, ["A", "B"], n=5)
        path = dmo.export_consensus_csv(out["consensus"],
                                        str(tmp_path / "consensus_targets.csv"))
        df = pd.read_csv(path, comment="#")
        assert "consensus_score" in df.columns and len(df) == len(out["consensus"])


class TestProgressHooks:
    """2026-11-28 UI 进度条钩子(§1.12):步进/逐置换回调只读不动结果——
    同输入同输出(可复现红线);CLI 不传即零开销,core 不解释回调语义。"""

    def test_diamond_step_progress_fires(self, ctx):
        calls = []
        res = dmo.diamond_expand(ctx, ["A", "B"], n=5,
                                 step_progress=lambda m, d, t: calls.append((m, d, t)))
        assert res["join_order"] == ["C", "F", "D"]  # 与无钩子逐字一致
        # 每入模一步上报一次:3 步扩张 → (diamond, 模块当前大小, n_target)
        assert calls == [("diamond", 3, 5), ("diamond", 4, 5), ("diamond", 5, 5)]

    def test_diamond_no_hook_same_result(self, ctx):
        assert (dmo.diamond_expand(ctx, ["A", "B"], n=5)["join_order"]
                == dmo.diamond_expand(ctx, ["A", "B"], n=5,
                                      step_progress=lambda *_: None)["join_order"])

    def test_consensus_passes_hook_through(self, ctx):
        calls = []
        dmo.consensus_targets(ctx, ["A", "B"], n=5, k=2,
                              step_progress=lambda m, d, t: calls.append((m, d, t)))
        assert calls and all(m == "diamond" for m, _d, _t in calls)  # 仅 DIAMOnD 步进

    def test_significance_perm_progress_fires(self, ctx):
        calls = []
        dmo.module_significance(ctx, ["A", "B"], n=3, n_perms=4,
                                perm_progress=lambda m, d, t: calls.append((m, d, t)))
        methods = [m for m, _d, _t in calls]
        assert methods == ["rwr"] * 4 + ["diamond"] * 4 + ["neighbor"] * 4
        assert all(d == i + 1 for i, (_m, d, _t) in enumerate(calls[:4]))
        assert {t for _m, _d, t in calls} == {4}
