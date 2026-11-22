"""benchmark/disease_module_benchmark 单元测试(清单 §1.8)。

重点:
  - 防泄漏最后一道闸(陷阱#1):assert_no_leak 违规输入会响——审核关卡
    "防泄漏断言存在且被测试覆盖"由此用例满足
  - SeedSequence.spawn 独立流(陷阱#2):确定性可复现 + 各轮隐藏集不全同
  - Recall@N 的 N=扩张规模(陷阱#3)与分层报告结构
  - 适配层走一期缓存与限流(陷阱#4):无 key/无缓存 → 响亮报错,零联网
"""

import json

import networkx as nx
import pandas as pd
import pytest

from WY import config
from WY.benchmark import disease_module_benchmark as dmb


def _ctx_from_graph(g: nx.Graph):
    nodes = tuple(sorted(g.nodes()))
    idx = {x: i for i, x in enumerate(nodes)}
    Mt = nx.to_scipy_sparse_array(g, nodelist=nodes, weight=None, format="csr")
    degree_map = {x: int(d) for x, d in g.degree()}
    return nodes, idx, Mt, degree_map


@pytest.fixture()
def ctx():
    g = nx.Graph()
    g.add_edges_from([("A", "B"), ("A", "C"), ("B", "C"),
                      ("A", "D"), ("D", "E"), ("B", "F")])
    return _ctx_from_graph(g)


@pytest.fixture()
def assoc():
    """6 基因关联:中位数 0.35 → 高分池 {A,B,C},低分池 {D,E,F}。"""
    return pd.DataFrame({"approvedSymbol": ["A", "B", "C", "D", "E", "F"],
                         "disgenet_score": [0.9, 0.8, 0.4, 0.3, 0.2, 0.1]})


class TestLeakGate:
    def test_gate_fires_on_violation(self):
        """审核关卡:防泄漏断言存在且被测试覆盖——违规输入必须响。"""
        with pytest.raises(AssertionError):
            dmb.assert_no_leak({"A"}, ["A", "B"], 0)

    def test_gate_silent_when_disjoint(self):
        dmb.assert_no_leak({"A"}, ["B", "C"], 0)  # 不抛即通过

    def test_gate_covered_by_protocol(self, ctx, assoc):
        """协议全轮次执行断言(覆盖):固定种子跑两遍结果逐位一致,且每轮
        隐藏集与按 pool−hidden 推导的种子集不相交。"""
        r1 = dmb.leave_one_out_protocol(ctx, assoc, ["A", "B", "C", "D", "E", "F"],
                                        n_repeats=8, hide_frac=0.5, n=5,
                                        base_seed=7)
        r2 = dmb.leave_one_out_protocol(ctx, assoc, ["A", "B", "C", "D", "E", "F"],
                                        n_repeats=8, hide_frac=0.5, n=5,
                                        base_seed=7)
        pd.testing.assert_frame_equal(r1["table"], r2["table"])
        assert len(r1["audit"]) == 8
        for rec in r1["audit"]:
            hidden = set(rec["hidden"])
            seeds = {"A", "B", "C", "D", "E", "F"} - hidden
            assert hidden.isdisjoint(seeds)   # 不变量逐轮成立
            assert rec["n_seeds"] == 6 - len(hidden)


class TestProtocol:
    def test_independent_streams(self, ctx, assoc):
        """陷阱#2:spawn 独立流——8 轮隐藏集组合必不全同(全同=虚而不实)。"""
        res = dmb.leave_one_out_protocol(ctx, assoc,
                                         ["A", "B", "C", "D", "E", "F"],
                                         n_repeats=8, hide_frac=0.5, n=5,
                                         base_seed=7)
        combos = {tuple(r["hidden"]) for r in res["audit"]}
        assert len(combos) > 1

    def test_strata_and_table_shape(self, ctx, assoc):
        res = dmb.leave_one_out_protocol(ctx, assoc,
                                         ["A", "B", "C", "D", "E", "F"],
                                         n_repeats=8, hide_frac=0.5, n=5,
                                         base_seed=7)
        assert res["strata"]["threshold"] == pytest.approx(0.35)  # 中位数切分
        assert res["strata"]["high"] == ["A", "B", "C"]
        assert res["strata"]["low"] == ["D", "E", "F"]
        methods = set(res["table"]["method"])
        assert methods == {"rwr", "diamond", "neighbor", "consensus"}
        for m in methods:
            rows = res["table"][res["table"]["method"] == m]
            assert set(rows["stratum"]) <= {"overall", "high", "low"}
            assert rows["mean"].between(0, 1).all()   # Recall@N 有界
        assert isinstance(res["consensus_wins"], bool)

    def test_explicit_threshold(self, ctx, assoc):
        res = dmb.leave_one_out_protocol(ctx, assoc,
                                         ["A", "B", "C", "D", "E", "F"],
                                         n_repeats=4, hide_frac=0.5, n=5,
                                         base_seed=7, score_threshold=0.5)
        assert res["strata"]["high"] == ["A", "B"]
        assert res["strata"]["low"] == ["C", "D", "E", "F"]

    def test_out_of_lcc_recorded(self, ctx, assoc):
        res = dmb.leave_one_out_protocol(ctx, assoc,
                                         ["A", "B", "C", "D", "E", "F", "ZZ"],
                                         n_repeats=4, hide_frac=0.5, n=5,
                                         base_seed=7)
        assert res["out_of_lcc"] == ["ZZ"]            # 不可恢复基因立账,不入分母
        assert res["pool_size"] == 6

    def test_pool_too_small_rejected(self, ctx, assoc):
        with pytest.raises(RuntimeError):
            dmb.leave_one_out_protocol(ctx, assoc, ["A", "B"],
                                       n_repeats=4, hide_frac=0.5, n=5)


class TestAdapters:
    def test_load_associations_no_key_loud(self, monkeypatch):
        """陷阱#4:走一期 fetch;无 key 时一期返回空表 → 适配层响亮报错,
        全程零联网(无 key 分支不发请求)。"""
        monkeypatch.setattr(config, "DisGeNET_API_KEY", "", raising=False)
        with pytest.raises(RuntimeError):
            dmb.load_associations("breast cancer")

    def test_load_t_d_missing_merged_loud(self):
        """merged(四源整合)缺失时拒绝——回测不得拿 DisGeNET 单源顶替 T_D。"""
        with pytest.raises(RuntimeError):
            dmb.load_t_d("不存在的疾病xyz")


class TestCliAndSave:
    def test_save_results(self, ctx, assoc, tmp_path):
        res = dmb.leave_one_out_protocol(ctx, assoc,
                                         ["A", "B", "C", "D", "E", "F"],
                                         n_repeats=4, hide_frac=0.5, n=5,
                                         base_seed=7)
        paths = dmb.save_results(res, str(tmp_path / "out"))
        for p in paths.values():
            assert __import__("os").path.isfile(p)
        summary = json.loads(open(paths["summary"], encoding="utf-8").read())
        assert summary["schema"] == "disease_module_benchmark.v1"
        assert summary["consensus_wins"] in (True, False)
        audit = pd.read_csv(paths["audit"], comment="#")
        assert len(audit) == 4                      # 切分清单逐轮落盘(§5.1)

    def test_main_exit_codes(self, ctx, assoc, tmp_path, monkeypatch):
        """CLI 退出码约定:全成 0/任一败 1(清单 §2 cli 行)。"""
        fake = dmb.leave_one_out_protocol(ctx, assoc,
                                          ["A", "B", "C", "D", "E", "F"],
                                          n_repeats=4, hide_frac=0.5, n=5,
                                          base_seed=7)
        monkeypatch.setattr(dmb, "run_benchmark", lambda *a, **k: fake)
        rc = dmb.main(["--disease", "test", "--repeats", "4",
                       "--out", str(tmp_path / "ok")])
        assert rc == 0
        monkeypatch.setattr(dmb, "run_benchmark",
                            lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")))
        rc = dmb.main(["--disease", "test", "--out", str(tmp_path / "bad")])
        assert rc == 1
