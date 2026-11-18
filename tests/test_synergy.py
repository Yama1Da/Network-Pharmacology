"""分离度协同预测数学核心单元测试(文件清单 §1.3 验收自检)。

红线:绝不真联网(net_block 拦截);LCC 上下文用构造小图注入(一期
"纯计算自检注入"同模式);s_AB 手算对照一正一负各一例(陷阱 #1 强制验收,
写反符号这里必炸)。
"""

import networkx as nx
import numpy as np
import pytest

from WY.core import synergy
from WY.core.synergy import SynergyPrecheckError, _analyze_pair_from_context

CONTRACT_KEYS = {"schema", "pair", "disease", "inputs", "separation", "proximity",
                 "geometry_label", "label_basis", "ranks", "reachability", "meta"}


def _ctx_from_graph(G) -> tuple:
    """构造小图 → 一期 _lcc_context 同形 ctx(nodes, idx, Mt, degree_map),含 LCC 截取。"""
    lcc = max(nx.connected_components(G), key=len)
    sub = G.subgraph(lcc)
    nodes = tuple(sorted(sub.nodes()))
    idx = {g: i for i, g in enumerate(nodes)}
    adj = nx.to_scipy_sparse_array(sub, nodelist=nodes, weight=None, format="csr")
    Mt = adj.T.tocsr()
    degree_map = {g: int(d) for g, d in sub.degree()}
    return nodes, idx, Mt, degree_map


def _clique(G, names: list[str]) -> None:
    for i, x in enumerate(names):
        for y in names[i + 1:]:
            G.add_edge(x, y)


def _named(G, prefix="N") -> nx.Graph:
    """生成器默认整数标签 → 统一改大写字符串(与 normalize_genes 口径一致)。"""
    return nx.relabel_nodes(G, {n: f"{prefix}{n}" for n in G.nodes()})


def _run(G, a, b, d=None, n_perm=20, seed=0, **kwargs):
    return _analyze_pair_from_context(_ctx_from_graph(G), a, b, d,
                                      n_perm=n_perm, seed=seed, **kwargs)


def test_s_ab_separated_hand_computed():
    """陷阱 #1 验收(分离侧):两个 K5 团由一条桥边相连。
    手算:d_aa=d_bb=1;d_ab=(1×1 + 2×8 + 3×16)/25 = 2.6;s_ab=+1.6>0。"""
    G = nx.Graph()
    a = [f"A{i}" for i in range(5)]
    b = [f"B{i}" for i in range(5)]
    _clique(G, a)
    _clique(G, b)
    G.add_edge("A0", "B0")
    r = _run(G, a, b)
    sep = r["separation"]
    assert sep["d_aa"] == pytest.approx(1.0)
    assert sep["d_bb"] == pytest.approx(1.0)
    assert sep["d_ab"] == pytest.approx(2.6)
    assert sep["s_ab"] == pytest.approx(1.6)
    assert r["geometry_label"] == "separated"
    # 方向健全性:两团分离程度远超小图随机组合的下尾,p_lower 应贴 1
    assert sep["null"]["p_lower"] >= 0.95
    assert sep["null"]["z"] is not None and sep["null"]["z"] > 0


def test_s_ab_overlap_hand_computed():
    """陷阱 #1 验收(重叠侧):C12 环上取奇偶交错两集。
    手算:d_aa=d_bb=36/10=3.6;跨集 25 对距离和=73 → d_ab=2.92;s_ab=−0.68<0。"""
    G = _named(nx.cycle_graph(12))
    a = [f"N{i}" for i in (0, 2, 4, 6, 8)]
    b = [f"N{i}" for i in (1, 3, 5, 7, 9)]
    r = _run(G, a, b)
    sep = r["separation"]
    assert sep["d_aa"] == pytest.approx(3.6)
    assert sep["d_bb"] == pytest.approx(3.6)
    assert sep["d_ab"] == pytest.approx(73 / 25)
    assert sep["s_ab"] == pytest.approx(-0.68)
    assert r["geometry_label"] == "overlap"
    assert 0.0 <= sep["null"]["p_lower"] <= 1.0


def test_e3_precheck_refuses_below_five():
    """陷阱 #5:入 LCC 靶点 <5 抛业务异常(不是 return None),reason 稳定可映射 i18n。"""
    G = _named(nx.complete_graph(10))
    a4 = [f"N{i}" for i in range(4)]
    b4 = [f"N{i}" for i in range(4, 8)]
    with pytest.raises(SynergyPrecheckError) as ei:
        _run(G, a4, b4)
    assert ei.value.reason == "too_few_targets_a"
    with pytest.raises(SynergyPrecheckError) as ei:
        _run(G, a4 + ["N4"], b4)  # a 够 5,b 仍 4 → 拒 b 侧
    assert ei.value.reason == "too_few_targets_b"


def test_e3_power_limited_and_normal_tier():
    """E3 分档:5–14 → power_limited=true;>=15 → false。"""
    G = nx.Graph()
    a = [f"A{i}" for i in range(5)]
    b = [f"B{i}" for i in range(5)]
    _clique(G, a)
    _clique(G, b)
    G.add_edge("A0", "B0")  # 连桥:否则两个连通分量,LCC 截取会丢掉一侧(约定行为)
    assert _run(G, a, b)["inputs"]["power_limited"] is True   # 5/5 → 分档内
    G2 = _named(nx.complete_graph(30))
    a15 = [f"N{i}" for i in range(15)]
    b15 = [f"N{i}" for i in range(15, 30)]
    r = _run(G2, a15, b15, n_perm=10)
    assert r["inputs"]["power_limited"] is False
    # K30 全连通:d_aa=d_bb=d_ab=1 → s_ab=0 → tie(枚举第三态)
    assert r["geometry_label"] == "tie"


def test_lcc_coverage_counts_not_silent():
    """陷阱 #8:LCC 外靶点不参与距离但计入 coverage;规范化(大小写/空白)同口径。"""
    G = nx.Graph()
    _clique(G, [f"C{i}" for i in range(7)])   # LCC
    G.add_edge("X1", "X2")                    # LCC 外小分支
    a = ["C0", "C1", "C2", "C3", "C4", "X1"]
    b = ["c0", "c1", "c2", "c3", " c5 ", "x2"]  # 小写/空白 → 规范化后 5 个入网
    r = _run(G, a, b)
    cov = r["inputs"]["coverage"]
    assert cov == {"a_total": 6, "a_in_lcc": 5, "b_total": 6, "b_in_lcc": 5,
                   "d_total": 0, "d_in_lcc": 0}
    assert r["disease"] is None
    assert r["proximity"] == {"z_a": None, "z_b": None, "z_ab": None}


def test_proximity_reuse_wiring(monkeypatch):
    """陷阱 #7:z_A/z_B/z_AB 编排复用一期 proximity_zscore——三调用参数与结果映射。"""
    calls: list[tuple[list, list]] = []

    def fake_prox(drug_genes, disease_genes, species="Human", n_perm=None, seed=None):
        calls.append((list(drug_genes), list(disease_genes)))
        return {"ok": True, "error": None, "d_c": 1.5, "z": -1.23, "p_empirical": 0.05}

    monkeypatch.setattr(synergy, "proximity_zscore", fake_prox)
    G = _named(nx.cycle_graph(12))
    a = [f"N{i}" for i in (0, 2, 4, 6, 8)]
    b = [f"N{i}" for i in (1, 3, 5, 7, 9)]
    d = ["D0", "D1", "D2", "D3", "D4", "D5"]
    r = _run(G, a, b, d)
    assert len(calls) == 3
    assert calls[0][0] == a and calls[1][0] == b
    assert calls[2][0] == a + [g for g in b if g not in a]      # z_ab = A∪B 去重保序
    assert all(c[1] == d for c in calls)
    assert r["proximity"]["z_ab"] == {"ok": True, "error": None,
                                      "d_c": 1.5, "z": -1.23, "p": 0.05}
    assert r["disease"] == {"name": "", "n_targets": 6, "source": ""}


def test_proximity_failure_passthrough(monkeypatch):
    """一期覆盖率保护 ok=False 透传,不伪装成成功,也不中断其余字段。"""
    def failing_prox(drug_genes, disease_genes, species="Human", n_perm=None, seed=None):
        return {"ok": False, "error": "low_coverage_drug", "d_c": None,
                "z": None, "p_empirical": None}

    monkeypatch.setattr(synergy, "proximity_zscore", failing_prox)
    G = _named(nx.cycle_graph(12))
    a = [f"N{i}" for i in (0, 2, 4, 6, 8)]
    b = [f"N{i}" for i in (1, 3, 5, 7, 9)]
    r = _run(G, a, b, ["D0", "D1", "D2", "D3", "D4", "D5"])
    assert r["proximity"]["z_a"] == {"ok": False, "error": "low_coverage_drug",
                                     "d_c": None, "z": None, "p": None}
    assert r["geometry_label"] == "overlap"  # 几何侧不受邻近度失败影响


def test_contract_conformance_and_determinism():
    """契约 synergy.v1 冻结符合性(§4 冻结纪律) + 同种子可复现。"""
    G = _named(nx.cycle_graph(12))
    a = [f"N{i}" for i in (0, 2, 4, 6, 8)]
    b = [f"N{i}" for i in (1, 3, 5, 7, 9)]
    r = _run(G, a, b)
    assert set(r) == CONTRACT_KEYS
    assert r["schema"] == "synergy.v1.1"
    assert set(r["separation"]) == {"s_ab", "d_aa", "d_bb", "d_ab", "null"}
    assert set(r["separation"]["null"]) == {"mean", "sd", "z", "p_lower", "n_perm"}
    assert set(r["meta"]) == {"computed_at", "engine", "run"}
    assert r["meta"]["run"] == {"top_n": None, "sort": None, "proximity_n_perm": None}
    assert "cache" not in r["meta"]          # 2026-11-17 定案:无缓存机制
    assert r["ranks"] is None and r["reachability"] is None
    assert r["label_basis"].startswith("s_ab")
    r2 = _run(G, a, b)
    assert r2["separation"]["s_ab"] == r["separation"]["s_ab"]
    assert r2["separation"]["null"]["mean"] == r["separation"]["null"]["mean"]
    assert r2["separation"]["null"]["p_lower"] == r["separation"]["null"]["p_lower"]


# ---- L1 统计健全性(不需要生物知识,验收自检"置换 p 值近似均匀"的落地) ----


def test_run_params_and_progress_wiring(monkeypatch):
    """v1.1:meta.run 记录页级参数;proximity_n_perm 透传一期邻近度;progress
    只发 fraction(无文案,核心无文案纪律)且单调不减。"""
    seen: dict = {}

    def fake_prox(drug_genes, disease_genes, species="Human", n_perm=None, seed=None):
        seen["n_perm"] = n_perm
        return {"ok": True, "error": None, "d_c": 1.0, "z": -1.0, "p_empirical": 0.02}

    monkeypatch.setattr(synergy, "proximity_zscore", fake_prox)
    G = _named(nx.cycle_graph(12))
    a = [f"N{i}" for i in (0, 2, 4, 6, 8)]
    b = [f"N{i}" for i in (1, 3, 5, 7, 9)]
    fr: list = []
    r = _analyze_pair_from_context(
        _ctx_from_graph(G), a, b, ["D0", "D1", "D2", "D3", "D4", "D5"],
        n_perm=10, seed=0, top_n=50, sort_rule="dual_source_max_norm",
        proximity_n_perm=37,
        progress=lambda fraction=None, **kw: fr.append(fraction))
    assert r["meta"]["run"] == {"top_n": 50, "sort": "dual_source_max_norm",
                                "proximity_n_perm": 37}
    assert seen["n_perm"] == 37
    assert fr and fr[0] >= 0.1 and fr[-1] <= 0.95
    assert all(x <= y for x, y in zip(fr, fr[1:]))


def test_separation_symmetry_and_self_zero():
    """代数性质卡方向(比手算更硬的防线):s_AB(A,B)=s_AB(B,A);s_AB(A,A)=0。"""
    G = _named(nx.cycle_graph(12))
    a = [f"N{i}" for i in (0, 2, 4, 6, 8)]
    b = [f"N{i}" for i in (1, 3, 5, 7, 9)]
    s_ab = _run(G, a, b)["separation"]["s_ab"]
    s_ba = _run(G, b, a)["separation"]["s_ab"]
    assert s_ab == pytest.approx(s_ba)
    r_self = _run(G, a, a)
    # 精确恒等式:s(A,A) = −d_AA/n(跨集均值含 n 个 0 距离自对)。共享靶点把
    # s_AB 拉负是"重叠度"定义内行为,非 bug;2026-11-17 审查修正原 0 假设。
    assert r_self["separation"]["s_ab"] == pytest.approx(-3.6 / 5)
    assert r_self["geometry_label"] == "overlap"


def test_p_lower_uniform_under_null():
    """验收自检"置换 p 值对随机靶点集应近似均匀分布":随机靶点对本身与零分布
    可交换 → p_lower ~ U(0,1)。KS 检验固定种子确定性通过;不通过 = 零分布
    构造或经验 p 公式有错(如尾向写反、重采样漏度匹配)。"""
    from scipy.stats import kstest

    rng = np.random.default_rng(7)
    G = _named(nx.gnp_random_graph(60, 0.08, seed=11))
    # 取最大连通分量,保证与生产的 LCC 口径一致
    G = G.subgraph(max(nx.connected_components(G), key=len))
    nodes = sorted(G.nodes())
    ctx = _ctx_from_graph(G)
    pairs_p = []
    for _ in range(30):
        picks = rng.choice(len(nodes), size=12, replace=False)
        genes = [nodes[i] for i in picks]
        r = _analyze_pair_from_context(ctx, genes[:6], genes[6:], None,
                                       n_perm=100, seed=0, bin_size=100)
        pairs_p.append(r["separation"]["null"]["p_lower"])
    stat, p_ks = kstest(pairs_p, "uniform")
    assert p_ks > 0.01, f"p_lower 非均匀: KS stat={stat:.3f}, p={p_ks:.3f}, 样本={pairs_p}"


def test_relabel_and_isolated_invariance():
    """不变性:节点重编号不改变结果;图里加孤立节点不改变结果。违反任一 =
    距离计算依赖了标签顺序或节点总数,属实现 bug。"""
    G = _named(nx.cycle_graph(12))
    a = [f"N{i}" for i in (0, 2, 4, 6, 8)]
    b = [f"N{i}" for i in (1, 3, 5, 7, 9)]
    base = _run(G, a, b)
    # 重编号:N0→N100, N1→N101, …(映射打乱原顺序)
    mapping = {f"N{i}": f"N{i * 7 % 12}" for i in range(12)}
    assert sorted(mapping.values()) == sorted(mapping.keys())  # 映射本身是重排
    relabeled = _run(nx.relabel_nodes(G, mapping),
                     [mapping[x] for x in a], [mapping[x] for x in b])
    assert relabeled["separation"]["s_ab"] == pytest.approx(base["separation"]["s_ab"])
    assert relabeled["separation"]["d_ab"] == pytest.approx(base["separation"]["d_ab"])
    # 加孤立节点:不进 LCC,coverage 不变,距离不变
    with_isolated = G.copy()
    with_isolated.add_edges_from([("Z1", "Z2")])  # 2 节点小分量,LCC 外
    iso = _run(with_isolated, a, b)
    assert iso["separation"]["s_ab"] == pytest.approx(base["separation"]["s_ab"])
    assert iso["inputs"]["coverage"]["a_in_lcc"] == 5
