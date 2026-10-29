"""core.analysis 的测试:RWR 查询集扩展的共享语义(网页 Build Report 与 CLI 共用)。

迁移自 selfcheck_network.py(44 项中的 1 项)。
"""
import WY.core.analysis as ana


def test_rwr_expansion_shared(monkeypatch):
    """RWR 查询集扩展共享语义(core.analysis.expand_targets_with_rwr):成功路径的
    平局排序/基因并集/前缀后缀、失败路径的原集降级与 error 透传、游走缓存 label 不带
    _rwr{N}、未启用或空基因集的原样返回。网页 Build Report 与 CLI run_pair 共用此
    唯一定义;旧 CLI 失败分支曾引用仅成功分支定义的变量直接 NameError,本测试锁住该路径。"""
    fake = {"ok": True, "cached": True,
            "top_nonseed": [{"gene": "EGFR", "score": 0.5},
                            {"gene": "BRAF", "score": 0.5},
                            {"gene": "AKT1", "score": 0.9}]}
    monkeypatch.setattr(ana, "run_rwr_cached", lambda genes, label=None, **k: fake)
    genes, prefix, info = ana.expand_targets_with_rwr(
        ["TP53", "KRAS"], 2, "metformin_Type_2_diabetes_mellitus")
    # 平局按基因名稳定排序:AKT1(0.9) 第一,0.5 的 BRAF/EGFR 按名字升序取 BRAF
    assert [r["gene"] for r in info["added"]] == ["AKT1", "BRAF"], info["added"]
    assert genes == ["AKT1", "BRAF", "KRAS", "TP53"], genes
    assert prefix == "metformin_Type_2_diabetes_mellitus_rwr2", prefix
    assert info["ok"] and info["cached"] and info["error"] is None, info
    assert info["attempted"] is True

    # 失败:按原基因集返回、prefix 仍带后缀(标记"请求了扩展")、error 透传
    monkeypatch.setattr(ana, "run_rwr_cached",
                        lambda genes, label=None, **k: {"ok": False, "error": "boom"})
    genes, prefix, info = ana.expand_targets_with_rwr(["TP53"], 30, "x_y")
    assert genes == ["TP53"] and info["ok"] is False, (genes, info)
    assert info["error"] == "boom" and info["attempted"] is True, info
    assert prefix == "x_y_rwr30", prefix

    # 游走缓存 label 不带 _rwr{N}(改 N 复用同一份游走)
    seen = {}
    monkeypatch.setattr(ana, "run_rwr_cached",
                        lambda genes, label=None, **k: seen.update(label=label)
                        or {"ok": False, "error": "probe"})
    ana.expand_targets_with_rwr(["TP53"], 5, "a_b")
    assert seen["label"] == "expand_a_b", seen

    # 未启用(0)或空基因集:不尝试游走,原样返回
    for g, n in ((["TP53"], 0), ([], 30)):
        genes, prefix, info = ana.expand_targets_with_rwr(g, n, "a_b")
        assert info["attempted"] is False and prefix == "a_b" and genes == list(g), (g, info)
