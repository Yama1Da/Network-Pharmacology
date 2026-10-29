"""core.opentargets 的测试:疾病编号搜索(词表优先级)、GraphQL errors 语义、
关联靶点解析与排序、进程内 memo 与失效。

阶段三第 2 批新写。打桩层级在 opentargets_post(模块内统一 POST 入口),不联网。
"""
import pandas as pd
import pytest

from WY import config
import WY.core.opentargets as ot


def _search_payload(hits):
    return {"data": {"search": {"hits": hits}}}


def _disease_payload(rows):
    return {"data": {"disease": {"id": "MONDO_0005148", "name": "type 2 diabetes",
                                 "associatedTargets": {"rows": rows}}}}


def _row(symbol, score, novelty=0.5, dt=None, organism_class=None):
    return {"score": score, "novelty": novelty,
            "datatypeScores": dt or [{"id": "genetic_association", "score": 0.7},
                                     {"id": "literature", "score": 0.9}],
            "target": {"id": f"ENSG_{symbol}", "approvedSymbol": symbol,
                       "approvedName": f"{symbol} full name", "biotype": "protein_coding",
                       "genomicLocation": {"chromosome": "17", "start": 1, "end": 2},
                       "targetClass": [{"label": "ENZYME"}]}}


def test_find_disease_id_prefers_mondo_then_umls(monkeypatch):
    """编号优先级:同名疾病命中多个词表时按 MONDO > UMLS > 其他 返回。"""
    monkeypatch.setattr(ot, "opentargets_post", lambda q, v, **k: _search_payload([
        {"id": "EFO_0000400", "entity": "disease", "name": "T2D"},
        {"id": "UMLS_C0011860", "entity": "disease", "name": "T2D"},
        {"id": "MONDO_0005148", "entity": "disease", "name": "T2D"},
    ]))
    assert ot.find_disease_id("type 2 diabetes") == "MONDO_0005148"


def test_find_disease_id_fallback_and_filters(monkeypatch):
    """无优先词表时取第一个疾病命中;非疾病实体与查询报错正确处理。"""
    monkeypatch.setattr(ot, "opentargets_post", lambda q, v, **k: _search_payload([
        {"id": "ENSG00000141510", "entity": "target", "name": "TP53"},   # 非疾病 → 滤掉
        {"id": "EFO_0000400", "entity": "disease", "name": "T2D"},
    ]))
    assert ot.find_disease_id("T2D") == "EFO_0000400"
    # GraphQL errors → None(查询报错不算网络失败,也不崩溃)
    monkeypatch.setattr(ot, "opentargets_post",
                        lambda q, v, **k: {"errors": [{"message": "boom"}]})
    assert ot.find_disease_id("T2D") is None
    # 网络彻底失败(None)→ None
    monkeypatch.setattr(ot, "opentargets_post", lambda q, v, **k: None)
    assert ot.find_disease_id("T2D") is None


def test_get_targets_full_flow(monkeypatch):
    """全链:搜索编号 → 关联靶点解析(分数取整/类别取一/证据分列)→ 降序截断 → 缓存+哨兵。"""
    def fake_post(query, variables, **k):
        if "search" in query:
            return _search_payload([{"id": "MONDO_0005148", "entity": "disease",
                                     "name": "T2D"}])
        rows = [_row("TP53", 0.5), _row("AKT1", 0.9),
                _row("EGFR", 0.7)]
        return _disease_payload(rows)

    monkeypatch.setattr(ot, "opentargets_post", fake_post)
    df, state = ot.get_disease_targets_stated("type 2 diabetes", need_rows=100)
    assert state == "ok"
    assert list(df["approvedSymbol"]) == ["AKT1", "EGFR", "TP53"]  # 按关联分降序
    top = df.iloc[0]
    assert top["association_score"] == 0.9 and top["novelty"] == 0.5
    assert top["genetic_association"] == 0.7 and top["literature"] == 0.9
    assert top["targetClass"] == "ENZYME" and top["chromosome"] == "17"
    # 缓存与哨兵已落盘(沙箱接管);行数不足抓取上限 → complete=True,后续不再重抓
    import os
    from WY.core.caching import load_json_file
    from WY.core.text_utils import safe_name
    stem = f"OpenTargets_{safe_name('type 2 diabetes')}.csv"
    assert os.path.isfile(os.path.join(config.DISEASE_CACHE_DIR, stem))
    meta = load_json_file(os.path.join(config.DISEASE_CACHE_DIR, "." + stem + ".meta"), {})
    assert meta.get("complete") is True


def test_unknown_disease_is_empty(monkeypatch):
    """查无此病 → (空表, empty),与请求失败区分。"""
    monkeypatch.setattr(ot, "opentargets_post",
                        lambda q, v, **k: _search_payload(
                            [{"id": "ENSG1", "entity": "target", "name": "X"}]))
    df, state = ot._get_disease_targets_impl("不存在的病")
    assert df.empty and state == "empty"


def test_disease_query_errors_is_failed(monkeypatch):
    """关联查询带 GraphQL errors → (空表, failed)。"""
    def fake_post(query, variables, **k):
        if "search" in query:
            return _search_payload([{"id": "MONDO_0005148", "entity": "disease",
                                     "name": "T2D"}])
        return {"errors": [{"message": "boom"}]}

    monkeypatch.setattr(ot, "opentargets_post", fake_post)
    df, state = ot._get_disease_targets_impl("type 2 diabetes")
    assert df.empty and state == "failed"


def test_memo_hit_and_invalidate(monkeypatch):
    """进程内 memo:成功结果被记住(第二次不发请求);失效指定疾病后重新抓取;
    失败结果不进 memo(下次自动重试)。"""
    monkeypatch.setattr(ot, "_success_memo", {})  # 隔离 memo,测试结束自动还原
    calls = {"n": 0}

    def fake_post(query, variables, **k):
        calls["n"] += 1
        if "search" in query:
            return _search_payload([{"id": "MONDO_0005148", "entity": "disease",
                                     "name": "T2D"}])
        return _disease_payload([_row("TP53", 0.5)])

    monkeypatch.setattr(ot, "opentargets_post", fake_post)
    df1, s1 = ot.get_disease_targets_stated("memo_test_disease")
    assert s1 == "ok" and calls["n"] == 2  # 搜索 + 关联查询
    df2, s2 = ot.get_disease_targets_stated("memo_test_disease")
    assert s2 == "ok" and calls["n"] == 2  # memo 命中,零新请求
    assert df2 is df1
    # 失效该疾病 → 重新联网。refresh 链路的真实顺序是先把缓存文件移走(改 .bak)再失效,
    # 否则重新抓取仍会命中文件缓存;测试按同样顺序删除缓存文件。
    import os
    from WY.core.text_utils import safe_name
    csv_file = os.path.join(config.DISEASE_CACHE_DIR,
                            f"OpenTargets_{safe_name('memo_test_disease')}.csv")
    os.remove(csv_file)
    ot.invalidate_disease_targets_cache("memo_test_disease")
    ot.get_disease_targets_stated("memo_test_disease")
    assert calls["n"] == 4
    # 只清该名称,其他疾病不受影响
    ot.get_disease_targets_stated("other_disease")
    ot.invalidate_disease_targets_cache("memo_test_disease")
    assert any(k[0] == "other_disease" for k in ot._success_memo)
