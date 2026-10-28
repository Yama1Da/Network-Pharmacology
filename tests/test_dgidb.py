"""core.dgidb 的测试:不打真实网络,用假 GraphQL 响应验证解析、去重与容错语义
(ok / empty / failed 三态分明,绝不把失败伪装成查无数据)。

阶段三第 2 批新写(仿原方案示例 B 的打桩手法)。
"""
import pandas as pd
import pytest

from WY import config
from WY.core.caching import save_csv_cache
from WY.core.text_utils import safe_name
import WY.core.dgidb as dgidb


class _FakeResp:
    """最小假响应:只需要 .json()。"""

    def __init__(self, payload):
        self._payload = payload

    def json(self):
        return self._payload


# 仿照 DGIdb GraphQL 真实返回结构捏造的假数据:两条交互 + 一条同基因同药重复
_GOOD_PAYLOAD = {
    "data": {"drugs": {"nodes": [{"interactions": [
        {"gene": {"name": "EGFR", "longName": "x", "conceptId": "C1"},
         "drug": {"name": "metformin", "conceptId": "C2", "approved": 1},
         "interactionScore": 7.5,
         "interactionClaims": [
             {"publications": [{"pmid": 123}, {"pmid": 456}],
              "source": {"sourceDbName": "DrugBank"}},
             {"publications": [], "source": {"sourceDbName": "DrugBank"}},
         ]},
        {"gene": {"name": "TP53", "longName": "y", "conceptId": "C3"},
         "drug": {"name": "metformin", "conceptId": "C2", "approved": 1},
         "interactionScore": 2.0, "interactionClaims": []},
        {"gene": {"name": "EGFR", "longName": "x", "conceptId": "C1"},
         "drug": {"name": "metformin", "conceptId": "C2", "approved": 1},
         "interactionScore": 7.5, "interactionClaims": []},  # 同基因同药重复
    ]}]}}}


def test_parses_interactions_dedup(monkeypatch):
    """正常响应:解析出基因/分数/来源/文献列,同基因同药去重,状态 ok。"""
    monkeypatch.setattr(dgidb, "request_with_retry",
                        lambda *a, **k: _FakeResp(_GOOD_PAYLOAD))
    df, state = dgidb.fetch_drug_targets_dgidb_stated("metformin")
    assert state == "ok"
    assert list(df["approvedSymbol"]) == ["EGFR", "TP53"]  # 重复已去重,保序
    row = df.iloc[0]
    assert row["interaction_score"] == 7.5
    assert row["interaction_sources"] == "DrugBank"        # 同源去重
    assert "123" in row["interaction_pmids"] and "456" in row["interaction_pmids"]


def test_errors_payload_is_failed(monkeypatch):
    """GraphQL 返回 errors → 状态 failed、空表(失败与查无数据绝不混淆)。"""
    monkeypatch.setattr(dgidb, "request_with_retry",
                        lambda *a, **k: _FakeResp({"errors": [{"message": "boom"}]}))
    df, state = dgidb.fetch_drug_targets_dgidb_stated("metformin")
    assert df.empty and state == "failed"


def test_missing_data_key_is_failed(monkeypatch):
    """响应缺 data 键(非法载荷)→ failed。"""
    monkeypatch.setattr(dgidb, "request_with_retry",
                        lambda *a, **k: _FakeResp({"something": "else"}))
    df, state = dgidb.fetch_drug_targets_dgidb_stated("metformin")
    assert df.empty and state == "failed"


def test_network_outage_is_failed(monkeypatch):
    """请求抛异常/重试耗尽返回 None → failed,不崩溃。"""
    monkeypatch.setattr(dgidb, "request_with_retry", lambda *a, **k: None)
    df, state = dgidb.fetch_drug_targets_dgidb_stated("metformin")
    assert df.empty and state == "failed"

    def boom(*a, **k):
        raise ConnectionError("simulated outage")

    monkeypatch.setattr(dgidb, "request_with_retry", boom)
    df, state = dgidb.fetch_drug_targets_dgidb_stated("metformin")
    assert df.empty and state == "failed"


def test_no_interactions_is_empty(monkeypatch):
    """请求成功但查无相互作用 → empty(与 failed 区分)。"""
    monkeypatch.setattr(dgidb, "request_with_retry",
                        lambda *a, **k: _FakeResp({"data": {"drugs": {"nodes": []}}}))
    df, state = dgidb.fetch_drug_targets_dgidb_stated("metformin")
    assert df.empty and state == "empty"


def test_cache_hit_skips_network(monkeypatch):
    """已有本地缓存:直接复用,不发任何网络请求。"""
    import os

    cached = pd.DataFrame({"approvedSymbol": ["TP53"], "drug_name": ["metformin"],
                           "interaction_score": [1.0]})
    # 缓存读取路径 = DRUG_CACHE_DIR/DGIdb_<safe_name>.csv(沙箱子目录),须按此路径预置
    cache_file = os.path.join(config.DRUG_CACHE_DIR,
                              f"DGIdb_{safe_name('metformin')}.csv")
    save_csv_cache(cached, cache_file)
    def must_not_call(*a, **k):
        raise AssertionError("network call made despite cache hit")
    monkeypatch.setattr(dgidb, "request_with_retry", must_not_call)
    df, state = dgidb.fetch_drug_targets_dgidb_stated("metformin")
    assert state == "ok" and list(df["approvedSymbol"]) == ["TP53"]
