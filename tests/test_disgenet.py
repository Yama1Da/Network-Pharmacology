"""core.disgenet 的测试:疾病名解析、配额限流、空密钥优雅降级、分页解析与去重。

阶段三第 2 批新写。HTTP 一律打桩(get_session 假会话),限流状态文件由
conftest 沙箱隔离,配额阈值经 monkeypatch 调小以离线复现限流分支。
"""
import os
import time

import pandas as pd
import pytest

from WY import config
import WY.core.disgenet as disgenet
import WY.core.opentargets as opentargets


class _FakeResp:
    """最小假响应:status_code/headers/json/raise_for_status。"""

    def __init__(self, payload, status=200, headers=None):
        self.status_code = status
        self.headers = headers or {}
        self._payload = payload

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"http {self.status_code}")


class _FakeSession:
    """假会话:记录收到的请求参数,按队列吐响应。"""

    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []

    def get(self, url, params=None, headers=None, timeout=None):
        self.calls.append({"url": url, "params": params, "headers": headers})
        return self.responses.pop(0)


# ---------- 疾病名 → DisGeNET 查询串(纯解析,零网络) ----------


def test_resolve_bare_cui_gets_umls_prefix():
    """裸 UMLS CUI(C+6~8位数字)自动补 UMLS_ 前缀。"""
    assert disgenet._resolve_disgenet_query("C0002395") == "UMLS_C0002395"


def test_resolve_prefixed_vocab_passthrough():
    """已带词表前缀(MONDO_/UMLS_/EFO_ 等)原样放行。"""
    assert disgenet._resolve_disgenet_query("MONDO_0005148") == "MONDO_0005148"
    assert disgenet._resolve_disgenet_query("UMLS_C0002395") == "UMLS_C0002395"


def test_resolve_free_text_via_opentargets(monkeypatch):
    """自由文本:借 Open Targets 解析成编号;解析失败时原样返回(行为不变)。"""
    monkeypatch.setattr(opentargets, "find_disease_id",
                        lambda name: "MONDO_0005148")
    assert disgenet._resolve_disgenet_query("type 2 diabetes") == "MONDO_0005148"
    monkeypatch.setattr(opentargets, "find_disease_id", lambda name: None)
    assert disgenet._resolve_disgenet_query("阿尔茨海默病") == "阿尔茨海默病"


# ---------- 主动限流(配额记账,状态文件已沙箱化) ----------


def test_rate_limit_slot_acquired(monkeypatch):
    """配额未满:占用一个名额并写入状态文件,直接放行。"""
    monkeypatch.setattr(config, "DISGENET_RATE_LIMIT_PER_MINUTE", 2)
    disgenet._acquire_disgenet_rate_limit_slot()
    disgenet._acquire_disgenet_rate_limit_slot()  # 第 2 个名额
    from WY.core.caching import load_json_file
    recorded = load_json_file(config.DISGENET_RATE_LIMIT_STATE_FILE, [])
    assert len(recorded) == 2


def test_rate_limit_daily_exhausted_raises(monkeypatch):
    """日配额用满:直接 RuntimeError(等待可能长达数小时,拒绝比等待合理)。"""
    monkeypatch.setattr(config, "DISGENET_RATE_LIMIT_PER_DAY", 2)
    monkeypatch.setattr(config, "DISGENET_RATE_LIMIT_PER_MONTH", 100)
    now = time.time()
    from WY.core.caching import save_json_file
    save_json_file(config.DISGENET_RATE_LIMIT_STATE_FILE,
                   [now - 10, now - 5])  # 两条都在今日窗口内
    with pytest.raises(RuntimeError, match="daily"):
        disgenet._acquire_disgenet_rate_limit_slot()


def test_rate_limit_minute_exceeded_aborts(monkeypatch):
    """分钟配额用满且需等待超过上限:报错放弃,不无限等(把 sleep 路径掐断在毫秒级)。"""
    monkeypatch.setattr(config, "DISGENET_RATE_LIMIT_PER_MINUTE", 1)
    monkeypatch.setattr(config, "DISGENET_RATE_LIMIT_MAX_WAIT", 0.0)
    from WY.core.caching import save_json_file
    save_json_file(config.DISGENET_RATE_LIMIT_STATE_FILE, [time.time()])
    with pytest.raises(RuntimeError, match="per-minute"):
        disgenet._acquire_disgenet_rate_limit_slot()


# ---------- 空密钥优雅降级(记日志返回空表,不崩溃不发请求) ----------


def test_empty_api_key_graceful(monkeypatch):
    """密钥为空:直接返回空表,连会话都不该碰(空密钥是配置错误,不是 0 靶点)。"""
    monkeypatch.setattr(config, "DisGeNET_API_KEY", "")
    def must_not_call():
        raise AssertionError("network session used despite empty API key")
    monkeypatch.setattr(disgenet, "get_session", must_not_call)
    df = disgenet.fetch_disease_targets_disgenet("type 2 diabetes")
    assert df.empty


def test_refresh_reports_no_api_key(monkeypatch):
    """管理页强制更新链路:密钥缺失 → 状态 failed:no_api_key(专用提示的依据)。"""
    monkeypatch.setattr(config, "DisGeNET_API_KEY", "")
    df, state = disgenet.fetch_disease_targets_refresh("type 2 diabetes")
    assert df.empty and state == "failed:no_api_key"


# ---------- 分页抓取:解析、跨页去重、提前收工 ----------


def _page_payload(rows):
    return {"status": "OK", "payload": rows}


def test_fetch_parses_rows_dedup_and_meta(monkeypatch):
    """正常抓取:解析字段、跨页同基因保留最高分、不满一页记 complete 哨兵、缓存落盘。"""
    monkeypatch.setattr(config, "DisGeNET_API_KEY", "test-key")
    monkeypatch.setattr(config, "DISGENET_PAGE_SIZE", 4)  # 3 行 < 4 → 单页取完,提前收工
    monkeypatch.setattr(config, "DISGENET_PAGES", 2)
    rows = [
        {"symbolOfGene": "TP53", "geneNcbiID": "7157",
         "geneEnsemblIDs": ["ENSG00000141510"], "diseaseName": "T2D",
         "diseaseUMLSCUI": "C0011860", "score": 0.3, "normalized_score": 0.5,
         "numPMIDs": 10, "yearInitial": 2000, "yearFinal": 2024},
        {"symbolOfGene": "AKT1", "geneNcbiID": "207",
         "geneEnsemblIDs": [], "diseaseName": "T2D",
         "diseaseUMLSCUI": "C0011860", "score": 0.9, "normalized_score": 0.7,
         "numPMIDs": 3, "yearInitial": 2005, "yearFinal": 2024},
        {"symbolOfGene": "TP53", "geneNcbiID": "7157",
         "geneEnsemblIDs": ["ENSG00000141510"], "diseaseName": "T2D",
         "diseaseUMLSCUI": "C0011860", "score": 0.8, "normalized_score": 0.6,
         "numPMIDs": 12, "yearInitial": 2001, "yearFinal": 2024},  # TP53 更高分记录
    ]
    session = _FakeSession([_FakeResp(_page_payload(rows))])
    monkeypatch.setattr(disgenet, "get_session", lambda: session)
    df = disgenet.fetch_disease_targets_disgenet("type 2 diabetes")
    # 请求头带密钥;单页取完后提前收工,不发第二页
    assert len(session.calls) == 1
    assert session.calls[0]["headers"]["Authorization"] == "test-key"
    assert session.calls[0]["params"]["page_number"] == 0
    # 跨页/页内同基因去重:TP53 保留高分 0.8(弃 0.3);全表按分数降序
    assert list(df["approvedSymbol"]) == ["AKT1", "TP53"]
    assert df.iloc[0]["disgenet_score"] == 0.9   # AKT1
    assert df.iloc[1]["disgenet_score"] == 0.8   # TP53 的高分记录
    assert df.iloc[0]["ensembl_id"] is None  # 空 Ensembl 列表 → None
    assert df.iloc[1]["ensembl_id"] == "ENSG00000141510"
    # 哨兵与缓存已落盘(路径由沙箱接管)
    from WY.core.caching import load_json_file
    from WY.core.text_utils import safe_name
    meta_file = os.path.join(
        config.DISEASE_CACHE_DIR,
        f".DisGeNET_{safe_name('type 2 diabetes')}.csv.meta")
    meta = load_json_file(meta_file, {})
    assert meta.get("complete") is True and meta.get("pages") == 1


def test_fetch_first_page_failure_is_empty(monkeypatch):
    """首页即失败(如状态非 OK):无保留行,返回空表。"""
    monkeypatch.setattr(config, "DisGeNET_API_KEY", "test-key")
    session = _FakeSession([_FakeResp({"status": "ERROR", "payload": "bad query"})])
    monkeypatch.setattr(disgenet, "get_session", lambda: session)
    df = disgenet.fetch_disease_targets_disgenet("type 2 diabetes")
    assert df.empty
