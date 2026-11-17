"""core.pk_mirror 的测试:不打真实网络,用假响应验证搬运语义——

- 逐字保存(单位/数值保持原始字符串,不做任何换算——分层红线,见 §0.3);
- licence/access 保留;非文本文件(png)记入 skipped 不丢弃;
- 断点续传:已存在的 sid 跳过(force=False),不发起该研究的任何请求;
- 单研究失败不断批(失败进 _meta.json);
- 状态读取(读目录 + _meta.json,不联网)。

打桩手法仿 tests/test_dgidb.py(monkeypatch 模块的 request_with_retry);
conftest 的沙箱夹具会把 CACHE_DIR 指到临时目录,故本测试不触真实 Static/。
"""
import json
import os

import WY.core.pk_mirror as pk_mirror
from WY import config


def _tsv(rows: list[list[str]]) -> str:
    return "\n".join("\t".join(r) for r in rows) + "\n"


# 仿 PK-DB 真实返回结构捏造的假数据(双层 data 分页信封别简化掉)。
_STUDY_SUMMARIES = {
    "PKDB000001": {"sid": "PKDB000001", "name": "StudyOne", "licence": "closed",
                   "access": "public", "date": "2024-01-01",
                   "substances": ["drug1"], "reference_pmid": "111",
                   "reference_title": "title one"},
    "PKDB000002": {"sid": "PKDB000002", "name": "StudyTwo", "licence": "open",
                   "access": "public", "date": "2024-02-02",
                   "substances": ["drug2"], "reference_pmid": "222",
                   "reference_title": "title two"},
}

_INDEX = {"current_page": 1, "last_page": 1, "next_page_url": None,
          "data": {"count": 2, "data": list(_STUDY_SUMMARIES.values())}}

# Tab 表即药代参数表(measurement_type=cmax 行),数值/单位保持原始字符串形态
_TAB_ROWS = [
    ["study", "group", "count", "intervention", "measurement_type", "substance",
     "tissue", "method", "mean", "sd", "unit", "comments"],
    ["PKDB000001", "g1", "5.0", "D1", "cmax", "drug1", "serum", "HPLC",
     "359.2", "98.3", "ng/ml", "NA"],
]

_INTERVENTIONS = [
    ["study", "name", "route", "substance", "value", "unit"],
    ["PKDB000001", "D1", "oral", "drug1", "3.0", "mg"],
]

_TSV_TEXTS = {
    # PK-DB 约定:点前缀 = 策展解析后的表(与真实 _datafiles URL 一致)
    ".PKDB000001_Tab1.tsv": _tsv(_TAB_ROWS),
    ".PKDB000001_TabInterventions.tsv": _tsv(_INTERVENTIONS),
    ".PKDB000002_Tab1.tsv": _tsv([_TAB_ROWS[0],
                                  ["PKDB000002", "g1", "4.0", "D2", "cmax",
                                   "drug2", "serum", "LCMS", "120.0", "30.0",
                                   "ng/ml", "NA"]]),
}


class _FakeResp:
    """比 conftest._FakeResp 多 text/content:镜像层要下载 TSV 文本。"""

    def __init__(self, payload=None, text=None):
        self._payload = payload
        self.text = text
        self.content = text.encode("utf-8") if text is not None else b""

    def json(self):
        if self._payload is None:
            raise ValueError("no json payload")
        return self._payload


# 每研究的文件清单(png 是非文本 → skipped;两研究 deliberately 不同构)
_STUDY_FILES = {
    "PKDB000001": [".PKDB000001_Tab1.tsv", ".PKDB000001_TabInterventions.tsv",
                   "PKDB000001_Fig1.png"],
    "PKDB000002": [".PKDB000002_Tab1.tsv", "PKDB000002_Fig1.png"],
}


def _mk_fake(calls: list, fail_sids: tuple = ()):
    """按 URL 分发假响应;fail_sids 中的研究在 _studies 步抛错(模拟远端故障)。"""

    def fake(url, **kwargs):
        calls.append(url)
        if "/pkdata/studies/" in url and url.endswith("/pkdata/studies/"):
            return _FakeResp(payload=_INDEX)
        if "/info_nodes/" in url:
            return _FakeResp(payload={"current_page": 1, "last_page": 1,
                                      "next_page_url": None,
                                      "data": {"count": 1, "data": [
                                          {"sid": "drug1", "ntype": "substance",
                                           "synonyms": ["metformin"]}]}})
        if "/_studies/" in url:
            sid = url.rstrip("/").rsplit("/", 1)[-1]
            if sid in fail_sids:
                raise OSError("boom")
            base = "https://pk-db.com/media/data/"
            return _FakeResp(payload={"files": [base + f for f in _STUDY_FILES[sid]]})
        if "/pkdata/studies/" in url:
            sid = url.rstrip("/").rsplit("/", 1)[-1]
            return _FakeResp(payload=_STUDY_SUMMARIES[sid])
        if "/media/data/" in url:
            name = url.rsplit("/", 1)[-1]
            return _FakeResp(text=_TSV_TEXTS[name])
        raise AssertionError(f"假响应未覆盖的 URL: {url}")

    return fake


def _load(sid: str) -> dict:
    path = os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR, f"{sid}.json")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def test_build_writes_verbatim_payload(monkeypatch):
    """搬运语义:数值/单位逐字保存字符串、licence/access 保留、png 记入 skipped、
    meta 落盘。全程无任何换算(分层红线)。"""
    calls: list = []
    monkeypatch.setattr(pk_mirror, "request_with_retry", lambda u, **k: _mk_fake(calls)(u, **k))
    meta = pk_mirror.build_mirror()
    assert meta["studies_built"] == 2 and meta["studies_failed"] == 0
    payload = _load("PKDB000001")
    assert payload["licence"] == "closed" and payload["access"] == "public"
    assert payload["reference_pmid"] == "111"
    # 逐字保存:cmax 行的数值与单位保持原始字符串,不做摩尔换算
    tab = next(t for t in payload["tsv"] if t["file"].endswith("Tab1.tsv"))
    row = tab["rows"][1]
    assert row[8] == "359.2" and row[10] == "ng/ml"
    # png 记入 skipped(不静默丢弃)
    assert any("Fig1.png" in s["file"] for s in payload["skipped"])
    # 物质映射落盘(归一化层对齐 ChEMBL 的桥)
    smeta = os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR, "substance_map.json")
    with open(smeta, encoding="utf-8") as f:
        assert json.load(f)["nodes"][0]["sid"] == "drug1"


def test_resume_skips_existing(monkeypatch):
    """断点续传:已存在的 sid 在 force=False 下跳过,不发起该研究的任何请求。"""
    calls: list = []
    monkeypatch.setattr(pk_mirror, "request_with_retry", lambda u, **k: _mk_fake(calls)(u, **k))
    pk_mirror.build_mirror()
    media_calls_first = sum(1 for u in calls if "/media/data/" in u)
    meta = pk_mirror.build_mirror()
    assert meta["studies_skipped_existing"] == 2 and meta["studies_built"] == 0
    media_calls_second = sum(1 for u in calls if "/media/data/" in u) - media_calls_first
    assert media_calls_second == 0   # 续传路径零下载


def test_failed_study_continues_batch(monkeypatch):
    """单研究远端故障:不断批,失败计入 _meta.json,其余研究照常落盘。"""
    calls: list = []
    monkeypatch.setattr(pk_mirror, "request_with_retry",
                        lambda u, **k: _mk_fake(calls, fail_sids=("PKDB000002",))(u, **k))
    meta = pk_mirror.build_mirror()
    assert meta["studies_failed"] == 1
    assert meta["failures"][0]["sid"] == "PKDB000002"
    assert os.path.exists(os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR,
                                       "PKDB000001.json"))


def test_status_reads_local_state(monkeypatch):
    """状态查询只读本地目录与 _meta.json(built/studies/ttl),不联网。"""
    calls: list = []
    monkeypatch.setattr(pk_mirror, "request_with_retry", lambda u, **k: _mk_fake(calls)(u, **k))
    before = pk_mirror.mirror_status()
    assert before["built"] is False
    pk_mirror.build_mirror()
    status = pk_mirror.mirror_status()
    assert status["built"] is True and status["studies"] == 2
    assert status["ttl_days"] == 180 and status["stale"] is False


def test_index_failure_aborts(monkeypatch):
    """研究清单拉取失败 → MirrorIndexError(整批终止,不产出半个镜像)。
    模拟真实行为:request_with_retry 重试耗尽后返回 None,不抛异常。"""
    def fail_index(url, **kwargs):
        if url.endswith("/pkdata/studies/"):
            return None
        raise AssertionError("清单失败后不应有后续请求")

    monkeypatch.setattr(pk_mirror, "request_with_retry", fail_index)
    try:
        pk_mirror.build_mirror()
        raise AssertionError("应抛 MirrorIndexError")
    except pk_mirror.MirrorIndexError:
        pass


def test_api_pk_mirror_status_mounted():
    """路由已挂载:GET /api/pk-mirror 在沙箱下返回未构建状态(不联网)。"""
    from fastapi.testclient import TestClient

    from WY.api import app
    client = TestClient(app)
    r = client.get("/api/pk-mirror")
    assert r.status_code == 200
    body = r.json()
    assert body["built"] is False and body["stale"] is True
