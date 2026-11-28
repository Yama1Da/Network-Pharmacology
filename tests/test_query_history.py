"""查询历史台账冒烟测试(2026-11-28 负责人提议,§1.14):
record/去重/损坏容错/首次种子/导入合并——文件落在 conftest 沙箱的
CACHE_DIR 下,测试间互不可见(每例重置)。
"""

import json
import os

import pytest

from WY.core import query_history as qh


@pytest.fixture(autouse=True)
def _fresh_ledger():
    """每例前清台账文件,保证种子/首录路径被真实走到。"""
    p = qh._path()
    if os.path.isfile(p):
        os.remove(p)
    yield
    if os.path.isfile(p):
        os.remove(p)


def test_record_dedupes_and_timestamps():
    # 键名大小写敏感,镜像 merged 缓存文件名键语义(两个大小写 = 两个缓存实体)
    qh.record("drug", "  Metformin ")
    qh.record("drug", "Metformin")  # 同名(去空格归一)只刷新 last_at
    qh.record("disease", "cancer")
    data = qh.load()
    assert len(data["drugs"]) == 1 and "Metformin" in data["drugs"]
    rec = data["drugs"]["Metformin"]
    assert rec["first_at"] <= rec["last_at"]
    assert set(data["diseases"]) == {"cancer"}
    # 无效入参静默不入账;不同大小写是不同缓存实体,如实分立
    qh.record("drug", "")
    qh.record("vehicle", "X")
    assert len(qh.load()["drugs"]) == 1


def test_names_sorted_and_unknown_kind_empty():
    qh.record("drug", "Zeta")
    qh.record("drug", "Alpha")
    assert qh.names("drug") == ["Alpha", "Zeta"]
    assert qh.names("vehicle") == []


def test_corrupt_file_treated_as_empty():
    qh.record("drug", "A")
    with open(qh._path(), "w", encoding="utf-8") as f:
        f.write("{broken json")
    # 损坏按空台账起底(以本次写入为准),不抛异常打断分析链路
    qh.record("drug", "B")
    assert set(qh.load()["drugs"]) == {"B"}


def test_merge_blob_union_never_deletes():
    qh.record("drug", "Local")
    before = qh.load()["drugs"]["Local"]
    blob = {"drugs": {"Local": {"first_at": "2000-01-01", "last_at": "2000-01-01"},
                      "Imported": {"first_at": "2001-01-01",
                                   "last_at": "2001-01-01"}},
            "diseases": {"Cancer": {"first_at": "2002-01-01",
                                    "last_at": "2002-01-01"}}}
    n = qh.merge_blob(blob)
    data = qh.load()
    assert n == 3 and set(data["drugs"]) == {"Local", "Imported"}
    assert set(data["diseases"]) == {"Cancer"}
    # 已有条目只刷 last_at,first_at/名字永不丢——台账只增不减
    assert data["drugs"]["Local"]["first_at"] == before["first_at"]


def test_export_blob_none_when_empty():
    assert qh.export_blob() is None
    qh.record("disease", "cancer")
    blob = qh.export_blob()
    assert blob is not None and set(blob["diseases"]) == {"cancer"}
