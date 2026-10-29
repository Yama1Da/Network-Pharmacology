"""core.interactome 的元数据/边表存取/版本缓存治理测试(纯本地,零网络)。

深水区冲线批。缓存目录由 conftest 沙箱接管;不构建真实背景网络。
"""
import json
import os
from datetime import datetime, timedelta

import pandas as pd
import pytest

from WY import config
from WY.core import interactome as itx


def _write_meta(species="Human", confidence=700, version="v12.0",
                *, with_files=True, built_at=None):
    """在沙箱缓存目录写一份完整的 _meta.json(可选连带边表/度表文件)。"""
    conf = int(confidence)
    meta = {"species": species, "confidence": conf, "string_version": version,
            "built_at": (built_at or datetime.now().isoformat(timespec="seconds")),
            "files": {"edges": f"{species}_{conf}_{version}_edges.parquet",
                      "degrees": f"{species}_{conf}_{version}_degrees.csv"}}
    os.makedirs(itx._cache_dir(), exist_ok=True)
    if with_files:
        edges = pd.DataFrame({"protein_a": ["A", "B"], "protein_b": ["B", "C"],
                              "score": [900, 800]})
        edges.to_parquet(os.path.join(itx._cache_dir(), meta["files"]["edges"]),
                         index=False)
        pd.DataFrame({"gene": ["A", "B", "C"], "degree": [2, 2, 1]}).to_csv(
            os.path.join(itx._cache_dir(), meta["files"]["degrees"]), index=False)
    path = os.path.join(itx._cache_dir(), itx._meta_path(species, conf, version))
    with open(path, "w", encoding="utf-8") as f:
        json.dump(meta, f)
    return meta


def test_is_built_requires_files(tmp_path, monkeypatch):
    """构建判定:meta 存在且边表/度表文件齐全才算 built;缺文件(登记在案但盘上没有)不算。"""
    monkeypatch.setattr(config, "INTERACTOME_CONFIDENCE", 700)
    monkeypatch.setattr(config, "STRING_VERSION", "v12.0")
    assert itx.is_built() is False                       # 无 meta
    _write_meta(with_files=False)                        # meta 在,文件缺
    assert itx.is_built() is False
    _write_meta()                                        # 文件齐全
    assert itx.is_built() is True
    assert itx.is_built(confidence=500) is False         # 其他置信度组合未构建


def test_require_meta_raises_when_missing(tmp_path, monkeypatch):
    """未构建时 require_meta 报错(绝不静默降级);已构建时返回元数据。"""
    monkeypatch.setattr(config, "INTERACTOME_CONFIDENCE", 700)
    monkeypatch.setattr(config, "STRING_VERSION", "v12.0")
    with pytest.raises(itx.InteractomeNotBuiltError):
        itx.require_meta()
    _write_meta()
    meta = itx.require_meta()
    assert meta["species"] == "Human" and meta["confidence"] == 700


def test_load_interactome_and_degrees(tmp_path, monkeypatch):
    """加载:图按边表构建(只读冻结)、度表列齐全;路径走 _load_graph_cached。"""
    monkeypatch.setattr(config, "INTERACTOME_CONFIDENCE", 700)
    monkeypatch.setattr(config, "STRING_VERSION", "v12.0")
    _write_meta()
    g = itx.load_interactome()
    assert g.number_of_nodes() == 3 and g.number_of_edges() == 2
    import networkx as nx
    assert nx.is_frozen(g), "背景网络必须只读,防后续计算误改"
    deg = itx.load_degree_table()
    assert list(deg.columns) == ["gene", "degree"] and len(deg) == 3


def test_interactome_age_and_stale(tmp_path, monkeypatch):
    """年龄/陈旧:按 built_at 计算;陈旧阈值是严格大于;未构建 age 为 None。"""
    monkeypatch.setattr(config, "INTERACTOME_CONFIDENCE", 700)
    monkeypatch.setattr(config, "STRING_VERSION", "v12.0")
    assert itx.interactome_age_days() is None and itx.is_stale() is False
    old = (datetime.now() - timedelta(days=400)).isoformat(timespec="seconds")
    _write_meta(built_at=old)
    assert itx.interactome_age_days() == 400
    assert itx.is_stale(days=365) is True
    assert itx.is_stale(days=10000) is False


def test_write_read_edges_auto_parquet(tmp_path):
    """边表自动格式:parquet 主路径落盘读回一致;CSV 旁路读回一致。"""
    df = pd.DataFrame({"protein_a": ["A"], "protein_b": ["B"], "score": [900]})
    path = itx._write_edges_auto(df, str(tmp_path / "edges"))
    assert path.endswith(".parquet")
    pd.testing.assert_frame_equal(itx._read_edges_auto(path), df)
    csv_path = str(tmp_path / "edges2.csv")
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    pd.testing.assert_frame_equal(itx._read_edges_auto(csv_path), df)


def test_clean_old_version_caches_keeps_current(tmp_path, monkeypatch):
    """版本清理:只删"其它版本/置信度"组合的三件套;当前组合与读不出的记录保守保留。"""
    monkeypatch.setattr(config, "INTERACTOME_CONFIDENCE", 700)
    monkeypatch.setattr(config, "STRING_VERSION", "v12.0")
    _write_meta(version="v12.0", confidence=700)          # 当前组合 → 保留
    stale = _write_meta(version="v11.5", confidence=700)  # 旧版本 → 删
    broken_dir = itx._cache_dir()
    with open(os.path.join(broken_dir, "Human_700_v11.0_meta.json"), "w") as f:
        f.write("{not json")                              # 读不出 → 保守保留
    n = itx.clean_old_version_caches()
    cur_meta = itx._meta_path("Human", 700, "v12.0")
    assert os.path.exists(os.path.join(broken_dir, cur_meta))
    assert not os.path.exists(os.path.join(broken_dir, stale["files"]["edges"]))
    assert not os.path.exists(os.path.join(
        broken_dir, itx._meta_path("Human", 700, "v11.5")))
    assert os.path.exists(os.path.join(broken_dir, "Human_700_v11.0_meta.json"))
    assert n >= 1
