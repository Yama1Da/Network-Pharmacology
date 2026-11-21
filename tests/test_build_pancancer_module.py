"""1.9 build_pancancer_module 单测：并集构建、验收断言、CancerMine 解析、
本地快照优先（Zenodo 不可达时的主路径）与离线语义。

全部打桩（net_block 下零外网）；不触碰真实 Static/ 与 WY/resources/。
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from WY.tools import build_pancancer_module as bpm


def _patch_sources(monkeypatch, ot_symbols, tsv):
    monkeypatch.setattr(bpm, "fetch_ot_targets",
                        lambda _disease_id=bpm.OT_DISEASE_ID: (list(ot_symbols), "cancer"))
    monkeypatch.setattr(bpm, "download_cancermine",
                        lambda url=None: (tsv, "stub://cancermine"))


def test_build_union_writes_parquet_and_meta(tmp_path, monkeypatch, cancermine_factory):
    tsv = cancermine_factory([
        ("TP53", "tumor_suppressor", 50),
        ("PIK3CA", "oncogene", 30),
        ("LOWCITE", "driver", 1),          # 引用数不足 → 过滤
        ("PASSENGER", "passenger", 99),    # 角色不符 → 过滤
    ])
    _patch_sources(monkeypatch, ["BRAF", "TP53", "EGFR"], tsv)
    snap = tmp_path / "snap.csv"

    meta = bpm.build(force=True, out=str(tmp_path), snapshot_path=str(snap))

    assert meta["ot_targets"] == 3
    assert meta["cancermine_genes"] == 2
    assert meta["union_genes"] == 4          # BRAF/TP53/EGFR + PIK3CA
    assert meta["overlap"] == 1              # TP53
    assert meta["markers_ok"] is True

    parquet = tmp_path / meta["parquet"]
    assert parquet.is_file()
    df = pd.read_parquet(parquet)
    assert list(df.columns) == ["gene", "sources"]
    assert set(df["gene"]) == {"BRAF", "TP53", "EGFR", "PIK3CA"}
    assert dict(zip(df["gene"], df["sources"]))["TP53"] == "OT;CANCERMINE"

    saved = json.loads((tmp_path / "meta.json").read_text(encoding="utf-8"))
    assert saved["ot_disease_id"] == bpm.OT_DISEASE_ID
    assert saved["cancermine_min_citations"] == 2


def test_snapshot_is_preferred_and_avoids_network(tmp_path, monkeypatch):
    """快照存在时，force 重建也不联网（Zenodo 被墙的主场景）。"""
    snap = tmp_path / "snap.csv"
    bpm.write_snapshot(str(snap),
                       [{"gene": "TP53", "roles": "tumor_suppressor", "citations": 9},
                        {"gene": "PIK3CA", "roles": "oncogene", "citations": 9}],
                       source="unit-test", min_citations=2)
    monkeypatch.setattr(bpm, "fetch_ot_targets",
                        lambda _disease_id=bpm.OT_DISEASE_ID: (["BRAF"], "cancer"))

    def _no_download(url=None):
        raise AssertionError("有快照时不应联网")

    monkeypatch.setattr(bpm, "download_cancermine", _no_download)
    meta = bpm.build(force=True, out=str(tmp_path), snapshot_path=str(snap))
    assert meta["cancermine_source"].startswith("snapshot:")
    assert meta["cancermine_genes"] == 2


def test_cancermine_file_wins_and_writes_snapshot(tmp_path, monkeypatch,
                                                  cancermine_factory):
    tsv_path = tmp_path / "cancermine_collated.tsv"
    tsv_path.write_text(cancermine_factory([
        ("BRAF", "driver", 10), ("TP53", "tumor_suppressor", 9),
        ("PIK3CA", "oncogene", 9)]), encoding="utf-8")
    snap = tmp_path / "snap.csv"
    monkeypatch.setattr(bpm, "fetch_ot_targets",
                        lambda _disease_id=bpm.OT_DISEASE_ID: (["BRAF"], "cancer"))

    meta = bpm.build(force=True, out=str(tmp_path), snapshot_path=str(snap),
                     cancermine_file=str(tsv_path), update_snapshot=True)
    assert meta["cancermine_source"].startswith("file:")
    assert "snapshot:" in meta["cancermine_source"]
    assert snap.is_file()
    assert {row["gene"] for row in bpm.read_snapshot(str(snap))} == {
        "BRAF", "TP53", "PIK3CA"}


def test_refresh_ignores_snapshot(tmp_path, monkeypatch, cancermine_factory):
    snap = tmp_path / "snap.csv"
    bpm.write_snapshot(str(snap), [{"gene": "TP53", "roles": "driver", "citations": 3}],
                       source="old", min_citations=2)
    _patch_sources(monkeypatch, ["BRAF"], cancermine_factory([
        ("PIK3CA", "oncogene", 5), ("TP53", "driver", 5)]))
    meta = bpm.build(force=True, refresh=True, out=str(tmp_path),
                     snapshot_path=str(snap), update_snapshot=True)
    assert not meta["cancermine_source"].startswith("snapshot:")
    assert meta["cancermine_genes"] == 2
    assert {row["gene"] for row in bpm.read_snapshot(str(snap))} == {"PIK3CA", "TP53"}


def test_read_snapshot_missing_returns_empty(tmp_path):
    assert bpm.read_snapshot(str(tmp_path / "nope.csv")) == []


def test_read_snapshot_rejects_bad_format(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("# c\nfoo,bar\n1,2\n", encoding="utf-8")
    with pytest.raises(bpm.PancancerError):
        bpm.read_snapshot(str(bad))


def test_build_fails_when_marker_missing(tmp_path, monkeypatch, cancermine_factory):
    _patch_sources(monkeypatch, ["BRAF", "EGFR"],
                   cancermine_factory([("TP53", "driver", 9)]))
    with pytest.raises(bpm.PancancerError) as exc:
        bpm.build(force=True, out=str(tmp_path),
                  snapshot_path=str(tmp_path / "snap.csv"))
    assert "PIK3CA" in str(exc.value)


def test_build_rejects_empty_union(tmp_path, monkeypatch, cancermine_factory):
    _patch_sources(monkeypatch, [], cancermine_factory([]))
    with pytest.raises(bpm.PancancerError):
        bpm.build(force=True, out=str(tmp_path),
                  snapshot_path=str(tmp_path / "snap.csv"))


def test_parse_filters_roles_and_citations(cancermine_factory):
    tsv = cancermine_factory([
        ("BRAF", "driver", 10),
        ("TP53", "tumor_suppressor", 8),
        ("NOTCH1", "oncogene", 0),          # 引用数不足
        ("FAKE", "biomarker", 100),         # 角色不符
    ])
    assert bpm.parse_cancermine(tsv, min_citations=2) == ["BRAF", "TP53"]


def test_parse_requires_gene_and_role_columns():
    with pytest.raises(bpm.PancancerError):
        bpm.parse_cancermine("foo\tbar\n1\t2\n")


def test_offline_without_any_source_fails_fast(tmp_path):
    with pytest.raises(bpm.PancancerError) as exc:
        bpm.build(offline=True, out=str(tmp_path),
                  snapshot_path=str(tmp_path / "snap.csv"))
    assert "离线" in str(exc.value)


def test_download_failure_mentions_manual_vpn_path(tmp_path, monkeypatch):
    monkeypatch.setattr(bpm, "resolve_cancermine_url", lambda: "https://zenodo/x")
    monkeypatch.setattr(bpm, "request_with_retry", lambda *a, **k: None)
    with pytest.raises(bpm.PancancerError) as exc:
        bpm.download_cancermine("https://zenodo/x")
    msg = str(exc.value)
    assert "--cancermine-file" in msg and "VPN" in msg


def test_resolve_cancermine_url_falls_back(monkeypatch):
    monkeypatch.setattr(bpm, "request_with_retry", lambda *a, **k: None)
    assert bpm.resolve_cancermine_url() == bpm.CANCERMINE_FALLBACK_URL


def test_resolve_cancermine_url_reads_zenodo_files(monkeypatch):
    class _Resp:
        def json(self):
            return {"hits": {"hits": [{"files": [
                {"key": bpm.CANCERMINE_FILE,
                 "links": {"self": "https://zenodo.org/api/files/x/cancermine_collated.tsv"}},
                {"key": "other.tsv", "links": {"self": "https://example/other"}},
            ]}]}}

    monkeypatch.setattr(bpm, "request_with_retry", lambda *a, **k: _Resp())
    assert bpm.resolve_cancermine_url().endswith(bpm.CANCERMINE_FILE)


def test_ot_empty_fetch_raises_and_skips_cache(tmp_path, monkeypatch):
    """OT 抓到空（HTTP 200 但疾病查空，2026-11-18 MONDO 迁移实证）必须当场抛错
    且不写缓存——空缓存落盘会让后续重建复用空臂（单臂降级无感知）。"""
    snap = tmp_path / "snap.csv"
    monkeypatch.setattr(bpm, "fetch_ot_targets",
                        lambda _disease_id=bpm.OT_DISEASE_ID: ([], "cancer"))

    def _no_download(url=None):
        raise AssertionError("OT 空结果应先抛错，不该走到 CancerMine")

    monkeypatch.setattr(bpm, "download_cancermine", _no_download)
    with pytest.raises(bpm.PancancerError) as exc:
        bpm.build(force=True, out=str(tmp_path), snapshot_path=str(snap))
    assert "一条靶点都没拿到" in str(exc.value)
    raw = Path(bpm.data_dir()) / "raw"
    assert not list(raw.glob("ot_*.json"))          # 缓存未落盘


def test_ot_empty_cache_is_treated_as_miss(tmp_path, monkeypatch, cancermine_factory):
    """历史坏缓存（空 symbols）不再被信任：视作未命中走重抓。"""
    raw = Path(bpm.data_dir()) / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    (raw / f"ot_{bpm.OT_DISEASE_ID}.json").write_text(
        json.dumps({"symbols": [], "disease_name": ""}), encoding="utf-8")
    _patch_sources(monkeypatch, ["BRAF", "TP53", "PIK3CA"],
                   cancermine_factory([("TP53", "driver", 9)]))
    meta = bpm.build(out=str(tmp_path), snapshot_path=str(tmp_path / "snap.csv"))
    assert meta["ot_targets"] == 3                  # 重抓成功，空缓存被绕过


def test_acceptance_gate_rejects_missing_arm(tmp_path, monkeypatch, cancermine_factory):
    """验收门 1：任一源零命中即拒绝（CancerMine 单臂 ~4000 条数量可观，
    不拦住则泛癌对照悄悄降级）。"""
    tsv = tmp_path / "filtered.tsv"
    tsv.write_text(cancermine_factory([("FAKE", "biomarker", 99)]), encoding="utf-8")
    monkeypatch.setattr(bpm, "fetch_ot_targets",
                        lambda _disease_id=bpm.OT_DISEASE_ID: (["BRAF"], "cancer"))
    with pytest.raises(bpm.PancancerError) as exc:
        bpm.build(force=True, out=str(tmp_path),
                  snapshot_path=str(tmp_path / "s.csv"), cancermine_file=str(tsv))
    assert "CancerMine" in str(exc.value) and "缺臂" in str(exc.value)
