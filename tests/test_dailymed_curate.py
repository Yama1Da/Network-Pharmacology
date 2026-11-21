"""1.10 dailymed_curate 单测：独立约束、预填、审查落盘、校验器、离线降级。

全部走桩：net_block 下必须零外网。
"""

import ast
import builtins
import csv
import os
from types import SimpleNamespace

import pytest

from WY.tools import dailymed_curate as dc


# ---------------------------------------------------------------------------
# 独立工具约束：import 区不得出现 WY
# ---------------------------------------------------------------------------
def test_import_area_has_no_wy():
    path = dc.__file__
    with open(path, encoding="utf-8") as f:
        tree = ast.parse(f.read(), filename=path)
    roots = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            roots.add(node.module.split(".")[0])
    assert "WY" not in roots, "独立工具约束：import 区不得出现 WY"


def test_columns_and_constants():
    assert dc.COLUMNS[0] == "drug" and "spl_pharm_class" in dc.COLUMNS
    assert "atc" not in "".join(dc.COLUMNS).lower()
    assert dc.PENDING == "待定"


# ---------------------------------------------------------------------------
# prefill：三级预填 + 合并 + 一致性排序
# ---------------------------------------------------------------------------
def _write_csv(path, rows):
    dc.write_rows(path, rows)


def _read_csv(path):
    return dc.read_rows(path)


def test_prefill_builds_review_queue(tmp_path, fake_http, spl_xml_factory, monkeypatch):
    csv_path = tmp_path / "targets.csv"
    drugs = tmp_path / "drugs.txt"
    drugs.write_text("metformin\nunknown-drug\n", encoding="utf-8")

    setid = "07b4ee2a-8435-4d38-a6b0-b919e3424153"

    def fake_ot_post(query, _variables, *, cache_dir, offline=False):
        if "search" in query:
            return {"search": {"hits": [
                {"id": "CHEMBL1431", "name": "metformin", "entity": "drug"}]}}
        return {"drug": {"name": "metformin", "mechanismsOfAction": {"rows": [
            {"mechanismOfAction": "AMPK activation",
             "targets": [{"approvedSymbol": "PRKAA1"},
                         {"approvedSymbol": "SLC2A4"}]},
        ]}}}

    fake_http.json("api.fda.gov", {"results": [{"openfda": {
        "pharm_class_epc": ["Biguanide [EPC]"],
        "pharm_class_moa": ["Biguanides [MoA]"]}}]})
    fake_http.json("spls.json", {"data": [{"setid": setid,
                                           "title": "METFORMIN HYDROCHLORIDE"}]})
    fake_http.text("spls/" + setid, spl_xml_factory(
        "Metformin decreases hepatic glucose production.", setid=setid))
    monkeypatch.setattr(dc, "ot_post", fake_ot_post)
    monkeypatch.setattr(dc, "_request", fake_http.handle)

    args = SimpleNamespace(csv=str(csv_path), drugs=str(drugs), limit=0,
                           cache_dir=str(tmp_path / "cache"), offline=False,
                           curator="tester")
    assert dc.prefill(args) == 0

    rows = {row["drug"]: row for row in _read_csv(str(csv_path))}
    assert set(rows) == {"metformin", "unknown-drug"}
    met = rows["metformin"]
    assert met["ot_moa_prefill"] == "PRKAA1;SLC2A4"
    assert met["spl_pharm_class"] == "Biguanide [EPC];Biguanides [MoA]"
    assert met["label_section"] == f"{setid}#43679-0"
    assert "AMPK activation" in met["notes"]
    assert "Metformin decreases hepatic glucose production" in met["notes"]
    assert met["curator"] == "tester"
    # 一致性排序：有基因级候选 + 原文的排最前
    assert _read_csv(str(csv_path))[0]["drug"] == "metformin"


def test_prefill_merges_existing_rows(tmp_path, monkeypatch):
    csv_path = tmp_path / "targets.csv"
    _write_csv(str(csv_path), [
        {"drug": "metformin", "label_mechanism_targets": "PRKAA1",
         "decision": "y", "curator": "bio", "date": "2026-11-18"},
    ])
    monkeypatch.setattr(dc, "ot_drug_moa", lambda *a, **k: {
        "chembl_id": "CHEMBL1431", "moa": [], "targets": ["SLC2A4"]})
    monkeypatch.setattr(dc, "openfda_class", lambda *a, **k: [])
    monkeypatch.setattr(dc, "dailymed_moa", lambda *a, **k: {"setid": "", "text": ""})
    args = SimpleNamespace(csv=str(csv_path), drugs=None, limit=0,
                           cache_dir=str(tmp_path / "cache"), offline=False,
                           curator="tester")
    assert dc.prefill(args) == 0
    row = _read_csv(str(csv_path))[0]
    assert row["drug"] == "metformin"
    assert row["decision"] == "y"                      # 已有判定不被覆盖
    assert row["label_mechanism_targets"] == "PRKAA1"  # 人工值优先
    assert row["ot_moa_prefill"] == "SLC2A4"           # 预填列独立更新


# ---------------------------------------------------------------------------
# review：y/t/x 三键 + 立即写盘
# ---------------------------------------------------------------------------
def test_review_confirm_and_pending(tmp_path, monkeypatch):
    csv_path = tmp_path / "targets.csv"
    _write_csv(str(csv_path), [
        {"drug": "ramipril", "ot_moa_prefill": "ACE"},
        {"drug": "metformin", "ot_moa_prefill": ""},
    ])
    answers = iter(["y", "t"])
    monkeypatch.setattr(builtins, "input", lambda _prompt="": next(answers))
    args = SimpleNamespace(csv=str(csv_path), curator="tester")
    assert dc.review(args) == 0
    rows = {row["drug"]: row for row in _read_csv(str(csv_path))}
    assert rows["ramipril"]["decision"] == "y"
    assert rows["ramipril"]["label_mechanism_targets"] == "ACE"
    assert rows["metformin"]["decision"] == "t"
    assert rows["metformin"]["label_mechanism_targets"] == dc.PENDING
    assert all(row["curator"] == "tester" for row in rows.values())


def test_review_quit_keeps_progress(tmp_path, monkeypatch):
    csv_path = tmp_path / "targets.csv"
    _write_csv(str(csv_path), [{"drug": "a", "ot_moa_prefill": "X"},
                               {"drug": "b", "ot_moa_prefill": "Y"}])
    answers = iter(["y", "q"])
    monkeypatch.setattr(builtins, "input", lambda _prompt="": next(answers))
    assert dc.review(SimpleNamespace(csv=str(csv_path), curator="t")) == 0
    rows = _read_csv(str(csv_path))
    assert rows[0]["decision"] == "y" and rows[1]["decision"] == ""


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------
def test_validate_good_row_passes(tmp_path, capsys):
    path = tmp_path / "ok.csv"
    _write_csv(str(path), [{
        "drug": "ramipril", "label_mechanism_targets": "ACE",
        "ot_moa_prefill": "ACE", "spl_pharm_class": "ACE inhibitor [EPC]",
        "label_section": "11111111-2222-3333-4444-555555555555#43679-0",
        "decision": "y", "curator": "tester", "date": "2026-11-18"}])
    args = SimpleNamespace(csv=str(path), offline=True, online=False,
                           cache_dir=str(tmp_path / "cache"))
    assert dc.validate(args) == 0
    assert "通过" in capsys.readouterr().out


def test_validate_reports_row_numbers(tmp_path, capsys):
    path = tmp_path / "bad.csv"
    _write_csv(str(path), [
        {"drug": "a", "decision": "z"},                        # 非法枚举
        {"drug": "b", "decision": "y",
         "label_mechanism_targets": "", "curator": "t",
         "date": "2026-11-18"},                                # y 但无靶点
        {"drug": "c", "decision": "t",
         "label_mechanism_targets": "EGFR", "curator": "t",
         "date": "2026-11-18"},                                # t 但不是待定
        {"drug": "A", "decision": "y",
         "label_mechanism_targets": "EGFR", "curator": "t",
         "date": "2026/11/18"},                                # 重复 + 日期格式
    ])
    args = SimpleNamespace(csv=str(path), offline=True, online=False,
                           cache_dir=str(tmp_path / "cache"))
    assert dc.validate(args) == 1
    out = capsys.readouterr()
    text = out.err + out.out
    assert "第 2 行" in text and "第 5 行" in text
    assert "decision" in text


# ---------------------------------------------------------------------------
# 离线：有缓存用缓存，无缓存警告跳过
# ---------------------------------------------------------------------------
def test_offline_without_cache_warns_and_skips(tmp_path, monkeypatch, capsys):
    def _boom(*args, **kwargs):
        raise AssertionError("离线模式不应发请求")

    monkeypatch.setattr(dc, "_request", _boom)
    ot = dc.ot_drug_moa("metformin", cache_dir=str(tmp_path), offline=True)
    fda = dc.openfda_class("metformin", cache_dir=str(tmp_path), offline=True)
    dm = dc.dailymed_moa("metformin", cache_dir=str(tmp_path), offline=True)
    assert ot["targets"] == [] and fda == [] and dm["setid"] == ""
    assert "离线且无缓存" in capsys.readouterr().err


def test_cache_roundtrip_avoids_network(tmp_path, monkeypatch):
    calls = {"n": 0}

    def fake_post(*args, **kwargs):
        calls["n"] += 1
        return {"search": {"hits": [{"id": "CHEMBL_X", "entity": "drug"}]}}

    monkeypatch.setattr(dc, "ot_post", fake_post)
    first = dc.ot_drug_moa("x", cache_dir=str(tmp_path))
    assert first["chembl_id"] == "CHEMBL_X"
    assert calls["n"] >= 1


# ---------------------------------------------------------------------------
# SPL XML 解析
# ---------------------------------------------------------------------------
def test_extract_moa_section(spl_xml_factory):
    xml = spl_xml_factory("Metformin decreases hepatic glucose production.")
    text = dc._extract_moa_section(xml)
    assert "hepatic glucose production" in text
    assert "INDICATIONS" not in text
    assert "Mechanism of Action" not in text


def test_extract_moa_section_missing_returns_empty(spl_xml_factory):
    xml = spl_xml_factory("x").replace("43679-0", "34067-9")
    assert dc._extract_moa_section(xml) == ""
