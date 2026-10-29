"""core.chembl 的测试:pChEMBL 排序、基因名提取、编号查找回退、批量分页抓取
与人类靶点过滤、查无此药/请求中止的三态语义。

阶段三第 2 批新写。打桩层级在 _chembl_get(带缓存的 GET 统一入口),不联网。
"""
import pandas as pd
import pytest

import WY.core.chembl as chembl


class _FakeResp:
    def __init__(self, payload):
        self._payload = payload

    def json(self):
        return self._payload


def test_sort_by_pchembl_desc_nan_last():
    """pChEMBL 降序排列,无活性值的排最后。"""
    df = pd.DataFrame({"approvedSymbol": ["A", "B", "C"],
                       "pchembl_value": [5.0, None, 9.0]})
    out = chembl._sort_drug_df_by_pchembl(df)
    assert list(out["approvedSymbol"]) == ["C", "A", "B"]
    assert pd.isna(out.iloc[2]["pchembl_value"])


def test_extract_gene_name_priority():
    """基因名提取:GENE_SYMBOL 同义词优先,退而 HGNC 交叉引用,都没有返回 None。"""
    t_symbol = {"target_components": [{
        "target_component_synonyms": [{"syn_type": "UNIPROT", "component_synonym": "P04637"},
                                      {"syn_type": "GENE_SYMBOL", "component_synonym": "TP53"}]}]}
    assert chembl._extract_gene_name(t_symbol) == "TP53"
    t_hgnc = {"target_components": [{
        "target_component_synonyms": [{"syn_type": "UNIPROT", "component_synonym": "P04637"}],
        "target_component_xrefs": [{"xref_src_db": "HGNC", "xref_name": "TP53"}]}]}
    assert chembl._extract_gene_name(t_hgnc) == "TP53"
    assert chembl._extract_gene_name({"target_components": []}) is None


def test_find_chembl_id_exact_then_fuzzy(monkeypatch):
    """编号查找:精确命中直接返回;精确未中转模糊;两者皆空返回 None。"""
    calls = []

    def fake_get(path, params, label):
        calls.append(params)
        if "pref_name__iexact" in params:
            return {"molecules": []}                     # 精确未中
        return {"molecules": [{"molecule_chembl_id": "CHEMBL1431"}]}  # 模糊命中

    monkeypatch.setattr(chembl, "_chembl_get", fake_get)
    assert chembl._find_chembl_id("metformin") == "CHEMBL1431"
    assert len(calls) == 2  # 先精确后模糊

    monkeypatch.setattr(chembl, "_chembl_get", lambda p, q, l: {"molecules": []})
    assert chembl._find_chembl_id("不存在的药") is None


def _install_chembl_world(monkeypatch, *, molecule_found=True, activities_fail=False):
    """搭一个完整的假 ChEMBL:分子、2 条人类活性(含重复)、1 条非人类靶点。"""
    def fake_get(path, params, label):
        if path == "molecule.json":
            if not molecule_found:
                return {"molecules": []}
            return {"molecules": [{"molecule_chembl_id": "CHEMBL1431"}]}
        if path == "activity.json":
            if activities_fail:
                return None
            return {"activities": [
                {"target_chembl_id": "T1", "assay_type": "B", "pchembl_value": 5.0},
                {"target_chembl_id": "T1", "assay_type": "B", "pchembl_value": 5.0},  # 重复
                {"target_chembl_id": "T2", "assay_type": "B", "pchembl_value": 8.0},
                {"target_chembl_id": None, "assay_type": "B", "pchembl_value": 1.0},  # 无靶点号
            ],
                "page_meta": {"limit": 1000}}
        if path == "target.json":
            return {"targets": [
                {"target_chembl_id": "T1", "pref_name": "TP53", "organism": "Homo sapiens",
                 "target_components": [{"target_component_synonyms": [
                     {"syn_type": "GENE_SYMBOL", "component_synonym": "TP53"}]}]},
                {"target_chembl_id": "T2", "pref_name": "Mouse thing",
                 "organism": "Mus musculus",  # 非人类 → 过滤
                 "target_components": [{"target_component_synonyms": [
                     {"syn_type": "GENE_SYMBOL", "component_synonym": "MOUSEGENE"}]}]},
            ]}
        raise AssertionError(f"unexpected path {path}")

    monkeypatch.setattr(chembl, "_chembl_get", fake_get)


def test_fetch_full_flow_human_only_sorted_dedup(monkeypatch, tmp_path):
    """全链:编号 → 活性翻页 → 批量靶点 → 人类过滤 + 去重 + pChEMBL 降序 + 缓存落盘。"""
    _install_chembl_world(monkeypatch)
    df, csv_path = chembl.fetch_drug_targets_chembl("metformin")
    assert csv_path is not None
    assert list(df["approvedSymbol"]) == ["TP53"]  # 非人类与无靶点号全被滤掉
    assert df.iloc[0]["pchembl_value"] == 5.0
    assert len(df) == 1  # 同 (基因,实验类型,pchembl) 去重
    import os
    assert os.path.exists(csv_path), "结果缓存应落盘"


def test_fetch_unknown_drug_is_empty_with_path(monkeypatch):
    """查无此药属正常业务:空表 + csv 路径(非 None),refresh 链路记 empty。"""
    _install_chembl_world(monkeypatch, molecule_found=False)
    df, csv_path = chembl.fetch_drug_targets_chembl("不存在的药")
    assert df.empty and csv_path is not None
    monkeypatch.setattr(chembl, "set_chembl_cache_ttl", lambda h: None)
    df2, state = chembl.fetch_drug_targets_refresh("不存在的药")
    assert df2.empty and state == "empty"


def test_fetch_abort_on_activity_failure(monkeypatch):
    """活性翻页失败 → 统一失败出口 (空表, None),refresh 链路记 failed。"""
    _install_chembl_world(monkeypatch, activities_fail=True)
    df, csv_path = chembl.fetch_drug_targets_chembl("metformin")
    assert df.empty and csv_path is None
    df2, state = chembl.fetch_drug_targets_refresh("metformin")
    assert df2.empty and state == "failed"
