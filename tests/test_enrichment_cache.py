"""core.enrichment 的缓存键/文件名/迁移/归一化测试(纯本地,零网络)。

深水区冲线批。熔断与编排语义已由搬家组覆盖,这里补缓存治理侧。
"""
import os

import pandas as pd
import pytest

from WY import config
import WY.core.enrichment as en


def test_enrichr_libs_mapping():
    """g:Profiler 数据源 → Enrichr 库名映射:未收录的数据源被过滤。"""
    libs = en._enrichr_libs(("GO:BP", "KEGG", "GO:CC"))
    assert all(isinstance(x, str) for x in libs)
    assert libs == [en.ENRICHR_LIB_MAP["GO:BP"], en.ENRICHR_LIB_MAP["KEGG"],
                    en.ENRICHR_LIB_MAP["GO:CC"]]
    assert en._enrichr_libs(("bogus",)) == []


def test_enrichment_cache_key_deterministic():
    """缓存键:同输入同键;基因顺序/重复不影响;参数变则键变。"""
    k1 = en._enrichment_cache_key(["TP53", "EGFR"], "hsapiens", ("GO:BP",), 0.05)
    k2 = en._enrichment_cache_key(["EGFR", "TP53", "TP53"], "hsapiens", ("GO:BP",), 0.05)
    k3 = en._enrichment_cache_key(["TP53"], "hsapiens", ("GO:BP",), 0.05)
    assert k1 == k2 and k1 != k3


def test_cache_filename_label_vs_bare():
    """文件名:有标识用同名覆盖式;无标识退回短哈希命名防互相覆盖。"""
    named = en.enrichment_cache_filename("gProfiler", "aspirin_cancer", "a" * 32)
    assert named == "gProfiler_aspirin_cancer.csv"
    bare = en.enrichment_cache_filename("gProfiler", None, "a" * 32)
    assert bare == "gProfiler_" + "a" * 12 + ".csv"


def test_parse_and_label_roundtrip():
    """文件名解析回环:新命名/旧哈希命名/无标识裸哈希/GSEA 侧别四种形态;
    无法识别的文件名 parse 返回 None、label 原样。"""
    info = en.parse_enrichment_filename("gProfiler_aspirin_cancer_12345678.csv")
    assert info == {"kind": "ora", "source": "gProfiler",
                    "label": "aspirin_cancer", "side": None}
    info2 = en.parse_enrichment_filename("GSEA_drug_metformin_T2D.csv")
    assert info2["kind"] == "gsea" and info2["side"] == "drug"
    assert en.parse_enrichment_filename("random_file.csv") is None
    assert en.enrichment_filename_label("random_file.csv") == ""  # 无法识别 → 空串


def test_migrate_legacy_cache_names(tmp_path, monkeypatch):
    """旧哈希命名迁移:同标识多变体保留最新、改写为新格式;新格式文件参与最新者胜;幂等。"""
    monkeypatch.setattr(config, "ENRICHMENT_CACHE_DIR", str(tmp_path))
    # 旧格式两个变体 + 一个新格式同名
    (tmp_path / "gProfiler_aspirin_cancer_11111111.csv").write_text("old", encoding="utf-8")
    (tmp_path / "gProfiler_aspirin_cancer_22222222.csv").write_text("newer", encoding="utf-8")
    (tmp_path / "gProfiler_aspirin_cancer.csv").write_text("newfmt", encoding="utf-8")
    (tmp_path / "unrelated.csv").write_text("keep", encoding="utf-8")
    n = en.migrate_legacy_cache_names()
    names = set(os.listdir(str(tmp_path)))
    assert "gProfiler_aspirin_cancer.csv" in names
    assert "gProfiler_aspirin_cancer_11111111.csv" not in names
    assert "gProfiler_aspirin_cancer_22222222.csv" not in names
    assert "unrelated.csv" in names
    assert n >= 1
    assert en.migrate_legacy_cache_names() == 0  # 幂等


def test_normalize_enrichment_gprofiler():
    """g:Profiler 结果归一化:列名映射齐全;缺必需列返回空表而非报错。"""
    res = pd.DataFrame([
        {"source": "KEGG", "term_id": "hsa04151", "term_name": "PI3K-Akt",
         "p_value": 0.001, "adjusted_p_value": 0.002, "intersection_size": 3,
         "query_size": 5, "genes": "A;B;C"},
    ])
    out = en._normalize_enrichment_gprofiler(res)
    assert list(out["term_name"]) == ["PI3K-Akt"]
    # 实现契约:source 必有;其余列按优先级取,genes 取自 intersection(缺失置空)
    assert list(out["source"]) == ["KEGG"]
    assert "term_id" in out.columns and "p_value" in out.columns
    # 缺列回退契约:source 与其他列一样走 _pick,缺 p 值列 → 全行被 dropna 剔除,
    # 返回空表而不是 KeyError(2026-10 修复:source 曾被直接引用)
    bad = en._normalize_enrichment_gprofiler(pd.DataFrame({"foo": [1]}))
    assert bad.empty
