"""core.targets 的测试:四源降噪过滤、双源合并、交集规则、并发编排与失败语义。

阶段三第 3 批新写(仿原方案示例 C:纯 DataFrame 为主,打桩仅用于并发编排)。
"""
import math
import os

import pandas as pd
import pytest

from WY import config
from WY.core.text_utils import safe_name
import WY.core.targets as tg


# ---------- 分数聚合 ----------


def test_max_score_map_uppercase_and_max():
    """基因符号统一大写、同基因多行取最大分;全为 None 的基因是 NaN(不是被剔除)。"""
    df = pd.DataFrame([
        {"approvedSymbol": "tp53", "score_pchembl": 5.0},
        {"approvedSymbol": "TP53", "score_pchembl": 9.0},
        {"approvedSymbol": "EGFR", "score_pchembl": None},
    ])
    m = tg._max_score_map(df, "score_pchembl")
    assert m["TP53"] == 9.0
    assert math.isnan(m["EGFR"])
    assert tg._max_score_map(pd.DataFrame(), "score_pchembl") == {}
    assert tg._max_score_map(df, "no_such_col") == {}


# ---------- 降噪过滤(阈值 monkeypatch 钉死,断言确定) ----------


def test_filter_dgidb_needs_score_and_pmid(monkeypatch):
    """DGIdb 降噪:交互分数达标 **且** 有 PMID 支撑才保留(两条都要)。"""
    monkeypatch.setattr(config, "DGIDB_MIN_SCORE", 0.3)
    df = pd.DataFrame([
        {"approvedSymbol": "TP53", "interaction_score": 0.9, "interaction_pmids": "123"},
        {"approvedSymbol": "EGFR", "interaction_score": 0.9, "interaction_pmids": ""},
        {"approvedSymbol": "KRAS", "interaction_score": 0.1, "interaction_pmids": "456"},
        {"approvedSymbol": "MTOR", "interaction_score": "x", "interaction_pmids": "789"},
    ])
    out = tg._filter_dgidb_targets(df)
    assert list(out["approvedSymbol"]) == ["TP53"]
    assert tg._filter_dgidb_targets(pd.DataFrame()).empty  # 空表原样


def test_filter_ot_score_and_biotype(monkeypatch):
    """Open Targets 降噪:总分达标 **且** 只留蛋白编码。"""
    monkeypatch.setattr(config, "OT_MIN_SCORE", 0.1)
    df = pd.DataFrame([
        {"approvedSymbol": "TP53", "association_score": 0.5, "biotype": "protein_coding"},
        {"approvedSymbol": "LINC1", "association_score": 0.5, "biotype": "lincRNA"},
        {"approvedSymbol": "KRAS", "association_score": 0.05, "biotype": "protein_coding"},
    ])
    assert list(tg._filter_ot_targets(df)["approvedSymbol"]) == ["TP53"]


def test_filter_chembl_blacklist_and_pchembl(monkeypatch):
    """ChEMBL 降噪:黑名单先剔;pChEMBL 达阈值保留;无定量基因须 DGIdb 佐证;
    仅有弱定量活性的基因剔除。"""
    monkeypatch.setattr(config, "CHEMBL_TARGET_BLACKLIST", ("ALB",))
    monkeypatch.setattr(config, "CHEMBL_PCHEMBL_THRESHOLD", 6.0)
    monkeypatch.setattr(config, "CHEMBL_CORROBORATION_REQUIRED", True)
    df = pd.DataFrame([
        {"approvedSymbol": "TP53", "pchembl_value": 7.0},   # 达阈值 → 留
        {"approvedSymbol": "EGFR", "pchembl_value": 5.0},   # 弱定量 → 剔
        {"approvedSymbol": "MTOR", "pchembl_value": None},  # 无定量 + 佐证 → 留
        {"approvedSymbol": "ALB", "pchembl_value": 9.0},    # 黑名单 → 剔
    ])
    out = tg._filter_chembl_targets(df, {"MTOR"})
    assert set(out["approvedSymbol"]) == {"TP53", "MTOR"}
    # 佐证开关关闭:无定量基因无条件保留
    monkeypatch.setattr(config, "CHEMBL_CORROBORATION_REQUIRED", False)
    out2 = tg._filter_chembl_targets(df, set())
    assert "MTOR" in set(out2["approvedSymbol"]) and "EGFR" not in set(out2["approvedSymbol"])
    # blacklist_only:跳过 pChEMBL 逻辑,只剔黑名单
    out3 = tg._filter_chembl_targets(df, set(), blacklist_only=True)
    assert set(out3["approvedSymbol"]) == {"TP53", "EGFR", "MTOR"}


def test_prepare_drug_sources_fallback(monkeypatch):
    """两源过滤后全空但原始数据非空 → ChEMBL 回退为仅黑名单过滤(避免整侧为空)。"""
    monkeypatch.setattr(config, "DGIDB_MIN_SCORE", 0.9)
    monkeypatch.setattr(config, "CHEMBL_PCHEMBL_THRESHOLD", 9.0)
    monkeypatch.setattr(config, "CHEMBL_TARGET_BLACKLIST", ())
    monkeypatch.setattr(config, "CHEMBL_CORROBORATION_REQUIRED", True)
    chmbl = pd.DataFrame([{"approvedSymbol": "TP53", "pchembl_value": 1.0}])   # 仅弱定量
    dgidb = pd.DataFrame([{"approvedSymbol": "TP53", "interaction_score": 0.1,
                           "interaction_pmids": "1"}])                          # 分数不足
    chmbl_f, dgidb_f = tg._prepare_drug_sources("metformin", chmbl, dgidb)
    assert dgidb_f.empty
    assert list(chmbl_f["approvedSymbol"]) == ["TP53"]  # 回退后保留(黑名单为空)


def test_prepare_disease_sources_truncates(monkeypatch):
    """疾病侧规模控制:OT 截前 OT_TOP_N、DisGeNET 截前 页数×页大小,OT 再过降噪。"""
    monkeypatch.setattr(config, "OT_TOP_N", 2)
    monkeypatch.setattr(config, "OT_MIN_SCORE", 0.1)
    monkeypatch.setattr(config, "DISGENET_PAGES", 1)
    monkeypatch.setattr(config, "DISGENET_PAGE_SIZE", 2)
    ot = pd.DataFrame([
        {"approvedSymbol": "A", "association_score": 0.9, "biotype": "protein_coding"},
        {"approvedSymbol": "B", "association_score": 0.8, "biotype": "protein_coding"},
        {"approvedSymbol": "C", "association_score": 0.7, "biotype": "protein_coding"},  # 被截掉
    ])
    dg = pd.DataFrame([{"approvedSymbol": g} for g in ("X", "Y", "Z")])
    ot_f, dg_f = tg._prepare_disease_sources("T2D", ot, dg)
    assert list(ot_f["approvedSymbol"]) == ["A", "B"]
    assert list(dg_f["approvedSymbol"]) == ["X", "Y"]


# ---------- 双源合并与交集 ----------


def test_merge_source_pair_union_and_flags(tmp_path):
    """按基因取并集:标记来源、source_count、各源分数最大值;缓存落盘。"""
    a = pd.DataFrame([
        {"approvedSymbol": "TP53", "assay_type": "F", "pchembl_value": 7.0},
        {"approvedSymbol": "EGFR", "assay_type": "F", "pchembl_value": 5.0},
    ])
    b = pd.DataFrame([{"approvedSymbol": "TP53", "interaction_score": 0.9},
                      {"approvedSymbol": "KRAS", "interaction_score": 0.4}])
    cache = str(tmp_path / "merged_x.csv")
    m = tg._merge_source_pair("metformin", a, b, cache_path=cache, **tg._DRUG_MERGE)
    assert list(m["approvedSymbol"]) == ["EGFR", "KRAS", "TP53"]  # 字母序并集
    tp = m.set_index("approvedSymbol").loc["TP53"]
    assert tp["source_chmbl"] == 1 and tp["source_dgidb"] == 1 and tp["source_count"] == 2
    assert tp["score_pchembl"] == 7.0 and tp["score_dgidb"] == 0.9
    assert os.path.exists(cache)
    # save=False:不落盘(某路失败时不缓存残缺结果)
    cache2 = str(tmp_path / "merged_y.csv")
    tg._merge_source_pair("metformin", a, b, cache_path=cache2, save=False, **tg._DRUG_MERGE)
    assert not os.path.exists(cache2)
    # 两源皆空 → 空表
    assert tg._merge_source_pair("x", pd.DataFrame(), pd.DataFrame(),
                                 cache_path=str(tmp_path / "z.csv"),
                                 **tg._DRUG_MERGE).empty


def test_drug_merge_v6_activity_details_roundtrip(tmp_path):
    """v6 药物 merged 保留逐条 ChEMBL 明细:位置对齐、确定性排序、CSV 往返后
    仍与原始记录等价,供可达性 B1 直接消费。"""
    a = pd.DataFrame([
        {"approvedSymbol": "tp53", "assay_type": "F", "pchembl_value": 7.0},
        {"approvedSymbol": "TP53", "assay_type": "B", "pchembl_value": 5.0},
        {"approvedSymbol": "TP53", "assay_type": "F", "pchembl_value": "6.0"},
        {"approvedSymbol": "TP53", "assay_type": "F", "pchembl_value": None},
        {"approvedSymbol": "EGFR", "assay_type": None, "pchembl_value": 9.0},
    ])
    b = pd.DataFrame([{"approvedSymbol": "TP53", "interaction_score": 0.9}])
    cache = str(tmp_path / "merged_drug_v6_probe.csv")
    merged = tg._merge_source_pair(
        "probe", a, b, cache_path=cache, **tg._DRUG_MERGE)
    tp = merged.set_index("approvedSymbol").loc["TP53"]
    assert tp["chmbl_assay_types"] == "B;F;F"
    assert tp["chmbl_pchembl_values"] == "5.0;6.0;7.0"

    # 输入顺序变化不得改变序列化结果(同缓存同输入可复现)
    shuffled = tg._merge_source_pair(
        "probe", a.sample(frac=1, random_state=7), b,
        cache_path=str(tmp_path / "shuffled.csv"), save=False, **tg._DRUG_MERGE)
    assert shuffled.set_index("approvedSymbol").loc["TP53", "chmbl_pchembl_values"] \
        == tp["chmbl_pchembl_values"]

    loaded = pd.read_csv(cache)
    from WY.core import reachability as rch

    b1 = rch._b1_activity(loaded).set_index("gene")
    assert b1.loc["TP53", "n_assays"] == 2
    assert b1.loc["TP53", "activity_M"] == pytest.approx(10 ** -6.5)
    assert b1.loc["TP53", "activity_alt_M"] == pytest.approx(10 ** -5.0)


def test_drug_merge_v6_rejects_missing_activity_columns(tmp_path):
    """v6 不兼容旧格式:A 源非空但缺 assay_type/pchembl_value 时拒绝写缓存。"""
    a = pd.DataFrame([{"approvedSymbol": "TP53", "pchembl_value": 7.0}])
    b = pd.DataFrame()
    with pytest.raises(RuntimeError, match="assay_type"):
        tg._merge_source_pair(
            "probe", a, b, cache_path=str(tmp_path / "bad.csv"),
            **tg._DRUG_MERGE)
    assert not (tmp_path / "bad.csv").exists()


def test_disease_merge_does_not_carry_drug_activity_columns(tmp_path):
    """明细列仅属药物侧;疾病 merged 列契约保持不变。"""
    ot = pd.DataFrame([{
        "approvedSymbol": "TP53", "association_score": 0.8,
        "biotype": "protein_coding",
    }])
    dg = pd.DataFrame([{"approvedSymbol": "TP53", "disgenet_score": 0.4}])
    merged = tg._merge_source_pair(
        "probe disease", ot, dg,
        cache_path=str(tmp_path / "merged_disease.csv"), **tg._DISEASE_MERGE)
    assert not any(c.startswith("chmbl_") for c in merged.columns)


def test_get_drug_targets_merged_writes_current_v6_and_hits_cache(monkeypatch):
    """完整入口:命中当前版本文件名、写出明细列;第二次调用只读 v6 缓存。"""
    monkeypatch.setattr(config, "CHEMBL_PCHEMBL_THRESHOLD", 6.0)
    monkeypatch.setattr(config, "CHEMBL_TARGET_BLACKLIST", ())
    monkeypatch.setattr(config, "CHEMBL_CORROBORATION_REQUIRED", True)
    monkeypatch.setattr(config, "DGIDB_MIN_SCORE", 0.3)
    calls = {"chembl": 0, "dgidb": 0}

    def fake_chembl(name):
        calls["chembl"] += 1
        return pd.DataFrame([
            {"approvedSymbol": "TP53", "assay_type": "B", "pchembl_value": 5.0},
            {"approvedSymbol": "TP53", "assay_type": "F", "pchembl_value": 6.0},
            {"approvedSymbol": "TP53", "assay_type": "F", "pchembl_value": 7.0},
        ]), "stub"

    def fake_dgidb(name):
        calls["dgidb"] += 1
        return pd.DataFrame([{
            "approvedSymbol": "TP53", "interaction_score": 0.9,
            "interaction_pmids": "123",
        }])

    monkeypatch.setattr(tg, "fetch_drug_targets_chembl", fake_chembl)
    monkeypatch.setattr(tg, "fetch_drug_targets_dgidb", fake_dgidb)
    first = tg.get_drug_targets_merged("v6-probe")
    cache = os.path.join(
        config.MERGED_CACHE_DIR,
        f"{config.MERGED_DRUG_PREFIX}{safe_name('v6-probe')}.csv")
    assert config.TARGET_FILTER_VERSION >= 6
    assert os.path.basename(cache).startswith(
        f"merged_drug_v{config.TARGET_FILTER_VERSION}_")
    assert os.path.exists(cache)
    assert first.iloc[0]["chmbl_assay_types"] == "B;F;F"
    assert first.iloc[0]["chmbl_pchembl_values"] == "5.0;6.0;7.0"
    assert calls == {"chembl": 1, "dgidb": 1}

    monkeypatch.setattr(tg, "fetch_drug_targets_chembl", lambda name: (
        (_ for _ in ()).throw(AssertionError("v6 cache miss"))))
    monkeypatch.setattr(tg, "fetch_drug_targets_dgidb", lambda name: (
        (_ for _ in ()).throw(AssertionError("v6 cache miss"))))
    second = tg.get_drug_targets_merged("v6-probe")
    assert list(second.columns) == list(first.columns)
    assert second.iloc[0]["chmbl_pchembl_values"] == "5.0;6.0;7.0"


def test_intersect_target_symbols():
    """交集规则:字母序排序、截断 top_n;任一侧为空/None → 空列表。"""
    drug = pd.DataFrame({"approvedSymbol": ["KRAS", "TP53", "EGFR"]})
    disease = pd.DataFrame({"approvedSymbol": ["BRAF", "EGFR", "TP53"]})
    assert tg.intersect_target_symbols(drug, disease, 10) == ["EGFR", "TP53"]
    assert tg.intersect_target_symbols(drug, disease, 1) == ["EGFR"]
    assert tg.intersect_target_symbols(drug, pd.DataFrame(), 10) == []
    assert tg.intersect_target_symbols(None, disease, 10) == []


# ---------- 并发编排 ----------


def _preload_merged_caches(drug: str, disease: str):
    from WY.core.caching import save_csv_cache
    save_csv_cache(pd.DataFrame({"approvedSymbol": ["TP53"]}),
                   os.path.join(config.MERGED_CACHE_DIR,
                                f"{config.MERGED_DRUG_PREFIX}{safe_name(drug)}.csv"),
                   skip_empty=False)
    save_csv_cache(pd.DataFrame({"approvedSymbol": ["TP53"]}),
                   os.path.join(config.MERGED_CACHE_DIR,
                                f"{config.MERGED_DISEASE_PREFIX}{safe_name(disease)}.csv"),
                   skip_empty=False)


def test_concurrent_both_cached_no_network(monkeypatch):
    """双侧 merged 缓存都命中:直接返回,四路数据源一个都不该被调用。"""
    _preload_merged_caches("drugX", "diseaseX")

    def must_not_call(*a, **k):
        raise AssertionError("network fetcher called despite merged cache hit")

    monkeypatch.setattr(tg, "fetch_drug_targets_chembl", must_not_call)
    monkeypatch.setattr(tg, "fetch_drug_targets_dgidb", must_not_call)
    monkeypatch.setattr(tg, "get_disease_targets", must_not_call)
    monkeypatch.setattr(tg, "fetch_disease_targets_disgenet", must_not_call)
    drug_df, disease_df, failed = tg.get_targets_concurrent("drugX", "diseaseX")
    assert failed == []
    assert not drug_df.empty and not disease_df.empty


def test_concurrent_failures_semantics(monkeypatch):
    """四路全挂:失败源列表如实上报(ChEMBL 中止=None 计失败、空密钥 DisGeNET 计失败),
    两侧仍返回空表而非崩溃,且残缺结果不写 merged 缓存。"""
    monkeypatch.setattr(config, "DisGeNET_API_KEY", "")
    monkeypatch.setattr(tg, "fetch_drug_targets_chembl",
                        lambda name: (pd.DataFrame(), None))       # 请求中止
    def boom(name):
        raise ConnectionError("simulated outage")
    monkeypatch.setattr(tg, "fetch_drug_targets_dgidb", boom)
    monkeypatch.setattr(tg, "get_disease_targets", boom)
    monkeypatch.setattr(tg, "fetch_disease_targets_disgenet",
                        lambda name, need_rows=None: pd.DataFrame())  # 空密钥空表
    ticks = []
    drug_df, disease_df, failed = tg.get_targets_concurrent(
        "drugY", "diseaseY", progress_cb=lambda done, total: ticks.append((done, total)))
    assert set(failed) == {"chembl", "dgidb", "opentargets", "disgenet"}
    assert isinstance(drug_df, pd.DataFrame) and isinstance(disease_df, pd.DataFrame)
    assert drug_df.empty and disease_df.empty
    # 残缺结果不落盘
    assert not os.path.exists(os.path.join(
        config.MERGED_CACHE_DIR,
        f"{config.MERGED_DRUG_PREFIX}{safe_name('drugY')}.csv"))
    # 进度回调至少被调过一次,且最终 done 数达到 total
    assert ticks and len(ticks[-1][0]) == ticks[-1][1] == 4  # done 为已完成的源 key 列表


def test_invalidate_merged_removes_only_target(tmp_path):
    """失效器:只删除对应 safe 名的 merged 缓存;不存在时静默无害。"""
    os.makedirs(config.MERGED_CACHE_DIR, exist_ok=True)
    keep = os.path.join(config.MERGED_CACHE_DIR,
                        f"{config.MERGED_DRUG_PREFIX}{safe_name('other')}.csv")
    open(keep, "w").close()
    victim = os.path.join(config.MERGED_CACHE_DIR,
                          f"{config.MERGED_DRUG_PREFIX}{safe_name('metformin')}.csv")
    open(victim, "w").close()
    tg.invalidate_merged_drug_targets("metformin")
    assert not os.path.exists(victim) and os.path.exists(keep)
    tg.invalidate_merged_drug_targets("never_existed")  # 无害
