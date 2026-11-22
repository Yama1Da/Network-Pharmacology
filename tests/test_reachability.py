"""core/reachability 单元测试(清单 §1.5):纯函数手算对照 + pkdb_raw 桩数据
沙箱端到端。零外网:conftest autouse 沙箱重定向 CACHE_DIR,桩数据按真实
API/TSV 结构缩样(分页信封双层 data 由镜像层剥好,这里落的是剥好后的形态)。

§5.6 三不变案例中本模块覆盖前两例(第三例置换零模型属 synergy 数学核心):
  二甲双胍: µM 暴露 × 经典靶点 mM 活性 → 不可达(单位换算与 pChEMBL 方向;
            若单位千倍错则 R>1 变可达,方向反则高 p 案例翻转——两案例互锁)
  雷米普利: nM 活性 vs 治疗 Cmax → 可达(端到端管线)
"""

import json
import os

import pandas as pd
import pytest

from WY import config
from WY.core import reachability as rch


# ---- 桩数据(真实结构缩样,§1.5-1) ----

def _substance_nodes() -> list[dict]:
    """物质映射节点:metformin/ramipril/caffeine/theophylline + 复合物质,
    字段形态对 2026-11-22 实抓节点核对(annotations/xrefs/parents/substance.mass)。"""
    def node(sid, mass, synonyms=(), chembl=None, chebi=None, parents=()):
        return {
            "sid": sid, "name": sid, "label": sid, "deprecated": False,
            "ntype": "substance", "dtype": "undefined",
            "synonyms": list(synonyms),
            "parents": [{"sid": p, "name": p, "label": p} for p in parents],
            "annotations": ([{"term": chebi, "collection": "chebi"}] if chebi else []),
            "xrefs": ([{"name": "chembl", "accession": chembl}] if chembl else []),
            "measurement_type": None,
            "substance": {"mass": mass, "charge": 0.0, "formula": None},
        }
    return [
        node("metformin", 129.167, synonyms=("glucophage",),
             chembl="CHEMBL901", chebi="CHEBI:6801"),
        node("ramipril", 416.5106, chembl="CHEMBL1135", chebi="CHEBI:8774"),
        node("caffeine", 194.19076, chembl="CHEMBL113", chebi="CHEBI:27732"),
        node("theophylline", 180.164, chembl="CHEMBL1234"),
        node("tp-plus-tb-plus-px-plus-caf", None,
             parents=("caf", "px", "tb", "tp")),
    ]


def _study(sid: str, tables: list[tuple[str, list[str], list[list]]]) -> dict:
    """研究 JSON(pk_mirror._fetch_study_payload 落盘形态):tsv 为逐字二维表。"""
    return {
        "sid": sid, "name": sid, "date": "2020-01-01",
        "licence": "CC0", "access": "open",
        "substances": [], "reference_pmid": "", "reference_title": "",
        "files": [],
        "tsv": [{"file": name, "url": f"https://pk-db.example/{name}",
                 "rows": [header] + rows} for name, header, rows in tables],
        "skipped": [], "fetched_at": "2026-11-22T00:00:00+00:00",
    }


_CMAX_HDR = ["study", "group", "intervention", "measurement_type", "substance",
             "tissue", "time", "time_unit", "mean", "se", "unit"]
_DOSE_HDR = ["study", "name", "substance", "measurement_type", "time",
             "time_unit", "route", "form", "value", "unit", "application"]


@pytest.fixture()
def sandbox_norm(monkeypatch):
    """在 conftest 自动沙箱 CACHE_DIR 下伪造最小 pkdb_raw(真实结构缩样):
    substance_map.json + 3 个研究(cmax/给药/fu 行),并挂 PK_NORM_VERSION=7
    ——走真实 _versioned_path 机制,不绕过版本化。返回 pkdb_raw 目录。"""
    monkeypatch.setattr(config, "PK_NORM_VERSION", 7, raising=False)
    raw = os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR)
    os.makedirs(raw, exist_ok=True)
    with open(os.path.join(raw, "substance_map.json"), "w", encoding="utf-8") as f:
        json.dump({"fetched_at": "2026-11-22T00:00:00+00:00",
                   "nodes": _substance_nodes()}, f, ensure_ascii=False)
    studies = [
        # metformin:1000 ng/ml → 1000×1e-6 g/L ÷129.167 = 7.7418e-6 M(≈7.74 µM)
        _study("PKDBTEST1", [
            (".Met_Tab1.tsv", _CMAX_HDR,
             [["PKDBTEST1", "control", "MET1000", "cmax", "metformin",
               "plasma", "NA", "NA", "1000", "na", "ng/ml"]]),
            (".Met_TabInterventions.tsv", _DOSE_HDR,
             [["PKDBTEST1", "MET1000", "metformin", "dosing", "0",
               "min", "po", "tablet", "1000", "mg", "single dose"]]),
        ]),
        # ramipril:52 ng/ml → 52e-6 g/L ÷416.5106 = 1.2485e-7 M(≈0.125 µM)
        _study("PKDBTEST2", [
            (".Ril_Tab1.tsv", _CMAX_HDR,
             [["PKDBTEST2", "all", "RIL10", "cmax", "ramipril",
               "plasma", "NA", "NA", "52", "na", "ng/ml"]]),
            (".Ril_TabInterventions.tsv", _DOSE_HDR,
             [["PKDBTEST2", "RIL10", "ramipril", "dosing", "0",
               "min", "po", "tablet", "10", "mg", "single dose"]]),
        ]),
        # fu:caffeine 直接小数;theophylline 百分数结合率 57.7 → fub=0.423;
        # 放射性剂量干预(示踪排除)+ 缺分子量物质(无法换算)各一
        _study("PKDBTEST3", [
            (".Caf_Tab1.tsv",
             ["study", "group", "measurement_type", "substance", "tissue",
              "unit", "mean"],
             [["PKDBTEST3", "control", "fraction_unbound", "caffeine",
               "plasma", "", "0.85"],
              ["PKDBTEST3", "g1", "plasma_binding", "theophylline",
               "serum", "percent", "57.7"]]),
            (".Caf_Tab2.tsv", _CMAX_HDR,
             [["PKDBTEST3", "control", "TR_C14", "cmax", "caffeine",
               "plasma", "NA", "NA", "10", "na", "ng/ml"],
              ["PKDBTEST3", "control", "UNK", "cmax", "mystery-substance",
               "plasma", "NA", "NA", "5", "na", "mg/l"]]),
            (".Caf_TabInterventions.tsv", _DOSE_HDR,
             [["PKDBTEST3", "TR_C14", "caffeine", "dosing", "0", "min",
               "iv", "solution", "1", "MBq", "single dose"]]),
        ]),
    ]
    for s in studies:
        with open(os.path.join(raw, f"{s['sid']}.json"), "w", encoding="utf-8") as f:
            json.dump(s, f, ensure_ascii=False)
    return raw


def _build_all():
    """沙箱内全量归一化(真实 _versioned_path:文件名带 _v7)。返回路径 dict。"""
    return rch.build_normalized_cache(force=True)


# ---- 契约(第 3 步) ----

class TestContract:
    def test_schema_and_enums(self):
        assert rch.REACHABILITY_SCHEMA == "reachability.v1"
        assert rch.REACH_TIERS == ("reachable", "partial", "unreachable", "unknown")
        assert rch.REACH_UNKNOWN_REASONS == ("no_exposure", "no_activity",
                                             "no_molweight")

    def test_exposure_confident(self):
        assert rch.exposure_confident("reachable", "partial") is True
        assert rch.exposure_confident("reachable", "unknown") is False
        assert rch.exposure_confident("unreachable", "unreachable") is False
        with pytest.raises(ValueError):
            rch.exposure_confident("bogus", "reachable")

    def test_versioned_path_gates_and_names(self, monkeypatch, sandbox_norm):
        path = rch._versioned_path("cmax_norm")
        assert path.endswith("cmax_norm_v7.parquet")
        assert os.path.join(config.CACHE_DIR, "reachability") in path
        monkeypatch.delattr(config, "PK_NORM_VERSION", raising=False)
        with pytest.raises(RuntimeError):
            rch._versioned_path("cmax_norm")  # 缺版本常量=响亮拒绝,非静默默认


# ---- substance_map 解析(第 4 步) ----

class TestSubstanceMap:
    def test_rows_parsing(self):
        df = rch._substance_rows(_substance_nodes()).set_index("sid")
        met = df.loc["metformin"]
        assert met["mol_weight_g_mol"] == pytest.approx(129.167)
        assert "CHEMBL901" in met["chembl_id"] and "CHEBI:6801" in met["chebi_ids"]
        assert "glucophage" in met["synonyms"]
        comp = df.loc["tp-plus-tb-plus-px-plus-caf"]
        assert comp["parent_sids"] == "caf;px;tb;tp"
        assert pd.isna(comp["mol_weight_g_mol"])  # 缺失=NaN(§0.4,下游 pd.isna)

    def test_rows_dedup_and_bad_mass(self):
        nodes = _substance_nodes() + [_substance_nodes()[0]]  # 重复 sid
        nodes[1]["substance"]["mass"] = "not-a-number"        # 非数值质量
        df = rch._substance_rows(nodes).set_index("sid")
        assert len(df) == 5  # 重复保留首条
        assert pd.isna(df.loc["ramipril", "mol_weight_g_mol"])

    def test_lookup_synonyms_and_mw(self):
        lookup = rch._substance_lookup(_substance_nodes())
        assert lookup["glucophage"][0] == "metformin"  # 同义词可查
        assert lookup["metformin"][1] == pytest.approx(129.167)


# ---- 单位规范与换算(第 5 步;§1.2-1 千倍错防线) ----

class TestUnits:
    def test_canon_unit(self):
        assert rch._canon_unit("NG / mL") == "ng/ml"
        assert rch._canon_unit("µmole/L") == "µmol/l"
        assert rch._canon_unit("umol/l") == "µmol/l"
        assert rch._canon_unit("microgr/ml") == "µg/ml"

    def test_parse_num(self):
        assert rch._parse_num("12.7") == 12.7
        assert rch._parse_num("NA") is None and rch._parse_num("") is None

    def test_hand_calc_factors(self):
        # 1000 ng/ml ÷ 129.167 = 7.7418 µM(手算对照,§5.6-二甲双胍)
        m = 1000 * rch._MASS_FACTOR["ng/ml"] / 129.167
        assert m == pytest.approx(7.7418e-6, rel=1e-4)
        assert 1000 * rch._MOLAR_FACTOR["nmol/l"] == pytest.approx(1e-6)
        assert 10 * rch._MOLAR_FACTOR["µm"] == pytest.approx(1e-5)


# ---- 归一化三件套端到端(第 5-6 步,桩数据沙箱) ----

class TestNormalizedCache:
    def test_three_parquets_written(self, sandbox_norm):
        paths = _build_all()
        for stem in ("substance_map", "cmax_norm", "fub_norm"):
            assert os.path.isfile(paths[stem]) and "_v7." in paths[stem]
        cmax = pd.read_parquet(paths["cmax_norm"])
        met = cmax[cmax["substance_raw"] == "metformin"].iloc[0]
        assert met["cmax_molar_M"] == pytest.approx(7.7418e-6, rel=1e-4)
        assert met["route"] == "po" and met["dose_value"] == 1000.0  # 给药关联

    def test_tracer_and_no_molweight(self, sandbox_norm):
        _build_all()
        agg, stats = rch._cmax_observations()
        assert stats["tracer_excluded"] >= 1      # MBq 剂量干预被结构性排除
        assert stats["no_molweight"] >= 1         # 无映射物质的 mg/l 行:无法换算,不硬算
        assert agg["cmax_molar_M"].notna().all()  # 聚合表只含可换算行

    def test_fub_norm_bounds_and_flags(self, sandbox_norm):
        paths = _build_all()
        fub = pd.read_parquet(paths["fub_norm"]).set_index("substance_sid")
        assert ((fub["fub_value"] > 0) & (fub["fub_value"] <= 1)).all()
        assert not fub["fub_value"].isna().any()  # 绝无 0/NaN(除零防线 §1.2-4)
        assert fub.loc["caffeine", "source"] == "pkdb"
        # 沙箱仅 1 条 fu 记录(真实镜像两条 0.85/0.86 → 中位数 0.855);桩=单行
        assert fub.loc["caffeine", "fub_value"] == pytest.approx(0.85)
        theo = fub.loc["theophylline"]
        assert theo["fub_value"] == pytest.approx(1 - 57.7 / 100)  # 百分数结合率换算
        assert theo["source"] == "pkdb" and theo["method"] == "plasma_binding"
        met = fub.loc["metformin"]
        assert met["source"] == "imputed" and bool(met["imputed_flag"])  # 可能高估标记
        assert met["fub_value"] == 1.0  # 保守上界

    def test_version_missing_rejects_build(self, sandbox_norm, monkeypatch):
        monkeypatch.delattr(config, "PK_NORM_VERSION", raising=False)
        with pytest.raises(RuntimeError):
            rch.build_normalized_cache(force=True)


# ---- B1 聚合与 R 值(第 7 步) ----

class TestB1AndR:
    def test_b1_cell_priority_and_fallback(self):
        df = pd.DataFrame({
            "approvedSymbol": ["GENE_F"] * 3 + ["GENE_B"] * 2 + ["GENE_N"]
                              + ["GENE_MIX"] * 2,
            "assay_type": ["F", "F", "B", "B", "B", "F", "F", "B"],
            "pchembl_value": [6, 7, 5, 8, 9, None, 6, 7],
        })
        b1 = rch._b1_activity(df).set_index("gene")
        assert b1.loc["GENE_F", "activity_M"] == pytest.approx(10 ** -6.5)
        assert not b1.loc["GENE_F", "fallback_flag"]
        assert b1.loc["GENE_F", "activity_alt_M"] == pytest.approx(1e-5)  # 生化对照列
        assert b1.loc["GENE_B", "activity_M"] == pytest.approx(10 ** -8.5)
        assert bool(b1.loc["GENE_B", "fallback_flag"])                   # 回退留标记
        assert b1.loc["GENE_MIX", "activity_M"] == pytest.approx(1e-6)   # 细胞优先
        assert "GENE_N" not in b1.index                                  # 无值缺席

    def test_r_value_and_guards(self):
        assert rch.r_value(7.64e-6, 1.0, 1e-3) == pytest.approx(7.64e-3)
        assert rch.r_value(1.25e-7, 1.0, 5e-9) == pytest.approx(25.0)
        assert rch.r_value(None, 1.0, 1e-3) is None   # 缺失显式 None,不硬算
        assert rch.r_value(1e-5, 1.0, 0.0) is None    # 活性≤0 防 inf

    def test_tier_boundaries_and_sensitivity(self):
        assert rch._tier(1.0) == "reachable"
        assert rch._tier(0.999999) == "partial"
        assert rch._tier(0.1) == "partial"
        assert rch._tier(0.099999) == "unreachable"
        assert rch._tier(2.0, reachable_min=3.0) == "partial"  # §3.2 敏感性传参


# ---- §5.6 两不变案例(端到端,互锁防单位/方向错) ----

class TestS56Cases:
    def test_case1_metformin_unreachable(self, sandbox_norm):
        """二甲双胍:µM 暴露 × 经典靶点 mM 活性 → 不可达。管线级:若单位
        千倍错(µM→mM)则 R>1 翻成可达,本案例即红(§5.6 验证什么:单位换算)。"""
        _build_all()
        targets = pd.DataFrame({"approvedSymbol": ["SLC22A1"],
                                "assay_type": ["F"], "pchembl_value": [3.0]})  # 1 mM
        res = rch.compute_reachability(targets, "metformin")
        assert res["r_primary"] < 0.1 and res["tier"] == "unreachable"
        assert res["exposure"]["cmax_molar_M"] == pytest.approx(7.7418e-6, rel=1e-3)

    def test_case2_ramipril_reachable(self, sandbox_norm):
        """雷米普利:nM 活性 vs 治疗 Cmax → 可达。pChEMBL 方向哨兵:若 10**(-p)
        写反,高 p(强效)案例活性剧增,R 翻转为不可达,本案例即红(§5.6)。"""
        _build_all()
        targets = pd.DataFrame({"approvedSymbol": ["ACE"],
                                "assay_type": ["F"], "pchembl_value": [8.3]})  # ≈5 nM
        res = rch.compute_reachability(targets, "ramipril")
        assert res["r_primary"] > 1 and res["tier"] == "reachable"
        assert res["exposure"]["cmax_molar_M"] == pytest.approx(1.2485e-7, rel=1e-3)

    def test_unknown_branches(self, sandbox_norm):
        """未知档显式:药名无物质映射→no_exposure;有暴露无活性→no_activity;
        活性缺失靶点 tier=unknown,永不并入不可达(§1.2-5)。"""
        _build_all()
        res = rch.compute_reachability(pd.DataFrame(
            {"approvedSymbol": ["X"], "assay_type": ["F"],
             "pchembl_value": [6.0]}), "不存在药名xyz")
        assert res["tier"] == "unknown" and res["unknown_reason"] == "no_exposure"
        res2 = rch.compute_reachability(pd.DataFrame(
            {"approvedSymbol": ["X"], "assay_type": ["F"],
             "pchembl_value": [None]}), "metformin")
        assert res2["tier"] == "unknown" and res2["unknown_reason"] == "no_activity"
        assert (res2["targets"]["tier"] == "unknown").all()

    def test_drug_tier_aggregation(self, sandbox_norm):
        """药级档=任一靶点达档取最高权档(合同冻结口径);r_primary=该档内最大 R。"""
        _build_all()
        targets = pd.DataFrame({
            "approvedSymbol": ["WEAK", "STRONG"],
            "assay_type": ["F", "F"], "pchembl_value": [3.0, 8.3],
        })
        res = rch.compute_reachability(targets, "ramipril")
        assert res["tier"] == "reachable" and res["r_primary"] > 1
        tiers = res["targets"].set_index("gene")["tier"]
        assert tiers["WEAK"] == "unreachable" and tiers["STRONG"] == "reachable"


# ---- 排序层(第 8 步;E1 键方向=强制审查点,逐字按 §12.1-E1) ----

class TestSorting:
    def test_weights(self):
        assert rch.reach_weights() == {"reachable": 1.0, "partial": 0.6,
                                       "unknown": 0.5, "unreachable": 0.2}
        assert rch.reach_weight("unknown") > rch.reach_weight("unreachable")
        with pytest.raises(ValueError):
            rch.reach_weight("bogus")
        with pytest.raises(RuntimeError):
            rch.reach_weight("reachable", {"reachable": 1.0})  # 残缺权重表拒绝

    def test_e1_sort_direction(self):
        df = pd.DataFrame({
            "drug_a": ["P1", "P2", "P3", "P4", "P5"],
            "drug_b": ["Q1", "Q2", "Q3", "Q4", "Q5"],
            "s_ab": [-1.0, -2.0, 1.5, 0.8, 0.0],
            "z_ab": [-5.0, -3.0, 0.5, 1.0, 2.0],
        })
        e1 = rch.apply_e1_sort(df)
        # 重叠型按 z_AB 降序(-3 > -5);互补型按分离强度降序;tie(s=0)归互补
        assert list(e1["drug_a"]) == ["P2", "P1", "P3", "P4", "P5"]
        assert list(e1["geometry_label"]) == ["overlap", "overlap",
                                              "separated", "separated", "separated"]
        assert list(rch.apply_e1_sort(df, top_k=2)["drug_a"]) == ["P2", "P1", "P3", "P4"]

    def test_e6_sort_and_confident(self):
        df = pd.DataFrame({
            "drug_a": ["R1", "R2", "R3"], "drug_b": ["S1", "S2", "S3"],
            "score": [1.0, 1.0, 2.0],
            "tier_a": ["reachable", "reachable", "partial"],
            "tier_b": ["reachable", "unknown", "partial"],
        })
        e6 = rch.apply_e6_sort(df, score_col="score",
                               tier_a_col="tier_a", tier_b_col="tier_b")
        assert list(e6["drug_a"]) == ["R3", "R1", "R2"]   # 1.2 > 1.0 > 0.5
        assert list(e6["exposure_confident"]) == [True, True, False]
        with pytest.raises(ValueError):
            bad = df.copy()
            bad.loc[0, "tier_b"] = "bogus"
            rch.apply_e6_sort(bad, score_col="score",
                              tier_a_col="tier_a", tier_b_col="tier_b")

    def test_exposure_compare(self):
        targets = pd.DataFrame([
            {"gene": "G1", "r_value": 2.0, "tier": "reachable"},
            {"gene": "G2", "r_value": None, "tier": "unknown"},
            {"gene": "G3", "r_value": 0.05, "tier": "unreachable"},
            {"gene": "G4", "r_value": 0.5, "tier": "partial"},
        ])
        cmp = rch.exposure_compare(["G1", "G2", "G3", "G4", "G5"],
                                   {"targets": targets}).set_index("gene")
        assert list(cmp.loc[["G1", "G4", "G2", "G5", "G3"], "new_rank"]) == [1, 2, 3, 4, 5]
        assert cmp.loc["G4", "delta"] == 2 and cmp.loc["G3", "delta"] == -2
        assert cmp.loc["G1", "delta"] == 0
        assert cmp.loc["G5", "tier"] == "unknown"  # baseline 外基因入未知档,不静默丢

    def test_tracer_detection(self):
        assert rch._is_tracer("caffeine", "TR1", "", "mbq")       # 放射性剂量单位
        assert rch._is_tracer("C14-labeled", "TR1", "", "mg")     # 物质名标记
        assert rch._is_tracer("x", "TR1", "scintillation counting", "mg")
        assert not rch._is_tracer("caffeine", "MET1000", "", "mg")  # 治疗剂量不误杀
