"""benchmark/drugcomb_benchmark 单元测试(清单 §1.6)。

手算对照:
  - 标签聚合(陷阱 B6):多数投票/平票剔除/<3 细胞系剔除/模糊带剔除/投票比例列
  - 防泄漏闸(第一审核项):药对跨切分两侧的违规输入必须响
  - 谱系映射(B5):CCLE 成功匹配 + 失败挂泛癌不剔除 + unmatched 报告
  - 分层报告:AUROC(sklearn)+ Top-K 回收手算
"""

import json

import pandas as pd
import pytest

from WY.benchmark import drugcomb_benchmark as dcb


@pytest.fixture()
def raw():
    """合成 DrugComb 导出(常见门户列名;块级多记录 + 模糊带 + 无效行)。"""
    rows = [
        # (aspirin, warfarin):4 细胞系,3 协同 1 拮抗 → 多数协同,frac 0.75
        ("Aspirin", "Warfarin", "MCF-7", 20.0),
        ("Aspirin", "Warfarin", "Hela", 15.0),
        ("Aspirin", "Warfarin", "A549", 12.0),
        ("Aspirin", "Warfarin", "A549", 18.0),  # 同格多 block → 中位数 15
        ("Aspirin", "Warfarin", "HEPG2", -12.0),
        # (aspirin, metformin):仅 2 细胞系 → <3 剔除
        ("Aspirin", "Metformin", "MCF-7", 20.0),
        ("Aspirin", "Metformin", "Hela", 18.0),
        # (aspirin, ibuprofen):4 细胞系 2v2 平票 → 剔除
        ("Aspirin", "Ibuprofen", "MCF-7", 12.0),
        ("Aspirin", "Ibuprofen", "Hela", 15.0),
        ("Aspirin", "Ibuprofen", "A549", -12.0),
        ("Aspirin", "Ibuprofen", "HEPG2", -15.0),
        # (aspirin, caffeine):全模糊带(|zip|<10) → 票前剔除
        ("Aspirin", "Caffeine", "MCF-7", 5.0),
        ("Aspirin", "Caffeine", "Hela", -5.0),
        # 无效 ZIP(NA)行 → 计数
        ("Aspirin", "Warfarin", "A549", None),
    ]
    return pd.DataFrame(rows, columns=["drug_row", "drug_col",
                                       "cell_line_name", "synergy_zip"])


class TestAggregate:
    def test_hand_calc(self, raw):
        out = dcb.aggregate_pair_labels(raw)
        labels = out["labels"].set_index(["drug_row", "drug_col"])
        stats = out["stats"]
        assert stats["rows_invalid_zip"] == 1
        assert stats["rows_fuzzy_band"] == 2
        assert stats["pairs_dropped_small"] == 1
        assert stats["pairs_dropped_tie"] == 1
        assert stats["pairs_kept"] == 1
        w = labels.loc[("aspirin", "warfarin")]
        assert w["vote"] == 1 and w["n_cell_lines"] == 4
        # 同格双 block(20/18)取中位数 15 → 4 票中 3 协同
        assert w["vote_fraction"] == pytest.approx(0.75)
        assert w["n_syn"] == 3

    def test_column_map_override(self, raw):
        renamed = raw.rename(columns={"synergy_zip": "zip_score",
                                      "cell_line_name": "cell"})
        out = dcb.aggregate_pair_labels(
            renamed, col_map={"zip": "zip_score", "cell_line": "cell"})
        assert out["stats"]["pairs_kept"] == 1

    def test_missing_column_loud(self, raw):
        with pytest.raises(RuntimeError):
            dcb.aggregate_pair_labels(raw.drop(columns=["synergy_zip"]))


class TestNormalize:
    def test_drug_and_cell_line(self):
        assert dcb.normalize_drug_name("Warfarin- Sodium") == "warfarinsodium"
        assert dcb.normalize_cell_line("MCF-7") == "MCF7"
        assert dcb.normalize_cell_line(" mcf7 ") == "MCF7"


class TestSplit:
    def _labels(self, n_pairs=10):
        rows = []
        for i in range(n_pairs):
            for j in range(4):
                rows.append({"drug_row": f"d{i}", "drug_col": "shared",
                             "cell_line": f"CL{j}", "vote": i % 2,
                             "vote_label": i % 2, "n_votes": 4,
                             "n_syn": 2, "vote_fraction": 0.5})
        return pd.DataFrame(rows)

    def test_pair_level_split_and_leak_gate(self):
        labels = self._labels()
        split = dcb.split_pairs(labels, test_frac=0.3, seed=0)
        assert len(split["pairs"]["train"]) == 7
        assert len(split["pairs"]["test"]) == 3
        # 药对粒度:同药对的全部记录只在一侧(陷阱#1,第一审核项)
        d = split["test"].merge(split["train"], on=["drug_row", "drug_col"],
                                how="inner")
        assert d.empty
        with pytest.raises(AssertionError):
            dcb.assert_pair_split_clean({"p1"}, {"p1", "p2"})  # 违规输入闸会响
        dcb.assert_pair_split_clean({"p1"}, {"p2"})            # 干净输入不抛

    def test_deterministic(self):
        a = dcb.split_pairs(self._labels(), seed=0)["pairs"]
        b = dcb.split_pairs(self._labels(), seed=0)["pairs"]
        assert a == b


class TestLineage:
    def test_map_and_pancancer_fallback(self):
        ccle = pd.DataFrame({"cell_line": ["MCF7", "HELA"],
                             "lineage": ["breast", "lung"]})
        res = dcb.map_lineage(["MCF-7", "Unknown-X"], ccle)
        assert res["lineage"]["MCF-7"] == "breast"        # 规范化后命中
        assert res["lineage"]["Unknown-X"] == "pancancer"  # B5:失败挂泛癌不剔除
        assert res["unmatched"] == ["Unknown-X"]           # 报告不静默丢


class TestReport:
    def _pred(self):
        return pd.DataFrame({
            "drug_row": ["d0", "d1", "d2", "d3"],
            "drug_col": ["shared"] * 4,
            "confidence": [10.0, 9.0, 8.0, 7.0],
            "vote_label": [1, 1, 0, 0],
            "geometry": [-1.0, 0.5, -0.2, 2.0],
            "lineage": ["breast", "lung", "pancancer", "pancancer"],
        })

    def test_auroc_and_topk(self):
        rep = dcb.stratified_report(self._pred(), score_col="confidence",
                                    geometry_col="geometry",
                                    lineage_col="lineage", top_k=2)
        overall = rep[rep["stratum"] == "overall"].iloc[0]
        assert overall["auroc"] == pytest.approx(1.0)   # 完美排序
        assert overall["topk_recall"] == pytest.approx(1.0)  # 前 2 全协同
        assert set(rep["stratum"]) == {"overall", "overlap", "separated",
                                       "lineage_mapped", "pancancer"}
        sep = rep[rep["stratum"] == "separated"].iloc[0]
        assert sep["auroc"] == pytest.approx(1.0)       # 互补侧 d3(0)>d1(1)? 否——
        # d1 协同分 9 与 d3 拮抗分 7:分层内 AUROC 仍应完美

    def test_single_class_stratum_reported_not_faked(self):
        df = self._pred()
        df.loc[df["geometry"] < 0, "vote_label"] = 1  # 重叠侧全正例
        rep = dcb.stratified_report(df, score_col="confidence",
                                    geometry_col="geometry")
        ov = rep[rep["stratum"] == "overlap"].iloc[0]
        assert pd.isna(ov["auroc"]) and ov["n_pairs"] == 2  # 单侧类如实缺报

    def test_topk_zero_positive_loud(self):
        df = self._pred()
        df["vote_label"] = 0
        with pytest.raises(RuntimeError):
            dcb.topk_recall(df, score_col="confidence", label_col="vote_label",
                            k=2)


class TestCli:
    def _pred(self):
        """CLI 预测表:只带分数与 geometry(真实管线形态,不含标签列)。"""
        return pd.DataFrame({
            "drug_row": ["d0", "d1", "d2", "d3"],
            "drug_col": ["shared"] * 4,
            "confidence": [10.0, 9.0, 8.0, 7.0],
        })

    def _write(self, tmp_path):
        preds = self._pred()
        raw_rows = []
        for _, r in preds.iterrows():
            for cl, z in (("MCF7", 20.0), ("HELA", 15.0), ("A549", 12.0)):
                raw_rows.append({"drug_row": r["drug_row"],
                                 "drug_col": r["drug_col"],
                                 "cell_line_name": cl,
                                 "synergy_zip": z})
        raw = pd.DataFrame(raw_rows)
        raw_path = tmp_path / "drugcomb.tsv"
        raw.to_csv(raw_path, sep="\t", index=False)
        pred_path = tmp_path / "pred.csv"
        preds.to_csv(pred_path, index=False)
        return str(raw_path), str(pred_path)

    def test_main_exit_codes(self, tmp_path):
        raw_path, pred_path = self._write(tmp_path)
        rc = dcb.main(["--labels", raw_path, "--predictions", pred_path,
                       "--out", str(tmp_path / "ok")])
        assert rc == 0
        summary = json.loads((tmp_path / "ok" / "summary.json")
                             .read_text(encoding="utf-8"))
        assert summary["schema"] == "drugcomb_benchmark.v1"
        # 预测与测试药对无交集 → 响亮失败,退出码 1
        preds_bad = pd.DataFrame({"drug_row": ["nope"], "drug_col": ["x"],
                                  "confidence": [1.0]})
        bad_path = tmp_path / "bad_pred.csv"
        preds_bad.to_csv(bad_path, index=False)
        rc = dcb.main(["--labels", raw_path, "--predictions", str(bad_path),
                       "--out", str(tmp_path / "bad")])
        assert rc == 1
