"""core.ppi 的 STRING 抓取与结果落盘测试(抓取打桩,缓存走沙箱)。

阶段三续:深水区补测第一批。
"""
import os

import pandas as pd
import pytest

from WY import config
from WY.core.text_utils import safe_name
import WY.core.ppi as ppi


class _FakeResp:
    def __init__(self, payload):
        self._payload = payload

    def json(self):
        return self._payload


def test_fetch_string_ppi_full_flow_and_cache(monkeypatch):
    """全链:基因→STRING ID 映射→网络边→列筛选改名;结果按缓存键落盘,二次调用命中缓存零请求。"""
    calls = {"n": 0}

    def fake_post(endpoint, params):
        calls["n"] += 1
        if endpoint == "get_string_ids":
            return _FakeResp([{"stringId": "9606.ENSP1"},
                              {"stringId": "9606.ENSP2"},
                              {"noId": True}])  # 无 stringId 的行被跳过
        return _FakeResp([
            {"preferredName_A": "TP53", "preferredName_B": "EGFR", "score": 900,
             "nscore": 0, "fscore": 0, "pscore": 0, "ascore": 0, "escore": 0,
             "dscore": 0, "tscore": 0, "extra_col": 1},
        ])

    monkeypatch.setattr(ppi, "_string_post", fake_post)
    df = ppi.fetch_string_ppi(["tp53", "EGFR"], score=0.7)
    assert calls["n"] == 2  # 两次请求:映射 + 网络
    assert list(df["protein_A"]) == ["TP53"] and list(df["protein_B"]) == ["EGFR"]
    assert df.iloc[0]["score"] == 900
    assert "extra_col" not in df.columns  # 只保留约定的证据分列

    df2 = ppi.fetch_string_ppi(["tp53", "EGFR"], score=0.7)
    assert calls["n"] == 2  # 二次命中缓存(同参数同基因列表指纹),零新请求
    pd.testing.assert_frame_equal(df, df2)


def test_fetch_string_ppi_unrecognized_genes(monkeypatch):
    """基因全部不被 STRING 识别:返回空表(不是报错)。"""
    monkeypatch.setattr(ppi, "_string_post",
                        lambda endpoint, params: _FakeResp([{"noId": True}]))
    assert ppi.fetch_string_ppi(["ZZZ1"]).empty


def test_fetch_string_ppi_none_response(monkeypatch):
    """请求重试耗尽(None):返回空表,不崩溃。"""
    monkeypatch.setattr(ppi, "_string_post", lambda endpoint, params: None)
    assert ppi.fetch_string_ppi(["TP53"]).empty


def test_base_disease_of_strips_rwr_suffix():
    """剥掉历史实验变体后缀(rwrNN,不区分大小写),干净疾病名用于下拉与 GEO 检索。"""
    assert ppi.base_disease_of("gastric cancer rwr72") == "gastric cancer"
    assert ppi.base_disease_of("T2D RWR3") == "T2D"
    assert ppi.base_disease_of("type 2 diabetes") == "type 2 diabetes"


def test_save_query_genes(tmp_path):
    """查询集清单:RWR 扩展表落盘为三列 CSV;未启用扩展(None/空表)不落盘。"""
    expansion = pd.DataFrame({"gene": ["TP53", "EGFR"],
                              "source": ["seed", "expanded"],
                              "rwr_score": [0.9, 0.4]})
    out = ppi.save_query_genes(str(tmp_path), "metformin_T2D", expansion)
    path = out["query_genes_csv"]
    assert os.path.exists(path)
    back = pd.read_csv(path)
    assert list(back["source"]) == ["seed", "expanded"]
    assert ppi.save_query_genes(str(tmp_path), "x", None) == {}
    assert ppi.save_query_genes(str(tmp_path), "x", pd.DataFrame()) == {}


def test_save_proximity_results(tmp_path):
    """邻近度落盘:metric/value 长表,只记非空指标;失败/未计算不落盘。"""
    result = {"ok": True, "d_c": 1.5, "z": -2.3, "p_empirical": 0.04,
              "d_c_excl": None,  # None 值不记
              "coverage": {"seed_in_lcc": 5, "seed_total": 5}}
    out = ppi.save_proximity_results(str(tmp_path), "metformin_T2D", result)
    path = out["proximity_csv"]
    back = pd.read_csv(path)
    metrics = dict(zip(back["metric"], back["value"]))
    assert metrics["z"] == -2.3 and metrics["coverage_seed_in_lcc"] == 5
    assert "d_c_excl" not in set(back["metric"])
    assert ppi.save_proximity_results(str(tmp_path), "x", None) == {}
    assert ppi.save_proximity_results(str(tmp_path), "x", {"ok": False}) == {}


def test_list_analyzed_diseases_and_drugs(monkeypatch):
    """疾病/药物下拉数据源:meta 精确 + 前缀反解 + rwr 变体剥离 + 大小写不敏感去重。"""
    dl = config.DOWNLOADS_DIR
    import json
    d1 = os.path.join(dl, "aspirin_cancer_data")           # meta 精确
    os.makedirs(d1, exist_ok=True)
    with open(os.path.join(d1, "analysis_meta.json"), "w", encoding="utf-8") as f:
        json.dump({"drug": "Aspirin", "disease": "cancer rwr50",
                   "string_version": "12.0", "analyzed_at": "2026-10-21"}, f)
    d2 = os.path.join(dl, "metformin_type_2_diabetes_data")  # 前缀反解
    os.makedirs(d2, exist_ok=True)
    d3 = os.path.join(dl, "manual_data")                     # 手动模式跳过
    os.makedirs(d3, exist_ok=True)

    diseases = ppi.list_analyzed_diseases()
    assert "cancer" in diseases                              # rwr 变体已剥离
    assert any("type 2 diabetes" == d for d in diseases)
    drugs = ppi.list_analyzed_drugs()
    assert "Aspirin" in drugs and "metformin" in drugs
