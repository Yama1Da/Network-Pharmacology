"""core.ppi 的小函数补测:蛋白球贴图、排序 CSV 发现、结果目录读取。

深水区冲线批 III(跨线补充)。
"""
import os

import pandas as pd

from WY import config
from WY.core.ppi import (
    _ball_sprite,
    list_ranked_csv_candidates,
    read_ppi_dataset,
    save_module_computation,
)


def test_ball_sprite_shape_and_cache():
    """贴图:rgb/alpha 均为 px×px,同色复用缓存,异色独立。"""
    rgb, alpha = _ball_sprite("#4C72B0", px=32)
    assert rgb.shape == (32, 32, 3) and alpha.shape == (32, 32)
    rgb2, _ = _ball_sprite("#4C72B0", px=32)
    assert rgb2 is rgb  # 同色命中缓存
    rgb3, _ = _ball_sprite("#DD8452", px=32)
    assert rgb3 is not rgb


def test_save_module_computation_and_read_ppi_dataset(tmp_path, monkeypatch):
    """模块计算落盘 + 结果目录回读:hub 集合来自 is_hub 列、富集基因集来自
    genes 列;缺列/缺文件按空集合兜底。"""
    from WY.core.ppi import save_module_results

    monkeypatch.setattr(config, "DOWNLOADS_DIR", str(tmp_path))
    directory = str(tmp_path / "metformin_T2D_data")
    os.makedirs(directory)
    computed = {"summary_df": pd.DataFrame({"module_id": [1], "size": [2]}),
                "long_df": pd.DataFrame({"module_id": [1], "genes": ["TP53"]})}
    out = save_module_computation(directory, "metformin_T2D", computed)
    assert out["saved_paths"]
    nodes = pd.DataFrame({"gene": ["TP53", "EGFR"], "is_hub": [True, False]})
    nodes.to_csv(os.path.join(directory, "metformin_T2D_nodes.csv"),
                 index=False, encoding="utf-8-sig")
    enr = pd.DataFrame({"genes": ["TP53;BRCA1"], "source": ["KEGG"]})
    enr.to_csv(os.path.join(directory, "metformin_T2D_enrichment.csv"),
               index=False, encoding="utf-8-sig")
    hub, genes = read_ppi_dataset("metformin_T2D_data")
    assert hub == {"TP53"} and genes == {"TP53", "BRCA1"}
    # 缺 is_hub 列的历史版本:按空集合兜底不崩溃
    bad_nodes = pd.DataFrame({"gene": ["TP53"]})
    bad_nodes.to_csv(os.path.join(directory, "metformin_T2D_nodes.csv"),
                     index=False, encoding="utf-8-sig")
    hub2, _ = read_ppi_dataset("metformin_T2D_data")
    assert hub2 == set()
    # 目录不存在:两个集合都为空
    hub3, genes3 = read_ppi_dataset("no_such_data")
    assert hub3 == set() and genes3 == set()


def test_list_ranked_csv_candidates(tmp_path, monkeypatch):
    """排序 CSV 发现:只认 *_query_genes.csv,label 由目录前缀/时间/版本拼成。"""
    import json

    dl = config.DOWNLOADS_DIR
    d1 = os.path.join(dl, "aspirin_cancer_data")
    os.makedirs(d1, exist_ok=True)
    with open(os.path.join(d1, "analysis_meta.json"), "w", encoding="utf-8") as f:
        json.dump({"drug": "Aspirin", "disease": "cancer", "string_version": "12.0",
                   "analyzed_at": "2026-10-21T10:00:00"}, f)
    pd.DataFrame({"gene": ["TP53"], "source": ["seed"],
                  "rwr_score": [1.0]}).to_csv(
        os.path.join(d1, "aspirin_cancer_query_genes.csv"), index=False)
    pd.DataFrame({"gene": ["TP53"]}).to_csv(
        os.path.join(d1, "aspirin_cancer_nodes.csv"), index=False)  # 非排序表,不收

    items = list_ranked_csv_candidates()
    hit = [i for i in items if i["drug"] == "Aspirin"]
    assert len(hit) == 1
    assert hit[0]["path"].endswith("aspirin_cancer_query_genes.csv")
    assert hit[0]["base_disease"] == "cancer"
    assert "STRING 12.0" in hit[0]["label"]
