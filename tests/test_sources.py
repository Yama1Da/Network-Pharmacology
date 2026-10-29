"""core.sources 的测试:药物/疾病管理规格的结构完整性,文件路径路由规则。

阶段三第 3 批新写(0% 覆盖小模块)。spec 是管理页与 API 层共用的唯一定义点。
"""
import os

from WY import config
from WY.core import sources


def test_drug_spec_structure():
    """药物规格:目录/前缀/抓取器/失效器/整合入口五件套齐备,来源前缀与 config 一致。"""
    spec = sources.drug_spec()
    assert spec["cache_dir"] == config.DRUG_CACHE_DIR
    assert spec["merged_dir"] == config.MERGED_CACHE_DIR
    assert spec["merged_prefix"] == config.MERGED_DRUG_PREFIX
    assert spec["source_prefixes"] == ("ChEMBL_", "DGIdb_")
    assert set(spec["fetch_source"]) == {"ChEMBL_", "DGIdb_"}
    assert callable(spec["fetch_merged"]) and callable(spec["invalidate_merged"])


def test_disease_spec_structure():
    """疾病规格:同构,来源前缀为 OpenTargets_/DisGeNET_。"""
    spec = sources.disease_spec()
    assert spec["cache_dir"] == config.DISEASE_CACHE_DIR
    assert spec["merged_prefix"] == config.MERGED_DISEASE_PREFIX
    assert spec["source_prefixes"] == ("OpenTargets_", "DisGeNET_")
    assert set(spec["fetch_source"]) == {"OpenTargets_", "DisGeNET_"}


def test_file_path_routing():
    """路径路由:整合前缀的文件在 merged 目录,来源文件在缓存目录(两页共用规则)。"""
    spec = {"merged_prefix": "merged_drug_v1_", "merged_dir": "/m", "cache_dir": "/c"}
    assert sources.file_path("merged_drug_v1_metformin.csv", spec) == \
        os.path.join("/m", "merged_drug_v1_metformin.csv")
    assert sources.file_path("ChEMBL_metformin.csv", spec) == \
        os.path.join("/c", "ChEMBL_metformin.csv")
