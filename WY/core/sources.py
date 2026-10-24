"""药物/疾病靶点来源管理的数据源规格(唯一定义点)。

原先在两个 Streamlit 管理页文件里各自定义 _SPEC,TS 前端阶段下沉到 core:
管理页与 API 层共用同一份 spec,改来源清单/目录只改这里。
抓取器本体在各数据源模块(opentargets/disgenet/chembl/dgidb 的 fetch_*_refresh),
更新编排(单文件/批量/进度)在 core.management。
"""

from WY.core.chembl import fetch_drug_targets_refresh
from WY.core.dgidb import fetch_drug_targets_refresh as fetch_dgidb_refresh
from WY.core.disgenet import fetch_disease_targets_refresh
from WY.core.opentargets import fetch_disease_targets_refresh as fetch_ot_refresh
from WY.core.targets import (get_disease_targets_merged, get_drug_targets_merged,
                             invalidate_merged_disease_targets,
                             invalidate_merged_drug_targets)


def drug_spec() -> dict:
    """药物靶点来源管理规格:来源文件目录、整合缓存目录与前缀、各来源抓取器、
    merged 缓存失效器、整合表构建入口。"""
    import WY.config as config

    return {
        "cache_dir": config.DRUG_CACHE_DIR,
        "merged_dir": config.MERGED_CACHE_DIR,
        "merged_prefix": config.MERGED_DRUG_PREFIX,
        "source_prefixes": ("ChEMBL_", "DGIdb_"),
        "fetch_source": {"ChEMBL_": fetch_drug_targets_refresh,
                         "DGIdb_": fetch_dgidb_refresh},
        "fetch_merged": get_drug_targets_merged,
        "invalidate_merged": invalidate_merged_drug_targets,
    }


def disease_spec() -> dict:
    """疾病靶点来源管理规格(字段含义同 drug_spec)。"""
    import WY.config as config

    return {
        "cache_dir": config.DISEASE_CACHE_DIR,
        "merged_dir": config.MERGED_CACHE_DIR,
        "merged_prefix": config.MERGED_DISEASE_PREFIX,
        "source_prefixes": ("OpenTargets_", "DisGeNET_"),
        "fetch_source": {"OpenTargets_": fetch_ot_refresh,
                         "DisGeNET_": fetch_disease_targets_refresh},
        "fetch_merged": get_disease_targets_merged,
        "invalidate_merged": invalidate_merged_disease_targets,
    }


def file_path(filename: str, spec: dict) -> str:
    """整合前缀的文件位于 merged 目录,其余位于来源缓存目录
    (原 core.management._mgmt_file_path,升为公开:API 层同样需要按 spec 解析路径)。"""
    import os

    if filename.startswith(spec["merged_prefix"]):
        return os.path.join(spec["merged_dir"], filename)
    return os.path.join(spec["cache_dir"], filename)
