"""药物/疾病靶点来源管理的数据源规格(唯一定义点)。

原先在两个 Streamlit 管理页文件里各自定义 _SPEC,TS 前端阶段下沉到 core:
管理页与 API 层共用同一份 spec,改来源清单/目录只改这里。
抓取器本体在各数据源模块(opentargets/disgenet/chembl/dgidb 的 fetch_*_refresh),
更新编排(单文件/批量/进度)在 core.management。
"""

from WY.core.chembl import CHEMBL_API_URL, fetch_drug_targets_refresh
from WY.core.dgidb import fetch_drug_targets_refresh as fetch_dgidb_refresh
from WY.core.disgenet import fetch_disease_targets_refresh
from WY.core.opentargets import (OPENTARGETS_API_URL,
                                 fetch_disease_targets_refresh as fetch_ot_refresh)
from WY.core.targets import (get_disease_targets_merged, get_drug_targets_merged,
                             invalidate_merged_disease_targets,
                             invalidate_merged_drug_targets)


def _hosts(*urls: str) -> tuple[str, ...]:
    """由 API 基址提取域名(强制刷新前定向清 HTTP 响应缓存用)。"""
    import urllib.parse

    return tuple(sorted({urllib.parse.urlsplit(u).netloc for u in urls}))


def drug_spec() -> dict:
    """药物靶点来源管理规格:来源文件目录、整合缓存目录与前缀、各来源抓取器、
    merged 缓存失效器、整合表构建入口。"""
    import WY.config as config

    return {
        "cache_dir": config.DRUG_CACHE_DIR,
        "merged_dir": config.MERGED_CACHE_DIR,
        "merged_prefix": config.MERGED_DRUG_PREFIX,
        "source_prefixes": ("ChEMBL_", "DGIdb_"),
        # 强制刷新前需定向清空的 HTTP 响应缓存域名(按来源前缀;整合文件取全部)
        "cache_hosts": {"ChEMBL_": _hosts(CHEMBL_API_URL),
                        "DGIdb_": _hosts(config.DGIDB_GRAPHQL_URL)},
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
        "cache_hosts": {"OpenTargets_": _hosts(OPENTARGETS_API_URL),
                        "DisGeNET_": _hosts(config.DisGeNET_URL)},
        "fetch_source": {"OpenTargets_": fetch_ot_refresh,
                         "DisGeNET_": fetch_disease_targets_refresh},
        "fetch_merged": get_disease_targets_merged,
        "invalidate_merged": invalidate_merged_disease_targets,
    }


def file_path(filename: str, spec: dict) -> str:
    """整合前缀的文件位于 merged 目录,其余位于来源缓存目录
    (原 core.management._mgmt_file_path,升为公开:API 层同样需要按 spec 解析路径)。
    旧过滤版本的整合缓存(merged_*_v<旧>_ 前缀,TARGET_FILTER_VERSION 升版遗留)
    同样按 merged 目录解析:它们实际就存放在那里,若按来源目录解析会变成
    "列表看得见、size 显示 0 B、删不掉"的僵尸条目。来源文件名不含 merged_ 前缀,不冲突。"""
    import os

    if filename.startswith(spec["merged_prefix"]) or filename.startswith("merged_"):
        return os.path.join(spec["merged_dir"], filename)
    return os.path.join(spec["cache_dir"], filename)


def cache_file_hosts(filename: str, spec: dict) -> tuple[str, ...]:
    """强制刷新某缓存文件前,需定向清空的 HTTP 响应缓存域名。
    来源文件只涉及该来源自己的域名;整合文件或未识别前缀取该类别全部来源域名
    (整合文件的更新链会先刷新两个来源)。
    spec 未声明 cache_hosts(测试自建的简化 spec)时返回空元组,
    由 clear_http_cache 按全量清空兜底,保证刷新必拿新数据。"""
    table = spec.get("cache_hosts") or {}
    for prefix, hosts in table.items():
        if filename.startswith(prefix):
            return hosts
    return tuple(sorted({h for hosts in table.values() for h in hosts}))
