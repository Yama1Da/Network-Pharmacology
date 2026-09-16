import os

import config
from core.chembl import fetch_drug_targets_chembl, set_chembl_cache_ttl
from core.dgidb import fetch_drug_targets_dgidb_stated
from core.management import refresh_cache_file
from core.targets import get_drug_targets_merged
from core.text_utils import safe_name
from ui.source_management import render_source_management_page

# 各来源文件的抓取逻辑:返回 (表, 状态 "ok"/"empty"/"failed"),供 refresh_cache_file 调用。
# "empty" 表示抓取成功但靶点数为 0——不算错误,由页面以警告提示检查查询物;请求失败等异常用 "failed"。
def _fetch_chembl(name: str, cache_hours):
    set_chembl_cache_ttl(cache_hours)
    df, csv_path = fetch_drug_targets_chembl(name)
    if csv_path is None:
        return df, "failed"  # 请求中止/异常(文件未写出),可可靠区分,保留为错误
    return df, ("empty" if df.empty else "ok")  # 查无此药 → empty


def _fetch_dgidb(name: str, cache_hours):
    # 带状态变体:请求失败(failed)与查无相互作用(empty)分开上报,不再把网络故障误报成 0 靶点
    df, status = fetch_drug_targets_dgidb_stated(name)
    return df, status


def _invalidate_merged_drug(drug_name: str) -> None:
    """更新来源文件后删除对应 merged 缓存，避免整合结果与来源不一致。"""
    # 传入的是原始药物名,文件名需走同样的 safe_name 规范化才能对上
    path = os.path.join(config.MERGED_CACHE_DIR,
                        f"{config.MERGED_DRUG_PREFIX}{safe_name(drug_name)}.csv")
    if os.path.exists(path):
        os.remove(path)


_SPEC = {
    "cache_dir": config.DRUG_CACHE_DIR,
    "merged_dir": config.MERGED_CACHE_DIR,
    "merged_prefix": config.MERGED_DRUG_PREFIX,
    "source_prefixes": ("ChEMBL_", "DGIdb_"),
    "fetch_source": {"ChEMBL_": _fetch_chembl, "DGIdb_": _fetch_dgidb},
    "fetch_merged": get_drug_targets_merged,
    "invalidate_merged": _invalidate_merged_drug,
}

# 渲染流程与疾病管理页共用(render_source_management_page);
# ChEMBL 已改为批量分页,单个文件更新通常几十秒内完成,短时间内重复更新只会刷新
# 时间戳,故全部/单文件更新的 TTL 固定 12h/1h
render_source_management_page(
    "drug", _SPEC,
    source_filters=[("ChEMBL", "flt_chembl"), ("DGIdb", "flt_dgidb")],
    cache_hours_all=12,
    cache_hours_single=1,
)
