import os

import config
import Utils

# 各来源文件的抓取逻辑:返回 (表, 是否成功),供 Utils.refresh_cache_file 调用
def _fetch_chembl(name: str, cache_hours):
    Utils.set_chembl_cache_ttl(cache_hours)
    df, csv_path = Utils.fetch_drug_targets_chembl(name)
    ok = not df.empty and bool(csv_path) and os.path.exists(csv_path)
    return df, ok


def _fetch_dgidb(name: str, cache_hours):
    df = Utils.fetch_drug_targets_dgidb(name)
    # 文件名经过规范化处理,按同样规则检查新文件是否生成
    ok = not df.empty and os.path.exists(
        os.path.join(config.DRUG_CACHE_DIR, f"DGIdb_{Utils.safe_name(name)}.csv"))
    return df, ok


def _invalidate_merged_drug(drug_name: str) -> None:
    """更新来源文件后删除对应 merged 缓存，避免整合结果与来源不一致。"""
    # 传入的是原始药物名,文件名需走同样的 safe_name 规范化才能对上
    path = os.path.join(config.MERGED_CACHE_DIR,
                        f"{config.MERGED_DRUG_PREFIX}{Utils.safe_name(drug_name)}.csv")
    if os.path.exists(path):
        os.remove(path)


_SPEC = {
    "cache_dir": config.DRUG_CACHE_DIR,
    "merged_dir": config.MERGED_CACHE_DIR,
    "merged_prefix": config.MERGED_DRUG_PREFIX,
    "source_prefixes": ("ChEMBL_", "DGIdb_"),
    "fetch_source": {"ChEMBL_": _fetch_chembl, "DGIdb_": _fetch_dgidb},
    "fetch_merged": Utils.get_drug_targets_merged,
    "invalidate_merged": _invalidate_merged_drug,
}

# 渲染流程与疾病管理页共用(Utils.render_source_management_page);
# ChEMBL 已改为批量分页,单个文件更新通常几十秒内完成,短时间内重复更新只会刷新
# 时间戳,故全部/单文件更新的 TTL 固定 12h/1h
Utils.render_source_management_page(
    "drug", _SPEC,
    source_filters=[("ChEMBL", "flt_chembl"), ("DGIdb", "flt_dgidb")],
    cache_hours_all=12,
    cache_hours_single=1,
)
