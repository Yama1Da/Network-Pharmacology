import os

import config
import Utils

# 各来源文件的抓取逻辑:返回 (表, 是否成功),供 Utils.refresh_cache_file 调用
def _fetch_opentargets(name: str, cache_hours):
    # 不使用缓存,强制重新请求 API
    Utils.get_disease_targets.clear()
    df = Utils.get_disease_targets(name)
    # 文件名经过规范化处理,按同样规则检查新文件是否生成
    ok = not df.empty and os.path.exists(
        os.path.join(config.DISEASE_CACHE_DIR, f"OpenTargets_{Utils.safe_name(name)}.csv"))
    return df, ok


def _fetch_disgenet(name: str, cache_hours):
    df = Utils.fetch_disease_targets_disgenet(name)
    ok = not df.empty and os.path.exists(
        os.path.join(config.DISEASE_CACHE_DIR, f"DisGeNET_{Utils.safe_name(name)}.csv"))
    return df, ok


def _invalidate_merged_disease(disease_name: str) -> None:
    """更新来源文件后删除对应 merged 缓存，避免整合结果与来源不一致。"""
    # 传入的是原始疾病名,文件名需走同样的 safe_name 规范化才能对上
    path = os.path.join(config.MERGED_CACHE_DIR,
                        f"{config.MERGED_DISEASE_PREFIX}{Utils.safe_name(disease_name)}.csv")
    if os.path.exists(path):
        os.remove(path)


_SPEC = {
    "cache_dir": config.DISEASE_CACHE_DIR,
    "merged_dir": config.MERGED_CACHE_DIR,
    "merged_prefix": config.MERGED_DISEASE_PREFIX,
    "source_prefixes": ("OpenTargets_", "DisGeNET_"),
    "fetch_source": {"OpenTargets_": _fetch_opentargets, "DisGeNET_": _fetch_disgenet},
    "fetch_merged": Utils.get_disease_targets_merged,
    "invalidate_merged": _invalidate_merged_disease,
}

# 渲染流程与药物管理页共用(Utils.render_source_management_page);
# 更新 TTL 走 refresh 系列默认值
Utils.render_source_management_page(
    "disease", _SPEC,
    source_filters=[("Open Targets", "flt_opentargets"), ("DisGeNET", "flt_disgenet")],
)
