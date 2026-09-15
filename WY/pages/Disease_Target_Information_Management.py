import os

import pandas as pd

import config
from core.disgenet import fetch_disease_targets_disgenet
from core.management import refresh_cache_file
from core.opentargets import get_disease_targets, invalidate_disease_targets_cache
from core.targets import get_disease_targets_merged
from core.text_utils import safe_name
from ui.source_management import render_source_management_page

# 各来源文件的抓取逻辑:返回 (表, 状态 "ok"/"empty"/"failed"),供 refresh_cache_file 调用。
# "empty" 表示抓取成功但靶点数为 0——不算错误,由页面以警告提示检查查询物;请求失败等异常用 "failed"。
def _fetch_opentargets(name: str, cache_hours):
    # 只失效该疾病名的进程内 memo:来源文件已被 refresh_cache_file 移到 .bak,
    # 失效后本次必然重新联网抓取。旧实现 .clear() 清空全部疾病的缓存,又在
    # 批量更新的线程池里并发调用,会把其他会话正在用的疾病缓存一起清掉。
    invalidate_disease_targets_cache(name)
    df = get_disease_targets(name)
    # 该源"查无此病"与"请求失败"都返回空表,无法区分,按 0 靶点处理
    return df, ("empty" if df.empty else "ok")


def _fetch_disgenet(name: str, cache_hours):
    if not config.DisGeNET_API_KEY:
        return pd.DataFrame(), "failed"  # 密钥缺失是配置错误,不属于"0 靶点"
    df = fetch_disease_targets_disgenet(name)
    # 该源"查无关联"与"请求失败"都返回空表,无法区分,按 0 靶点处理
    return df, ("empty" if df.empty else "ok")


def _invalidate_merged_disease(disease_name: str) -> None:
    """更新来源文件后删除对应 merged 缓存，避免整合结果与来源不一致。"""
    # 传入的是原始疾病名,文件名需走同样的 safe_name 规范化才能对上
    path = os.path.join(config.MERGED_CACHE_DIR,
                        f"{config.MERGED_DISEASE_PREFIX}{safe_name(disease_name)}.csv")
    if os.path.exists(path):
        os.remove(path)


_SPEC = {
    "cache_dir": config.DISEASE_CACHE_DIR,
    "merged_dir": config.MERGED_CACHE_DIR,
    "merged_prefix": config.MERGED_DISEASE_PREFIX,
    "source_prefixes": ("OpenTargets_", "DisGeNET_"),
    "fetch_source": {"OpenTargets_": _fetch_opentargets, "DisGeNET_": _fetch_disgenet},
    "fetch_merged": get_disease_targets_merged,
    "invalidate_merged": _invalidate_merged_disease,
}

# 渲染流程与药物管理页共用(render_source_management_page);
# 更新 TTL 走 refresh 系列默认值
render_source_management_page(
    "disease", _SPEC,
    source_filters=[("Open Targets", "flt_opentargets"), ("DisGeNET", "flt_disgenet")],
)
