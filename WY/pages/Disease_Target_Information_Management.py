import config
from core.disgenet import fetch_disease_targets_refresh as fetch_disgenet_refresh
from core.management import refresh_cache_file
from core.opentargets import fetch_disease_targets_refresh as fetch_opentargets_refresh
from core.targets import (get_disease_targets_merged,
                          invalidate_merged_disease_targets)
from ui.source_management import render_source_management_page

# 各来源文件的抓取逻辑与 merged 缓存失效已下沉 core:
#   opentargets.fetch_disease_targets_refresh / disgenet.fetch_disease_targets_refresh /
#   targets.invalidate_merged_disease_targets,经 _SPEC 交给 core.management.refresh_cache_file 调用。
# 抓取器返回 (表, 状态 "ok"/"empty"/"failed"):状态 "empty" 表示抓取成功但靶点数为 0
# ——不算错误,由页面以警告提示检查查询物;请求失败等异常用 "failed"
# (DisGeNET 密钥缺失为 "failed:no_api_key",管理页据此显示专用提示)。

_SPEC = {
    "cache_dir": config.DISEASE_CACHE_DIR,
    "merged_dir": config.MERGED_CACHE_DIR,
    "merged_prefix": config.MERGED_DISEASE_PREFIX,
    "source_prefixes": ("OpenTargets_", "DisGeNET_"),
    "fetch_source": {"OpenTargets_": fetch_opentargets_refresh,
                     "DisGeNET_": fetch_disgenet_refresh},
    "fetch_merged": get_disease_targets_merged,
    "invalidate_merged": invalidate_merged_disease_targets,
}

# 渲染流程与药物管理页共用(render_source_management_page);
# 更新 TTL 走 refresh 系列默认值
render_source_management_page(
    "disease", _SPEC,
    source_filters=[("Open Targets", "flt_opentargets"), ("DisGeNET", "flt_disgenet")],
)
