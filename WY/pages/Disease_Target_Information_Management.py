from WY.core import sources
from WY.ui.source_management import render_source_management_page

# 数据源规格(抓取器/缓存失效/目录)唯一定义在 core.sources,与 TS 前端 API 层共用;
# 状态语义("ok"/"empty"/"failed"、"failed:no_api_key")见各抓取器 docstring。
_SPEC = sources.disease_spec()

# 渲染流程与药物管理页共用(render_source_management_page);
# 更新 TTL 走 refresh 系列默认值
render_source_management_page(
    "disease", _SPEC,
    source_filters=[("Open Targets", "flt_opentargets"), ("DisGeNET", "flt_disgenet")],
)
