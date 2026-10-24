from WY.core import sources
from WY.ui.source_management import render_source_management_page

# 数据源规格(抓取器/缓存失效/目录)唯一定义在 core.sources,与 TS 前端 API 层共用;
# 状态语义("ok"/"empty"/"failed")见各抓取器 docstring。
_SPEC = sources.drug_spec()

# 渲染流程与疾病管理页共用(render_source_management_page);
# ChEMBL 已改为批量分页,单个文件更新通常几十秒内完成,短时间内重复更新只会刷新
# 时间戳,故全部/单文件更新的 TTL 固定 12h/1h
render_source_management_page(
    "drug", _SPEC,
    source_filters=[("ChEMBL", "flt_chembl"), ("DGIdb", "flt_dgidb")],
    cache_hours_all=12,
    cache_hours_single=1,
)
