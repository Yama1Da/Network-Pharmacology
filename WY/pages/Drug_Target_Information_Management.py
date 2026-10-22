import config
from core.chembl import fetch_drug_targets_refresh as fetch_chembl_refresh
from core.dgidb import fetch_drug_targets_refresh as fetch_dgidb_refresh
from core.management import refresh_cache_file
from core.targets import get_drug_targets_merged, invalidate_merged_drug_targets
from ui.source_management import render_source_management_page

# 各来源文件的抓取逻辑与 merged 缓存失效已下沉 core:
#   chembl.fetch_drug_targets_refresh / dgidb.fetch_drug_targets_refresh /
#   targets.invalidate_merged_drug_targets,经 _SPEC 交给 core.management.refresh_cache_file 调用。
# 抓取器返回 (表, 状态 "ok"/"empty"/"failed"):状态 "empty" 表示抓取成功但靶点数为 0
# ——不算错误,由页面以警告提示检查查询物;请求失败等异常用 "failed"。

_SPEC = {
    "cache_dir": config.DRUG_CACHE_DIR,
    "merged_dir": config.MERGED_CACHE_DIR,
    "merged_prefix": config.MERGED_DRUG_PREFIX,
    "source_prefixes": ("ChEMBL_", "DGIdb_"),
    "fetch_source": {"ChEMBL_": fetch_chembl_refresh,
                     "DGIdb_": fetch_dgidb_refresh},
    "fetch_merged": get_drug_targets_merged,
    "invalidate_merged": invalidate_merged_drug_targets,
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
