"""core 业务域包：数据抓取、分析与缓存管理的领域/编排层（不依赖任何 UI 框架）。

模块与依赖方向（上层引下层，禁止反向）：
    management / network_algorithms / ppi / enrichment / docking / targets / files
        → opentargets / disgenet / chembl / dgidb
        → http_client / caching / text_utils / logging_setup
        → config（路径常量统一见 config.py：CACHE_DIR / DRUG_CACHE_DIR 等）

界面层请直接 import 本包与 ui/（界面组件见 ui/widgets.py 与
ui/source_management.py），没有集中门面。"""
