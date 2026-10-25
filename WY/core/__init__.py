"""core 业务域包：数据抓取、分析与缓存管理的领域/编排层（不依赖任何 UI 框架）。

模块与依赖方向（上层引下层，禁止反向）：
    analysis / management / network_algorithms / ppi / enrichment / docking / targets / files
        → opentargets / disgenet / chembl / dgidb
        → http_client / caching / text_utils / logging_setup
        → config（路径常量统一见 config.py：CACHE_DIR / DRUG_CACHE_DIR 等）

界面层已迁 TS+Vue(WY/web),经 WY/api.py(FastAPI)访问本包;
没有集中门面,API 端点即门面。"""
