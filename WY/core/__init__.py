"""core 业务域包：按数据源 / 功能拆分的 Utils 实现。

模块与依赖方向（上层引下层，禁止反向）：
    ppi / enrichment / docking / targets / files
        → opentargets / disgenet / chembl / dgidb
        → http_client / caching / text_utils / logging_setup / paths
        → config

对外入口是 WY/Utils.py 门面（re-export），页面代码无需直接 import 本包。
"""
