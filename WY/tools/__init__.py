"""WY.tools:独立命令行工具入口(互作组构建 / i18n 词典导出 / GEO-ML 筛选)。

与 benchmark/ 同为"人工运行的入口层":不进覆盖率分母,不与 web/cli 构成导入关系;
各脚本直接运行场景下的条件 bootstrap(sys.path 兜底)见 test_architecture 白名单。
"""
