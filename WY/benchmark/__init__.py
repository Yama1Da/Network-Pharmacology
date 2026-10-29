"""WY.benchmark:离线基准验证入口(LIT-PCBA / 恢复型 / 规模扩展 / 跨数据集)。

与 tools/ 同为"人工运行的入口层":不进覆盖率分母,不与 web/cli 构成导入关系;
指标计算保持纯函数(可直接被 pytest 复用),各脚本直接运行场景下的条件
bootstrap(sys.path 兜底)见 test_architecture 白名单。
"""
