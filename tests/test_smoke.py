"""第一个测试:确认测试环境本身是通的。"""
from WY import config


def test_config_importable():
    """config 能导入,且一个已知常量取值正确。"""
    assert config.STRING_QUERY_LIMIT == 1000
