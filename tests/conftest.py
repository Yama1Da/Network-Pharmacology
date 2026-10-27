"""pytest 公共配置:让测试能找到项目代码。

pytest 约定:tests/ 下的 conftest.py 会在收集测试前自动执行,
不需要任何人 import 它,这里做的事对全部测试生效。
"""
import os
import sys

# 把仓库根(而非 WY/ 目录)加进 Python 的模块搜索路径。
# 全仓已统一为 WY. 绝对导入(cli.py 与自检脚本同款),测试里写
# from WY.core.xxx import ...、from WY import cli,而不是裸名 import config。
# 条件引导:可编辑安装(pip install -e .)已提供 WY 可见性时不动 sys.path,
# 与 selfcheck_architecture 的「sys.path 零容忍」约定保持同一立场。
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
try:
    import WY  # noqa: F401
except ImportError:
    sys.path.insert(0, REPO_ROOT)
