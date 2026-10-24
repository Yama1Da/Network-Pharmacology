"""WY 包标记:使 WY 成为常规 Python 包,供 pyproject.toml 的 [project.scripts]
入口点(npap = "WY.cli:main")以可预测的方式发现与安装。
保持空实现:包内各模块一律以 WY. 绝对导入(可编辑安装提供可见性);
App.py 与 scripts/selfcheck_architecture.py 保留条件 bootstrap 供未安装环境兜底。"""
