"""WY 包标记:使 WY 成为常规 Python 包,供 pyproject.toml 的 [project.scripts]
入口点(npap = "WY.cli:main")以可预测的方式发现与安装。
保持空实现:包内各模块(config/core/i18n/ui...)一律按"脚本目录在 sys.path"的
顶层模块方式导入,与就地运行(streamlit run App.py)的既有约定保持一致。"""
