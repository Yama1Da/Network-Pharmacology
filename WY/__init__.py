"""WY 包标记:使 WY 成为常规 Python 包,供 pyproject.toml 的 [project.scripts]
入口点(npap = "WY.cli:main")以可预测的方式发现与安装。
包内各模块一律以 WY. 绝对导入(可编辑安装提供可见性);
tools/ 与 benchmark/ 下各独立入口脚本保留条件 bootstrap 供未安装环境兜底。"""

# matplotlib 后端全局钉 Agg(必须在任何 pyplot 导入前生效):本产品只在服务端
# 出文件不出窗,而分析任务跑在 jobs 工作线程——GUI 后端在非主线程建画布会刷
# "Starting a Matplotlib GUI outside of the main thread" 警告(enrichment/GSEA/
# ppi_plots 均曾触发);Agg 无 GUI,工作线程安全,一劳永逸(2026-11-27)。
import os

os.environ.setdefault("MPLBACKEND", "Agg")
