"""pytest 公共夹具:路径引导、绘图后端、临时缓存目录。

夹具(fixture)= pytest 替你"自动备料 + 自动收拾"的机制:
测试函数把夹具名写进参数表,pytest 就在测试开始前把东西准备好、
结束后自动回收。
"""
import os
import sys

import pytest

# ---- 让 from WY import config / from WY.core.* 生效(与业务代码同款绝对导入)----
# 可编辑安装(pip install -e .)已提供 WY 可见性时不动 sys.path,
# 与 tests/test_architecture.py 的「sys.path 零容忍」约定保持同一立场。
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
try:
    import WY  # noqa: F401
except ImportError:
    sys.path.insert(0, REPO_ROOT)

# ---- 无界面绘图后端:必须在任何测试 import pyplot 之前钉死 ----
import matplotlib  # noqa: E402

matplotlib.use("Agg", force=True)

# 会写盘的 config 属性全清单(逐一核对过 config.py)。子目录在 config 导入期
# 已从 CACHE_DIR/DOWNLOADS_DIR 派生拼接好,所以必须逐个替换,只换 CACHE_DIR 无效:
# - 缓存目录:十类数据源/算法缓存
# - 对接缓存:四个子目录 + UniProt 基因映射缓存文件
# - 下载/日志:MD 结果、星形文件、GEO 结果、日志目录
# - 单文件:HTTP 缓存、DisGeNET 限流状态、UI 设置、搜索历史
_SANDBOX_ATTRS = (
    "CACHE_DIR",
    "PPI_CACHE_DIR", "ENRICHMENT_CACHE_DIR", "DRUG_CACHE_DIR",
    "DISEASE_CACHE_DIR", "MERGED_CACHE_DIR", "GSEA_GMT_DIR",
    "INTERACTOME_CACHE_DIR", "PROXIMITY_CACHE_DIR", "RWR_CACHE_DIR",
    "DOCKING_CACHE_DIR", "DOCKING_RECEPTOR_CACHE_DIR",
    "DOCKING_LIGAND_CACHE_DIR", "DOCKING_RAW_CACHE_DIR",
    "UNIPROT_LOOKUP_CACHE",
    "DOWNLOADS_DIR", "MD_RESULT_DIR", "MD_STAR_DIR", "GEO_RESULT_DIR",
    "LOG_DIR",
    "HTTP_CACHE_FILE", "DISGENET_RATE_LIMIT_STATE_FILE",
    "UI_SETTINGS_FILE", "SEARCH_HISTORY_FILE",
)


@pytest.fixture(autouse=True)
def _sandbox_caches(tmp_path, monkeypatch):
    """把所有会写盘的缓存/输出路径统一指到本次测试的临时目录。

    autouse=True:不用声明就自动作用于每个测试;tmp_path 是 pytest 自带夹具,
    每个测试独享一个临时文件夹,测试结束自动删除——测试之间互不干扰,
    也绝不触碰真实 Static/ 下的缓存与下载数据。

    有意不进沙箱的属性:
    - CLI_RESULTS_DIR:默认空串即「不写 CLI 结果」,强行指向临时目录反而改变
      语义;需要它的测试自行 patch;
    - VINA_DIR/OBABEL_DIR 及各 EXE:外部工具发现结果,不写盘;
    - GEO_SERIES_CACHE_SUBDIR:相对子目录名,运行时由 CACHE_DIR 拼接,
      换 CACHE_DIR 即已生效。
    """
    from WY import config

    for attr in _SANDBOX_ATTRS:
        monkeypatch.setattr(config, attr, str(tmp_path / attr.lower()))
