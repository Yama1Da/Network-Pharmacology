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

# 互作组三目录仅在「双条件齐全」时豁免沙箱:NPAP_TEST_REAL_INTERACTOME=1
# 且当前测试带 @pytest.mark.real_interactome 标记。单靠环境变量不够——
# 缓存治理类测试(test_interactome_cache.py)默认假设自己在空沙箱里演练,
# 让它们见到真实目录会覆盖/删除真文件(2026-10 曾因此毁掉真实 v12.0 三件套,
# 靠 tools/build_interactome.py 重建才恢复)。
# 标记的 4 条集成测试在无数据的机器上 is_built() 为假 → 自动跳过,开关无副作用。
# 豁免期间 *_real_and_cached 会把 selfcheck_* 前缀的缓存文件写进真实缓存
# 目录(带前缀自清理,与自检脚本时代同款)。
_REAL_INTERACTOME_ATTRS = (
    "INTERACTOME_CACHE_DIR", "PROXIMITY_CACHE_DIR", "RWR_CACHE_DIR",
)


@pytest.fixture(autouse=True)
def _sandbox_caches(request, tmp_path, monkeypatch):
    """把所有会写盘的缓存/输出路径统一指到本次测试的临时目录。

    autouse=True:不用声明就自动作用于每个测试;tmp_path 是 pytest 自带夹具,
    每个测试独享一个临时文件夹,测试结束自动删除——测试之间互不干扰,
    也绝不触碰真实 Static/ 下的缓存与下载数据。

    有意不进沙箱的属性:
    - CLI_RESULTS_DIR:默认空串即「不写 CLI 结果」,强行指向临时目录反而改变
      语义;需要它的测试自行 patch;
    - VINA_DIR/OBABEL_DIR 及各 EXE:外部工具发现结果,不写盘;
    - GEO_SERIES_CACHE_SUBDIR:相对子目录名,运行时由 CACHE_DIR 拼接,
      换 CACHE_DIR 即已生效;
    - 互作组三目录:仅在 NPAP_TEST_REAL_INTERACTOME=1 且测试带
      @pytest.mark.real_interactome 时豁免(双条件,缺一不可)。
    """
    from WY import config

    attrs = _SANDBOX_ATTRS
    if os.environ.get("NPAP_TEST_REAL_INTERACTOME") == "1" \
            and request.node.get_closest_marker("real_interactome"):
        attrs = tuple(a for a in _SANDBOX_ATTRS if a not in _REAL_INTERACTOME_ATTRS)
        # 可见性防线:碰真实数据必须留在终端里。print 连同 sys.__stdout__
        # 都会被 pytest 的 fd 级捕获吞掉,只有 terminal reporter 不受捕获影响
        reporter = request.config.pluginmanager.getplugin("terminalreporter")
        if reporter is not None:
            reporter.write_line(
                f"[real-interactome] 沙箱已豁免: {request.node.nodeid} "
                f"正接触真实目录 {', '.join(_REAL_INTERACTOME_ATTRS)}",
                yellow=True)
    for attr in attrs:
        monkeypatch.setattr(config, attr, str(tmp_path / attr.lower()))


# ---- ML 测试降耗开关:默认瘦身跑,设 NPAP_ML_TEST_FULL=1 恢复生产参数 ----
# Boruta 影子检验每次迭代都要训一棵 200 树随机森林,是 ML 套件约八成耗时的来源。
# 日常 CI 只验证行为语义(信号召回/种子一致/LOOCV 切换/分母回归),不需要
# 12/5 轮的完整统计动力学;改统计参数(_BORUTA_ITER_*/算法超参)或发版前,
# 手动跑一次全参档:NPAP_ML_TEST_FULL=1 pytest
# 降档语义(有意为之,知悉再用):
# - 全量档 12→6:二项确认门仍可达(需 6/6 命中,p=0.016<0.05),确认分支仍被走到;
# - bootstrap 档 5→3:确认门在 3 轮下数学上不可达(p 最小 0.125>0.05),
#   rf_boruta 退化为纯 RF 增益(走产品代码里 confirmed 为空的分支,同样是定义行为)。
@pytest.fixture(autouse=True)
def _ml_slim_boruta(monkeypatch):
    if os.environ.get("NPAP_ML_TEST_FULL") == "1":
        yield
        return
    import WY.core.ml_screening as ml

    monkeypatch.setattr(ml, "_BORUTA_ITER_FULL", 6)
    monkeypatch.setattr(ml, "_BORUTA_ITER_BOOT", 3)
    yield
