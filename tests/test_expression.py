"""core.expression 的测试:BH-FDR 手算对照、探针识别回归、log2 启发式回归、
RNA-seq DESeq2 路径。

迁移自 selfcheck_ml.py(17 项中的 4 项)。
"""
import numpy as np
import pandas as pd

import WY.core.expression as expr


def test_bh_fdr():
    """BH 手算对照:BH 1995 经典 15 个 p 值的 FDR。"""
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205,
                  0.212, 0.216, 0.222, 0.251, 0.269, 0.275, 0.340])
    q = expr.benjamini_hochberg(p)
    # 手算:BH 标准公式 p(i)*n/rank,再取反向累积最小值
    n = len(p)
    order = np.argsort(p)
    expected = np.empty(n)
    raw = p[order] * n / (np.arange(n) + 1)
    expected[order] = np.clip(np.minimum.accumulate(raw[::-1])[::-1], 0, 1)
    assert np.allclose(q, expected, atol=1e-10), f"BH 不匹配:\n{q}\n{expected}"
    # 边界:单基因
    assert abs(expr.benjamini_hochberg(np.array([0.01]))[0] - 0.01) < 1e-10


def test_probe_identifier():
    """BUG-1 回归:Ensembl/Illumina/Agilent 探针不再被误判为 symbol 级。"""
    assert expr.looks_symbol_level(["ENSG00000141510", "ENSG00000139618"]) is False, \
        "Ensembl ID 被误判为 symbol 级"
    assert expr.looks_symbol_level(["ILMN_123456", "ILMN_234567"]) is False, \
        "Illumina 探针被误判为 symbol 级"
    assert expr.looks_symbol_level(["A_23_P123456", "A_24_P234567"]) is False, \
        "Agilent 探针被误判为 symbol 级"
    assert expr.looks_symbol_level(["TP53", "TNF", "VEGFA", "AKT1"]) is True, \
        "真 symbol 被误拒"
    assert expr.looks_symbol_level(["1007_s_at", "1053_at"]) is False, \
        "Affy 探针被误判为 symbol 级"


def test_log2_heuristic():
    """BUG-3 回归:大尺度含负值应做 log2;z-score/已 log 不变换。"""
    rng = np.random.default_rng(0)
    # 大尺度含负值(均值 500,标准差 200)→ log2
    _, b1 = expr.preprocess(pd.DataFrame(rng.normal(500, 200, (18, 50))))
    assert "log2" in b1, f"大尺度含负值未做 log2: {b1}"
    # z-score(均值 0 含负值)→ 不变换
    _, b2 = expr.preprocess(pd.DataFrame(rng.normal(0, 1, (18, 50))))
    assert "不变换" in b2, f"z-score 被误变换: {b2}"
    # 已 log(均值 8 无负值)→ 不变换
    _, b3 = expr.preprocess(pd.DataFrame(rng.normal(8, 2, (18, 50))))
    assert "不变换" in b3, f"已 log 被误变换: {b3}"
    # force_log 覆盖
    _, b4 = expr.preprocess(pd.DataFrame(rng.normal(30, 15, (18, 50))), force_log=True)
    assert "log2" in b4, f"force_log 未生效: {b4}"


def test_rnaseq_deg_matrix():
    """二期 RNA-seq 路径:合成 NB counts → DESeq2 一次拟合产出 DEG + ML 矩阵。"""
    rng = np.random.default_rng(21)
    genes = [f"G{i:02d}" for i in range(20)]
    ctrl = rng.poisson(50, size=(20, 3))
    case = rng.poisson(50, size=(20, 3))
    case[:6] = rng.poisson(200, size=(6, 3))
    counts = pd.DataFrame(np.hstack([ctrl, case]), index=genes,
                          columns=[f"C{i}" for i in range(3)] + [f"T{i}" for i in range(3)])
    labels = pd.Series([0] * 3 + [1] * 3, index=list(counts.columns))
    deg, ml_matrix = expr.rnaseq_deg_and_matrix(counts, labels)
    assert len(set(genes[:6]) & set(deg.head(8)["gene"])) >= 5, "DESeq2 信号召回不足"
    assert set(ml_matrix.columns) == set(counts.index), "ML 矩阵列应为基因"
    assert set(ml_matrix.index) == set(labels.index), "ML 矩阵行应为样本"
    assert ml_matrix.to_numpy().min() >= 0, "log2 标准化矩阵不应有负值"
    # 非 counts 拒收
    try:
        expr.rnaseq_deg_and_matrix(counts.astype(float) / 3, labels)
        raise AssertionError("非 counts 输入未拒收")
    except RuntimeError:
        pass
