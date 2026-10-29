"""core.expression 的组间差异检验测试:Welch t + BH 校正的输出契约与坏数据防线。

阶段三第 5 批新写(自检只测过 BH 手算小例,这里补完整管线)。
"""
import numpy as np
import pandas as pd

import WY.core.expression as expr


def _synthetic(n_ctrl=3, n_case=3, n_genes=8, n_planted=4, seed=7):
    """合成矩阵:前 n_planted 个基因在 case 组抬升,末尾补一个全 NaN 基因。"""
    rng = np.random.default_rng(seed)
    cols = [f"C{i}" for i in range(n_ctrl)] + [f"T{i}" for i in range(n_case)]
    X = rng.normal(0, 1, (n_genes, n_ctrl + n_case))
    X[:n_planted, n_ctrl:] += 3.0
    genes = [f"G{i:02d}" for i in range(n_genes)]
    matrix = pd.DataFrame(X, index=genes, columns=cols)
    labels = pd.Series([0] * n_ctrl + [1] * n_case, index=cols)
    return matrix, labels, genes[:n_planted]


def test_differential_expression_contract():
    """输出契约:列齐全、按 p 升序、方向标注正确、FDR 单调不小于 p。"""
    matrix, labels, planted = _synthetic()
    res = expr.differential_expression(matrix, labels)
    assert list(res.columns) == ["gene", "log2FC", "t", "p", "FDR", "direction"]
    assert res["p"].is_monotonic_increasing
    # 植入基因(强信号)应占居前若干名
    assert set(planted) >= set(res.head(len(planted))["gene"])
    up_rows = res[res["gene"].isin(planted)]
    assert (up_rows["direction"] == "up").all()
    assert (up_rows["log2FC"] > 0).all()
    # FDR 是 p 的校正,单侧不低于原 p(数值容差内)
    assert (res["FDR"] + 1e-12 >= res["p"]).all()


def test_differential_expression_drops_all_nan_genes():
    """两组中全为 NaN 的基因被过滤,不产出 NaN 的 p 值行。"""
    matrix, labels, _ = _synthetic()
    matrix.loc["G00", :] = np.nan  # 某基因整行 NaN
    res = expr.differential_expression(matrix, labels)
    assert "G00" not in set(res["gene"])
    assert res["p"].notna().all()


def test_differential_expression_down_direction():
    """反向植入(case 组下调):方向列如实标 down、log2FC 为负。"""
    matrix, labels, _ = _synthetic(seed=11)
    matrix.loc["G07", ["T0", "T1", "T2"]] -= 5.0
    res = expr.differential_expression(matrix, labels)
    row = res[res["gene"] == "G07"].iloc[0]
    assert row["direction"] == "down" and row["log2FC"] < 0
