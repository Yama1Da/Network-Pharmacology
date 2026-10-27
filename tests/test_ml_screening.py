"""core.ml_screening 的测试:六算法共识的信号召回/种子一致性/退化对照/LOOCV
自动切换/方向防御,以及 bootstrap 分母回归。

迁移自 selfcheck_ml.py(17 项中的 6 项)。纯本地合成数据,不联网。
"""
import numpy as np
import pandas as pd

import WY.core.ml_screening as ml


def _make_synthetic(n_ctrl=10, n_case=10, n_genes=200, n_planted=10, seed=42):
    """合成数据:前 n_planted 个基因在 case 组抬升 shift(区分信号),其余为噪声。"""
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    X = rng.normal(0, 1, (n_ctrl + n_case, n_genes))
    X[n_ctrl:, :n_planted] += 2.0
    X = pd.DataFrame(X, columns=genes)
    y = pd.Series([0] * n_ctrl + [1] * n_case)
    return X, y, genes, genes[:n_planted]


def test_signal_recall():
    """合成数据信号:共识 top15 应召回 ≥8 个植入基因,方向正确。"""
    X, y, genes, planted = _make_synthetic()
    res = ml.run_ml_screening(X, y, bootstrap_rounds=20)
    top15 = res["consensus"].head(15).index.tolist()
    recall = len(set(planted) & set(top15))
    assert recall >= 8, f"信号召回不足: {recall}/15 (期望 ≥8)"
    cs = res["consensus"]["consensus_score"]
    assert cs[planted].mean() > cs[[g for g in genes if g not in planted]].median(), \
        "方向错误:植入基因平均分应高于噪声中位数"


def test_seed_consistency():
    """种子一致性:同种子两次运行逐基因一致。"""
    X1, y, _, _ = _make_synthetic(seed=42)
    X2, _, _, _ = _make_synthetic(seed=42)
    r1 = ml.run_ml_screening(X1, y, bootstrap_rounds=10, seed=0)
    r2 = ml.run_ml_screening(X2, y, bootstrap_rounds=10, seed=0)
    cs1 = r1["consensus"]["consensus_score"].round(6).to_dict()
    cs2 = r2["consensus"]["consensus_score"].round(6).to_dict()
    assert cs1 == cs2, "同种子两次运行结果不一致"


def test_label_shuffle():
    """打乱标签:AUC 应回落 ~0.5,植入基因不再富集于 top15。"""
    X, y, _, planted = _make_synthetic(seed=99)
    rng = np.random.default_rng(77)
    y_shuf = pd.Series(rng.permutation(y.to_numpy()))
    res = ml.run_ml_screening(X, y_shuf, bootstrap_rounds=10)
    auc = res["meta"]["cv_auc"]
    assert 0.20 < auc < 0.80, f"打乱后 AUC={auc} 异常(期望 ~0.5)"
    top15 = res["consensus"].head(15).index.tolist()
    n_planted_top = len(set(planted) & set(top15))
    assert n_planted_top <= 4, f"打乱后植入基因仍有 {n_planted_top} 个在 top15(不应富集)"


def test_loocv_auto():
    """LOOCV 自动切换:最小类 ≤6 时 cv_mode 应为 loocv。"""
    X, y, _, _ = _make_synthetic(n_ctrl=6, n_case=6, n_genes=50, seed=42)
    res = ml.run_ml_screening(X, y, bootstrap_rounds=5)
    assert res["meta"]["cv_mode"] == "loocv", f"期望 loocv,实际 {res['meta']['cv_mode']}"


def test_direction_defense():
    """方向防御:传入转置矩阵应报错(行数≠标签数 或 行>列)。"""
    X, y, _, _ = _make_synthetic(n_genes=200, seed=42)
    X_t = X.T  # 故意传反:基因×样本 (200行 × 20列)
    try:
        ml.run_ml_screening(X_t, y, bootstrap_rounds=1)
        raise AssertionError("方向防御未触发(传反矩阵未报错)")
    except RuntimeError as e:
        assert "方向" in str(e) or "转置" in str(e), f"报错信息不含方向提示: {e}"


def test_bootstrap_denominator(monkeypatch):
    """分母回归(2026-10 审查修复):失败轮不计入分母,成功轮必选基因频率应为 1.0。"""
    X, y, genes, _ = _make_synthetic(n_ctrl=6, n_case=6, n_genes=30, seed=1)
    calls = {"n": 0}

    def fake_select(Xr, yr, seed, boruta_iter=5):
        calls["n"] += 1
        if calls["n"] <= 2:
            raise RuntimeError("模拟选择器失败")
        # 确定性重要性 = 列序:每个成功轮都选中末 k 个基因
        return {alg: pd.Series(np.arange(len(Xr.columns), dtype=float), index=Xr.columns)
                for alg in ("lasso", "enet", "ridge", "rf_boruta", "svm_rfe", "xgboost")}

    monkeypatch.setattr(ml, "_select_all", fake_select)
    freq, n_done = ml.bootstrap_stability(X, y, rounds=10, seed=0)
    assert n_done == 8, f"成功轮数应为 8(10 轮 - 2 次失败),实际 {n_done}"
    k = ml._top_k(30)
    assert k == 10, f"_top_k(30) 应为 10,实际 {k}"
    always = genes[-k:]  # 末 k 个基因在每个成功轮都被选中
    got = freq.loc[always].to_numpy()
    assert (got == 1.0).all(), \
        f"成功轮必选基因频率应为 1.0(分母=成功轮数),实际最大 {got.max()}"
    assert (freq.to_numpy() <= 1.0).all(), "频率超出 [0,1]"
