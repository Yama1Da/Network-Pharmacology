"""core.ml_screening 的测试:六算法共识的信号召回/种子一致性/退化对照/LOOCV
自动切换/方向防御,以及 bootstrap 分母回归。

迁移自 selfcheck_ml.py(17 项中的 6 项)。纯本地合成数据,不联网。

耗时口径:数据规模与 bootstrap 轮数按「足以支撑各断言」的下限取值;
Boruta 迭代与 RF 树数由 conftest._ml_slim_params 统一降档(NPAP_ML_TEST_FULL=1 恢复全参)。
本文件的断言全部是行为语义,与这些规模参数无关——改统计参数后请跑一次全参档。
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
    """合成数据信号:共识 top15 应召回 ≥8 个植入基因,方向正确。

    100 基因中植入 10(shift=2.0):信号强度足够,共识召回与 200 基因口径一致
    (2026-10 瘦身时全参档对照验证过,阈值未动)。
    cv_repeats=1:本测试断言的是共识召回,与 CV 重复精度无关,不付重复成本。"""
    X, y, genes, planted = _make_synthetic(n_genes=100)
    res = ml.run_ml_screening(X, y, bootstrap_rounds=10, cv_repeats=1)
    top15 = res["consensus"].head(15).index.tolist()
    recall = len(set(planted) & set(top15))
    assert recall >= 8, f"信号召回不足: {recall}/15 (期望 ≥8)"
    cs = res["consensus"]["consensus_score"]
    assert cs[planted].mean() > cs[[g for g in genes if g not in planted]].median(), \
        "方向错误:植入基因平均分应高于噪声中位数"


def test_seed_consistency():
    """种子一致性:同种子两次运行逐基因一致。

    确定性只依赖固定 seed 贯穿全部算法,与规模/轮数无关——50 基因 3 轮足矣。
    cv_repeats=1:共识表不经过 CV,重复只添成本。"""
    X1, y, _, _ = _make_synthetic(n_genes=50, seed=42)
    X2, _, _, _ = _make_synthetic(n_genes=50, seed=42)
    r1 = ml.run_ml_screening(X1, y, bootstrap_rounds=3, seed=0, cv_repeats=1)
    r2 = ml.run_ml_screening(X2, y, bootstrap_rounds=3, seed=0, cv_repeats=1)
    cs1 = r1["consensus"]["consensus_score"].round(6).to_dict()
    cs2 = r2["consensus"]["consensus_score"].round(6).to_dict()
    assert cs1 == cs2, "同种子两次运行结果不一致"


def test_label_shuffle():
    """打乱标签:AUC 应回落 ~0.5,植入基因不显著富集于 top15。

    规模与轮数保持原口径(200 基因 / 10 轮):全参档对照实测,该断言对这两者
    敏感——100 基因或 5 轮时,打乱后植入基因随机混入 top15 的方差在生产 Boruta
    参数下会击穿阈值,故不参与瘦身,只吃 conftest 的 Boruta 降档。
    阈值按档位区分(2026-11 实测,种子 99/77):瘦身档 ≤4;全参档生产 Boruta
    的秩方差更大,同样种子下植入基因可随机混入 5 个 → 放宽到 ≤6。两档阈值都
    远高于随机期望(200 基因含 10 植入,top15 期望命中 ~0.75),仅作回归护栏。"""
    import os

    X, y, _, planted = _make_synthetic(seed=99)
    rng = np.random.default_rng(77)
    y_shuf = pd.Series(rng.permutation(y.to_numpy()))
    # cv_repeats=1:AUC 窗口断言(0.2~0.8)对重复精度不敏感,不付 ×5 成本
    res = ml.run_ml_screening(X, y_shuf, bootstrap_rounds=10, cv_repeats=1)
    auc = res["meta"]["cv_auc"]
    assert 0.20 < auc < 0.80, f"打乱后 AUC={auc} 异常(期望 ~0.5)"
    top15 = res["consensus"].head(15).index.tolist()
    n_planted_top = len(set(planted) & set(top15))
    threshold = 6 if os.environ.get("NPAP_ML_TEST_FULL") == "1" else 4
    assert n_planted_top <= threshold, \
        f"打乱后植入基因仍有 {n_planted_top} 个在 top15(阈值 {threshold},不应富集)"


def test_loocv_auto():
    """LOOCV 自动切换:最小类 ≤6 时 cv_mode 应为 loocv。

    4+4 样本仍在切换阈值内(≤6),LOOCV 折数 12→8 进一步省时。
    cv_repeats=1:断言的只是模式切换;LOOCV 本就强制归 1,显式声明免歧义。"""
    X, y, _, _ = _make_synthetic(n_ctrl=4, n_case=4, n_genes=50, seed=42)
    res = ml.run_ml_screening(X, y, bootstrap_rounds=3, cv_repeats=1)
    assert res["meta"]["cv_mode"] == "loocv", f"期望 loocv,实际 {res['meta']['cv_mode']}"


def test_direction_defense():
    """方向防御:传入转置矩阵应报错(行数≠标签数 或 行>列)。"""
    X, y, _, _ = _make_synthetic(n_genes=100, seed=42)
    X_t = X.T  # 故意传反:基因×样本 (100行 × 20列)
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


def test_cv_caliber_meta():
    """CV 口径记录:cv_matrix 提供时 meta.cv_pool="fold",缺省为 "pool"(兼容旧路径)。

    40 基因窄池(模拟全量预筛)+ 200 基因全矩阵:共识跑在窄池上,CV 折内从
    全矩阵重筛(cv_top_n=15,折内 FDR 应把植入基因排进池子),植入信号保持可分。
    cv_repeats=1:断言的是口径字段,不是重复精度。
    bootstrap_rounds=2:共识稳定性不在断言内(cv_auc 来自折内 CV),不付 5 轮成本。"""
    X_full, y, _, _ = _make_synthetic(n_genes=200, seed=42)
    X_pool = X_full.iloc[:, :40]
    res_pool = ml.run_ml_screening(X_pool, y, bootstrap_rounds=2, cv_repeats=1)
    assert res_pool["meta"]["cv_pool"] == "pool"
    assert res_pool["meta"]["cv_top_n"] is None
    res_fold = ml.run_ml_screening(X_pool, y, bootstrap_rounds=2, cv_repeats=1,
                                   cv_matrix=X_full, cv_top_n=15)
    assert res_fold["meta"]["cv_pool"] == "fold"
    assert res_fold["meta"]["cv_top_n"] == 15
    assert res_fold["meta"]["cv_auc"] >= 0.8, \
        f"折内重筛后植入信号应仍可分,AUC={res_fold['meta']['cv_auc']}"


def test_cv_matrix_label_shuffle():
    """折内重筛口径下打乱标签:AUC 应回落 ~0.5(无泄漏口径的退化对照)。

    cv_repeats=1:窗口断言不依赖重复精度(与 test_label_shuffle 同理)。
    bootstrap_rounds=1:AUC 来自折内 CV,与 bootstrap 共识无关,不付
    200 基因 × 10 轮的纯展示成本。"""
    X_full, y, _, _ = _make_synthetic(seed=99)
    rng = np.random.default_rng(77)
    y_shuf = pd.Series(rng.permutation(y.to_numpy()))
    res = ml.run_ml_screening(X_full, y_shuf, bootstrap_rounds=1, cv_repeats=1,
                              cv_matrix=X_full)
    auc = res["meta"]["cv_auc"]
    assert 0.20 < auc < 0.80, f"折内口径打乱后 AUC={auc} 异常(期望 ~0.5)"


def test_cv_repeats_meta():
    """重复 CV:5-fold 下 repeats=3 时 cv_repeats=3、cv_mean/cv_sd 落盘;
    LOOCV(最小类 ≤6)划分确定性,强制归 1 且 sd=None。

    本测试是重复成本唯一该付的地方:规模压到断言下限(40 基因 / 2 轮 bootstrap)。"""
    X, y, _, _ = _make_synthetic(n_ctrl=10, n_case=10, n_genes=40, seed=42)
    res = ml.run_ml_screening(X, y, bootstrap_rounds=2, cv_repeats=3)
    m = res["meta"]
    assert m["cv_mode"].endswith("-fold"), f"10v10 应走分层,实际 {m['cv_mode']}"
    assert m["cv_repeats"] == 3
    assert m["cv_mean"] is not None and m["cv_sd"] is not None
    assert 0.8 <= m["cv_mean"] <= 1.0, f"植入信号下 cv_mean 应高,实际 {m['cv_mean']}"
    Xs, ys, _, _ = _make_synthetic(n_ctrl=4, n_case=4, n_genes=40, seed=42)
    res_loo = ml.run_ml_screening(Xs, ys, bootstrap_rounds=2, cv_repeats=3)
    m2 = res_loo["meta"]
    assert m2["cv_mode"] == "loocv" and m2["cv_repeats"] == 1 and m2["cv_sd"] is None
    assert m2["cv_mean"] == m2["cv_auc"]


def test_permutation_test():
    """置换检验:n_perm=3 跑通,p∈(0,1],obs 与主 CV 同量级;合成信号下 p 应小。

    规模压到断言下限(40 基因 / 6v6,4 次完整 CV 已是本测试的固有成本)。"""
    X, y, _, _ = _make_synthetic(n_ctrl=6, n_case=6, n_genes=40, seed=42)
    perm = ml.permutation_test(X, y, n_perm=3, seed=0)
    assert set(perm) == {"obs_auc", "aucs", "p_value", "n_perm"}
    assert perm["n_perm"] == 3 and len(perm["aucs"]) == 3
    assert 0 < perm["p_value"] <= 1
    assert 0.5 <= perm["obs_auc"] <= 1.0
    # 观测 AUC 应明显高于置换分布均值(合成信号可分,随机标签不可分)
    perm_mean = float(np.mean([a for a in perm["aucs"] if a is not None]))
    assert perm["obs_auc"] > perm_mean, \
        f"obs {perm['obs_auc']} 应高于置换均值 {perm_mean:.3f}"
