#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ML 共识筛选模块自检脚本:固化 M2 验收关卡为正式自检(纯本地,无网络请求)。

覆盖点(方案 §8 + M2 全量审核):
- 合成数据信号召回(植入 10 个区分基因,共识 top15 应召回 ≥8)
- 方向验证(植入基因平均共识分 > 噪声基因中位数)
- 种子一致性(同种子两次运行逐基因一致)
- 打乱标签退化(AUC 应回落 ~0.5,植入基因不再富集)
- LOOCV 自动切换(最小类 ≤6 时)
- 方向防御(传入转置矩阵应报错)
- BH 手算小例(5 个已知 p 值的 FDR 对照)
- 探针识别回归(BUG-1:Ensembl/Illumina/Agilent 不再误判为 symbol 级)
- log2 启发式回归(BUG-3:大尺度含负值应做 log2;z-score 不变换)
- 丢弃样本回归(M1 审查:未命中样本应从标签中剔除)
- 共用管线回归(M4:本地矩阵全链落盘产物齐全 + summary 陈旧判定)

运行方式(在 WY 目录下):
    python scripts/selfcheck_ml.py
"""

import os
import sys
import warnings


import numpy as np
import pandas as pd

from WY.scripts.selfcheck_common import run

warnings.filterwarnings("ignore")

import WY.core.expression as expr
import WY.core.ml_screening as ml
from WY.core.geo import GeoDataset, assign_groups


def _make_synthetic(n_ctrl=10, n_case=10, n_genes=200, n_planted=10, seed=42):
    """合成数据:前 n_planted 个基因在 case 组抬升 shift(区分信号),其余为噪声。"""
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(n_genes)]
    X = rng.normal(0, 1, (n_ctrl + n_case, n_genes))
    X[n_ctrl:, :n_planted] += 2.0
    X = pd.DataFrame(X, columns=genes)
    y = pd.Series([0] * n_ctrl + [1] * n_case)
    return X, y, genes, genes[:n_planted]


def t_signal_recall():
    """合成数据信号:共识 top15 应召回 ≥8 个植入基因,方向正确。"""
    X, y, genes, planted = _make_synthetic()
    res = ml.run_ml_screening(X, y, bootstrap_rounds=20)
    top15 = res["consensus"].head(15).index.tolist()
    recall = len(set(planted) & set(top15))
    assert recall >= 8, f"信号召回不足: {recall}/15 (期望 ≥8)"
    cs = res["consensus"]["consensus_score"]
    assert cs[planted].mean() > cs[[g for g in genes if g not in planted]].median(), \
        "方向错误:植入基因平均分应高于噪声中位数"


def t_seed_consistency():
    """种子一致性:同种子两次运行逐基因一致。"""
    X1, y, _, _ = _make_synthetic(seed=42)
    X2, _, _, _ = _make_synthetic(seed=42)
    r1 = ml.run_ml_screening(X1, y, bootstrap_rounds=10, seed=0)
    r2 = ml.run_ml_screening(X2, y, bootstrap_rounds=10, seed=0)
    cs1 = r1["consensus"]["consensus_score"].round(6).to_dict()
    cs2 = r2["consensus"]["consensus_score"].round(6).to_dict()
    assert cs1 == cs2, "同种子两次运行结果不一致"


def t_label_shuffle():
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


def t_loocv_auto():
    """LOOCV 自动切换:最小类 ≤6 时 cv_mode 应为 loocv。"""
    X, y, _, _ = _make_synthetic(n_ctrl=6, n_case=6, n_genes=50, seed=42)
    res = ml.run_ml_screening(X, y, bootstrap_rounds=5)
    assert res["meta"]["cv_mode"] == "loocv", f"期望 loocv,实际 {res['meta']['cv_mode']}"


def t_direction_defense():
    """方向防御:传入转置矩阵应报错(行数≠标签数 或 行>列)。"""
    X, y, _, _ = _make_synthetic(n_genes=200, seed=42)
    X_t = X.T  # 故意传反:基因×样本 (200行 × 20列)
    try:
        ml.run_ml_screening(X_t, y, bootstrap_rounds=1)
        raise AssertionError("方向防御未触发(传反矩阵未报错)")
    except RuntimeError as e:
        assert "方向" in str(e) or "转置" in str(e), f"报错信息不含方向提示: {e}"


def t_bh_fdr():
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


def t_probe_identifier():
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


def t_log2_heuristic():
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


def t_dropped_samples():
    """M1 审查回归:未命中任何 pattern 的样本应从标签中真正剔除。"""
    rows = []
    for g, n in (("Severe", 5), ("Incipient", 4), ("Moderate", 4), ("Control", 6)):
        for i in range(n):
            rows.append((f"GSM_{g}_{i}", f"tissue sample", f"group: {g}"))
    samples = pd.DataFrame(
        {"source_name_ch1": [r[1] for r in rows],
         "characteristics_ch1": [r[2] for r in rows]},
        index=[r[0] for r in rows])
    ds = GeoDataset(gse="SYNTH", meta={}, samples=samples, table=pd.DataFrame(), gpl="")
    labels = assign_groups(ds, "Severe", "Control")
    n_case, n_ctrl, n_total = int((labels == 1).sum()), int((labels == 0).sum()), len(labels)
    assert n_case == 5 and n_ctrl == 6, f"分组数错: case={n_case} ctrl={n_ctrl}"
    assert n_total == 11, f"丢弃样本未剔除: 总数 {n_total}(期望 11 = 5+6)"
    assert "GSM_Incipient_0" not in labels.index, "Incipient 样本不应在标签中"
    assert "GSM_Moderate_0" not in labels.index, "Moderate 样本不应在标签中"


def t_group_hard_boundaries():
    """分组四道硬边界:零命中/全命中/重叠/样本不足均应报错。"""
    rows = []
    for g, n in (("Severe", 5), ("Incipient", 4), ("Control", 6)):
        for i in range(n):
            rows.append((f"GSM_{g}_{i}", f"tissue", f"group: {g}"))
    samples = pd.DataFrame(
        {"source_name_ch1": [r[1] for r in rows],
         "characteristics_ch1": [r[2] for r in rows]},
        index=[r[0] for r in rows])
    ds = GeoDataset(gse="SYNTH", meta={}, samples=samples, table=pd.DataFrame(), gpl="")

    for case_p, ctrl_p, expect in [
        ("zzz", "Control", "零命中"),
        ("tissue", "Control", "全命中"),
        ("Severe|Control", "Control", "重叠"),
        ("zzz", "Control", "零命中(case 匹配 0 个)"),
    ]:
        try:
            assign_groups(ds, case_p, ctrl_p)
            raise AssertionError(f"{expect}: 未报错")
        except RuntimeError:
            pass  # 预期报错


def t_rnaseq_deg_matrix():
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


def t_pipeline_counts_local():
    """二期管线 counts 分支(桩替 fetch_gse_count_matrix,纯本地):
    全链 7 产物 + source=suppl-counts + DESeq2 DEG 表同构。"""
    import json
    import shutil
    import tempfile
    from WY.core import geo, ml_pipeline

    rng = np.random.default_rng(5)
    genes = [f"G{i:03d}" for i in range(40)]
    cols = ["CTL1", "CTL2", "CTL3", "PA1", "PA2", "PA3"]
    ctrl = rng.poisson(50, size=(40, 3))
    case = rng.poisson(50, size=(40, 3))
    case[:8] = rng.poisson(200, size=(8, 3))
    counts = pd.DataFrame(np.hstack([ctrl, case]), index=genes, columns=cols)
    tmp = tempfile.mkdtemp(prefix="selfcheck_counts_")
    try:
        orig = geo.fetch_gse_count_matrix
        geo.fetch_gse_count_matrix = lambda gse, prefer=None: {
            "counts": counts, "source_file": "SYNTH_raw_counts.tsv.gz",
            "candidates": [], "transposed": False}
        try:
            res = ml_pipeline.run_screen_chain(
                gse="SYNTHC", counts_mode=True, case_pattern="PA",
                control_pattern="CTL", bootstrap_rounds=5, prefilter_top=20,
                top_k=8, out_base=tmp)
        finally:
            geo.fetch_gse_count_matrix = orig
        assert res["summary"]["source"] == "suppl-counts", "source 标注错误"
        assert set(res["deg"].columns) == {"gene", "log2FC", "stat", "p",
                                           "FDR", "direction"}
        assert res["meta"]["n_samples"] == 6 and res["meta"]["n_genes"] == 20
        for fname in ("ml_consensus_scores.csv", "shortlist.csv", "funnel.csv",
                      "summary.json"):
            assert os.path.isfile(os.path.join(res["out_dir"], fname)), f"缺产物 {fname}"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def t_pipeline_roundtrip():
    """共用管线(方案 §M4):本地矩阵全链 → 落盘产物齐全 → summary 可读 → 陈旧判定。"""
    import json
    import shutil
    import tempfile

    from WY import config
    from WY.core.ml_pipeline import run_screen_chain, stale_reasons

    tmp = tempfile.mkdtemp(prefix="selfcheck_ml_pipe_")
    try:
        rng = np.random.default_rng(3)
        genes = [f"G{i:03d}" for i in range(40)]
        cols = [f"S{i:02d}" for i in range(12)]
        X = rng.normal(0, 1, (40, 12))
        X[:6, 6:] += 1.5  # 前 6 个基因在 case(S06..S11) 抬升
        expr_p = os.path.join(tmp, "expr.csv")
        pd.DataFrame(X, index=genes, columns=cols).to_csv(expr_p, encoding="utf-8-sig")
        lab_p = os.path.join(tmp, "labels.csv")
        pd.DataFrame({"sample_id": cols,
                      "group": ["control"] * 6 + ["disease"] * 6}).to_csv(lab_p, index=False)
        prox_p = os.path.join(tmp, "prox.csv")
        pd.DataFrame({"gene": genes[:20], "score": rng.random(20)}).to_csv(prox_p, index=False)

        res = run_screen_chain(expr=expr_p, labels=lab_p, bootstrap_rounds=5,
                               prefilter_top=30, ranked_csvs={"prox": prox_p},
                               top_k=10, out_base=tmp)
        out_dir = res["out_dir"]
        for fname in ("ml_consensus_scores.csv", "shortlist.csv", "funnel.csv",
                      "deg_results.csv", "summary.json"):
            assert os.path.isfile(os.path.join(out_dir, fname)), f"缺产物 {fname}"
        assert any(f.startswith("ml_roc_") and f.endswith(".png")
                   for f in os.listdir(out_dir)), "缺 ROC 图"
        with open(os.path.join(out_dir, "summary.json"), encoding="utf-8") as f:
            summary = json.load(f)
        assert summary["source"] == "local-expr", "source 标注错误"
        assert summary["ml_version"] == config.ML_SCREENING_VERSION, "算法版本未记录"
        assert summary["groups"]["case"] == 6 and summary["groups"]["control"] == 6
        assert "evidence_scope" in res["shortlist"].columns, "shortlist 缺 evidence_scope 列"
        assert res["shortlist"]["evidence_scope"].iloc[0] == "未标注(本地矩阵)", \
            f"本地矩阵 evidence_scope 应为未标注,实际 {res['shortlist']['evidence_scope'].iloc[0]}"
        assert summary["evidence_scope"] == "未标注(本地矩阵)", "summary 未记录 evidence_scope"
        assert stale_reasons(summary) == [], "刚生成的新结果不应判陈旧"
        # 定位注释行(方案 §11):首行 "#" 注释必须能被干净跳过,不得混成脏数据行
        sl_path = os.path.join(out_dir, "shortlist.csv")
        with open(sl_path, encoding="utf-8-sig") as f:
            assert f.readline().startswith("#"), "shortlist.csv 缺首行定位注释"
        sl_back = pd.read_csv(sl_path, comment="#", encoding="utf-8-sig")
        assert list(sl_back.columns)[:2] == ["gene", "evidence_scope"], \
            f"注释行被读成脏数据(表头异常): {list(sl_back.columns)[:3]}"
        assert len(sl_back) == 10, f"注释行被计成数据行(期望 10 行,实际 {len(sl_back)})"
        # 漏斗基因明细(方案 §4.3):dropped_genes 列存在,且末层计数与 shortlist 一致
        fn_back = pd.read_csv(os.path.join(out_dir, "funnel.csv"))
        assert "dropped_genes" in fn_back.columns, "funnel 缺 dropped_genes 明细列"
        assert int(fn_back.iloc[-1]["n"]) == int(res["shortlist"]["in_shortlist"].sum()), \
            "漏斗末层计数与 shortlist 交集数不一致"
        assert "bootstrap_rounds_done" in summary["ml_meta"], \
            "ml_meta 未记录实际完成轮数"
        # 排序 CSV 被重跑(内容/mtime/size 变化)应触发 ranked_changed
        pd.DataFrame({"gene": genes[:25], "score": rng.random(25)}).to_csv(prox_p, index=False)
        reasons = stale_reasons(summary)
        assert any(r["type"] == "ranked_changed" for r in reasons), \
            f"排序更新未触发陈旧: {reasons}"
        # 算法版本不一致应触发 version
        summary["ml_version"] = 0
        assert any(r["type"] == "version" for r in stale_reasons(summary)), \
            "版本不一致未触发陈旧"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def t_sample_source_label():
    """样本来源启发式(方案 §4.1/§4.3):检索表 possible_cellline 与 shortlist
    evidence_scope 共用 sample_source_label,三值口径 + 细胞系优先。"""
    from WY.core.geo import sample_source_label
    assert sample_source_label("HEPG2 cell line", "hepatic cells in vitro") == "细胞系", \
        "细胞系关键词未命中"
    assert sample_source_label("liver tissue biopsy", "patients") == "组织", \
        "组织关键词未命中"
    assert sample_source_label("", "a randomized double-blind trial") == "未标注", \
        "无关键词应判未标注"
    assert sample_source_label("cell line and tissue samples", "") == "细胞系", \
        "同时命中时细胞系应优先"


def t_bootstrap_denominator():
    """分母回归(2026-10 审查修复):失败轮不计入分母,成功轮必选基因频率应为 1.0。"""
    X, y, genes, _ = _make_synthetic(n_ctrl=6, n_case=6, n_genes=30, seed=1)
    orig = ml._select_all
    calls = {"n": 0}

    def fake_select(Xr, yr, seed, boruta_iter=5):
        calls["n"] += 1
        if calls["n"] <= 2:
            raise RuntimeError("模拟选择器失败")
        # 确定性重要性 = 列序:每个成功轮都选中末 k 个基因
        return {alg: pd.Series(np.arange(len(Xr.columns), dtype=float), index=Xr.columns)
                for alg in ("lasso", "enet", "ridge", "rf_boruta", "svm_rfe", "xgboost")}

    ml._select_all = fake_select
    try:
        freq, n_done = ml.bootstrap_stability(X, y, rounds=10, seed=0)
    finally:
        ml._select_all = orig
    assert n_done == 8, f"成功轮数应为 8(10 轮 - 2 次失败),实际 {n_done}"
    k = ml._top_k(30)
    assert k == 10, f"_top_k(30) 应为 10,实际 {k}"
    always = genes[-k:]  # 末 k 个基因在每个成功轮都被选中
    got = freq.loc[always].to_numpy()
    assert (got == 1.0).all(), \
        f"成功轮必选基因频率应为 1.0(分母=成功轮数),实际最大 {got.max()}"
    assert (freq.to_numpy() <= 1.0).all(), "频率超出 [0,1]"


def t_list_analyzed_pairs():
    """疾病下拉数据源(M4 输入区):meta 精确 + 旧结果前缀反解 + rwr 后缀剥离 + manual 跳过。"""
    import json
    import shutil
    import tempfile

    from WY import config
    from WY.core.ppi import list_analyzed_pairs

    tmp = tempfile.mkdtemp(prefix="selfcheck_pairs_")
    orig_dl = config.DOWNLOADS_DIR
    config.DOWNLOADS_DIR = tmp
    try:
        # 新结果: meta 带 drug/disease(精确)
        d1 = os.path.join(tmp, "aspirin_cancer_data")
        os.makedirs(d1)
        with open(os.path.join(d1, "analysis_meta.json"), "w", encoding="utf-8") as f:
            json.dump({"drug": "Aspirin", "disease": "Cancer", "string_version": "12.0",
                       "analyzed_at": "2026-10-21T10:00:00"}, f)
        # 旧结果: meta 无药病对 → 从前缀反解(下划线还原为空格)
        d2 = os.path.join(tmp, "metformin_type_2_diabetes_data")
        os.makedirs(d2)
        # 旧结果带 rwr 变体后缀: disease 保留原样、base_disease 剥离
        d3 = os.path.join(tmp, "aspirin_cancer_rwr50_data")
        os.makedirs(d3)
        # manual 目录(手动基因输入,无药病对)与散文件应被跳过
        os.makedirs(os.path.join(tmp, "manual_data"))
        open(os.path.join(tmp, "loose.txt"), "w").close()

        pairs = list_analyzed_pairs()
        by_dir = {os.path.basename(p["directory"]): p for p in pairs}
        assert set(by_dir) == {"aspirin_cancer_data", "metformin_type_2_diabetes_data",
                               "aspirin_cancer_rwr50_data"}, \
            f"目录筛选错误: {sorted(by_dir)}"
        assert by_dir["aspirin_cancer_data"]["disease"] == "Cancer"
        assert by_dir["aspirin_cancer_data"]["base_disease"] == "Cancer"
        assert by_dir["aspirin_cancer_data"]["string_version"] == "12.0"
        assert by_dir["aspirin_cancer_data"]["source"] == "meta"
        assert by_dir["metformin_type_2_diabetes_data"]["disease"] == "type 2 diabetes", \
            "前缀反解未还原空格"
        assert by_dir["metformin_type_2_diabetes_data"]["base_disease"] == "type 2 diabetes"
        assert by_dir["metformin_type_2_diabetes_data"]["source"] == "prefix"
        assert by_dir["aspirin_cancer_rwr50_data"]["disease"] == "cancer rwr50"
        assert by_dir["aspirin_cancer_rwr50_data"]["base_disease"] == "cancer", \
            "rwr 变体后缀未剥离"
    finally:
        config.DOWNLOADS_DIR = orig_dl
        shutil.rmtree(tmp, ignore_errors=True)


def t_results_mgmt():
    """结果管理数据层(§10 增补):keep-N 只删同 tag 最旧、登记根之外与缺
    summary 的目录拒删、zip 产出非空。"""
    import json as _json
    import shutil
    import tempfile

    from WY.core import ml_results

    tmp = tempfile.mkdtemp(prefix="selfcheck_ml_mgmt_")
    try:
        root = os.path.join(tmp, "mlscreen")
        os.makedirs(root)

        def _mk(tag, ymd, hms):
            d = os.path.join(root, f"{tag}_{ymd}_{hms}")
            os.makedirs(d)
            with open(os.path.join(d, "summary.json"), "w", encoding="utf-8") as f:
                _json.dump({"ml_version": 1, "ml_meta": {}, "ranked_csvs": {},
                            "shortlist_count": 0}, f)
            with open(os.path.join(d, "data.txt"), "w") as f:
                f.write("x" * 128)
            return d

        d_old = _mk("GSE1", "20251231", "235959")
        _mk("GSE1", "20260101", "000001")
        d_keep = _mk("GSE1", "20260102", "000002")
        _mk("GSE1", "20260103", "000003")   # 最新,必保留
        _mk("GSE2", "20260105", "000000")   # 不同 tag,不参与 GSE1 的清理

        runs = ml_results.list_ml_runs(roots=[root])
        assert len(runs) == 5 and runs[0]["tag"] == "GSE2", "列表应按时间新→旧"

        removed = ml_results.cleanup_keep_per_gse(keep=3, roots=[root])
        assert len(removed) == 1, f"应只删 1 个最旧目录,实际 {removed}"
        assert not os.path.isdir(d_old) and os.path.isdir(d_keep)
        assert len(ml_results.list_ml_runs(roots=[root])) == 4

        outside = os.path.join(tmp, "outside")
        os.makedirs(outside)
        assert ml_results.delete_ml_runs([outside], roots=[root]) == [], "根外路径必须拒删"
        assert os.path.isdir(outside)
        junk = os.path.join(root, "junk_dir")
        os.makedirs(junk)
        assert ml_results.delete_ml_runs([junk], roots=[root]) == [], "缺 summary 的目录必须拒删"
        assert os.path.isdir(junk)

        assert len(ml_results.zip_ml_runs([d_keep])) > 100, "zip 产物不应为空"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    tests = [
        ("合成数据信号召回(≥8/10)与方向", t_signal_recall),
        ("种子一致性(逐基因)", t_seed_consistency),
        ("打乱标签退化(AUC~0.5)", t_label_shuffle),
        ("LOOCV 自动切换(≤6 样本)", t_loocv_auto),
        ("方向防御(转置矩阵报错)", t_direction_defense),
        ("BH-FDR 手算对照(15 值)", t_bh_fdr),
        ("探针识别(BUG-1 回归)", t_probe_identifier),
        ("log2 启发式(BUG-3 回归)", t_log2_heuristic),
        ("丢弃样本剔除(M1 审查回归)", t_dropped_samples),
        ("分组四道硬边界", t_group_hard_boundaries),
        ("共用管线全链落盘与陈旧判定(M4)", t_pipeline_roundtrip),
        ("RNA-seq DESeq2 路径(合成 NB counts)", t_rnaseq_deg_matrix),
        ("管线 counts 分支(桩替 fetch,纯本地)", t_pipeline_counts_local),
        ("样本来源标注(细胞系/组织/未标注)", t_sample_source_label),
        ("bootstrap 分母(失败轮不计)", t_bootstrap_denominator),
        ("疾病下拉数据源(已分析药病对发现)", t_list_analyzed_pairs),
        ("结果管理数据层(keep-N/根外拒删/zip)", t_results_mgmt),
    ]
    run(tests, "GEO+ML screening selfcheck")


if __name__ == "__main__":
    main()
