"""共用管线(core.ml_pipeline)、结果管理(core.ml_results)与药病对发现
(core.ppi.list_analyzed_pairs)的测试。

迁移自 selfcheck_ml.py(17 项中的 4 项)。落盘产物由 conftest 沙箱与
tmp_path 接管,不触碰真实数据目录。
"""
import json
import os

import numpy as np
import pandas as pd

from WY import config
from WY.core.ml_pipeline import run_screen_chain, stale_reasons


def test_pipeline_counts_local(monkeypatch, tmp_path):
    """二期管线 counts 分支(桩替 fetch_gse_count_matrix,纯本地):
    全链 7 产物 + source=suppl-counts + DESeq2 DEG 表同构。

    轮数 5→3(2026-10 瘦身):分支持的全链产物结构不受轮数影响。"""
    from WY.core import geo, ml_pipeline

    rng = np.random.default_rng(5)
    genes = [f"G{i:03d}" for i in range(40)]
    cols = ["CTL1", "CTL2", "CTL3", "PA1", "PA2", "PA3"]
    ctrl = rng.poisson(50, size=(40, 3))
    case = rng.poisson(50, size=(40, 3))
    case[:8] = rng.poisson(200, size=(8, 3))
    counts = pd.DataFrame(np.hstack([ctrl, case]), index=genes, columns=cols)
    monkeypatch.setattr(geo, "fetch_gse_count_matrix",
                        lambda gse, prefer=None: {
                            "counts": counts, "source_file": "SYNTH_raw_counts.tsv.gz",
                            "candidates": [], "transposed": False})
    res = ml_pipeline.run_screen_chain(
        gse="SYNTHC", counts_mode=True, case_pattern="PA",
        control_pattern="CTL", bootstrap_rounds=3, prefilter_top=20,
        top_k=8, out_base=str(tmp_path))
    assert res["summary"]["source"] == "suppl-counts", "source 标注错误"
    assert set(res["deg"].columns) == {"gene", "log2FC", "stat", "p",
                                       "FDR", "direction"}
    assert res["meta"]["n_samples"] == 6 and res["meta"]["n_genes"] == 20
    for fname in ("ml_consensus_scores.csv", "shortlist.csv", "funnel.csv",
                  "summary.json"):
        assert os.path.isfile(os.path.join(res["out_dir"], fname)), f"缺产物 {fname}"


def test_pipeline_roundtrip(tmp_path):
    """共用管线(方案 §M4):本地矩阵全链 → 落盘产物齐全 → summary 可读 → 陈旧判定。

    规模按断言所需下限取值(2026-10 瘦身):30 基因 × 4+4 样本(仍走 LOOCV 折内
    全选择路径),产物结构/注释行/漏斗/陈旧判定的断言与规模无关。"""
    tmp = str(tmp_path)
    rng = np.random.default_rng(3)
    genes = [f"G{i:03d}" for i in range(30)]
    cols = [f"S{i:02d}" for i in range(8)]
    X = rng.normal(0, 1, (30, 8))
    X[:6, 4:] += 1.5  # 前 6 个基因在 case(S04..S07) 抬升
    expr_p = os.path.join(tmp, "expr.csv")
    pd.DataFrame(X, index=genes, columns=cols).to_csv(expr_p, encoding="utf-8-sig")
    lab_p = os.path.join(tmp, "labels.csv")
    pd.DataFrame({"sample_id": cols,
                  "group": ["control"] * 4 + ["disease"] * 4}).to_csv(lab_p, index=False)
    prox_p = os.path.join(tmp, "prox.csv")
    pd.DataFrame({"gene": genes[:15], "score": rng.random(15)}).to_csv(prox_p, index=False)

    res = run_screen_chain(expr=expr_p, labels=lab_p, bootstrap_rounds=3,
                           prefilter_top=20, ranked_csvs={"prox": prox_p},
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
    assert summary["groups"]["case"] == 4 and summary["groups"]["control"] == 4
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


def test_list_analyzed_pairs():
    """疾病下拉数据源(M4 输入区):meta 精确 + 旧结果前缀反解 + rwr 后缀剥离 + manual 跳过。
    conftest 沙箱已把 DOWNLOADS_DIR 指到本测试的临时目录,直接在里面造结果目录。"""
    from WY.core.ppi import list_analyzed_pairs

    tmp = config.DOWNLOADS_DIR
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


def test_results_mgmt(tmp_path):
    """结果管理数据层(§10 增补):keep-N 只删同 tag 最旧、登记根之外与缺
    summary 的目录拒删、zip 产出非空。"""
    from WY.core import ml_results

    # zip 的双重校验(P1-B)针对登记结果根(_roots()),自造根须对齐:
    # 沙箱里 default_results_root() = tmp_path/results,登记根即 results/mlscreen
    root = os.path.join(str(tmp_path), "results", "mlscreen")
    os.makedirs(root, exist_ok=True)

    def _mk(tag, ymd, hms):
        d = os.path.join(root, f"{tag}_{ymd}_{hms}")
        os.makedirs(d)
        with open(os.path.join(d, "summary.json"), "w", encoding="utf-8") as f:
            json.dump({"ml_version": 1, "ml_meta": {}, "ranked_csvs": {},
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

    outside = os.path.join(str(tmp_path), "outside")
    os.makedirs(outside)
    assert ml_results.delete_ml_runs([outside], roots=[root]) == [], "根外路径必须拒删"
    assert os.path.isdir(outside)
    junk = os.path.join(root, "junk_dir")
    os.makedirs(junk)
    assert ml_results.delete_ml_runs([junk], roots=[root]) == [], "缺 summary 的目录必须拒删"
    assert os.path.isdir(junk)

    assert len(ml_results.zip_ml_runs([d_keep])) > 100, "zip 产物不应为空"


# ---------------- load_local_expr 健壮化(2026-11 GSE253502 实跑教训固化) ----------------

def _write_local(tmp_path, expr_df, labels_rows):
    expr_p = tmp_path / "expr.csv"
    lab_p = tmp_path / "labels.csv"
    expr_df.to_csv(expr_p, encoding="utf-8-sig")
    pd.DataFrame(labels_rows, columns=["sample_id", "group"]).to_csv(
        lab_p, index=False, encoding="utf-8-sig")
    return str(expr_p), str(lab_p)


def test_load_local_expr_drops_comment_cols(tmp_path):
    """作者分析列(FDR/PValue/log2FC 等精确匹配)被剔除并报告,不污染样本轴。"""
    from WY.core.ml_pipeline import load_local_expr

    rng = np.random.default_rng(7)
    genes = [f"G{i:03d}" for i in range(20)]
    cols = ["S1", "S2", "S3", "S4", "S5", "S6"]
    df = pd.DataFrame(rng.normal(0, 1, (20, 6)), index=genes, columns=cols)
    # 伪装成样本列的作者分析列:数值型,不会被 to_numeric 丢弃
    df["FDR"] = 0.01
    df["log2FC"] = 1.5
    df["PValue"] = 0.03
    labels = [(c, "control" if j < 3 else "case") for j, c in enumerate(cols)]
    expr_p, lab_p = _write_local(tmp_path, df, labels)
    notes: list[str] = []

    matrix, y, _tag, info = load_local_expr(
        expr_p, lab_p, case_label="case", progress=lambda s, m: notes.append(m))
    assert list(matrix.columns) == cols, "注释列应被剔除,只剩真实样本列"
    assert info["n_case"] == 3 and info["n_control"] == 3
    assert any("剔除疑似作者分析列" in n for n in notes), notes


def test_load_local_expr_duplicate_symbols_mean(tmp_path):
    """重复基因符号按组均值聚合(替代旧 keep=first 的任意性)。"""
    from WY.core.ml_pipeline import load_local_expr

    cols = ["C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"]  # 4v4,过每组 ≥3 下限
    expr = pd.DataFrame([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
                         [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0],
                         [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]],
                        index=["T1", "T1", "T2"], columns=cols)
    labels = [(c, "control" if c.startswith("C") else "case") for c in cols]
    expr_p, lab_p = _write_local(tmp_path, expr, labels)

    matrix, _y, _tag, _info = load_local_expr(expr_p, lab_p, case_label="case")
    assert list(matrix.index) == ["T1", "T2"], "重复符号应聚合为一行"
    # T1 两行均值 = [(1+10)/2, (2+20)/2, ...]
    assert matrix.loc["T1"].tolist() == [5.5, 11.0, 16.5, 22.0, 27.5, 33.0, 38.5, 44.0]


def test_load_local_expr_ensembl_error_guides_map_ids(tmp_path):
    """Ensembl 级矩阵:map_ids 关 → 响亮报错并给出 --map-ids 引导。"""
    from WY.core.ml_pipeline import load_local_expr

    cols = ["C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4"]
    expr = pd.DataFrame(np.random.default_rng(1).normal(0, 1, (20, 8)),
                        index=[f"ENSG{i:08d}" for i in range(20)], columns=cols)
    labels = [(c, "control" if c.startswith("C") else "case") for c in cols]
    expr_p, lab_p = _write_local(tmp_path, expr, labels)

    try:
        load_local_expr(expr_p, lab_p, case_label="case")
        raise AssertionError("Ensembl 级矩阵未开 map_ids 应报错")
    except RuntimeError as e:
        assert "--map-ids" in str(e) and "mygene" in str(e)


# ---------- counts 模式符号映射重复列回归(P0-C,2026-11-08 审查) ----------


def test_pipeline_counts_duplicate_symbol_mapping(monkeypatch, tmp_path):
    """回归测试:counts 模式两个 Ensembl ID 映射到同一符号时整链不得崩溃。

    缺陷:counts 分支 rename 后只对 deg 去重,X 未去重——mygene 多对一映射
    (常见)使 X 出现重复列,ML 阶段整链崩溃(审查记录旧栈:bootstrap 的
    X.columns.get_loc(g) 对重复列返回布尔数组抛 ValueError;当前依赖版本由
    sklearn/narwhals 列名校验先拦截,抛 DuplicateError——同一根因),且
    DESeq2 已算完,结果全盘作废。
    复现要点:两行不同 ENS ID 强信号 counts,桩替 map_ids_to_symbol 使二者
    映射到同一符号;修复前 run_screen_chain 在 ML 阶段崩溃,修复后整链跑通,
    重复符号按均值聚合(与 load_local_expr 本地矩阵既有口径一致)。
    """
    from WY.core import geo, ml_pipeline
    from WY.core import expression as expr_mod

    dup_a, dup_b = "ENSG000002DUPA", "ENSG000002DUPB"
    rng = np.random.default_rng(7)
    genes = [dup_a, dup_b] + [f"ENSG{i:08d}" for i in range(38)]
    cols = ["CTL1", "CTL2", "CTL3", "PA1", "PA2", "PA3"]
    ctrl = rng.poisson(50, size=(40, 3))
    case = rng.poisson(50, size=(40, 3))
    case[:20] = rng.poisson(200, size=(20, 3))  # 前 20 行强信号,必进预过滤池
    case[:2] = rng.poisson(400, size=(2, 3))    # 重复对最响,保证二者同入池
    counts = pd.DataFrame(np.hstack([ctrl, case]), index=genes, columns=cols)
    monkeypatch.setattr(geo, "fetch_gse_count_matrix",
                        lambda gse, prefer=None: {
                            "counts": counts, "source_file": "SYNTH_raw_counts.tsv.gz",
                            "candidates": [], "transposed": False})
    # 桩替 mygene 批量映射:两个 ENS ID 映射到同一符号,其余未命中(None 保留原 ID)
    monkeypatch.setattr(expr_mod, "map_ids_to_symbol",
                        lambda ids: {g: ("TP53" if g in (dup_a, dup_b) else None)
                                     for g in ids})

    res = ml_pipeline.run_screen_chain(
        gse="SYNTHDUP", counts_mode=True, case_pattern="PA",
        control_pattern="CTL", bootstrap_rounds=3, prefilter_top=20,
        top_k=8, out_base=str(tmp_path))

    cons = res["consensus"]
    assert not cons.index.duplicated().any(), "共识表基因必须唯一(X 重复列已聚合)"
    assert (cons.index == "TP53").sum() == 1, "重复符号应聚合为一条(均值口径)"
    assert (res["deg"]["gene"] == "TP53").sum() == 1, "deg 重复符号 keep=first"
    for fname in ("ml_consensus_scores.csv", "shortlist.csv", "summary.json"):
        assert os.path.isfile(os.path.join(res["out_dir"], fname)), f"缺产物 {fname}"
