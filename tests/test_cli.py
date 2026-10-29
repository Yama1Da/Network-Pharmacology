"""WY.cli 的测试:参数解析与落盘纪律(不打网络,core 层一律打桩)。

阶段三第 4 批新写(仿原方案示例 D,import 按现行 WY. 绝对导入约定)。
打桩打在"使用者的命名空间"里(cli 按名字导入的 core 函数,
替换 cli 命名空间里的那个名字)。
"""
import argparse
import json
import os

import pandas as pd
import pytest

from WY import cli


# ---------- 通用小工具 ----------


def test_slug_basic_and_truncation():
    """目录名:各段 safe_name 后双下划线连接;超长截断并带短哈希防撞。"""
    assert cli._slug("metformin", "Type 2 diabetes") == "metformin__Type_2_diabetes"
    long = cli._slug("x" * 300, "y" * 300)
    assert len(long) <= cli._SLUG_MAX_LEN and "~" in long


def test_resolve_out_dir(tmp_path, monkeypatch):
    """--out 三级解析:显式参数优先;无参数回退 ./results;--out 撞上已存在文件 → None。"""
    out = cli._resolve_out_dir(str(tmp_path / "made"))
    assert os.path.isdir(out)
    monkeypatch.setattr(config_CLI(), "CLI_RESULTS_DIR", "")
    fallback = cli._resolve_out_dir(None)  # 兜底 ./results(仓库根,已存在)
    assert os.path.isdir(fallback)
    existing_file = tmp_path / "occupied.txt"
    existing_file.write_text("x", encoding="utf-8")
    assert cli._resolve_out_dir(str(existing_file)) is None


def config_CLI():
    from WY import config
    return config


def test_save_csv_discipline(tmp_path, monkeypatch):
    """落盘纪律:utf-8-sig(Excel 不乱码)、空表不落盘、被占用(PermissionError)不崩溃。"""
    path = str(tmp_path / "t.csv")
    assert cli._save_csv(pd.DataFrame({"g": ["TP53"]}), path) == path
    assert "TP53" in open(path, encoding="utf-8-sig").read()
    assert cli._save_csv(pd.DataFrame(), path) is None
    assert cli._save_csv(None, path) is None
    monkeypatch.setattr(pd.DataFrame, "to_csv",
                        lambda *a, **k: (_ for _ in ()).throw(PermissionError("locked")))
    assert cli._save_csv(pd.DataFrame({"g": [1]}), path) is None


def test_read_genes_formats(tmp_path):
    """基因清单三种来源:字面列表(逗号/分号/空白分隔,去重保序)、CSV 取基因列、
    文本文件(# 注释忽略)。目录 → 空列表。"""
    assert cli._read_genes("TP53, EGFR;EGFR KRAS") == ["TP53", "EGFR", "KRAS"]
    csv = tmp_path / "common_targets.csv"
    pd.DataFrame({"gene": ["TP53", "EGFR"]}).to_csv(csv, index=False)
    assert cli._read_genes(str(csv)) == ["TP53", "EGFR"]
    txt = tmp_path / "seeds.txt"
    txt.write_text("TP53  # 种子\nAKT1\n", encoding="utf-8")
    assert cli._read_genes(str(txt)) == ["TP53", "AKT1"]
    assert cli._read_genes(str(tmp_path)) == []


def test_parse_sources_whitelist():
    """--sources:大小写不敏感、去重保序;白名单外的值报错返回 None(拼错不静默走空)。"""
    valid = cli._parse_sources("kegg,go:bp,KEGG")
    assert valid == ("KEGG", "GO:BP")
    assert cli._parse_sources("bogus_source") is None


def test_pvalue_validation():
    """--pvalue:必须在 (0,1] 内,越界/非数值按用法错误拦下(在进网络请求之前)。"""
    assert cli._pvalue("0.05") == 0.05
    assert cli._pvalue("1") == 1.0
    for bad in ("1.5", "0", "-1", "abc"):
        with pytest.raises(argparse.ArgumentTypeError):
            cli._pvalue(bad)


# ---------- 子命令:targets ----------


def test_targets_drug_only_writes_csv(tmp_path, monkeypatch):
    """npap targets(仅药物):打桩数据源整合 → 退出码 0 → 靶点表落在 --out 目录。"""
    monkeypatch.setattr(cli, "get_drug_targets_merged",
                        lambda drug: pd.DataFrame({"approvedSymbol": ["TP53", "EGFR"],
                                                   "score_dgidb": [0.9, 0.8]}))
    rc = cli.main(["targets", "--drug", "metformin", "--out", str(tmp_path)])
    assert rc == 0
    out_csv = tmp_path / "drug_targets.csv"
    assert out_csv.exists()
    assert "EGFR" in out_csv.read_text(encoding="utf-8-sig")


def test_targets_empty_side_is_failure(tmp_path, monkeypatch):
    """明确请求的一侧查无靶点按失败计(防拼写错误静默通过)。"""
    monkeypatch.setattr(cli, "get_drug_targets_merged", lambda drug: pd.DataFrame())
    rc = cli.main(["targets", "--drug", "不存在", "--out", str(tmp_path)])
    assert rc == 1


def test_targets_with_disease(tmp_path, monkeypatch):
    """带疾病:走并发编排入口,两张表都落盘;失败源列表非空时按失败退出。"""
    drug_df = pd.DataFrame({"approvedSymbol": ["TP53"]})
    disease_df = pd.DataFrame({"approvedSymbol": ["EGFR"]})
    monkeypatch.setattr(cli, "get_targets_concurrent",
                        lambda d, dis, progress_cb=None: (drug_df, disease_df, []))
    rc = cli.main(["targets", "--drug", "d1", "--disease", "d2", "--out", str(tmp_path)])
    assert rc == 0
    assert (tmp_path / "drug_targets.csv").exists()
    assert (tmp_path / "disease_targets.csv").exists()
    monkeypatch.setattr(cli, "get_targets_concurrent",
                        lambda d, dis, progress_cb=None: (drug_df, disease_df, ["dgidb"]))
    rc2 = cli.main(["targets", "--drug", "d1", "--disease", "d2", "--out",
                    str(tmp_path / "again")])
    assert rc2 == 1


# ---------- 子命令:ppi / enrich ----------


def _stub_ppi_chain(monkeypatch, edges):
    """PPI 链打桩:抓边打桩,图构建/指标/MCC 用真实小图跑(纯本地核心)。"""
    monkeypatch.setattr(cli, "fetch_string_ppi", lambda genes, **k: edges)
    import networkx as nx

    def fake_build(edges_df, genes):
        g = nx.Graph()
        g.add_edges_from([("A", "B"), ("B", "C")])
        return g

    monkeypatch.setattr(cli, "build_ppi_graph", fake_build)
    monkeypatch.setattr(cli, "compute_ppi_metrics",
                        lambda g: pd.DataFrame({"gene": ["A", "B", "C"],
                                                "degree": [1, 2, 1],
                                                "betweenness": [0.0, 0.6, 0.0],
                                                "closeness": [0.4, 0.5, 0.4]}))
    monkeypatch.setattr(cli, "compute_mcc_scores",
                        lambda g: {"A": 1.0, "B": 2.0, "C": 1.0})


def test_ppi_command_happy_path(tmp_path, monkeypatch):
    """npap ppi:三张结果表落盘,退出码 0。"""
    _stub_ppi_chain(monkeypatch, pd.DataFrame({"protein1": ["A"], "protein2": ["B"]}))
    rc = cli.main(["ppi", "--genes", "A,B,C", "--out", str(tmp_path)])
    assert rc == 0
    for name in ("ppi_edges.csv", "ppi_metrics.csv", "core_targets.csv"):
        assert (tmp_path / name).exists(), name


def test_ppi_no_string_edges_is_failure(tmp_path, monkeypatch):
    """STRING 没返回任何互作(基因不被识别/无达标边):明确失败退出。"""
    monkeypatch.setattr(cli, "fetch_string_ppi", lambda genes, **k: None)
    rc = cli.main(["ppi", "--genes", "ZZZ1,ZZZ2", "--out", str(tmp_path)])
    assert rc == 1


def test_ppi_no_genes_parsed_is_failure(tmp_path):
    """--genes 解析不出任何基因:进网络请求前就失败退出。"""
    rc = cli.main(["ppi", "--genes", "   ", "--out", str(tmp_path)])
    assert rc == 1


def test_enrich_command_semantics(tmp_path, monkeypatch):
    """npap enrich 三态:有显著项 → 0 + 落盘;真无显著项(可信)→ 0;服务不可信 → 1。"""
    terms = pd.DataFrame([{"source": "KEGG", "term_name": "x", "p_value": 0.01}])
    monkeypatch.setattr(cli, "run_enrichment",
                        lambda genes, **k: (terms, "gProfiler"))
    rc = cli.main(["enrich", "--genes", "TP53", "--out", str(tmp_path)])
    assert rc == 0 and (tmp_path / "enrichment.csv").exists()

    def no_sig(genes, **k):
        status = k["status"]
        status.update(ok=True, reason="no_significant")
        return pd.DataFrame(), "g:Profiler"

    monkeypatch.setattr(cli, "run_enrichment", no_sig)
    rc2 = cli.main(["enrich", "--genes", "TP53", "--out", str(tmp_path / "b")])
    assert rc2 == 0  # 真无显著项属正常

    def down(genes, **k):
        status = k["status"]
        status.update(ok=False, reason="both_sources_down")
        return pd.DataFrame(), None

    monkeypatch.setattr(cli, "run_enrichment", down)
    rc3 = cli.main(["enrich", "--genes", "TP53", "--out", str(tmp_path / "c")])
    assert rc3 == 1  # 服务不可用必须非零退出暴露给脚本


def test_enrich_bad_sources_rejected(tmp_path):
    """--sources 拼错:用法层直接失败,不发任何请求。"""
    rc = cli.main(["enrich", "--genes", "TP53", "--sources", "bogus",
                   "--out", str(tmp_path)])
    assert rc == 1


# ---------- 子命令:proximity / rwr(依赖背景互作组) ----------


def test_proximity_and_rwr_skip_when_interactome_missing(tmp_path, monkeypatch):
    """背景互作组未构建:警告并跳过,CLI 绝不自动触发数 GB 下载。"""
    monkeypatch.setattr(cli.interactome, "is_built", lambda *a, **k: False)
    rc1 = cli.main(["proximity", "--drug", "d", "--disease", "x",
                    "--out", str(tmp_path)])
    rc2 = cli.main(["rwr", "--seeds", "TP53", "--out", str(tmp_path)])
    assert rc1 == 1 and rc2 == 1


# ---------- 子命令:analyze / batch ----------


def _write_summary(pair_dir, summary):
    """与真实 run_pair 同职责:把 summary.json 落盘到结果目录。"""
    os.makedirs(pair_dir, exist_ok=True)
    with open(os.path.join(pair_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f)
    return summary


def _fake_summary(ok=True, failed_step=None):
    steps = {"targets": {"status": "ok", "detail": "x"},
             "ppi": {"status": "ok", "detail": "x"},
             "enrichment": {"status": "ok", "detail": "x"},
             "proximity": {"status": "ok", "detail": "x"},
             "rwr": {"status": "ok", "detail": "x"}}
    if failed_step:
        steps[failed_step] = {"status": "failed", "detail": "boom"}
    return {"ok": ok, "steps": steps, "elapsed_seconds": 1.0}


def test_analyze_command_writes_summary(tmp_path, monkeypatch):
    """npap analyze:run_pair 打桩 → 每对一个子目录 + summary.json,全 ok 退出码 0。"""
    monkeypatch.setattr(cli, "run_pair",
                        lambda drug, disease, pair_dir:
                        _write_summary(pair_dir, _fake_summary(True)))
    rc = cli.main(["analyze", "--drug", "metformin", "--disease", "T2D",
                   "--out", str(tmp_path)])
    assert rc == 0
    pair_dir = tmp_path / "metformin__T2D"
    summary_file = pair_dir / "summary.json"
    assert summary_file.exists()
    assert json.loads(summary_file.read_text(encoding="utf-8"))["ok"] is True


def test_analyze_failed_steps_nonzero_exit(tmp_path, monkeypatch):
    """任一步骤失败:summary 如实记录,退出码 1(某步失败不阻断其余步骤)。"""
    monkeypatch.setattr(cli, "run_pair",
                        lambda drug, disease, pair_dir: _fake_summary(
                            ok=False, failed_step="ppi"))
    rc = cli.main(["analyze", "--drug", "d", "--disease", "x", "--out", str(tmp_path)])
    assert rc == 1


def test_analyze_keyboard_interrupt_exit_130(tmp_path, monkeypatch):
    """Ctrl+C:退出码 130(batch 复用同约定),不按普通崩溃处理。"""
    def interrupted(drug, disease, pair_dir):
        raise KeyboardInterrupt
    monkeypatch.setattr(cli, "run_pair", interrupted)
    rc = cli.main(["analyze", "--drug", "d", "--disease", "x", "--out", str(tmp_path)])
    assert rc == 130


def test_read_pairs_validation(tmp_path):
    """pairs.csv 解析:表头大小写不敏感、多余列忽略、空缺行记非法(行号从 2 计)、缺列抛错。"""
    csv = tmp_path / "pairs.csv"
    pd.DataFrame({"No.": [1, 2, 3], "Drug": ["a", "b", ""],
                  "Disease": ["x", "", "z"]}).to_csv(csv, index=False)
    pairs, bad = cli._read_pairs(str(csv))
    assert pairs == [("a", "x")]
    assert bad == [(3, "empty drug or disease"), (4, "empty drug or disease")]
    bad_csv = tmp_path / "bad.csv"
    pd.DataFrame({"foo": [1]}).to_csv(bad_csv, index=False)
    with pytest.raises(ValueError, match="drug"):
        cli._read_pairs(str(bad_csv))


def test_batch_happy_and_summary(tmp_path, monkeypatch):
    """npap batch:逐对跑完出 batch_summary.csv,全成功退出码 0。"""
    monkeypatch.setattr(cli, "run_pair",
                        lambda drug, disease, pair_dir:
                        _write_summary(pair_dir, _fake_summary(True)))
    csv = tmp_path / "pairs.csv"
    pd.DataFrame({"drug": ["a", "b"], "disease": ["x", "y"]}).to_csv(csv, index=False)
    rc = cli.main(["batch", "--pairs", str(csv), "--out", str(tmp_path)])
    assert rc == 0
    summary = pd.read_csv(tmp_path / "batch_summary.csv", encoding="utf-8-sig")
    assert len(summary) == 2 and set(summary["status"]) == {"ok"}


def test_batch_partial_failure_and_invalid_rows(tmp_path, monkeypatch):
    """批量:单对失败不中断;非法行按失败入汇总对账;有失败退出码 1。"""
    calls = {"n": 0}

    def flaky(drug, disease, pair_dir):
        calls["n"] += 1
        if drug == "b":
            raise RuntimeError("crashed")
        return _fake_summary(True)

    monkeypatch.setattr(cli, "run_pair", flaky)
    csv = tmp_path / "pairs.csv"
    pd.DataFrame({"drug": ["a", "b", ""], "disease": ["x", "y", "z"]}).to_csv(
        csv, index=False)
    rc = cli.main(["batch", "--pairs", str(csv), "--out", str(tmp_path)])
    assert rc == 1 and calls["n"] == 2
    summary = pd.read_csv(tmp_path / "batch_summary.csv", encoding="utf-8-sig")
    statuses = list(summary["status"])
    assert statuses.count("ok") == 1 and statuses.count("failed") == 2
    assert "invalid_row" in set(summary["failed_steps"])


def test_batch_resume_skips_successful(tmp_path, monkeypatch):
    """--resume:已有成功 summary.json 的对直接跳过(状态沿用落盘记录),不重跑。"""
    calls = {"n": 0}

    def counting(drug, disease, pair_dir):
        calls["n"] += 1
        return _write_summary(pair_dir, _fake_summary(True))

    monkeypatch.setattr(cli, "run_pair", counting)
    csv = tmp_path / "pairs.csv"
    pd.DataFrame({"drug": ["a"], "disease": ["x"]}).to_csv(csv, index=False)
    out1 = tmp_path / "run1"
    assert cli.main(["batch", "--pairs", str(csv), "--out", str(out1)]) == 0
    assert calls["n"] == 1
    rc = cli.main(["batch", "--pairs", str(csv), "--out", str(out1), "--resume"])
    assert rc == 0 and calls["n"] == 1  # 第二轮被 resume 跳过
