"""WY.cli 的成功路径补测:proximity/rwr 全绿路径、run_pair 真骨架端到端、
batch 失败后续跑、serve 端口预检。core 重活全部打桩。

深水区冲线批 II。
"""
import json
import os

import pandas as pd
import pytest

from WY import config
from WY import cli


def _prox_result(ok=True, **over):
    base = {"ok": ok, "d_c": 1.5, "z": -2.3, "p_empirical": 0.04, "n_excl": 5,
            "d_c_excl": 1.6, "z_excl": -2.0, "p_excl": 0.05, "cached": False,
            "n_perm": 100}
    base.update(over)
    return base


def _rwr_result(ok=True):
    return {"ok": ok, "converged": True, "n_iter": 12, "total_mass": 1.0,
            "seed_mass": 0.4, "coverage": {"seed_in_lcc": 5, "seed_total": 5},
            "seeds": [{"gene": "TP53", "score": 0.5}],
            "top_nonseed": [{"gene": "EGFR", "score": 0.3},
                            {"gene": "KRAS", "score": 0.2}],
            "error": None}


def test_proximity_command_success(tmp_path, monkeypatch):
    """npap proximity 全绿:两侧靶点表 + 邻近度计算(打桩)→ proximity.json 落盘,退出码 0。"""
    monkeypatch.setattr(cli.interactome, "is_built", lambda *a, **k: True)
    drug_df = pd.DataFrame({"approvedSymbol": ["TP53", "EGFR"]})
    disease_df = pd.DataFrame({"approvedSymbol": ["TP53", "BRCA1"]})
    monkeypatch.setattr(cli, "get_drug_targets_merged", lambda drug: drug_df)
    monkeypatch.setattr(cli, "get_disease_targets_merged", lambda disease: disease_df)
    monkeypatch.setattr(cli, "run_proximity_cached",
                        lambda drug_syms, disease_syms, label=None: _prox_result())
    rc = cli.main(["proximity", "--drug", "metformin", "--disease", "T2D",
                   "--out", str(tmp_path)])
    assert rc == 0
    out = json.loads((tmp_path / "proximity.json").read_text(encoding="utf-8"))
    assert out["ok"] is True and out["z"] == -2.3


def test_rwr_command_success(tmp_path, monkeypatch):
    """npap rwr 全绿:种子排序表落盘(gene,score,is_seed),退出码 0。"""
    monkeypatch.setattr(cli.interactome, "is_built", lambda *a, **k: True)
    monkeypatch.setattr(cli, "run_rwr_cached",
                        lambda seeds, label=None, annotate=None: _rwr_result())
    rc = cli.main(["rwr", "--seeds", "TP53,EGFR", "--out", str(tmp_path)])
    assert rc == 0
    out = pd.read_csv(tmp_path / "rwr.csv")
    assert list(out.columns) == ["gene", "score", "is_seed"]
    assert out.iloc[0]["gene"] == "TP53" and out.iloc[0]["is_seed"] == 1
    # 失败路径:游走不收敛/失败 → 非零退出
    monkeypatch.setattr(cli, "run_rwr_cached",
                        lambda seeds, label=None: {"ok": False, "error": "no converge"})
    rc2 = cli.main(["rwr", "--seeds", "TP53", "--out", str(tmp_path / "bad")])
    assert rc2 == 1


def test_save_json_locked_is_none(tmp_path, monkeypatch):
    """JSON 落盘被占用(PermissionError):警告并返回 None,不崩溃。"""
    from WY.core.caching import save_json_file
    monkeypatch.setattr(cli, "save_json_file",
                        lambda *a, **k: (_ for _ in ()).throw(PermissionError("locked")))
    assert cli._save_json({"a": 1}, str(tmp_path / "x.json")) is None


def test_run_pair_full_skeleton(tmp_path, monkeypatch):
    """run_pair 真骨架(全部 core 重活打桩):靶点→交集→PPI→富集 ok;
    邻近度/RWR 成功;summary.json 落盘且全步 ok。"""
    drug_df = pd.DataFrame({"approvedSymbol": ["TP53", "EGFR"]})
    disease_df = pd.DataFrame({"approvedSymbol": ["TP53", "BRCA1"]})
    monkeypatch.setattr(cli, "get_targets_concurrent",
                        lambda d, dis, progress_cb=None: (drug_df, disease_df, []))
    edges = pd.DataFrame([{"protein_A": "TP53", "protein_B": "EGFR", "score": 900}])
    monkeypatch.setattr(cli, "fetch_string_ppi", lambda genes, **k: edges)
    monkeypatch.setattr(cli, "run_enrichment",
                        lambda genes, **k: (pd.DataFrame(
                            [{"source": "KEGG", "term_name": "x", "p_value": 0.01}]),
                            "gProfiler"))
    monkeypatch.setattr(cli.interactome, "is_built", lambda *a, **k: True)
    monkeypatch.setattr(cli, "run_proximity_cached",
                        lambda drug_syms, disease_syms, label=None: _prox_result())
    monkeypatch.setattr(cli, "run_rwr_cached",
                        lambda seeds, label=None, annotate=None: _rwr_result())
    monkeypatch.setattr(config, "RWR_EXPAND_TOP_N", 0)
    out_dir = str(tmp_path / "pair")
    summary = cli.run_pair("metformin", "T2D", out_dir)
    assert summary["ok"] is True, summary["steps"]
    assert summary["steps"]["targets"]["status"] == "ok"
    assert summary["steps"]["proximity"]["status"] == "ok"
    assert summary["steps"]["rwr"]["status"] == "ok"
    assert (tmp_path / "pair" / "summary.json").exists()
    assert (tmp_path / "pair" / "common_targets.csv").exists()


def test_run_pair_empty_intersection_skips_rest(tmp_path, monkeypatch):
    """交集为空:ppi/enrichment/proximity/rwr 标记 skipped,summary.ok 为 False,
    不发任何后续请求。"""
    monkeypatch.setattr(cli, "get_targets_concurrent",
                        lambda d, dis, progress_cb=None: (
                            pd.DataFrame(), pd.DataFrame(), []))
    monkeypatch.setattr(config, "RWR_EXPAND_TOP_N", 0)

    def must_not_call(*a, **k):
        raise AssertionError("downstream step ran despite empty intersection")

    monkeypatch.setattr(cli, "fetch_string_ppi", must_not_call)
    monkeypatch.setattr(cli, "run_enrichment", must_not_call)
    summary = cli.run_pair("a", "b", str(tmp_path / "pair"))
    assert summary["ok"] is False
    assert all(summary["steps"][n]["status"] == "skipped"
               for n in ("ppi", "enrichment", "proximity", "rwr"))


def test_batch_retries_failed_previous_run(tmp_path, monkeypatch):
    """--resume:此前失败(非 ok summary)的对会被重跑;成功过的才跳过。"""
    calls = {"n": 0}
    pair_dir = tmp_path / "a__x"
    pair_dir.mkdir()
    (pair_dir / "summary.json").write_text(
        json.dumps({"ok": False, "steps": {}, "elapsed_seconds": 1.0}),
        encoding="utf-8")

    def run(drug, disease, out):
        calls["n"] += 1
        os.makedirs(out, exist_ok=True)
        summary = {"ok": True, "steps": {"targets": {"status": "ok", "detail": "d"}},
                   "elapsed_seconds": 2.0}
        (os.path.join(out, "summary.json"))
        with open(os.path.join(out, "summary.json"), "w", encoding="utf-8") as f:
            json.dump(summary, f)
        return summary

    monkeypatch.setattr(cli, "run_pair", run)
    csv = tmp_path / "pairs.csv"
    pd.DataFrame({"drug": ["a"], "disease": ["x"]}).to_csv(csv, index=False)
    rc = cli.main(["batch", "--pairs", str(csv), "--out", str(tmp_path), "--resume"])
    assert rc == 0 and calls["n"] == 1  # 失败过的被重试
