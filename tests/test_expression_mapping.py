"""core.expression 的探针→基因符号映射测试(注释文件打桩,纯本地)。

阶段三续:深水区补测第一批。GPL 注释 gz 由测试合成,_fetch_gpl_annotation 打桩。
"""
import gzip

import numpy as np
import pandas as pd
import pytest

import WY.core.expression as expr


def _write_annot(path, rows):
    """合成 GPL 注释 gz(platform_table 段,ID + Gene symbol 两列)。"""
    lines = ["!platform_table_begin", "ID\tGene symbol"]
    lines += [f"{probe}\t{sym}" for probe, sym in rows]
    lines.append("!platform_table_end")
    with open(path, "wb") as f:
        f.write(gzip.compress(("\n".join(lines) + "\n").encode("utf-8")))
    return path


def _probe_matrix():
    """4 探针 × 4 样本(行=探针):ILMN_1/2 映射同一基因 TP53(考察 IQR 择优),
    ILMN_3→EGFR,ILMN_4→注释缺失。"""
    return pd.DataFrame.from_dict(
        {
            "ILMN_1": [10.0, 11.0, 12.0, 13.0],   # IQR 大 → 胜出
            "ILMN_2": [10.1, 10.2, 10.3, 10.4],   # 同基因,几乎恒定 → 被弃
            "ILMN_3": [5.0, 6.0, 7.0, 8.0],
            "ILMN_4": [1.0, 2.0, 3.0, 4.0],       # 注释为 --- → 不映射
        },
        orient="index", columns=["S1", "S2", "S3", "S4"])


def test_map_probe_to_symbol_iqr_best_and_skip_missing(monkeypatch, tmp_path):
    """映射:注释缺失行(-- -)跳过;多探针同基因取 IQR 最大者;覆盖率达标正常输出。"""
    annot = _write_annot(str(tmp_path / "annot.gz"),
                         [("ILMN_1", "TP53"), ("ILMN_2", "TP53"),
                          ("ILMN_3", "EGFR"), ("ILMN_4", "---")])
    monkeypatch.setattr(expr, "_fetch_gpl_annotation", lambda gpl: annot)
    out = expr.map_probe_to_symbol(_probe_matrix(), "GPL570")
    assert list(out.index) == ["EGFR", "TP53"]  # 排序后输出
    # TP53 行 = IQR 更大的 ILMN_1 的数据
    assert list(out.loc["TP53"]) == [10.0, 11.0, 12.0, 13.0]


def test_map_probe_coverage_below_half_raises(monkeypatch, tmp_path):
    """映射覆盖率 < 50%:直接报错(平台注释不完整,结果不可信),不猜测。"""
    annot = _write_annot(str(tmp_path / "annot.gz"),
                         [("ILMN_1", "TP53"),
                          ("UNKN_a", "X"), ("UNKN_b", "Y"), ("UNKN_c", "Z")])
    monkeypatch.setattr(expr, "_fetch_gpl_annotation", lambda gpl: annot)
    with pytest.raises(RuntimeError, match="覆盖率"):
        expr.map_probe_to_symbol(_probe_matrix(), "GPL570")


def test_map_symbol_level_passthrough(monkeypatch):
    """行标识本就是基因符号级:原样返回,不碰注释、不需要 GPL。"""
    matrix = pd.DataFrame({"s1": [1.0]}, index=["TP53"])
    monkeypatch.setattr(expr, "_fetch_gpl_annotation",
                        lambda gpl: (_ for _ in ()).throw(AssertionError("annotation fetched")))
    out = expr.map_probe_to_symbol(matrix, "")
    assert out is matrix or out.equals(matrix)


def test_map_probe_without_gpl_raises(monkeypatch):
    """探针级数据缺平台信息:明确报错指引换 symbol 级数据集。"""
    monkeypatch.setattr(expr, "_fetch_gpl_annotation",
                        lambda gpl: (_ for _ in ()).throw(AssertionError("annotation fetched")))
    with pytest.raises(RuntimeError, match="GPL"):
        expr.map_probe_to_symbol(_probe_matrix(), "")


def test_looks_count_matrix():
    """raw counts 判定:非负整数真;含小数/负值假。"""
    assert expr.looks_count_matrix(pd.DataFrame({"a": [1, 2], "b": [0, 5]})) is True
    assert expr.looks_count_matrix(pd.DataFrame({"a": [1.5, 2], "b": [0, 5]})) is False
    assert expr.looks_count_matrix(pd.DataFrame({"a": [-1, 2], "b": [0, 5]})) is False
