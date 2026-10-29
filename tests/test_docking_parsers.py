"""core.docking 的纯解析与归档函数测试(PDB/盒子/构象切分/分数解析/元数据/版本控制)。

深水区补测第二批;全部离线,文件走 tmp_path 与 conftest 沙箱。
"""
import os

import pandas as pd
import pytest

from WY import config
from WY.core.text_utils import remember_name_mapping, safe_name
import WY.core.docking as dock


def _atom(rec, serial, name, resn, chain, resi, x, y, z, elem="C"):
    return (f"{rec:<6s}{serial:5d}  {name:<3s} {resn:>3s} {chain}{resi:>4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n")


def _write_pdb(path, lines):
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(lines))
    return path


# ---------- 指纹与防覆盖 ----------


def test_docking_fingerprint_and_unique_path(tmp_path):
    """指纹 = 内容 MD5 前 12 位,文件缺失返回空串;防覆盖:已存在时改带时间戳的新名。"""
    f1 = tmp_path / "a.pdb"
    f1.write_text("ATOM", encoding="utf-8")
    fp = dock._docking_fingerprint(str(f1))
    assert len(fp) == 12 and fp == dock._docking_fingerprint(str(f1))
    assert dock._docking_fingerprint(str(tmp_path / "nope")) == ""
    assert dock._unique_path(str(tmp_path / "free.csv")) == str(tmp_path / "free.csv")
    existing = tmp_path / "taken.csv"
    existing.write_text("x", encoding="utf-8")
    new = dock._unique_path(str(existing))
    assert new != str(existing) and not os.path.exists(new)


# ---------- PDB 解析:共晶配体 / 盒内残基 ----------


def test_parse_hetatm_and_pdb_coords(tmp_path):
    """共晶配体识别:HETATM 保留、水(HOH)剔除;坐标按蛋白/非水其他分子分列。"""
    pdb = _write_pdb(str(tmp_path / "x.pdb"), [
        _atom("ATOM", 1, "CA", "ALA", "A", 10, 0.0, 0.0, 0.0),
        _atom("HETATM", 2, "O", "HOH", "A", 999, 50.0, 50.0, 50.0),   # 水 → 剔
        _atom("HETATM", 3, "FE", "HEM", "A", 200, 10.0, 10.0, 10.0),  # 辅因子 → 留
    ])
    het = dock.parse_hetatm(pdb)
    assert [(h[0], h[3]) for h in het] == [("HEM", (10.0, 10.0, 10.0))]
    protein, others = dock._pdb_coords(pdb)
    assert len(protein) == 1 and len(others) == 1  # 水两边都不在


def test_residues_in_box(tmp_path):
    """盒内残基清单:按 链+残基号:残基名 去重排序,盒外残基不进。"""
    pdb = _write_pdb(str(tmp_path / "x.pdb"), [
        _atom("ATOM", 1, "CA", "ALA", "A", 10, 1.0, 1.0, 1.0),     # 盒内
        _atom("ATOM", 2, "CA", "GLY", "A", 11, 2.0, 2.0, 2.0),     # 盒内
        _atom("ATOM", 3, "CA", "VAL", "B", 12, 99.0, 99.0, 99.0),  # 盒外
    ])
    box = {"center": (1.5, 1.5, 1.5), "size": (4.0, 4.0, 4.0)}
    assert dock.residues_in_box(pdb, box) == ["A10:ALA", "A11:GLY"]


def test_infer_box_from_cofactor_prefers_ligand(tmp_path):
    """盒子推断:优先用共晶辅因子定位;无辅因子时退回全蛋白范围;完全无坐标报错。"""
    pdb = _write_pdb(str(tmp_path / "x.pdb"), [
        _atom("ATOM", 1, "CA", "ALA", "A", 10, 100.0, 100.0, 100.0),
        _atom("HETATM", 2, "FE", "HEM", "A", 200, 10.0, 20.0, 30.0),
    ])
    box = dock.infer_box_from_cofactor(pdb, padding=8.0)
    assert box["center"] == (10.0, 20.0, 30.0)          # 只按辅因子算中心
    assert box["size"] == (16.0, 16.0, 16.0)            # 单点 + 2*padding
    empty = _write_pdb(str(tmp_path / "e.pdb"), ["END\n"])
    with pytest.raises(RuntimeError, match="no coordinates"):
        dock.infer_box_from_cofactor(empty)


# ---------- vina 输出:构象切分与分数解析 ----------


def _vina_output(affinities, skip_remark_for=None):
    lines = []
    for i, aff in enumerate(affinities, 1):
        lines.append(f"MODEL {i}")
        lines.append(f"REMARK VINA RESULT:    {aff:>9.3f}      0.000      0.000")
        if skip_remark_for != i:
            lines.append("REMARK x = 1")
        lines.append("ENDMDL")
    return "\n".join(lines) + "\n"


def test_iter_and_split_vina_models(tmp_path):
    """多模型切分:按 MODEL/ENDMDL 逐块,内置拆分器产出与构象数一致的独立文件。"""
    out = tmp_path / "job_out.pdbqt"
    out.write_text(_vina_output([-7.9, -6.5]), encoding="utf-8")
    blocks = list(dock._iter_vina_models(str(out)))
    assert len(blocks) == 2 and blocks[0][0].startswith("MODEL")
    out_dir = tmp_path / "split"
    out_dir.mkdir()
    files = dock._split_models_python(str(out), str(out_dir))
    assert len(files) == 2 and all(os.path.exists(f) for f in files)


def test_parse_vina_scores_sorted_and_tolerant(tmp_path):
    """分数解析:按结合能升序(越负越强);数值残缺的构象跳过;无打分行 → 空表。"""
    text = []
    for i, aff in enumerate([-6.5, -7.9, None], 1):
        text.append(f"MODEL {i}")
        if aff is not None:
            text.append(f"REMARK VINA RESULT:    {aff:>9.3f}      0.000      0.000")
        else:
            text.append("REMARK VINA RESULT:    garbage")  # 数值残缺 → 跳过
        text.append("ENDMDL")
    out = tmp_path / "out.pdbqt"
    out.write_text("\n".join(text) + "\n", encoding="utf-8")
    df = dock.parse_vina_scores(str(out))
    # mode 是原始模型号,行序按结合能升序(越负越强)
    assert list(df["mode"]) == [2, 1]
    assert list(df["affinity"]) == [-7.9, -6.5]
    # 完全无打分行 → 空表
    empty = tmp_path / "empty.pdbqt"
    empty.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    assert dock.parse_vina_scores(str(empty)).empty


# ---------- 元数据 / 收藏 / 版本控制 ----------


def test_write_docking_metadata(tmp_path):
    """元数据:# 键值对 + 时间戳落盘,返回实际路径;同名自动换新文件名。"""
    p1 = dock.write_docking_metadata(str(tmp_path / "meta.txt"),
                                     {"drug": "metformin", "gene": "TP53"})
    text = open(p1, encoding="utf-8").read()
    assert "# drug: metformin" in text and "# gene: TP53" in text
    assert "# timestamp:" in text
    p2 = dock.write_docking_metadata(p1, {"drug": "again"})
    assert p2 != p1, "同名应自动换新文件名"


def test_star_docking_result_copies_to_star_dir(tmp_path, monkeypatch):
    """收藏:整个成果目录复制进 MD_STAR_DIR(沙箱),原目录不动。"""
    src = tmp_path / "result_v1"
    src.mkdir()
    (src / "out.pdbqt").write_text("data", encoding="utf-8")
    star_dir = tmp_path / "star"
    monkeypatch.setattr(config, "MD_STAR_DIR", str(star_dir))
    dest = dock.star_docking_result(str(src))
    assert dest.startswith(str(star_dir))
    assert os.path.exists(os.path.join(dest, "out.pdbqt"))
    assert (src / "out.pdbqt").exists()  # 原目录不受影响


def test_prune_docking_results_version_control(tmp_path, monkeypatch):
    """版本控制:同前缀超过保留上限删最旧;只在结果根目录内动手;不合命名约定的目录不管。"""
    md_root = tmp_path / "MDresult"
    monkeypatch.setattr(config, "MD_RESULT_DIR", str(md_root))
    md_root.mkdir()
    parent = md_root / "metformin_T2D_TP53"
    parent.mkdir()
    for ts in ("20260101_000001", "20260102_000002", "20260103_000003"):
        (parent / f"metformin_T2D_TP53_{ts}").mkdir()
    removed = dock.prune_docking_results(
        str(parent / "metformin_T2D_TP53_20260103_000003"), max_versions=2)
    assert len(removed) == 1 and removed[0].endswith("20260101_000001")
    assert (parent / "metformin_T2D_TP53_20260102_000002").exists()
    # 目录名无时间戳段 → 不做版本控制
    plain = md_root / "no_ts_dir"
    plain.mkdir()
    assert dock.prune_docking_results(str(plain), max_versions=1) == []


def test_split_result_prefix_uses_name_map(tmp_path):
    """前缀反拆:优先借名称映射还原含空格的药物名;无映射退回旧拆法(首个下划线切段)。"""
    remember_name_mapping("Salicylic acid")          # 有损 → 记映射
    remember_name_mapping("type 2 diabetes")
    drug, disease = dock.split_result_prefix(
        f"{safe_name('Salicylic acid')}_{safe_name('type 2 diabetes')}")
    assert drug == "Salicylic acid" and disease == "type 2 diabetes"
    drug2, disease2 = dock.split_result_prefix("aspirin_cancer")
    assert drug2 == "aspirin" and disease2 == "cancer"  # 无映射退回旧拆法
