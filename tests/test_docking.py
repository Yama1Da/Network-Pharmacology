"""core.docking 的测试:盒子收拢、受体 PDB 清洗、UniProt 缓存、AlphaFold 兜底、
RCSB 候选检索、组成防呆标注、口袋残基/相互作用几何判定、3D HTML、前缀消毒、缓存清理。

迁移自 selfcheck_docking.py(18 项中的 14 项)。网络请求一律打桩;
缓存路径由 conftest 沙箱指到每个测试的临时目录——自检脚本里的
_tmp_uniprot_cache / 手工临时目录换装由此整个省去。
"""
import json
import os

import requests

from WY import config
import WY.core.docking as dock


def _json_nonempty(path):
    try:
        with open(path, encoding="utf-8") as f:
            return bool(json.load(f))
    except Exception:
        return False


def _atom_line(rec: str, serial: int, name: str, resn: str, chain: str,
               resi: int, x: float, y: float, z: float, elem: str) -> str:
    """按 PDB/pdbqt 列格式拼一条原子记录(elem 处 pdbqt 传 AutoDock 原子类型)。"""
    return (f"{rec:<6s}{serial:5d}  {name:<3s} {resn:>3s} {chain}{resi:>4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n")


def _write_tmp(path: str, lines: list[str]) -> str:
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(lines))
    return path


def test_clamp_box():
    """盒子尺寸收拢:只有超限的边被收拢,中心必须保持不变,未超限的盒子原样返回。"""
    box = {"center": (1.0, 2.0, 3.0), "size": (30.0, 35.0, 60.0)}
    under = dock.clamp_box(box, 40.0)
    assert under["size"] == (30.0, 35.0, 40.0), under  # 只有超限的边被收拢
    assert under["center"] == (1.0, 2.0, 3.0), "center must be preserved"
    small = dock.clamp_box({"center": (0, 0, 0), "size": (20, 20, 20)}, 40.0)
    assert small["size"] == (20.0, 20.0, 20.0), "under-limit box must be unchanged"


def test_strip_af_style_pdb(tmp_path):
    """AlphaFold 模型的 PDB 格式特殊,清洗后应只保留原子记录并去掉重复原子。"""
    src = os.path.join(str(tmp_path), "af_like.pdb")
    out = src + ".clean.pdb"
    with open(src, "w", encoding="utf-8") as f:
        f.write("REMARK   1 ALPHAFOLD MONOMER V2.0 PREDICTION\n")
        f.write("TITLE    ALPHAFOLD MONOMER PREDICTION FOR TEST\n")
        f.write("ATOM      1  N   MET A   1      63.063 -13.528  -5.589  1.00 42.72           N\n")
        f.write("ATOM      2  CA  MET A   1      62.454 -12.831  -4.439  1.00 42.72           C\n")
        f.write("ATOM      3  CA  MET A   1      62.454 -12.831  -4.439  1.00 42.72 B\n")  # altloc B → 去重
        f.write("HETATM    4  O   HOH A 999      10.000  10.000  10.000  1.00 20.00           O\n")
        f.write("END\n")
    dock._strip_receptor_pdb(src, out)
    lines = open(out, encoding="utf-8").read().splitlines()
    assert all(ln.startswith(("ATOM", "TER", "END")) for ln in lines), lines
    assert sum(1 for ln in lines if ln.startswith("ATOM")) == 2, "altloc B must be dropped"


def test_ebi_gene_names():
    """EBI 条目的正式名与同义词都应收入基因名集合;空条目返回空集。"""
    entry = {"gene": [{"name": {"value": "ANO9"},
                       "synonyms": [{"value": "PIG5"}, {"value": "TP53I5"}]}]}
    names = dock._gene_names_from_ebi_entry(entry)
    assert names == {"ANO9", "PIG5", "TP53I5"}, names
    assert dock._gene_names_from_ebi_entry({}) == set()


def test_uniprot_negative_cache(monkeypatch):
    """数据库里查无此基因时记入缓存,第二次查询不再重复发请求。"""
    calls = {"n": 0}

    def fake_ebi(gene, taxon):
        calls["n"] += 1
        return ""  # 确定性"查无"

    monkeypatch.setattr(dock, "_lookup_uniprot_ebi", fake_ebi)
    assert dock.fetch_uniprot_accession("NOPE_X", 9606) == ""
    assert dock.fetch_uniprot_accession("NOPE_X", 9606) == ""
    assert calls["n"] == 1, f"negative cache failed: {calls['n']} calls"


def test_uniprot_network_failure_not_cached(monkeypatch):
    """网络全部失败时应当报错且不写缓存,避免把临时故障当成"查无此基因"。"""
    def boom(gene, taxon):
        raise ConnectionError("simulated network outage")

    monkeypatch.setattr(dock, "_lookup_uniprot_ebi", boom)
    monkeypatch.setattr(dock, "_lookup_uniprot_rest", boom)
    try:
        dock.fetch_uniprot_accession("TP53", 9606)
        raise AssertionError("all sources down must raise RuntimeError")
    except RuntimeError:
        pass
    tmp = config.UNIPROT_LOOKUP_CACHE
    assert not os.path.exists(tmp) or not _json_nonempty(tmp), "failure must not be cached"


def test_af_pdb_skip_existing(monkeypatch):
    """AlphaFold 模型文件已存在时直接复用,不再发网络请求。"""
    os.makedirs(config.DOCKING_RECEPTOR_CACHE_DIR, exist_ok=True)
    pre = os.path.join(config.DOCKING_RECEPTOR_CACHE_DIR,
                       f"AF_P99999_{config.ALPHAFOLD_MODEL_VERSION}.pdb")
    with open(pre, "w", encoding="utf-8") as f:
        f.write("REMARK x\n" + "ATOM " * 300)
    monkeypatch.setattr(
        requests, "get",
        lambda *a, **k: (_ for _ in ()).throw(AssertionError("network call made")))
    path = dock.fetch_alphafold_pdb("P99999")
    assert path == pre


def test_candidates_af_fallback_row(monkeypatch):
    """实验结构查不到且允许 AlphaFold 时,追加预测结构作候选;关闭时不追加。"""
    # 单基因入口是批量版的薄包装,检索统一走 _rcsb_search_entry_ids(带 HTTP 缓存会话),
    # 打桩层级在检索函数而不是 requests.post
    monkeypatch.setattr(dock, "fetch_uniprot_accession", lambda gene, taxon=9606: "P99999")
    monkeypatch.setattr(dock, "_rcsb_search_entry_ids", lambda gene, max_hits: [])
    df = dock.fetch_receptor_pdb_candidates("FAKEGENE", allow_alphafold=True)
    assert len(df) == 1 and df.iloc[0]["pdb_id"] == "AFDB:P99999", df
    assert df.iloc[0]["source"] == "alphafold"
    df2 = dock.fetch_receptor_pdb_candidates("FAKEGENE", allow_alphafold=False)
    assert df2.empty and "source" in df2.columns


def test_candidates_search_failure_friendly_error(monkeypatch):
    """检索服务不可达时,单基因入口抛出普通用户可读的中文错误(而不是静默按无结构处理)。"""
    def _boom(gene, max_hits):
        raise requests.exceptions.ConnectionError("simulated outage")

    monkeypatch.setattr(dock, "_rcsb_search_entry_ids", _boom)
    try:
        dock.fetch_receptor_pdb_candidates("TP53")
        raise AssertionError("expected RuntimeError")
    except RuntimeError as e:
        msg = str(e)
        assert "RCSB 搜索失败（TP53）" in msg, msg
        assert "无法连接 RCSB 服务器" in msg, msg
    # 批量模式保持原语义:单基因失败不中断整批,按无结构处理
    monkeypatch.setattr(dock, "_rcsb_search_entry_ids", _boom)
    results = dock.fetch_receptor_candidates_batch(["TP53", "EGFR"])
    assert results["TP53"].empty and results["EGFR"].empty


def test_receptor_composition_flag(monkeypatch):
    """受体检索防呆标注(AGT 三重撞名事故的防再犯措施):纯核酸→阻断级,蛋白-核酸复合物→警示级。"""
    f = dock.composition_flag
    assert f(0, 2, 0) == "dna_only"          # 事故同款:纯 DNA 双链
    assert f(0, 0, 1) == "dna_only"          # 纯 RNA 同理阻断
    assert f(1, 1, 0) == "has_nucleic"       # 蛋白-DNA 复合物:警示不阻断
    assert f(3, None, None) == ""            # 正常蛋白结构
    assert f(1, 0, 0) == ""                  # 无核酸
    assert f(None, None, None) == ""         # 计数全缺失(旧条目):按未知,不拦截
    assert f(None, 2, 0) == "has_nucleic"    # 蛋白数未知但确含核酸:保守警示
    # 组装链:批量结果表的 flag 列由元数据计数生成
    monkeypatch.setattr(dock, "_rcsb_search_entry_ids",
                        lambda gene, max_hits: ["2MMQ", "2WXW"])
    monkeypatch.setattr(dock, "_rcsb_graphql_entries", lambda ids: {
        "2MMQ": {"resolution": None, "method": "SOLUTION NMR", "title": "AGT FAPY duplex",
                 "n_protein": 0, "n_dna": 2, "n_rna": 0},
        "2WXW": {"resolution": 3.3, "method": "X-RAY DIFFRACTION",
                 "title": "Crystal structure of human angiotensinogen",
                 "n_protein": 1, "n_dna": 0, "n_rna": 0},
    })
    df = dock.fetch_receptor_pdb_candidates("AGT")
    assert df.set_index("pdb_id").loc["2MMQ", "flag"] == "dna_only", df
    assert df.set_index("pdb_id").loc["2WXW", "flag"] == "", df


def test_residues_near_ligand_cutoff(tmp_path):
    """口袋残基判定半径应可调:同一体系 4.0 Å 判不进、5.0 Å 判得进。"""
    lig = _write_tmp(os.path.join(str(tmp_path), "lig.pdbqt"), [
        _atom_line("ATOM", 1, "N1", "LIG", "B", 1, 0.0, 0.0, 0.0, "N")])
    # CA 在 4.5 Å 处:默认 4 Å 不算口袋残基,放宽到 5 Å 才算
    rec = _write_tmp(os.path.join(str(tmp_path), "rec.pdb"), [
        _atom_line("ATOM", 1, "CA", "ALA", "A", 10, 4.5, 0.0, 0.0, "C"),
        _atom_line("ATOM", 2, "CB", "ALA", "A", 10, 5.6, 0.9, 0.0, "C")])
    assert dock.residues_near_ligand(rec, lig) == [], "4.5 Å residue must be out at 4.0 Å"
    near5 = dock.residues_near_ligand(rec, lig, cutoff=5.0)
    assert len(near5) == 1 and near5[0]["resn"] == "ALA", near5
    assert near5[0]["chain"] == "A" and near5[0]["resi"] == "10", near5


def test_interactions_near_ligand(tmp_path):
    """几何判据检测氢键/盐桥/疏水:类型、去重(盐桥优先)、排序与坐标须正确。"""
    # 配体三个原子彼此相距 5 Å 以上,互不干扰:N1(原点) O1(5,0,0) C1(0,6,0)
    lig = _write_tmp(os.path.join(str(tmp_path), "lig.pdbqt"), [
        _atom_line("ATOM", 1, "N1", "LIG", "B", 1, 0.0, 0.0, 0.0, "N"),
        _atom_line("ATOM", 2, "O1", "LIG", "B", 1, 5.0, 0.0, 0.0, "OA"),
        _atom_line("ATOM", 3, "C1", "LIG", "B", 1, 0.0, 6.0, 0.0, "A")])
    rec = _write_tmp(os.path.join(str(tmp_path), "rec.pdb"), [
        # LYS NZ:碱性 N,距 O1 3.0 Å → 盐桥(不再重复记氢键);距 N1 2.0 Å → 氢键
        _atom_line("ATOM", 1, "NZ", "LYS", "A", 10, 2.0, 0.0, 0.0, "N"),
        # GLU OE1:酸性 O,距 O1 2.5 Å → 氢键(酸 O↔配体 O 不算盐桥)
        _atom_line("ATOM", 2, "OE1", "GLU", "A", 11, 5.0, 2.5, 0.0, "O"),
        # ALA CB:距 C1 2.2 Å → 疏水接触
        _atom_line("ATOM", 3, "CB", "ALA", "A", 12, 0.0, 3.8, 0.0, "C"),
        # 远处残基:不应产生任何相互作用
        _atom_line("ATOM", 4, "CA", "VAL", "A", 13, 30.0, 30.0, 30.0, "C")])
    its = dock.interactions_near_ligand(rec, lig)
    kinds = [it["kind"] for it in its]
    assert kinds == ["salt", "hbond", "hbond", "hydro"], kinds  # 盐桥优先、同类型按距离
    salt = its[0]
    assert (salt["resn"], salt["res_atom"], salt["lig_atom"], salt["dist"]) \
        == ("LYS", "NZ", "O1", 3.0), salt
    hb_lig_n = next(it for it in its if it["lig_atom"] == "N1")
    assert (hb_lig_n["resn"], hb_lig_n["res_atom"], hb_lig_n["dist"]) \
        == ("LYS", "NZ", 2.0), hb_lig_n
    hb_glu = next(it for it in its if it["resn"] == "GLU")
    assert (hb_glu["kind"], hb_glu["lig_atom"], hb_glu["dist"]) == ("hbond", "O1", 2.5), hb_glu
    hydro = its[3]
    assert (hydro["resn"], hydro["res_atom"], hydro["lig_atom"], hydro["dist"]) \
        == ("ALA", "CB", "C1", 2.2), hydro
    # 画虚线用的两端坐标须与原子坐标一致
    assert salt["pos1"] == (2.0, 0.0, 0.0) and salt["pos2"] == (5.0, 0.0, 0.0), salt


def test_interactive_3d_html():
    """交互式 3D HTML:坐标须以合法 JSON 内嵌,离线/在线两种模式标记正确。"""
    receptor = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    ligand = "ATOM      1  N   LIG B   1       1.000   1.000   1.000  1.00  0.00           N\n"
    # 离线模式:传入库源码,HTML 直接内嵌且标记为离线可用
    html, offline = dock.build_interactive_3d_html(
        receptor, ligand, caption="gene · pose 1", lib_js="var $3Dmol = {};")
    assert offline and "$3Dmol.createViewer" in html
    assert 'var receptor = "ATOM' in html, "receptor must be embedded as a JSON string"
    # 内嵌的坐标文本必须是合法 JSON(转义正确,浏览器端才能还原)
    embedded = html.split("var receptor = ", 1)[1].split(";\n", 1)[0]
    assert json.loads(embedded) == receptor, "embedded receptor must round-trip"
    # 在线模式:不传库源码,应改用 CDN 引用并明确标记非离线
    html2, offline2 = dock.build_interactive_3d_html(receptor, ligand)
    assert not offline2 and 'src="https://3dmol.org' in html2
    # 3Dmol.js 许可标注:内嵌/CDN 两种模式措辞不同,但版权行与完整许可文本都必须随文件提供
    assert "内嵌 3Dmol.js (BSD-3-Clause)" in html, html[html.find("#lic"):html.find("#lic") + 200]
    assert "本页面通过 CDN 加载 3Dmol.js (BSD-3-Clause)" in html2
    for h in (html, html2):
        assert "Copyright (c) 2014, University of Pittsburgh" in h, "attribution missing"
        assert "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS" in h, "license text missing"
    # 默认导出:无口袋残基、无相互作用、无共晶配体参照,浏览器端不渲染对应元素
    assert "var pocket = []" in html and "var het = null" in html
    assert "var interactions = []" in html, "default interactions must be empty"
    assert "Focus pocket" in html, "focus button must exist"
    # 聚焦口袋须把周围蛋白卡通虚化(半透明),双击复位还原
    assert "applyProteinStyles" in html and "opacity: 0.3" in html, "focus must fade protein cartoon"
    # 传入口袋残基、相互作用与共晶配体:三类数据须以合法 JSON 嵌入(浏览器端据此还原样式与标注)
    pocket = [{"chain": "A", "resi": "1", "resn": "ALA", "pos": (0.0, 0.0, 0.0)}]
    inter = [{"kind": "salt", "chain": "A", "resi": "1", "resn": "ALA",
              "res_atom": "CA", "lig_atom": "N", "dist": 3.0,
              "pos1": (0.0, 0.0, 0.0), "pos2": (1.0, 1.0, 1.0)}]
    html3, _ = dock.build_interactive_3d_html(
        receptor, ligand, pocket=pocket, interactions=inter, focus_label="聚焦口袋",
        het_pdb="HETATM    1  O   LIG B   2       2.000   1.000   1.000  1.00  0.00           O\n")
    assert 'var pocket = [{"chain": "A"' in html3, "pocket must be embedded as JSON"
    assert 'var interactions = [{"kind": "salt"' in html3, "interactions must be embedded as JSON"
    assert 'var het = "HETATM' in html3, "co-crystal ligand text must be embedded as JSON"
    assert "聚焦口袋" in html3, "focus button label must be localized"
    # 相互作用用虚线圆柱绘制(1px 线白底下看不清),且配色与页面图例一致
    assert "addCylinder" in html3 and 'hbond: "#0f62fe"' in html3, html3[html3.find("lineColors"):][:200]


def test_docking_result_prefix():
    """归档前缀消毒:空段剔除(不出现 unnamed)、路径穿越/分隔符净化、正常名与主分析页命名对齐。"""
    f = dock.docking_result_prefix
    # 正常名:与主分析页(Data_Network)的前缀命名规则一致
    assert f("metformin", "Type 2 diabetes mellitus", "EGFR") == \
        "metformin_Type_2_diabetes_mellitus_EGFR"
    # 手动基因模式 drug/disease 为空:空段剔除,不出现 unnamed、不产生前导下划线
    assert f("", "", "TP53") == "TP53"
    assert f("  ", "", "TP53") == "TP53"
    # 路径穿越与分隔符:目录名不含 / 等字符,不逃出结果目录(旧实现 ".." 可逃逸一级)
    assert f("../evil", "a/b", "TP53") == "evil_a_b_TP53"
    assert f("..", "", "TP53") == "unnamed_TP53"
    # 含点/撇号的名字:统一归下划线,与分析页产物命名对齐
    assert f("St. John's Wort", "", "TP53") == "St_John_s_Wort_TP53"
    # 中文输入保留
    assert f("二甲双胍", "2型糖尿病", "EGFR") == "二甲双胍_2型糖尿病_EGFR"
    # 全部为纯特殊字符段:safe_name 各自兜底为 unnamed
    assert f(".", "..", "!") == "unnamed_unnamed_unnamed"
    # 全空:回退固定名,避免空文件名
    assert f("", "", "") == "docking"


def test_clean_docking_cache_exempts_permanent():
    """缓存清理只扫 raw/ligand/receptor 三个子目录:过期文件按各自 TTL 删除,
    缓存根目录下的长期缓存(UniProt 基因映射 / 3Dmol.js 库)即使 mtime 很老也不动。"""
    import time

    cache_dir = config.DOCKING_CACHE_DIR
    raw_dir = config.DOCKING_RAW_CACHE_DIR
    lig_dir = config.DOCKING_LIGAND_CACHE_DIR
    rec_dir = config.DOCKING_RECEPTOR_CACHE_DIR
    # conftest 沙箱把四个目录指为临时目录下的平级路径,此处按需建目录
    for d in (cache_dir, raw_dir, lig_dir, rec_dir):
        os.makedirs(d)

    def _aged(path: str, days: float) -> None:
        with open(path, "w", encoding="utf-8") as f:
            f.write("x")
        old = time.time() - days * 86400 - 60
        os.utime(path, (old, old))

    old_raw = os.path.join(raw_dir, "vina_out_old.pdbqt")
    old_lig = os.path.join(lig_dir, "ligand_old.pdbqt")
    new_lig = os.path.join(lig_dir, "ligand_new.pdbqt")
    new_rec = os.path.join(rec_dir, "receptor_fresh.pdbqt")
    uniprot_map = os.path.join(cache_dir, "uniprot_gene_map.json")
    threedmol = os.path.join(cache_dir, "3Dmol-min.js")
    _aged(old_raw, 3)        # raw 超过 1 天 → 删
    _aged(old_lig, 9)        # ligand 超过 7 天 → 删
    _aged(new_lig, 3)        # ligand 未过期 → 留
    _aged(new_rec, 0)        # receptor 新文件 → 留
    _aged(uniprot_map, 30)   # 根目录长期缓存(故意做旧) → 必须豁免
    _aged(threedmol, 30)     # 同上
    removed = dock.clean_docking_cache()
    assert removed == 2, removed
    assert not os.path.exists(old_raw) and not os.path.exists(old_lig)
    assert os.path.exists(new_lig) and os.path.exists(new_rec)
    assert os.path.exists(uniprot_map), "UniProt 映射缓存不得被 TTL 清理"
    assert os.path.exists(threedmol), "3Dmol.js 库缓存不得被 TTL 清理"
    # ttl_days=0 的"全清"同样只清三个子目录,根目录长期缓存仍豁免
    removed0 = dock.clean_docking_cache(ttl_days=0)
    assert removed0 == 2, removed0  # new_lig / new_rec 此时全部过期
    assert os.path.exists(uniprot_map) and os.path.exists(threedmol)
