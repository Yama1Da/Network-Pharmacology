"""分子对接的结构文件处理与相互作用分析。

PDB/pdbqt 原子记录解析、共晶配体与口袋残基定位、受体 PDB 清洗、
obabel 格式转换封装,以及配体-受体相互作用(氢键/盐桥/疏水接触)检测。
供 core.docking(对接执行链)与 api 对接页共用。
"""

import os
from typing import Iterator, NamedTuple


# PDB 结构中的水分子残基名,不参与结合位点定位与共晶配体识别
_WATER_RESNAMES = ("HOH", "H2O", "WAT", "DOD")


def _strip_receptor_pdb(pdb_path: str, out_path: str) -> None:
    """清洗受体 PDB:只保留蛋白原子,去掉水、配体、离子和重复构象,便于后续对接。"""
    seen = set()
    lines = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "ATOM":
                if line[16] not in (" ", "A"):
                    continue  # 跳过 altloc 冗余构象
                key = (line[21], line[22:26].strip(), line[26], line[17:20],
                       line[12:16].strip())  # 键含插入码(line[26]):82 与 82A 是不同残基
                if key in seen:
                    continue
                seen.add(key)
                lines.append(line.rstrip("\n"))
            elif rec in ("TER", "END"):
                lines.append(line.rstrip("\n"))
    if not lines:
        raise RuntimeError(f"no ATOM records in receptor PDB: {pdb_path}")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


class _AtomRecord(NamedTuple):
    """PDB 原子记录的公共字段(仅产出坐标解析成功的行)。"""
    rec: str      # 记录类型:"ATOM" / "HETATM"
    resn: str     # 残基名
    chain: str    # 链 ID(单字符,原样保留)
    resi: str     # 残基编号
    x: float
    y: float
    z: float
    altloc: str   # 冗余构象标志(单字符)
    atom: str     # 原子名
    elem: str     # 元素(PDB 元素列; pdbqt 为 AutoDock 原子类型,如 OA/NA/A)


def _guess_element(atom_name: str) -> str:
    """元素列缺失时从原子名推断元素;仅个别两字符元素需特殊区分,相互作用判定够用。"""
    name = atom_name.strip()
    two = name[:2].upper()
    if len(name) >= 2 and two in ("CL", "BR", "ZN", "FE", "MG", "MN"):
        return two
    return name[:1].upper() if name[:1].isalpha() else ""


def _iter_atom_records(pdb_path: str) -> Iterator[_AtomRecord]:
    """逐行流式读取 PDB 的 ATOM/HETATM 记录并解析坐标,逐条产出公共字段;
    非原子记录与坐标残缺的坏行直接跳过(配体 pdbqt 文件同样按此规则解析)。
    """
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            # 元素列(PDB)/原子类型列(pdbqt)都在第 77 列之后,取该区域最后一个字母词
            tail = line[76:].split()
            elem = tail[-1] if tail and tail[-1].isalpha() else _guess_element(line[12:16].strip())
            yield _AtomRecord(rec, line[17:20].strip(), line[21], line[22:26].strip(),
                              x, y, z, line[16], line[12:16].strip(), elem)


def _pdb_coords(pdb_path: str) -> tuple[list, list]:
    """解析 PDB 文件中的原子坐标(单位 Å),返回 (蛋白原子坐标列表, 非水其他分子坐标列表)。"""
    protein, het = [], []
    for r in _iter_atom_records(pdb_path):
        if r.resn in _WATER_RESNAMES:
            continue  # 水分子不参与结合位点定位
        (het if r.rec == "HETATM" else protein).append((r.x, r.y, r.z))
    return protein, het


def parse_hetatm(pdb_path: str) -> list[tuple[str, str, str, tuple]]:
    """找出 PDB 中的共晶配体/辅因子(水除外),返回 (残基名, 链, 残基编号, (x, y, z)) 列表,供对接盒推断与 3D 预览使用。"""
    het = []
    for r in _iter_atom_records(pdb_path):
        if r.rec == "HETATM" and r.resn not in _WATER_RESNAMES:
            het.append((r.resn, r.chain, r.resi, (r.x, r.y, r.z)))
    return het


def residues_in_box(pdb_path: str, box: dict) -> list[str]:
    """列出落在对接盒内的蛋白残基,用于检查盒子是否覆盖了结合口袋。"""
    c, s = box["center"], box["size"]
    lo = (c[0] - s[0] / 2, c[1] - s[1] / 2, c[2] - s[2] / 2)
    hi = (c[0] + s[0] / 2, c[1] + s[1] / 2, c[2] + s[2] / 2)
    res = set()
    for r in _iter_atom_records(pdb_path):
        if r.rec != "ATOM":
            continue
        if lo[0] <= r.x <= hi[0] and lo[1] <= r.y <= hi[1] and lo[2] <= r.z <= hi[2]:
            res.add(f"{r.chain}{r.resi}:{r.resn}")
    return sorted(res)


def residues_near_ligand(receptor_pdb: str, ligand_pdbqt: str, cutoff: float = 4.0) -> list[dict]:
    """找出对接构象(配体)周围 cutoff Å 以内的结合位点残基,只统计受体 ATOM 记录(水与共晶配体不参与)。
    返回按链与残基编号排序、去重后的残基列表,每项含链号 chain、残基编号 resi、残基名 resn 及标注名称用的代表原子坐标 pos(CA 原子优先)。
    """
    import re

    lig = [(r.x, r.y, r.z) for r in _iter_atom_records(ligand_pdbqt)]
    if not lig:
        return []

    cutoff2 = cutoff * cutoff
    found: dict[tuple, dict] = {}
    for r in _iter_atom_records(receptor_pdb):
        if r.rec != "ATOM" or r.altloc not in ("", " ", "A"):
            continue  # 非蛋白记录与 altloc 冗余构象不参与位点判定
        xyz = (r.x, r.y, r.z)
        key = (r.chain.strip(), r.resi)
        info = found.get(key)
        if info is None:
            info = {"resn": r.resn, "ca": None, "pos": xyz, "near": False}
            found[key] = info
        if info["ca"] is None and r.atom == "CA":
            info["ca"] = xyz
        if not info["near"]:
            for lx, ly, lz in lig:
                dx, dy, dz = xyz[0] - lx, xyz[1] - ly, xyz[2] - lz
                if (abs(dx) <= cutoff and abs(dy) <= cutoff and abs(dz) <= cutoff
                        and dx * dx + dy * dy + dz * dz <= cutoff2):
                    info["near"] = True
                    break
    residues = [
        {"chain": chain, "resi": resi, "resn": info["resn"],
         "pos": info["ca"] or info["pos"]}
        for (chain, resi), info in found.items() if info["near"]
    ]
    residues.sort(key=lambda r: (r["chain"], int(re.match(r"\d+", r["resi"]).group())
                                 if re.match(r"\d+", r["resi"]) else 0, r["resi"]))
    return residues


# 相互作用几何判据(自实现,参数引自文献):
# 氢键 N/O/S 重原子间距 ≤3.5 Å(Baker & Hubbard 1984;IUPAC 2011 建议)
# 盐桥 反向带电基团间距 ≤4.0 Å(Barlow & Thornton 1983)
# 疏水接触 C–C 间距 ≤4.0 Å(常用接触半径)
# 无氢坐标下按重原子距离判定,不区分供体/受体方向;π-π 堆积不检测。
_HBOND_CUTOFF = 3.5
_SALT_CUTOFF = 4.0
_HYDROPHOBIC_CUTOFF = 4.0

# 标准氨基酸侧链带电基团原子(盐桥判定用)
_SALT_ACIDIC = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
_SALT_BASIC = {"LYS": ("NZ",), "ARG": ("NH1", "NH2"), "HIS": ("ND1", "NE2")}

# pdbqt 原子类型(后缀 A 是芳香/极性标记)与 PDB 元素符号 → 相互作用判定用元素
_INTERACTION_ELEMS = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
                      "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S"}


def _interaction_elem(elem: str) -> str:
    """把元素列(PDB 元素或 pdbqt 原子类型)归一成 N/O/S/C,其余元素(氢/金属/卤素)返回空串。"""
    return _INTERACTION_ELEMS.get(elem.upper(), "")


def interactions_near_ligand(receptor_pdb: str, ligand_pdbqt: str) -> list[dict]:
    """检测对接构象与受体之间的氢键、盐桥与疏水接触,返回按类型与距离排序的原子对列表。

    判据与文献出处见函数上方常量说明;水、共晶配体(HETATM)与 altloc 冗余构象不参与。
    不解析 pdbqt 部分电荷,盐桥按"受体酸性 O↔配体 N、受体碱性 N↔配体 O"近似判定;
    同一原子对同时满足多种类型时按 盐桥>氢键>疏水 只记一种。
    每项含类型 kind(hbond/salt/hydro)、残基 chain/resi/resn、原子名 res_atom/lig_atom、
    距离 dist(Å)与画虚线用的两端坐标 pos1/pos2。
    """
    import re

    def _resi_num(resi: str) -> int:
        m = re.match(r"\d+", resi)
        return int(m.group()) if m else 0

    lig_atoms = [(r, _interaction_elem(r.elem)) for r in _iter_atom_records(ligand_pdbqt)]
    lig_atoms = [(r, e) for r, e in lig_atoms if e]
    rec_atoms = [(r, _interaction_elem(r.elem)) for r in _iter_atom_records(receptor_pdb)
                 if r.rec == "ATOM" and r.altloc in ("", " ", "A")]
    rec_atoms = [(r, e) for r, e in rec_atoms if e]

    def _entry(kind: str, r, l, d2: float) -> dict:
        return {"kind": kind, "chain": r.chain, "resi": r.resi, "resn": r.resn,
                "res_atom": r.atom, "lig_atom": l.atom, "dist": round(d2 ** 0.5, 2),
                "pos1": (r.x, r.y, r.z), "pos2": (l.x, l.y, l.z)}

    def _pair_d2(r, l) -> float:
        dx, dy, dz = r.x - l.x, r.y - l.y, r.z - l.z
        return dx * dx + dy * dy + dz * dz

    out: list[dict] = []
    # 1) 盐桥:酸性 O↔配体 N,碱性 N↔配体 O(近似,不区分质子化状态)
    salted: set = set()  # (chain, resi, res_atom, 配体原子序号) 已归入盐桥的原子对
    for r, _ in rec_atoms:
        if r.atom in _SALT_ACIDIC.get(r.resn, ()):
            targets = ("N",)
        elif r.atom in _SALT_BASIC.get(r.resn, ()):
            targets = ("O",)
        else:
            continue
        for i, (l, le) in enumerate(lig_atoms):
            if le not in targets:
                continue
            d2 = _pair_d2(r, l)
            if d2 <= _SALT_CUTOFF ** 2:
                salted.add((r.chain, r.resi, r.atom, i))
                out.append(_entry("salt", r, l, d2))
    # 2) 氢键:双方均为 N/O/S 重原子且距离足够近;已记为盐桥的原子对不重复记
    polar = ("N", "O", "S")
    for r, e in rec_atoms:
        if e not in polar:
            continue
        for i, (l, le) in enumerate(lig_atoms):
            if le not in polar or (r.chain, r.resi, r.atom, i) in salted:
                continue
            d2 = _pair_d2(r, l)
            if d2 <= _HBOND_CUTOFF ** 2:
                out.append(_entry("hbond", r, l, d2))
    # 3) 疏水接触:C–C
    for r, e in rec_atoms:
        if e != "C":
            continue
        for i, (l, le) in enumerate(lig_atoms):
            if le != "C":
                continue
            d2 = _pair_d2(r, l)
            if d2 <= _HYDROPHOBIC_CUTOFF ** 2:
                out.append(_entry("hydro", r, l, d2))

    rank = {"salt": 0, "hbond": 1, "hydro": 2}
    out.sort(key=lambda it: (rank[it["kind"]], it["dist"], it["chain"],
                             _resi_num(it["resi"]), it["resi"]))
    return out


def cif_atom_site_to_pdb(cif_path: str, out_path: str) -> int:
    """把 RCSB mmCIF 的 _atom_site 循环转写成 PDB 定长原子行,返回写出的原子行数。

    供 receptor 兜底下载用:部分条目(如 8VHE)RCSB 不提供旧版 .pdb,只有 mmCIF。
    不走 obabel——它的 mmCIF 读取会把配体/离子混写成 ATOM(污染受体)且丢失
    HETATM 语义;这里按字段逐行转写,ATOM/HETATM、altloc、插入码、链号原样保留,
    坐标列位与原生 PDB 完全一致,下游清洗/盒子推断/相互作用链路无需区分来源。
    条目一旦发布内容冻结,输入只假设为 RCSB 导出的规范 mmCIF。
    """
    import re as _re

    def _tokens(line: str) -> list[str]:
        out, i, n = [], 0, len(line)
        while i < n:
            ch = line[i]
            if ch in " \t":
                i += 1
            elif ch in "'\"":
                j = line.find(ch, i + 1)
                j = n if j < 0 else j
                out.append(line[i + 1:j])
                i = j + 1
            else:
                j = i
                while j < n and line[j] not in " \t":
                    j += 1
                out.append(line[i:j])
                i = j
        return out

    with open(cif_path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    tags: list[str] = []
    rows: list[list[str]] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            i += 1
            block: list[str] = []
            while i < len(lines) and lines[i].strip().startswith("_"):
                block.append(lines[i].strip().split()[0])
                i += 1
            if any(t.startswith("_atom_site.") for t in block):
                tags = block
                while i < len(lines):
                    s = lines[i].strip()
                    if not s or s.startswith("#") or s == "loop_" or s.startswith("_"):
                        break
                    rows.append(_tokens(s))
                    i += 1
                break
        else:
            i += 1
    if not tags:
        raise RuntimeError(f"mmCIF 中未找到 _atom_site 数据:{cif_path}")

    def idx(*names: str) -> int:
        for n in names:
            if f"_atom_site.{n}" in tags:
                return tags.index(f"_atom_site.{n}")
        return -1

    i_group = idx("group_PDB")
    i_name = idx("auth_atom_id", "label_atom_id")
    i_alt = idx("label_alt_id")
    i_comp = idx("auth_comp_id", "label_comp_id")
    i_chain = idx("auth_asym_id", "label_asym_id")
    i_seq = idx("auth_seq_id", "label_seq_id")
    i_ins = idx("pdbx_PDB_ins_code")
    i_x, i_y, i_z = idx("Cartn_x"), idx("Cartn_y"), idx("Cartn_z")
    i_occ = idx("occupancy")
    i_b = idx("B_iso_or_equiv")
    i_elem = idx("type_symbol")
    if min(i_group, i_name, i_comp, i_chain, i_seq, i_x, i_y, i_z) < 0:
        raise RuntimeError(f"mmCIF _atom_site 缺少必要字段:{cif_path}")

    def val(row, j: int, default: str) -> str:
        v = row[j] if 0 <= j < len(row) else "."
        return default if v in (".", "?") else v

    out_lines: list[str] = []
    serial = 0
    for row in rows:
        try:
            x, y, z = float(row[i_x]), float(row[i_y]), float(row[i_z])
        except (ValueError, IndexError):
            continue  # 坐标残缺的坏行直接跳过
        serial += 1
        record = val(row, i_group, "ATOM")
        record = record if record in ("ATOM", "HETATM") else "ATOM"
        name = val(row, i_name, "")[:4]
        alt = val(row, i_alt, " ")[:1]
        comp = val(row, i_comp, "UNK")[:3]
        chain = val(row, i_chain, " ")[:1]
        m = _re.match(r"-?\d+", val(row, i_seq, "1"))
        seqn = int(m.group()) if m else 1
        icode = val(row, i_ins, " ")[:1]
        occ = float(val(row, i_occ, "1.00"))
        b = float(val(row, i_b, "0.00"))
        elem = val(row, i_elem, name[:1].upper())[:2]
        # 列位:record1-6 serial7-11 name13-16 alt17 resn18-20 chain22 resi23-26
        # icode27 x31-38 y39-46 z47-54 occ55-60 b61-66 elem77-78
        out_lines.append(
            f"{record:<6s}{serial % 100000:>5d} {name:<4s}{alt:1s}{comp:>3s} {chain:1s}"
            f"{seqn:>4d}{icode:1s}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{b:>6.2f}"
            f"          {elem:>2s}")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines) + ("\n" if out_lines else ""))
    return len(out_lines)


def _run_obabel(source: str, out_path: str, *flags: str) -> None:
    """以子进程调用 obabel 命令行(与 Vina 同类的外部程序,不依赖其 Python 绑定)。

    obabel 对读不到输入等情况只打错误日志、仍返回退出码 0 并留下空输出文件,
    因此除退出码外还须校验输出文件存在且非空。找不到可执行文件、转换失败或超时报错。
    """
    import subprocess

    from WY import config

    if not os.path.isfile(config.OBABEL_EXE):
        raise RuntimeError(
            f"未找到 Open Babel 可执行文件（obabel）：{config.OBABEL_EXE or '未探测到'}。"
            "请设置环境变量 OBABEL_DIR 指向 obabel 所在目录，或将 obabel 加入 PATH。")
    try:
        proc = subprocess.run(
            [config.OBABEL_EXE, source, "-O", out_path, *flags], capture_output=True,
            text=True, encoding="utf-8", errors="replace", timeout=config.OBABEL_TIMEOUT)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"obabel 转换超时（>{config.OBABEL_TIMEOUT} 秒）") from exc
    if proc.returncode != 0 or not os.path.isfile(out_path) or os.path.getsize(out_path) == 0:
        lines = (proc.stderr or proc.stdout or "").strip().splitlines()
        detail = lines[-1] if lines else f"exit code {proc.returncode}"
        raise RuntimeError(f"obabel 转换失败（{source}）：{detail}")
