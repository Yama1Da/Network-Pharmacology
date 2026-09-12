"""分子对接模块(AutoDock Vina)。

流程:把受体蛋白和配体小分子转成对接所需格式,确定结合口袋范围,
运行 Vina 对接并解析打分。中间文件会缓存一段时间以便复用,对接结果只归档不缓存。
"""

import hashlib
import logging
import os
from datetime import datetime
from typing import Iterator, NamedTuple

import pandas as pd
import requests

from .caching import load_json_file, save_json_file, save_text_atomic
from i18n import t

# PDB 结构中的水分子残基名,不参与结合位点定位与共晶配体识别
_WATER_RESNAMES = ("HOH", "H2O", "WAT", "DOD")


def _docking_fingerprint(path: str) -> str:
    """计算文件内容的指纹编号,用于判断该文件之前是否处理过、能否直接复用结果。"""
    try:
        with open(path, "rb") as f:
            return hashlib.md5(f.read()).hexdigest()[:12]
    except OSError:
        return ""


def _unique_path(path: str) -> str:
    """防止覆盖旧结果:目标路径已存在时改用带时间戳的新文件名。"""
    if not os.path.exists(path):
        return path
    stem, ext = os.path.splitext(path)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = f"{stem}_{ts}{ext}"
    n = 1
    while os.path.exists(candidate):
        candidate = f"{stem}_{ts}_{n}{ext}"
        n += 1
    return candidate


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


def download_pdb(pdb_id: str) -> str:
    """按 PDB 编号下载受体结构文件;之前下载过则直接复用本地文件。返回本地路径。"""
    import config

    pdb_id = pdb_id.strip().lower()
    os.makedirs(config.DOCKING_RAW_CACHE_DIR, exist_ok=True)
    path = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{pdb_id}.pdb")
    if os.path.exists(path):
        return path
    url = config.PDB_DOWNLOAD_URL.format(pdb_id=pdb_id)
    resp = requests.get(url, timeout=config.API_TIMEOUT)
    resp.raise_for_status()
    # 原子写:同一 PDB 被两个会话并发处理时,半截文件不会被"存在即复用"判定永久信任
    save_text_atomic(path, resp.text)
    return path


# -------------------- 配体链 --------------------


def resolve_cid(query: str) -> int | None:
    """把药物名换成 PubChem 化合物编号(CID);输入本身就是编号时直接使用。"""
    import pubchempy

    query = query.strip()
    if query.isdigit():
        return int(query)
    try:
        comps = pubchempy.get_compounds(query, "name")
        return comps[0].cid if comps else None
    except Exception:
        return None


def smiles_to_sdf_3d(smiles: str, out_path: str, seed: int = 0, n_confs: int = 20) -> str:
    """由 SMILES 生成配体 3D 结构:生成多个构象,取能量最低的一个存为 SDF 文件。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(t("dock.err_invalid_smiles", s=smiles))
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if not cids:
        raise RuntimeError(t("dock.err_rdkit_3d"))
    best_cid, best_e = None, None
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, confId=cid)
        if ff is None:
            continue
        ff.Minimize()
        e = ff.CalcEnergy()
        if best_e is None or e < best_e:
            best_e, best_cid = e, cid
    Chem.MolToMolFile(mol, out_path, confId=int(best_cid if best_cid is not None else cids[0]))
    return out_path


def fetch_ligand_sdf_3d(cid: int, out_dir: str | None = None) -> str:
    """从 PubChem 下载配体的现成 3D 结构(SDF);无此化合物或无 3D 结构时返回空字符串。
    同一化合物下载过一次后直接复用缓存。
    """
    import config
    import pubchempy

    out_dir = out_dir or config.DOCKING_RAW_CACHE_DIR
    os.makedirs(out_dir, exist_ok=True)
    # 同一化合物固定用同一个缓存文件名:已下载过就直接复用,不重复下载。
    path = os.path.join(out_dir, f"ligand_{cid}.sdf")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        logging.info("Reusing cached 3D SDF for CID %d: %s", cid, path)
        return path
    sdf_text = pubchempy.get_sdf(cid, namespace="cid", record_type="3d")
    if not sdf_text:
        logging.warning("pubchempy: no 3D SDF for CID %s", cid)
        return ""
    save_text_atomic(path, sdf_text)
    logging.info("Fetched 3D SDF for CID %d -> %s", cid, path)
    return path


def prepare_ligand_pdbqt(sdf_path: str, out_path: str | None = None) -> str:
    """把配体 3D SDF 转成对接所需的 pdbqt 格式,转换结果可复用;读取失败或缺少 3D 坐标时直接报错。"""
    import config
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem

    mol = None
    try:
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False)  # 标准方式读取,保证键序和电荷正确
    except Exception:
        logging.warning("ligand SDF sanitize failed, fallback to raw parse: %s", sdf_path)
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"rdkit failed to read ligand SDF: {sdf_path}")
    if mol.GetNumConformers() == 0:
        raise RuntimeError(f"ligand SDF has no 3D conformer: {sdf_path}（需 record_type=3d）")
    mol = Chem.AddHs(mol)

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
    if not is_ok:
        raise RuntimeError(f"meeko failed to write ligand pdbqt: {error_msg}")

    if not out_path:
        fingerprint = _docking_fingerprint(sdf_path)
        out_path = os.path.join(config.DOCKING_LIGAND_CACHE_DIR, f"ligand_{fingerprint}.pdbqt")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(pdbqt_string)
    logging.info("Prepared ligand pdbqt: %s", out_path)
    return out_path


# -------------------- 受体链 --------------------


# RCSB 批量检索:逐基因搜索(召回与旧版一致)按 4 路并发执行,元数据阶段
# 合并成少数几次 GraphQL 批量请求(实测把每基因最多 20 次元数据请求压到
# 全局约 1 次/100 个结构)。曾尝试把多个基因合并进一次 OR 搜索,实测相关度
# 排序会把部分基因的结果挤出截断上限(如 TP53 完全丢失),故不采用。
_RCSB_GRAPHQL_CHUNK = 100
_RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
_RCSB_BATCH_WORKERS = 4
_RCSB_HEADERS = {"User-Agent": "NetworkPharmacologyPlatform/1.0"}

_RCSB_ENTRY_QUERY = """
query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl { method }
    struct { title }
    rcsb_entry_info { resolution_combined }
  }
}
"""


def _rcsb_search_entry_ids(gene_symbol: str, max_hits: int) -> list[str]:
    """在 RCSB 中用基因符号精确检索结构,返回 entry 编号列表(最多 max_hits 条)。
    与旧的逐基因检索完全同一查询,保证召回一致;失败向上抛错。
    """
    import config
    from .http_client import request_with_retry

    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": f'"{gene_symbol.strip().upper()}"'},
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    resp = request_with_retry(
        _RCSB_SEARCH_URL, json=query,
        max_retries=config.RCSB_SEARCH_MAX_RETRIES, timeout=config.API_TIMEOUT,
        label=f"RCSB search for {gene_symbol}", throttle_key="rcsb_search",
        headers=_RCSB_HEADERS, raise_on_exhausted=True, use_cache=True,
        no_retry_status=(204,),
    )
    # RCSB 官方行为:查询成功但 0 命中时返回 204 No Content(空响应体),
    # 不能按 JSON 解析——按"无实验结构"处理,返回空列表
    if resp is None or resp.status_code == 204 or not resp.content.strip():
        logging.info("RCSB search: no experimental structures for %s", gene_symbol)
        return []
    return [hit.get("identifier", "") for hit in resp.json().get("result_set", [])
            if hit.get("identifier")][:max_hits]


def _rcsb_graphql_entries(entry_ids: list[str]) -> dict[str, dict]:
    """按 entry 编号批量取分辨率/实验方法/标题元数据;返回 {entry号: 信息字典}。
    一次请求即可覆盖最多 100 个 entry;GraphQL 层面的软错误记录日志后按缺数据处理,不中断整体。
    """
    import config
    from .http_client import request_with_retry

    out: dict[str, dict] = {}
    for start in range(0, len(entry_ids), _RCSB_GRAPHQL_CHUNK):
        chunk = entry_ids[start:start + _RCSB_GRAPHQL_CHUNK]
        resp = request_with_retry(
            _RCSB_GRAPHQL_URL, json={"query": _RCSB_ENTRY_QUERY,
                                     "variables": {"ids": chunk}},
            max_retries=config.RCSB_SEARCH_MAX_RETRIES, timeout=config.API_TIMEOUT,
            label=f"RCSB GraphQL batch ({len(chunk)} entries)", throttle_key="rcsb_data",
            headers=_RCSB_HEADERS, raise_on_exhausted=True, use_cache=True,
        )
        payload = resp.json()
        if payload.get("errors"):
            logging.warning("RCSB GraphQL errors on %d entries: %s",
                            len(chunk), payload["errors"][:1])
        for e in payload.get("data", {}).get("entries") or []:
            eid = e and e.get("rcsb_id")
            if not eid:
                continue
            res = (e.get("rcsb_entry_info") or {}).get("resolution_combined")
            res = res[0] if isinstance(res, list) and res else res
            out[eid] = {
                "resolution": res,
                "method": ", ".join(m.get("method", "") for m in e.get("exptl") or []),
                "title": (e.get("struct") or {}).get("title", "") or "",
            }
    return out


def _friendly_rcsb_error(gene_symbol: str, exc: Exception) -> RuntimeError:
    """把 RCSB 检索的底层网络异常翻译成普通用户可读的失败提示(原始异常已在日志中)。"""
    if isinstance(exc, requests.exceptions.Timeout):
        hint = "RCSB 检索请求超时，请稍后重试"
    elif isinstance(exc, requests.exceptions.ConnectionError):
        hint = "无法连接 RCSB 服务器，请检查网络连接后重试"
    elif isinstance(exc, requests.exceptions.HTTPError):
        code = getattr(exc.response, "status_code", None)
        if code == 429:
            hint = "RCSB 检索请求过于频繁被限流，请等待约 1 分钟后再试"
        elif code and code >= 500:
            hint = "RCSB 服务器暂时故障，请稍后重试"
        else:
            hint = f"RCSB 检索返回异常状态（HTTP {code}），请稍后重试"
    else:
        hint = "RCSB 检索失败，请稍后重试"
    return RuntimeError(f"RCSB 搜索失败（{gene_symbol}）：{hint}")


def fetch_receptor_candidates_batch(genes: list[str], max_hits: int = 20,
                                    allow_alphafold: bool = False, taxon: int = 9606,
                                    progress_cb=None,
                                    raise_on_search_failure: bool = False) -> dict[str, pd.DataFrame]:
    """批量查询多个基因的 PDB 候选结构,返回 {基因符号: 候选表};progress_cb(0~1) 可选,用于页面进度条。
    搜索与旧的逐基因检索同一查询保证召回一致、按 4 路并发,元数据合并为 GraphQL 批量请求使总请求数与基因数解耦;查无实验结构的基因在 allow_alphafold 开启时追加 AlphaFold 预测结构兜底。
    raise_on_search_failure=True 时检索服务不可达直接抛出底层异常(单基因入口用,
    由包装层翻译成可读错误);默认 False 按无结构处理,不让整批报废。
    """
    import config

    uniq = sorted({g.strip().upper() for g in genes if g and g.strip()})
    empty = pd.DataFrame(columns=["pdb_id", "resolution", "method", "title", "source"])
    results: dict[str, pd.DataFrame] = {g: empty.copy() for g in uniq}
    if not uniq:
        return results

    # 1) 并发逐基因搜索(与旧版同一查询,召回一致)
    from concurrent.futures import ThreadPoolExecutor

    gene_ids: dict[str, list[str]] = {}

    def _search(g: str) -> tuple[str, list[str]]:
        try:
            return g, _rcsb_search_entry_ids(g, max_hits)
        except Exception:
            if raise_on_search_failure:
                raise
            # 单基因检索失败只影响该基因:记日志、按无结构处理,不让整批报废;
            # 无结构基因仍可被后续的 AlphaFold 兜底接住
            logging.exception("RCSB batch search failed for gene %s; treating as no-hit", g)
            return g, []

    with ThreadPoolExecutor(max_workers=_RCSB_BATCH_WORKERS,
                            thread_name_prefix="rcsb-search") as pool:
        for i, (g, ids) in enumerate(pool.map(_search, uniq)):
            gene_ids[g] = ids
            if progress_cb:
                progress_cb((i + 1) / len(uniq) * 0.7)

    all_ids = sorted({eid for ids in gene_ids.values() for eid in ids})
    meta: dict[str, dict] = {}
    gql_chunks = [all_ids[i:i + _RCSB_GRAPHQL_CHUNK]
                  for i in range(0, len(all_ids), _RCSB_GRAPHQL_CHUNK)]
    for idx, chunk in enumerate(gql_chunks):
        meta.update(_rcsb_graphql_entries(chunk))
        if progress_cb:
            progress_cb(0.7 + (idx + 1) / max(len(gql_chunks), 1) * 0.3)

    # 3) 按基因组装结果表(与旧版字段完全一致;个别 entry 元数据缺失按缺信息展示)
    for g, ids in gene_ids.items():
        rows = [{"pdb_id": eid,
                 "resolution": meta.get(eid, {}).get("resolution"),
                 "method": meta.get(eid, {}).get("method", ""),
                 "title": meta.get(eid, {}).get("title", ""),
                 "source": "experiment"}
                for eid in ids]
        if rows:
            results[g] = pd.DataFrame(rows)

    if allow_alphafold and config.ALPHAFOLD_ENABLED:
        for g in uniq:
            if not results[g].empty:
                continue
            try:
                acc = fetch_uniprot_accession(g, taxon=taxon)
                if acc:
                    results[g] = pd.DataFrame([{
                        "pdb_id": f"AFDB:{acc}",
                        "resolution": None,
                        "method": f"AlphaFold predicted ({config.ALPHAFOLD_MODEL_VERSION})",
                        "title": f"AlphaFold {config.ALPHAFOLD_MODEL_VERSION} monomer model for {g}",
                        "source": "alphafold",
                    }])
            except Exception:
                # 预测结构查询失败不影响空表正常返回
                logging.exception("AlphaFold fallback lookup failed for gene %s", g)
    return results


def fetch_receptor_pdb_candidates(gene_symbol: str, max_hits: int = 50,
                                  allow_alphafold: bool = False,
                                  taxon: int = 9606) -> pd.DataFrame:
    """按基因符号在 PDB 数据库(RCSB)中查找该靶蛋白的候选结构(单基因入口,即批量版的单基因形式)。
    查无实验结构时返回空表;开启 AlphaFold 选项且查无实验结构时,追加一条 AlphaFold 预测结构作兜底。
    与批量版的差异:检索服务不可达(网络断开/限流/服务器故障)时抛出带中文提示的 RuntimeError
    供页面展示,而不是按"无结构"静默处理——单基因场景没有其它候选来源,静默会让用户误以为
    该蛋白确实没有可用结构。
    """
    gene = gene_symbol.strip().upper()
    try:
        results = fetch_receptor_candidates_batch(
            [gene], max_hits=max_hits, allow_alphafold=allow_alphafold, taxon=taxon,
            raise_on_search_failure=True)
    except Exception as exc:
        logging.exception("RCSB search failed for gene %s", gene_symbol)
        raise _friendly_rcsb_error(gene_symbol, exc)
    return results[gene]


def prepare_receptor_pdbqt(pdb_path: str, out_path: str | None = None) -> str:
    """把受体蛋白 PDB 转成对接所需的 pdbqt 格式:去水去杂原子、加氢(pH 7.4)、分配 AutoDock 原子类型,结果可复用。
    失败时报错,成功返回 pdbqt 文件路径。
    """
    import config
    from openbabel import openbabel as ob

    os.makedirs(config.DOCKING_RAW_CACHE_DIR, exist_ok=True)
    cleaned = os.path.join(
        config.DOCKING_RAW_CACHE_DIR,
        f"receptor_clean_{os.path.splitext(os.path.basename(pdb_path))[0]}.pdb",
    )
    try:
        _strip_receptor_pdb(pdb_path, cleaned)
        source = cleaned
    except Exception:
        logging.exception("receptor PDB cleaning failed, fallback to raw PDB")
        source = pdb_path

    if not out_path:
        fingerprint = _docking_fingerprint(source)
        out_path = os.path.join(config.DOCKING_RECEPTOR_CACHE_DIR, f"receptor_{fingerprint}.pdbqt")
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdb", "pdbqt")
    conv.SetOptions("r", ob.OBConversion.OUTOPTIONS)  # 受体模式:加氢、去水、AutoDock 原子类型
    mol = ob.OBMol()
    if not conv.ReadFile(mol, source):
        raise RuntimeError(f"openbabel failed to read receptor PDB: {source}")
    if not conv.WriteFile(mol, out_path):
        raise RuntimeError(f"openbabel failed to write receptor pdbqt: {out_path}")
    logging.info("Prepared receptor pdbqt: %s", out_path)
    return out_path


# -------------------- 对接盒 --------------------


def infer_box_from_cofactor(pdb_path: str, padding: float = 8.0) -> dict:
    """根据 PDB 中共晶配体/辅因子的位置推断对接盒(结合口袋范围);没有配体/辅因子时用整个蛋白的范围,外扩默认 8 Å。
    结果仅供参考,使用前需人工确认口袋位置。返回口袋中心和尺寸(单位 Å)。
    """
    protein, het = _pdb_coords(pdb_path)
    if not protein and not het:
        raise RuntimeError(f"no coordinates found in PDB: {pdb_path}")
    source = het if het else protein  # 优先配体/辅因子定位口袋
    xs, ys, zs = zip(*source)
    center = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    size = (
        max(xs) - min(xs) + 2 * padding,
        max(ys) - min(ys) + 2 * padding,
        max(zs) - min(zs) + 2 * padding,
    )
    return {
        "center": tuple(round(v, 4) for v in center),
        "size": tuple(round(v, 4) for v in size),
    }


# -------------------- 执行链 --------------------


def run_vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    box: dict,
    out_path: str | None = None,
) -> "subprocess.CompletedProcess":
    """调用 AutoDock Vina,在指定口袋范围内对受体和配体执行分子对接;超时或运行失败时报错,运行日志(含打分)保留在返回结果中。"""
    import subprocess
    import config

    if not os.path.isfile(config.VINA_EXE):
        raise RuntimeError(
            f"未找到 AutoDock Vina 可执行文件：{config.VINA_EXE}。"
            "请设置环境变量 VINA_DIR 指向 vina 所在目录，或将 vina 加入 PATH。")

    out_path = out_path or _unique_path(os.path.join(config.MD_RESULT_DIR, "vina_out.pdbqt"))
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    cx, cy, cz = box["center"]
    sx, sy, sz = box["size"]
    cmd = [
        config.VINA_EXE,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str(config.EXHAUSTIVENESS),
        "--num_modes", str(config.NUM_MODES),
        "--cpu", str(config.DOCKING_CPU),
        "--seed", str(config.DOCKING_SEED),
        "--out", out_path,
    ]
    logging.info("Running vina: %s", " ".join(cmd))
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=config.DOCKING_TIMEOUT,
            encoding="utf-8", errors="replace",
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"vina timed out after {config.DOCKING_TIMEOUT}s")
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "")[-500:]
        raise RuntimeError(f"vina failed (rc={proc.returncode}): {tail}")
    return proc


def _iter_vina_models(out_pdbqt: str) -> Iterator[list[str]]:
    """按 MODEL/ENDMDL 把 vina 多模型输出切成逐个构象的行块(含 MODEL/ENDMDL 行)。"""
    current: list[str] | None = None
    with open(out_pdbqt, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                current = [line]
            elif line.startswith("ENDMDL") and current is not None:
                current.append(line)
                yield current
                current = None
            elif current is not None:
                current.append(line)


def _split_models_python(out_pdbqt: str, out_dir: str) -> list[str]:
    """内置方式:按构象把 vina 多模型输出逐个拆成独立文件(拆分工具不可用时的备用方案)。"""
    stem = os.path.splitext(os.path.basename(out_pdbqt))[0]
    files = []
    for i, block in enumerate(_iter_vina_models(out_pdbqt), 1):
        path = os.path.join(out_dir, f"{stem}_{i}.pdbqt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("".join(block))
        files.append(path)
    return files


def split_vina(out_pdbqt: str, out_dir: str | None = None) -> list[str]:
    """把 vina 输出按结合模式拆成一个个独立的构象(pose)文件,返回各 pose 文件路径。
    拆分工具不可用时自动改用内置方式,不影响流程;两种方式都失败才报错。
    """
    import subprocess
    import glob
    import config

    out_dir = out_dir or os.path.dirname(out_pdbqt)
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(out_pdbqt))[0]

    split_exe = config.VINA_SPLIT_EXE
    if os.path.isfile(split_exe) and os.path.getsize(split_exe) > 0:
        cmd = [split_exe, "--input", out_pdbqt]
        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, cwd=out_dir,
                timeout=config.DOCKING_TIMEOUT, encoding="utf-8", errors="replace",
            )
            if proc.returncode == 0:
                files = sorted(glob.glob(os.path.join(out_dir, f"{stem}*")))
                if files:
                    return files
        except OSError:
            logging.warning("vina_split unavailable (%s), fallback to python split", split_exe)

    files = _split_models_python(out_pdbqt, out_dir)
    if not files:
        raise RuntimeError(f"failed to split vina output: {out_pdbqt}")
    return files


def extract_model(out_pdbqt: str, mode: int, out_path: str) -> str:
    """从 Vina 输出中取出第 mode 个对接构象,单独保存为文件。"""
    n = 0
    for i, block in enumerate(_iter_vina_models(out_pdbqt), 1):
        n = i
        if i == mode:
            with open(out_path, "w", encoding="utf-8") as f:
                f.write("".join(block))
            return out_path
    raise RuntimeError(t("dock.err_mode_range", mode=mode, n=n))


def parse_vina_scores(out_pdbqt: str) -> pd.DataFrame:
    """从 vina 输出中读取各结合模式的打分(结合能)与 RMSD;无可解析结果时返回空表。
    按结合能从低到高排序(数值越负结合越强),第 1 名即最佳构象。
    """
    modes = []
    for block in _iter_vina_models(out_pdbqt):
        info: dict = {}
        for line in block:
            if line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                try:
                    info["affinity"] = float(parts[3])
                    info["rmsd_lb"] = float(parts[4])
                    info["rmsd_ub"] = float(parts[5]) if len(parts) > 5 else None
                except (IndexError, ValueError):
                    # 个别构象的打分行数值残缺时跳过该行,不让整个解析中止
                    logging.warning("Malformed VINA RESULT line skipped: %r", line.strip())
                    info.clear()
        if info:  # 与原逻辑一致:没有打分行的构象不计入
            modes.append(info)
    df = pd.DataFrame(modes)
    if df.empty:
        return pd.DataFrame(columns=["mode", "affinity", "rmsd_lb", "rmsd_ub"])
    df.insert(0, "mode", range(1, len(df) + 1))
    return df.sort_values("affinity").reset_index(drop=True)


# -------------------- 元数据 / 归档 / 清理 --------------------


def write_docking_metadata(meta_path: str, fields: dict) -> str:
    """把本次对接的关键信息(药物、靶点、参数等)写成 # 开头的文本记录,便于追溯;同名时自动换新文件名,返回实际写入的文件路径。"""
    lines = [f"# {key}: {value}" for key, value in fields.items()]
    lines.append(f"# timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    path = _unique_path(meta_path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    logging.info("Wrote docking metadata: %s", path)
    return path


def star_docking_result(result_dir: str) -> str:
    """把本次对接成果目录完整复制到收藏区长期保存,不会被自动清理;原目录仍按版本控制自动管理,返回收藏副本路径。"""
    import shutil
    import config

    os.makedirs(config.MD_STAR_DIR, exist_ok=True)
    dest = _unique_path(os.path.join(config.MD_STAR_DIR, os.path.basename(result_dir)))
    shutil.copytree(result_dir, dest)
    logging.info("Starred docking result: %s -> %s", result_dir, dest)
    return dest


def prune_docking_results(result_dir: str, max_versions: int | None = None) -> list[str]:
    """成果版本控制:同一药物-疾病-靶点组合的结果目录超过保留上限时删除最旧的;返回本次删除的目录列表,无需删除时返回空列表。"""
    import re
    import shutil
    import glob
    import config

    max_versions = max_versions or config.DOCKING_MAX_VERSIONS
    parent = os.path.dirname(result_dir)
    basename = os.path.basename(result_dir)
    prefix = re.sub(r"_\d{8}_\d{6}$", "", basename)  # 剥离末尾时间戳段
    if prefix == basename:
        return []  # 目录名不符合命名约定，不做版本控制
    siblings = sorted(
        d for d in glob.glob(os.path.join(parent, prefix + "_*"))
        if os.path.isdir(d)
    )
    # 删除前确认目录确在结果目录内,防止误删其他位置的文件夹
    md_root = os.path.realpath(config.MD_RESULT_DIR)
    siblings = [d for d in siblings
                if os.path.realpath(d).startswith(md_root + os.sep)]
    removed = []
    while len(siblings) > max_versions:
        oldest = siblings.pop(0)
        shutil.rmtree(oldest, ignore_errors=True)
        removed.append(oldest)
        logging.info("Pruned old docking result: %s", oldest)
    return removed


def clean_docking_cache(ttl_days: int | None = None) -> int:
    """清理对接过程中的缓存中间文件:原始文件保留 1 天,加工产物保留 7 天;传 ttl_days 时统一按该天数清理(0 可全清),返回删除的文件数。
    只扫 raw/ligand/receptor 三个已知子目录;缓存根目录下的文件(UniProt 基因映射、
    3Dmol.js 库等按设计长期有效的缓存)一律不动,否则会被 TTL 误删后反复重下、浪费配额。"""
    import time
    import config

    if ttl_days is not None:
        cache_ttl = raw_ttl = ttl_days
    else:
        cache_ttl = config.DOCKING_CACHE_TTL_DAYS
        raw_ttl = config.DOCKING_RAW_TTL_DAYS
    now = time.time()
    removed = 0
    scan_dirs = (
        (config.DOCKING_RAW_CACHE_DIR, raw_ttl),
        (config.DOCKING_LIGAND_CACHE_DIR, cache_ttl),
        (config.DOCKING_RECEPTOR_CACHE_DIR, cache_ttl),
    )
    for scan_dir, ttl in scan_dirs:
        if not os.path.isdir(scan_dir):
            continue
        cutoff = now - ttl * 86400
        for root, _, files in os.walk(scan_dir):
            for name in files:
                path = os.path.join(root, name)
                try:
                    if os.path.getmtime(path) < cutoff:
                        os.remove(path)
                        removed += 1
                except OSError:
                    logging.warning("clean_docking_cache: cannot remove %s", path)
    if removed:
        logging.info("Docking cache cleaned: %d files removed", removed)
    return removed


# -------------------- AlphaFold 预测结构兜底 --------------------
# 靶点在 PDB 中没有实验结构时,改用 AlphaFold 预测模型作受体。
# 注意:预测模型没有共晶配体,口袋只能粗略估计,必须人工确认;打分仅作初筛参考。


def _gene_names_from_ebi_entry(entry: dict) -> set:
    """提取条目中记录的基因名(含别名),统一转为大写。"""
    names: set = set()
    for g in entry.get("gene") or []:
        name = (g.get("name") or {}).get("value")
        if name:
            names.add(str(name).upper())
        for s in g.get("synonyms") or []:
            v = (s or {}).get("value")
            if v:
                names.add(str(v).upper())
    return names


def _lookup_uniprot_ebi(gene_symbol: str, taxon: int) -> str:
    """在 EBI 数据库按基因符号查找已校验(reviewed)蛋白的条目号。
    接口按名称模糊匹配,易混入同源基因,拿到结果后需再按基因名精确筛选;只取主条目,跳过异构体(isoform)。
    """
    import config
    from .http_client import throttle

    throttle("ebi_proteins")
    resp = requests.get(
        config.EBI_PROTEINS_API_URL,
        params={"gene": gene_symbol, "taxid": taxon, "reviewed": "true", "format": "json"},
        timeout=config.API_TIMEOUT, headers=_RCSB_HEADERS,
    )
    resp.raise_for_status()
    for entry in resp.json():
        acc = str(entry.get("accession", ""))
        if "-" in acc:  # isoform（如 P04637-2）不作 AF 兜底
            continue
        if gene_symbol.strip().upper() in _gene_names_from_ebi_entry(entry):
            return acc
    return ""


def _lookup_uniprot_rest(gene_symbol: str, taxon: int) -> str:
    """在 UniProt 官方接口按基因符号精确查找已校验蛋白的条目号;查询失败时报错。"""
    import config
    from .http_client import throttle

    throttle("uniprot")
    query = f"gene_exact:{gene_symbol} AND organism_id:{taxon} AND reviewed:true"
    resp = requests.get(
        config.UNIPROT_SEARCH_URL,
        params={"query": query, "format": "tsv", "fields": "accession,id", "size": "10"},
        timeout=config.API_TIMEOUT, headers=_RCSB_HEADERS,
    )
    resp.raise_for_status()
    for line in resp.text.strip().splitlines()[1:]:  # 首行为表头
        acc = line.split("\t")[0].strip()
        if acc and "-" not in acc:  # 跳过 isoform
            return acc
    return ""


def fetch_uniprot_accession(gene_symbol: str, taxon: int = 9606) -> str:
    """基因符号 → UniProt 蛋白条目号,依次尝试多个数据源;全失败报错,确认查无条目返回空字符串。
    查询结果永久缓存,包括"确认查无"的情况;但网络失败不写入缓存,下次会重查。
    """
    import config

    gene = gene_symbol.strip().upper()
    cache_file = config.UNIPROT_LOOKUP_CACHE
    db = load_json_file(cache_file, {}) or {}
    key = f"{gene}_{taxon}"
    if key in db:
        return str(db[key])

    lookup = {"ebi": _lookup_uniprot_ebi, "uniprot": _lookup_uniprot_rest}
    acc, last_exc, definitive = "", None, False
    for src in config.UNIPROT_LOOKUP_SOURCES:
        fn = lookup.get(src)
        if fn is None:
            continue
        try:
            acc = fn(gene, taxon)
            definitive = True  # 请求有了确定结论(无论是否命中)
            if acc:
                break
        except Exception as exc:
            last_exc = exc
            logging.exception("UniProt lookup via %s failed for %s", src, gene)
    if not acc and not definitive and last_exc is not None:
        raise RuntimeError(f"UniProt lookup failed for {gene} (all sources): {last_exc}")

    db[key] = acc
    try:
        save_json_file(cache_file, db, sort_keys=True)
    except Exception:
        logging.exception("Failed to write UniProt lookup cache: %s", cache_file)
    if acc:
        logging.info("UniProt lookup: %s (taxon %s) -> %s", gene, taxon, acc)
    else:
        logging.warning("UniProt lookup: %s (taxon %s) has no reviewed entry", gene, taxon)
    return acc


def fetch_alphafold_pdb(accession: str) -> str:
    """下载 AlphaFold 预测结构 PDB 文件,下载过的直接复用缓存。
    下载内容异常(文件过小或不含原子坐标)时报错,不会保存坏文件。
    """
    import requests

    import config
    from .http_client import throttle

    acc = accession.strip()
    version = config.ALPHAFOLD_MODEL_VERSION
    os.makedirs(config.DOCKING_RECEPTOR_CACHE_DIR, exist_ok=True)
    path = os.path.join(config.DOCKING_RECEPTOR_CACHE_DIR, f"AF_{acc}_{version}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 1024:
        logging.info("Reusing cached AlphaFold model: %s", path)
        return path
    throttle("alphafold")
    url = config.ALPHAFOLD_FILE_URL_TEMPLATE.format(accession=acc, version=version)
    resp = requests.get(url, timeout=config.API_TIMEOUT,
                        headers=_RCSB_HEADERS)
    resp.raise_for_status()
    if len(resp.text) < 1024 or "\nATOM" not in resp.text:
        # 目录在函数开头已创建;这里只报错,不落盘
        raise RuntimeError(f"AlphaFold PDB unexpected content for {acc}: {url}")
    save_text_atomic(path, resp.text)
    logging.info("Fetched AlphaFold model %s -> %s", url, path)
    return path


def fetch_alphafold_receptor(gene_symbol: str, taxon: int = 9606) -> tuple:
    """由基因符号获取 AlphaFold 预测结构文件;查不到对应蛋白条目时报错。"""
    acc = fetch_uniprot_accession(gene_symbol, taxon=taxon)
    if not acc:
        raise RuntimeError(
            f"no reviewed UniProt entry for gene {gene_symbol} (taxon {taxon}) — "
            "cannot fetch AlphaFold model")
    return fetch_alphafold_pdb(acc), acc


def clamp_box(box: dict, max_size: float) -> dict:
    """对接盒边长超过上限时,以中心为准把尺寸压缩到上限。
    这只保证对接能运行,不保证口袋位置正确,需人工在 3D 预览中确认。
    """
    size = tuple(round(min(float(v), float(max_size)), 4) for v in box["size"])
    return {"center": tuple(round(float(v), 4) for v in box["center"]), "size": size}


# ---------- 交互式 3D 视图导出 ----------

# 3Dmol.js 库的下载地址:官方源在前,jsdelivr 镜像兜底
_3DMOL_JS_URLS = (
    "https://3dmol.org/build/3Dmol-min.js",
    "https://cdn.jsdelivr.net/npm/3dmol@2/build/3Dmol-min.js",
)

# 3Dmol.js 的许可文本(BSD-3-Clause,取自官方仓库 LICENSE):导出 HTML 内嵌或引用其库时
# 随文件一并提供,满足许可证"再分发须保留版权声明与许可文本"的要求。
_3DMOL_LICENSE_TEXT = """3Dmol.js incorporates code from GLmol, Three.js, and jQuery and is licensed under a BSD-3-Clause license.

Copyright (c) 2014, University of Pittsburgh
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."""

_3D_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
  body { margin: 0; font-family: "Segoe UI", Arial, sans-serif; }
  #cap { padding: 8px 14px; background: #f5f7fa; border-bottom: 1px solid #ddd;
         font-size: 13px; color: #333; }
  #lic { display: block; margin-top: 2px; font-size: 11px; color: #777; }
  #bar { padding: 4px 14px; border-bottom: 1px solid #eee; }
  #focus { cursor: pointer; }
  #viewer { width: 100vw; height: calc(100vh - 88px); position: relative; }
</style>
__LIB__
</head>
<body>
<div id="cap">__CAPTION__<span id="lic">__LICENSE_NOTE__</span></div>
<div id="bar"><button id="focus">__FOCUS_BTN__</button></div>
<div id="viewer"></div>
<script>
(function () {
  var receptor = __RECEPTOR__;
  var ligand = __LIGAND__;
  var het = __HET__;
  var pocket = __POCKET__;
  var interactions = __INTERACTIONS__;
  var viewer = $3Dmol.createViewer("viewer", { backgroundColor: "white" });
  viewer.addModel(receptor, "pdb");
  viewer.addModel(ligand, "pdb");
  viewer.setStyle({ model: 1 }, { stick: { colorscheme: "greenCarbon" } });
  if (het) {
    viewer.addModel(het, "pdb");
    viewer.setStyle({ model: 2 }, { stick: { colorscheme: "greyCarbon", opacity: 0.5, radius: 0.15 } });
  }
  // 口袋残基按链分组选择(无链号时不限定 chain)
  var byChain = {};
  (pocket || []).forEach(function (r) {
    (byChain[r.chain] = byChain[r.chain] || []).push(r.resi);
  });
  // 受体卡通样式统一由此函数设置:faded=true 时卡通整体半透明(虚化周围蛋白),
  // 突出口袋棍状、相互作用虚线与距离标注;虚化只是显示开关,信息不删减,双击复位即可还原。
  // 口袋残基的卡通与棍状在一条语句里显式声明完整样式,不依赖多次 setStyle 的合并行为
  function applyProteinStyles(faded) {
    var cartoon = faded ? { color: "spectrum", opacity: 0.3 } : { color: "spectrum" };
    viewer.setStyle({ model: 0 }, { cartoon: cartoon });
    Object.keys(byChain).forEach(function (ch) {
      var sel = { model: 0, resi: byChain[ch] };
      if (ch) { sel.chain = ch; }
      viewer.setStyle(sel, { cartoon: cartoon,
                             stick: { colorscheme: "orangeCarbon", radius: 0.2 } });
    });
    viewer.render();
  }
  applyProteinStyles(false);
  if (pocket && pocket.length) {
    pocket.forEach(function (r) {
      viewer.addLabel(r.resn + "-" + r.resi, {
        position: { x: r.pos[0], y: r.pos[1], z: r.pos[2] },
        fontSize: 11, fontColor: "black",
        showBackground: true, backgroundColor: "white",
        backgroundOpacity: 0.7
      });
    });
  }
  // 相互作用虚线圆柱:蓝=氢键 紫=盐桥 灰=疏水(用虚线圆柱而非 1px 线,白底下清晰);
  // 氢键与盐桥在中点标注距离
  var lineColors = { hbond: "#0f62fe", salt: "#c026d3", hydro: "#808080" };
  (interactions || []).forEach(function (it) {
    viewer.addCylinder({
      start: { x: it.pos1[0], y: it.pos1[1], z: it.pos1[2] },
      end: { x: it.pos2[0], y: it.pos2[1], z: it.pos2[2] },
      radius: 0.06, dashed: true, dashLength: 0.3, gapLength: 0.22,
      color: lineColors[it.kind] || "#777777"
    });
    if (it.kind !== "hydro") {
      viewer.addLabel(it.dist.toFixed(1) + "\\u00C5", {
        position: { x: (it.pos1[0] + it.pos2[0]) / 2,
                    y: (it.pos1[1] + it.pos2[1]) / 2,
                    z: (it.pos1[2] + it.pos2[2]) / 2 },
        fontSize: 10, fontColor: "#333333",
        showBackground: true, backgroundColor: "white",
        backgroundOpacity: 0.6
      });
    }
  });
  viewer.zoomTo();
  viewer.render();
  // 聚焦口袋:以配体为中心、拉远到能框住口袋残基,并把周围蛋白卡通虚化;
  // 双击复位为整体视图并取消虚化
  document.getElementById("focus").onclick = function () {
    applyProteinStyles(true);
    viewer.zoomTo({ model: 1 });
    viewer.zoom(0.55, 400);
    viewer.render();
  };
  document.getElementById("viewer").ondblclick = function () {
    applyProteinStyles(false);
    viewer.zoomTo();
    viewer.render();
  };
})();
</script>
<!--
__3DMOL_LICENSE__
-->
</body>
</html>
"""


def _load_3dmol_js() -> str | None:
    """读取 3Dmol.js 库源码:本地有缓存直接用,没有则下载一次并长期缓存。
    两个下载地址都失败时返回 None,HTML 将改用在线引用(打开时需联网)。
    """
    import config

    cache_path = os.path.join(config.DOCKING_CACHE_DIR, "3Dmol-min.js")
    if os.path.exists(cache_path) and os.path.getsize(cache_path) > 100_000:
        try:
            with open(cache_path, encoding="utf-8") as f:
                return f.read()
        except OSError:
            logging.exception("cached 3Dmol.js unreadable, re-downloading")
    for url in _3DMOL_JS_URLS:
        try:
            resp = requests.get(url, timeout=config.API_TIMEOUT,
                                headers=_RCSB_HEADERS)
            resp.raise_for_status()
            if "$3Dmol" not in resp.text or "</script" in resp.text:
                raise RuntimeError(f"unexpected content from {url}")
            os.makedirs(os.path.dirname(cache_path) or ".", exist_ok=True)
            # 原子写:下载中途被打断不会留下半截 JS 被上面"大小达标即信任"的判定永久复用
            save_text_atomic(cache_path, resp.text)
            logging.info("Cached 3Dmol.js (%d KB) from %s", len(resp.text) // 1024, url)
            return resp.text
        except Exception:
            logging.exception("3Dmol.js download failed: %s", url)
    return None


def pose_pdbqt_to_pdb(pose_pdbqt: str, out_path: str) -> str:
    """把对接构象从 pdbqt 转成通用 PDB 格式(供各查看器显示);失败时抛 RuntimeError。"""
    from openbabel import openbabel as ob

    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdbqt", "pdb")
    mol = ob.OBMol()
    if not (conv.ReadFile(mol, pose_pdbqt) and conv.WriteFile(mol, out_path)):
        raise RuntimeError(f"pose pdbqt -> pdb conversion failed: {pose_pdbqt}")
    return out_path


def build_interactive_3d_html(receptor_pdb: str, ligand_pdb: str, caption: str = "",
                              lib_js: str | None = None,
                              pocket: list[dict] | None = None,
                              het_pdb: str | None = None,
                              interactions: list[dict] | None = None,
                              focus_label: str = "Focus pocket") -> tuple[str, bool]:
    """生成自包含的交互式 3D 视图 HTML,样式与页面内预览一致(受体卡通彩虹色、配体绿色棍状、结合位点残基橙色棍状并标注名称、相互作用虚线、共晶配体灰色半透明参照),返回 (html 文本, 是否离线可用)。
    pocket 传入 residues_near_ligand 的结果时高亮结合位点残基;interactions 传入 interactions_near_ligand 的结果时画相互作用虚线并标注距离;
    het_pdb 传入共晶配体的 HETATM 文本时以半透明灰棍叠加作天然口袋参照;focus_label 为「聚焦口袋」按钮文案;
    lib_js 传入 3Dmol.js 源码时内嵌进文件双击即可离线打开,传 None 改用在线引用(打开时需联网)。
    """
    import html as _html
    import json

    if lib_js and "</script" not in lib_js:
        lib_tag = f"<script>\n{lib_js}\n</script>"
        offline_ok = True
    else:
        lib_tag = f'<script src="{_3DMOL_JS_URLS[0]}"></script>'
        offline_ok = False
    caption = _html.escape(caption)
    # 3Dmol.js 许可标注:内嵌时"内嵌",CDN 引用时"加载";完整许可文本以注释随文件提供
    note = ("内嵌 3Dmol.js (BSD-3-Clause)" if offline_ok
            else "本页面通过 CDN 加载 3Dmol.js (BSD-3-Clause)")
    out = (_3D_HTML_TEMPLATE
           .replace("__3DMOL_LICENSE__", _3DMOL_LICENSE_TEXT)
           .replace("__LICENSE_NOTE__",
                    f"{note} · Copyright (c) 2014, University of Pittsburgh · https://3dmol.org")
           .replace("__LIB__", lib_tag)
           .replace("__TITLE__", caption or "3D view")
           .replace("__CAPTION__", caption)
           .replace("__FOCUS_BTN__", _html.escape(focus_label))
           .replace("__INTERACTIONS__", json.dumps(interactions or []))
           .replace("__POCKET__", json.dumps(pocket or []))
           .replace("__HET__", json.dumps(het_pdb))
           .replace("__RECEPTOR__", json.dumps(receptor_pdb))
           .replace("__LIGAND__", json.dumps(ligand_pdb)))
    return out, offline_ok


def cofactor_hetatm_pdb(pdb_path: str, receptor_text: str | None = None) -> str:
    """提取受体 PDB 中的非水 HETATM(共晶配体/辅因子)行,拼成 PDB 文本,供 3D 渲染叠加天然口袋参照;
    无共晶配体时返回空字符串。receptor_text 传入已读取的 PDB 全文时不再重复读盘。"""
    if receptor_text is None:
        with open(pdb_path, encoding="utf-8", errors="replace") as f:
            receptor_text = f.read()
    return "\n".join(ln for ln in receptor_text.splitlines()
                     if ln.startswith("HETATM")
                     and ln[17:20].strip() not in _WATER_RESNAMES)


def interactive_3d_html(pdb_path: str, pose_pdbqt: str, caption: str = "",
                        pocket: list[dict] | None = None,
                        receptor_text: str | None = None,
                        pose_pdb: str | None = None,
                        interactions: list[dict] | None = None,
                        focus_label: str = "Focus pocket") -> tuple[str, bool]:
    """由受体 PDB 与所选构象 pdbqt 生成可保存的交互式 3D 视图 HTML。
    pocket 传入口袋残基列表时一并嵌入;interactions 传入相互作用列表时一并嵌入;
    受体含共晶配体(非水 HETATM)时以半透明灰棍叠加,过滤规则与页面预览一致。
    receptor_text 传入已读取的受体 PDB 全文、pose_pdb 传入已转换好的构象 PDB 路径时直接复用,
    供页面在同一次刷新里让预览与保存 HTML 共用一次转换、一次读盘;不传时保持原有的自足行为。
    """
    if pose_pdb is None:
        pose_pdb = os.path.join(os.path.dirname(pose_pdbqt), "pose_view.pdb")
        pose_pdbqt_to_pdb(pose_pdbqt, pose_pdb)
    if receptor_text is None:
        with open(pdb_path, encoding="utf-8", errors="replace") as f:
            receptor_text = f.read()
    with open(pose_pdb, encoding="utf-8") as f:
        ligand_pdb = f.read()
    het_pdb = cofactor_hetatm_pdb(pdb_path, receptor_text=receptor_text)
    return build_interactive_3d_html(receptor_text, ligand_pdb, caption,
                                     lib_js=_load_3dmol_js(),
                                     pocket=pocket, het_pdb=het_pdb or None,
                                     interactions=interactions, focus_label=focus_label)


def save_interactive_3d_html(pdb_path: str, pose_pdbqt: str, out_path: str,
                             caption: str = "",
                             pocket: list[dict] | None = None,
                             interactions: list[dict] | None = None,
                             focus_label: str = "Focus pocket") -> tuple[str, bool]:
    """把交互式 3D 视图保存为单个 HTML 文件(同名覆盖,便于同一构象反复保存);返回 (文件路径, 是否离线可用)。"""
    html_text, offline_ok = interactive_3d_html(pdb_path, pose_pdbqt, caption,
                                                pocket=pocket, interactions=interactions,
                                                focus_label=focus_label)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html_text)
    logging.info("Saved interactive 3D view: %s (offline=%s)", out_path, offline_ok)
    return out_path, offline_ok
