"""分子对接模块(AutoDock Vina)。

流程:把受体蛋白和配体小分子转成对接所需格式,确定结合口袋范围,
运行 Vina 对接并解析打分。中间文件会缓存一段时间以便复用,对接结果只归档不缓存。
"""

import hashlib
import logging
import os
from datetime import datetime
from typing import Iterator

import pandas as pd
import requests

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
                key = (line[21], line[22:26].strip(), line[17:20], line[12:16].strip())
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


def _pdb_coords(pdb_path: str) -> tuple[list, list]:
    """解析 PDB 文件中的原子坐标(单位 Å)。

    返回 (蛋白原子坐标列表, 非水其他分子坐标列表),坐标为 (x, y, z)。
    """
    protein, het = [], []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if line[17:20].strip() in _WATER_RESNAMES:
                continue  # 水分子不参与结合位点定位
            try:
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            (het if rec == "HETATM" else protein).append(coord)
    return protein, het


def parse_hetatm(pdb_path: str) -> list[tuple[str, str, str, tuple]]:
    """找出 PDB 结构中的共结晶配体/辅因子(水分子除外)及其位置。

    返回 (残基名, 链, 残基编号, (x, y, z)) 列表,供对接盒推断与 3D 预览使用。
    """
    het = []
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            resn = line[17:20].strip()
            if resn in _WATER_RESNAMES:
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            het.append((resn, line[21], line[22:26].strip(), xyz))
    return het


def residues_in_box(pdb_path: str, box: dict) -> list[str]:
    """列出落在对接盒内的蛋白残基,用于检查盒子是否覆盖了结合口袋。"""
    c, s = box["center"], box["size"]
    lo = (c[0] - s[0] / 2, c[1] - s[1] / 2, c[2] - s[2] / 2)
    hi = (c[0] + s[0] / 2, c[1] + s[1] / 2, c[2] + s[2] / 2)
    res = set()
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            if lo[0] <= xyz[0] <= hi[0] and lo[1] <= xyz[1] <= hi[1] and lo[2] <= xyz[2] <= hi[2]:
                res.add(f"{line[21]}{line[22:26].strip()}:{line[17:20].strip()}")
    return sorted(res)


def residues_near_ligand(receptor_pdb: str, ligand_pdbqt: str, cutoff: float = 4.0) -> list[dict]:
    """找出对接构象(配体)周围 cutoff Å 以内的蛋白残基,即结合位点残基。

    只统计受体 PDB 中的蛋白原子(ATOM 记录),水与共晶配体不参与;
    返回按链与残基编号排序、去重后的残基列表,每项含链号 chain、
    残基编号 resi、残基名 resn,以及标注名称用的代表原子坐标 pos(CA 原子优先)。
    """
    import re

    lig: list[tuple] = []
    with open(ligand_pdbqt, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            try:
                lig.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    if not lig:
        return []

    cutoff2 = cutoff * cutoff
    found: dict[tuple, dict] = {}
    with open(receptor_pdb, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM") or line[16:17] not in ("", " ", "A"):
                continue  # 非蛋白记录与 altloc 冗余构象不参与位点判定
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            key = (line[21:22].strip(), line[22:26].strip())
            info = found.get(key)
            if info is None:
                info = {"resn": line[17:20].strip(), "ca": None, "pos": xyz, "near": False}
                found[key] = info
            if info["ca"] is None and line[12:16].strip() == "CA":
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
    with open(path, "w", encoding="utf-8") as f:
        f.write(resp.text)
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
    """从 PubChem 下载配体的现成 3D 结构(SDF 文件)。

    cid 为 PubChem 化合物编号;无此化合物或无 3D 结构时返回空字符串。
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
    with open(path, "w", encoding="utf-8") as f:
        f.write(sdf_text)
    logging.info("Fetched 3D SDF for CID %d -> %s", cid, path)
    return path


def prepare_ligand_pdbqt(sdf_path: str, out_path: str | None = None) -> str:
    """把配体 3D SDF 转成对接所需的 pdbqt 格式,转换结果可复用。

    读取失败或缺少 3D 坐标时直接报错。返回生成的 pdbqt 文件路径。
    """
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
    mol = Chem.AddHs(mol)  # 补全氢原子

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
    )
    return [hit.get("identifier", "") for hit in resp.json().get("result_set", [])
            if hit.get("identifier")][:max_hits]


def _rcsb_graphql_entries(entry_ids: list[str]) -> dict[str, dict]:
    """按 entry 编号批量取分辨率/实验方法/标题元数据;返回 {entry号: 信息字典}。

    一次请求即可覆盖最多 100 个 entry;GraphQL 层面的软错误记录日志后
    按缺数据处理,不中断整体。
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


def fetch_receptor_candidates_batch(genes: list[str], max_hits: int = 20,
                                    allow_alphafold: bool = False, taxon: int = 9606,
                                    progress_cb=None) -> dict[str, pd.DataFrame]:
    """批量查询多个基因的 PDB 候选结构,返回 {基因符号: 候选表}。

    搜索阶段与旧的逐基因检索完全同一查询(保证召回一致),但按 4 路并发执行;
    元数据阶段合并为 GraphQL 批量请求,总请求数与基因数解耦。
    查无实验结构的基因在 allow_alphafold 开启时追加 AlphaFold 预测结构兜底。
    progress_cb(0~1) 可选,用于页面进度条。
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
        return g, _rcsb_search_entry_ids(g, max_hits)

    with ThreadPoolExecutor(max_workers=_RCSB_BATCH_WORKERS,
                            thread_name_prefix="rcsb-search") as pool:
        for i, (g, ids) in enumerate(pool.map(_search, uniq)):
            gene_ids[g] = ids
            if progress_cb:
                progress_cb((i + 1) / len(uniq) * 0.7)

    # 2) 所有命中 entry 去重后,合并为 GraphQL 批量元数据请求
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
                 "title": meta.get(eid, {}).get("title", "")}
                for eid in ids]
        if rows:
            results[g] = pd.DataFrame(rows)

    # 4) 无实验结构且开启开关时,逐基因追加 AlphaFold 预测结构兜底
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
    """按基因符号在 PDB 数据库(RCSB)中查找该靶蛋白的候选结构。

    返回表格含结构编号、分辨率、实验方法、标题等信息;查无结果时返回空表。
    开启 AlphaFold 选项且查无实验结构时,追加一条 AlphaFold 预测结构作兜底
    (预测结构查询失败不影响空表正常返回)。
    """
    import config
    from .http_client import request_with_retry

    # 1) 在 RCSB 中用基因符号精确检索结构
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": f'"{gene_symbol.strip().upper()}"'},
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    # 查询有频率限制,失败会自动重试几次,全部失败才报错。
    # 注意:检索成功但没有命中结构不算失败。
    try:
        resp = request_with_retry(
            "https://search.rcsb.org/rcsbsearch/v2/query", json=query,
            max_retries=config.RCSB_SEARCH_MAX_RETRIES, timeout=config.API_TIMEOUT,
            label=f"RCSB search for {gene_symbol}", throttle_key="rcsb_search",
            headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"},
            raise_on_exhausted=True,
        )
    except Exception as exc:
        # 搜索失败要报错供页面提示;查无结构属正常情况
        logging.exception("RCSB search failed for gene %s", gene_symbol)
        raise RuntimeError(f"RCSB 搜索失败（{gene_symbol}）：{exc}")
    hits = resp.json().get("result_set", [])[:max_hits]

    # 2) 补充每个结构的分辨率、实验方法、标题等信息
    def _fetch_meta(entry_id: str) -> tuple[str, float | None, str, str]:
        resp = request_with_retry(
            f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}", method="get",
            max_retries=2, timeout=15, label="RCSB meta", throttle_key="rcsb_data",
            headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"}, retry_delay=1,
        )
        if resp is not None:
            try:
                d = resp.json()
                res = d.get("rcsb_entry_info", {}).get("resolution_combined")
                res = res[0] if isinstance(res, list) and res else res
                method = ", ".join(m.get("method", "") for m in (d.get("exptl") or [])) or ""
                title = (d.get("struct") or {}).get("title", "") or ""
                return entry_id, res, method, title
            except Exception:
                pass  # 元数据解析失败按缺信息处理,不影响主流程
        return entry_id, None, "", ""

    ids = [hit.get("identifier", "") for hit in hits if hit.get("identifier")]
    meta = {}
    if ids:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=config.RCSB_DATA_MAX_WORKERS) as pool:
            for eid, res, method, title in pool.map(_fetch_meta, ids):
                meta[eid] = (res, method, title)

    rows = []
    for eid in ids:
        res, method, title = meta.get(eid, (None, "", ""))
        rows.append({"pdb_id": eid, "resolution": res, "method": method,
                     "title": title, "source": "experiment"})
    df = pd.DataFrame(rows)

    # 3) 无实验结构且开启开关时,追加 AlphaFold 预测结构
    if allow_alphafold and config.ALPHAFOLD_ENABLED and df.empty:
        try:
            acc = fetch_uniprot_accession(gene_symbol, taxon=taxon)
            if acc:
                df = pd.DataFrame([{
                    "pdb_id": f"AFDB:{acc}",
                    "resolution": None,
                    "method": f"AlphaFold predicted ({config.ALPHAFOLD_MODEL_VERSION})",
                    "title": f"AlphaFold {config.ALPHAFOLD_MODEL_VERSION} monomer model "
                             f"for {gene_symbol.strip().upper()}",
                    "source": "alphafold",
                }])
        except Exception:
            # 预测结构查询失败不影响空表正常返回
            logging.exception("AlphaFold fallback lookup failed for gene %s", gene_symbol)
    if df.empty:
        # 空表也保留统一的列名,方便后续处理
        df = pd.DataFrame(columns=["pdb_id", "resolution", "method", "title", "source"])
    return df


def prepare_receptor_pdbqt(pdb_path: str, out_path: str | None = None) -> str:
    """把受体蛋白 PDB 转成对接所需的 pdbqt 格式。

    处理包括:去水去杂原子、加氢(pH 7.4)、分配 AutoDock 原子类型;
    同一文件转换结果可复用。失败时报错,成功返回 pdbqt 文件路径。
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
    """根据 PDB 中共晶配体/辅因子的位置推断对接盒(结合口袋范围)。

    有配体/辅因子时以其位置圈定口袋,没有时用整个蛋白的范围;
    外扩距离默认 8 Å。结果仅供参考,使用前需人工确认口袋位置。
    返回口袋中心和尺寸(单位 Å)。
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
    """调用 AutoDock Vina,在指定口袋范围内对受体和配体执行分子对接。

    超时或运行失败时报错;运行日志(含打分)保留在返回结果中。
    """
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
    """把 vina 输出按结合模式拆成一个个独立的构象(pose)文件。

    拆分工具不可用时自动改用内置方式,不影响流程;两种方式都失败才报错。
    返回各 pose 文件路径。
    """
    import subprocess
    import glob
    import config

    out_dir = out_dir or os.path.dirname(out_pdbqt)
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(out_pdbqt))[0]

    # 优先使用 vina_split 工具拆分
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

    # 备用:内置方式拆分
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
    """从 vina 输出中读取各结合模式的打分(结合能)与 RMSD。

    按结合能从低到高排序(数值越负结合越强),第 1 名即最佳构象;
    无可解析结果时返回空表。
    """
    modes = []
    for block in _iter_vina_models(out_pdbqt):
        info: dict = {}
        for line in block:
            if line.startswith("REMARK VINA RESULT"):
                parts = line.split()
                if len(parts) >= 5:
                    info["affinity"] = float(parts[3])
                    info["rmsd_lb"] = float(parts[4])
                    info["rmsd_ub"] = float(parts[5]) if len(parts) > 5 else None
        if info:  # 与原逻辑一致:没有打分行的构象不计入
            modes.append(info)
    df = pd.DataFrame(modes)
    if df.empty:
        return pd.DataFrame(columns=["mode", "affinity", "rmsd_lb", "rmsd_ub"])
    df.insert(0, "mode", range(1, len(df) + 1))
    return df.sort_values("affinity").reset_index(drop=True)


# -------------------- 元数据 / 归档 / 清理 --------------------


def write_docking_metadata(meta_path: str, fields: dict) -> str:
    """把本次对接的关键信息(药物、靶点、参数等)写成 # 开头的文本记录,便于追溯。

    同名时自动换新文件名,不覆盖旧记录;返回实际写入的文件路径。
    """
    lines = [f"# {key}: {value}" for key, value in fields.items()]
    lines.append(f"# timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    path = _unique_path(meta_path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    logging.info("Wrote docking metadata: %s", path)
    return path


def star_docking_result(result_dir: str) -> str:
    """把本次对接成果目录完整复制到收藏区长期保存,不会被自动清理。

    原目录仍按版本控制自动管理,收藏副本独立存在;返回收藏副本路径。
    """
    import shutil
    import config

    os.makedirs(config.MD_STAR_DIR, exist_ok=True)
    dest = _unique_path(os.path.join(config.MD_STAR_DIR, os.path.basename(result_dir)))
    shutil.copytree(result_dir, dest)
    logging.info("Starred docking result: %s -> %s", result_dir, dest)
    return dest


def prune_docking_results(result_dir: str, max_versions: int | None = None) -> list[str]:
    """成果版本控制:同一药物-疾病-靶点组合的结果目录超过保留上限时,删除最旧的。

    返回本次删除的目录列表;无需删除时返回空列表。
    """
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
    """清理对接过程中的缓存中间文件:原始文件保留 1 天,加工产物保留 7 天。

    传入 ttl_days 时所有目录统一按该天数清理(如传 0 可全清);返回删除的文件数。
    """
    import time
    import config

    if ttl_days is not None:
        cache_ttl = raw_ttl = ttl_days
    else:
        cache_ttl = config.DOCKING_CACHE_TTL_DAYS
        raw_ttl = config.DOCKING_RAW_TTL_DAYS
    removed = 0
    if not os.path.isdir(config.DOCKING_CACHE_DIR):
        return 0
    raw_dir = config.DOCKING_RAW_CACHE_DIR
    for root, _, files in os.walk(config.DOCKING_CACHE_DIR):
        for name in files:
            path = os.path.join(root, name)
            ttl = raw_ttl if root == raw_dir else cache_ttl
            cutoff = time.time() - ttl * 86400
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

    接口按名称模糊匹配,易混入同源基因,拿到结果后需再按基因名精确筛选;
    只取主条目,跳过异构体(isoform)。
    """
    import config
    from .http_client import throttle

    throttle("ebi_proteins")
    resp = requests.get(
        config.EBI_PROTEINS_API_URL,
        params={"gene": gene_symbol, "taxid": taxon, "reviewed": "true", "format": "json"},
        timeout=config.API_TIMEOUT, headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"},
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
        timeout=config.API_TIMEOUT, headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"},
    )
    resp.raise_for_status()
    for line in resp.text.strip().splitlines()[1:]:  # 首行为表头
        acc = line.split("\t")[0].strip()
        if acc and "-" not in acc:  # 跳过 isoform
            return acc
    return ""


def fetch_uniprot_accession(gene_symbol: str, taxon: int = 9606) -> str:
    """
    基因符号 → UniProt 蛋白条目号,依次尝试多个数据源。

    查询结果会永久缓存,包括"确认查无"的情况;但网络失败不写入缓存,下次会重查。
    所有数据源都失败时报错;确认查无条目返回空字符串。
    """
    import json

    import config

    gene = gene_symbol.strip().upper()
    cache_file = config.UNIPROT_LOOKUP_CACHE
    os.makedirs(os.path.dirname(cache_file) or ".", exist_ok=True)
    db: dict = {}
    try:
        with open(cache_file, "r", encoding="utf-8") as f:
            db = json.load(f)
    except Exception:
        db = {}
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
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(db, f, ensure_ascii=False, sort_keys=True)
    except Exception:
        logging.exception("Failed to write UniProt lookup cache: %s", cache_file)
    if acc:
        logging.info("UniProt lookup: %s (taxon %s) -> %s", gene, taxon, acc)
    else:
        logging.warning("UniProt lookup: %s (taxon %s) has no reviewed entry", gene, taxon)
    return acc


def fetch_alphafold_pdb(accession: str) -> str:
    """
    下载 AlphaFold 预测结构 PDB 文件,下载过的直接复用缓存。

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
                        headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"})
    resp.raise_for_status()
    if len(resp.text) < 1024 or "\nATOM" not in resp.text:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        raise RuntimeError(f"AlphaFold PDB unexpected content for {acc}: {url}")
    with open(path, "w", encoding="utf-8") as f:
        f.write(resp.text)
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
    """
    对接盒边长超过上限时,以中心为准把尺寸压缩到上限。

    这只保证对接能运行,不保证口袋位置正确,需人工在 3D 预览中确认。
    """
    size = tuple(round(min(float(v), float(max_size)), 4) for v in box["size"])
    return {"center": tuple(round(float(v), 4) for v in box["center"]), "size": size}
