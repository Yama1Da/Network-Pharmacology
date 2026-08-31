"""分子对接模块（AutoDock Vina v1.2.7）。

配体 pdbqt 用 meeko（确定性）；受体 pdbqt 用 openbabel -xr（去水/加氢 pH7.4/AD typing）；
配体 3D SDF 统一 pubchempy.get_sdf(record_type="3d")。
缓存策略：受体/配体 pdbqt 与对接盒是确定性中间物 → 安全缓存（TTL 见
config.DOCKING_CACHE_TTL_DAYS）；vina 输出（best pose）是随机搜索产物 → 只归档不缓存。
"""

import hashlib
import logging
import os
from datetime import datetime

import pandas as pd
import requests


def _docking_fingerprint(path: str) -> str:
    """文件内容 MD5 指纹（前 12 位），用于 pdbqt 缓存命中判断。

    结构变化（不同 PDB 版本 / 不同 SDF 内容）指纹即变化，保证缓存不串货。
    """
    try:
        with open(path, "rb") as f:
            return hashlib.md5(f.read()).hexdigest()[:12]
    except OSError:
        return ""


def _unique_path(path: str) -> str:
    """防覆盖：目标路径已存在时追加时间戳（必要时再加序号），绝不覆写固定名。

    对接成果/元数据文件名都应经此函数生成，避免 open(..., "w") 直接覆盖旧结果。
    """
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
    """文本层清洗受体 PDB：仅保留蛋白 ATOM/TER/END。

    移除 HETATM（水、配体、离子——apo 受体对接不需要）并去重 altloc
    （只保留主构象，即 altLoc 为空或 'A'）。清洗失败抛 RuntimeError。
    """
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
    """解析 PDB 原子坐标（单位 Å）。

    Returns
    -------
    tuple[list, list]
        (蛋白 ATOM 坐标列表, 非水 HETATM 坐标列表)；每项为 (x, y, z) float 元组。
    """
    protein, het = [], []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            if line[17:20].strip() in ("HOH", "H2O", "WAT", "DOD"):
                continue  # 水分子不参与结合位点定位
            try:
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            (het if rec == "HETATM" else protein).append(coord)
    return protein, het


# -------------------- 配体链 --------------------


def fetch_ligand_sdf_3d(cid: int, out_dir: str | None = None) -> str:
    """通过 pubchempy 获取含 3D 坐标的配体 SDF 并落盘。

    Parameters
    ----------
    cid : int
        PubChem Compound ID（配体来源方案 A 的现成 3D 结构）。
    out_dir : str | None
        输出目录，默认 config.DOCKING_RAW_CACHE_DIR（原始文件，随 TTL 清理）。

    Returns
    -------
    str
        SDF 文件路径；PubChem 无此 CID / 无 3D 结构时返回空串。
    """
    import config
    import pubchempy

    out_dir = out_dir or config.DOCKING_RAW_CACHE_DIR
    os.makedirs(out_dir, exist_ok=True)
    # 确定性中间物：同 CID 同名幂等覆盖（不使用 _unique_path，避免缓存堆积）。
    # 已存在且非空时直接复用，不重复下载；过期由 clean_docking_cache 按 TTL 清理。
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
    """用 meeko 将配体 3D SDF 转为对接用 pdbqt（确定性中间物，可缓存）。

    Parameters
    ----------
    sdf_path : str
        含 3D 坐标的 SDF（meeko 要求分子带 conformer）。
    out_path : str | None
        输出路径；默认写入 config.DOCKING_LIGAND_CACHE_DIR，文件名带内容指纹。

    Returns
    -------
    str
        生成的 pdbqt 文件路径。

    Raises
    ------
    RuntimeError
        rdkit 读 SDF 失败 / meeko 原子类型分配或写出失败时抛出。
    """
    import config
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem

    mol = None
    try:
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False)  # 默认 sanitize，保证键序/形式电荷正确
    except Exception:
        logging.warning("ligand SDF sanitize failed, fallback to raw parse: %s", sdf_path)
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"rdkit failed to read ligand SDF: {sdf_path}")
    if mol.GetNumConformers() == 0:
        raise RuntimeError(f"ligand SDF has no 3D conformer: {sdf_path}（需 record_type=3d）")
    mol = Chem.AddHs(mol)  # 补显式氢，保证原子类型分配完整

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


def fetch_receptor_pdb_candidates(gene_symbol: str, max_hits: int = 50) -> pd.DataFrame:
    """按基因符号查询 RCSB PDB 候选蛋白结构。

    实现说明：RCSB Search API 已不再对 gene/organism 类 attribute 开放 text 检索
    （"search is not enabled"），故改用 full-text 精确短语匹配基因符号，再经
    RCSB Data API 并发补充每个条目的分辨率/方法/标题。

    Returns
    -------
    pd.DataFrame
        列：pdb_id / resolution / method / title；查询失败或无结果返回空表。
    """
    import config

    # 1) RCSB Search API：full-text 精确短语匹配基因符号
    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": f'"{gene_symbol.strip().upper()}"'},
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    try:
        resp = requests.post(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            json=query, timeout=config.API_TIMEOUT,
            headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"},
        )
        resp.raise_for_status()
        hits = resp.json().get("result_set", [])[:max_hits]
    except Exception as exc:
        # 搜索失败必须向上抛（带原因），供页面提示；200 空 result_set 才是"无结构"
        logging.exception("RCSB search failed for gene %s", gene_symbol)
        raise RuntimeError(f"RCSB 搜索失败（{gene_symbol}）：{exc}") from exc

    # 2) RCSB Data API：并发补充元数据（分辨率/方法/标题）
    def _fetch_meta(entry_id: str) -> tuple[str, float | None, str, str]:
        try:
            d = requests.get(
                f"https://data.rcsb.org/rest/v1/core/entry/{entry_id}",
                timeout=15, headers={"User-Agent": "NetworkPharmacologyPlatform/1.0"},
            ).json()
            res = d.get("rcsb_entry_info", {}).get("resolution_combined")
            res = res[0] if isinstance(res, list) and res else res
            method = ", ".join(m.get("method", "") for m in (d.get("exptl") or [])) or ""
            title = (d.get("struct") or {}).get("title", "") or ""
            return entry_id, res, method, title
        except Exception:
            return entry_id, None, "", ""

    ids = [hit.get("identifier", "") for hit in hits if hit.get("identifier")]
    meta = {}
    if ids:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=8) as pool:
            for eid, res, method, title in pool.map(_fetch_meta, ids):
                meta[eid] = (res, method, title)

    rows = []
    for eid in ids:
        res, method, title = meta.get(eid, (None, "", ""))
        rows.append({"pdb_id": eid, "resolution": res, "method": method, "title": title})
    return pd.DataFrame(rows)


def prepare_receptor_pdbqt(pdb_path: str, out_path: str | None = None) -> str:
    """将受体 PDB 转为对接用 pdbqt（确定性中间物，可缓存）。

    处理链：文本层清洗（去水/去 HETATM 配体/altloc 去重）→
    openbabel 受体模式（加氢 pH7.4、AutoDock 原子类型分配）。

    Parameters
    ----------
    pdb_path : str
        受体原始 PDB 文件路径（如 RCSB 下载的 {pdb_id}.pdb）。
    out_path : str | None
        输出路径；默认写入 config.DOCKING_RECEPTOR_CACHE_DIR，文件名带内容指纹。

    Returns
    -------
    str
        生成的 pdbqt 文件路径。

    Raises
    ------
    RuntimeError
        openbabel 读 PDB 或写 pdbqt 失败时抛出。
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
    conv.SetOptions("r", ob.OBConversion.OUTOPTIONS)  # 受体模式：加氢 / 去水 / AD typing
    mol = ob.OBMol()
    if not conv.ReadFile(mol, source):
        raise RuntimeError(f"openbabel failed to read receptor PDB: {source}")
    if not conv.WriteFile(mol, out_path):
        raise RuntimeError(f"openbabel failed to write receptor pdbqt: {out_path}")
    logging.info("Prepared receptor pdbqt: %s", out_path)
    return out_path


# -------------------- 对接盒 --------------------


def infer_box_from_cofactor(pdb_path: str, padding: float = 8.0) -> dict:
    """依据 PDB 中的配体/辅因子（HETATM，非水）坐标推断对接盒。

    PDB 含 HETATM 时以其质心/包围盒定位结合位点；无 HETATM（apo 结构）
    时退化为蛋白整体包围盒。返回值为科研决策点参考，仍需用户人工确认。

    Parameters
    ----------
    pdb_path : str
        PDB 文件路径。
    padding : float
        结合位点包围盒外扩距离（Å），默认 8.0（覆盖口袋残基）。

    Returns
    -------
    dict
        {"center": (x, y, z), "size": (sx, sy, sz)}，单位 Å，保留 4 位小数。
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
    """subprocess 调用 AutoDock Vina 执行对接（随机搜索）。

    Parameters
    ----------
    receptor_pdbqt / ligand_pdbqt : str
        预处理后的 pdbqt 路径。
    box : dict
        {"center": (x,y,z), "size": (sx,sy,sz)}，由 infer_box_from_cofactor 得到。
    out_path : str | None
        vina 输出 pdbqt 路径；默认写入 config.MD_RESULT_DIR 且自动防覆盖。

    Returns
    -------
    subprocess.CompletedProcess
        包含 stdout/stderr（vina 打分 log 在 stdout，供前端展示）。

    Raises
    ------
    RuntimeError
        vina 超时（config.DOCKING_TIMEOUT）或非零退出时抛出。
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


def _split_models_python(out_pdbqt: str, out_dir: str) -> list[str]:
    """纯 Python 按 MODEL/ENDMDL 拆分 vina 多模型输出（vina_split 不可用时的降级）。

    文件名格式：{输入名}_{序号}.pdbqt（vina_split 的替代，不依赖外部 exe）。
    """
    stem = os.path.splitext(os.path.basename(out_pdbqt))[0]
    blocks = []
    current = None
    with open(out_pdbqt, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                current = [line]
            elif line.startswith("ENDMDL") and current is not None:
                current.append(line)
                blocks.append(current)
                current = None
            elif current is not None:
                current.append(line)
    files = []
    for i, block in enumerate(blocks, 1):
        path = os.path.join(out_dir, f"{stem}_{i}.pdbqt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("".join(block))
        files.append(path)
    return files


def split_vina(out_pdbqt: str, out_dir: str | None = None) -> list[str]:
    """将 vina 多模型输出拆分为独立 pose 文件（每个结合模式一个 pdbqt）。

    优先调用 vina_split exe；若 exe 缺失/损坏（如 0 字节）或执行失败，
    自动降级为纯 Python 按 MODEL 拆分，保证流程不被外部工具阻塞。

    Parameters
    ----------
    out_pdbqt : str
        run_vina 的输出文件。
    out_dir : str | None
        拆分产物输出目录；默认与 out_pdbqt 同目录。

    Returns
    -------
    list[str]
        拆分出的 pose 文件路径（按文件名排序）。

    Raises
    ------
    RuntimeError
        vina_split 与纯 Python 拆分均失败时抛出。
    """
    import subprocess
    import glob
    import config

    out_dir = out_dir or os.path.dirname(out_pdbqt)
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(out_pdbqt))[0]

    # 优先：vina_split exe（须是有效 PE 文件，非 0 字节）
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

    # 降级：纯 Python 拆分
    files = _split_models_python(out_pdbqt, out_dir)
    if not files:
        raise RuntimeError(f"failed to split vina output: {out_pdbqt}")
    return files


def parse_vina_scores(out_pdbqt: str) -> pd.DataFrame:
    """解析 vina 输出 pdbqt 中各结合模式的打分与 RMSD。

    Vina 在每个 MODEL 块写入 "REMARK VINA RESULT <affinity> <rmsd_lb> <rmsd_ub>"，
    本函数按此提取并按亲和力升序排序（mode 1 即最佳 pose）。

    Returns
    -------
    pd.DataFrame
        列：mode / affinity / rmsd_lb / rmsd_ub；无可解析模式时返回空表。
    """
    modes = []
    current = None
    with open(out_pdbqt, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                current = {}
            elif line.startswith("REMARK VINA RESULT") and current is not None:
                parts = line.split()
                if len(parts) >= 5:
                    current["affinity"] = float(parts[3])
                    current["rmsd_lb"] = float(parts[4])
                    current["rmsd_ub"] = float(parts[5]) if len(parts) > 5 else None
            elif line.startswith("ENDMDL") and current:
                modes.append(current)
                current = None
    df = pd.DataFrame(modes)
    if df.empty:
        return pd.DataFrame(columns=["mode", "affinity", "rmsd_lb", "rmsd_ub"])
    df.insert(0, "mode", range(1, len(df) + 1))
    return df.sort_values("affinity").reset_index(drop=True)


# -------------------- 元数据 / 归档 / 清理 --------------------


def write_docking_metadata(meta_path: str, fields: dict) -> str:
    """写对接元数据文本（# 注释行，对齐富集/PPI 缓存的 # 头风格）。

    Parameters
    ----------
    meta_path : str
        目标路径（自动防覆盖，同名时追加时间戳）。
    fields : dict
        键值对，按序写为 "# key: value"；时间戳自动追加到末尾。

    Returns
    -------
    str
        实际写入的文件路径。
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
    """将成果目录完整复制到收藏区（MD_STAR_DIR/star），长期保留（仅可手动删除）。

    收藏副本与原成果互不影响：原成果仍受 prune_docking_results 版本控制，
    收藏副本独立存在于 star/ 下，不受自动删除影响。

    Parameters
    ----------
    result_dir : str
        本次对接成果目录（MD_RESULT_DIR 下，时间戳命名）。

    Returns
    -------
    str
        收藏区中的目录路径（自动防覆盖）。
    """
    import shutil
    import config

    os.makedirs(config.MD_STAR_DIR, exist_ok=True)
    dest = _unique_path(os.path.join(config.MD_STAR_DIR, os.path.basename(result_dir)))
    shutil.copytree(result_dir, dest)
    logging.info("Starred docking result: %s -> %s", result_dir, dest)
    return dest


def prune_docking_results(result_dir: str, max_versions: int | None = None) -> list[str]:
    """成果版本控制：同条件（drug_disease_gene）目录超过上限时，直接删除最旧版本。

    目录名格式约定：{prefix}_{YYYYMMDD_HHMMSS}，prefix 如 aspirin_cancer_BRCA1。

    Parameters
    ----------
    result_dir : str
        本次对接成果目录（MD_RESULT_DIR 下）。
    max_versions : int | None
        最多保留版本数，默认 config.DOCKING_MAX_VERSIONS（3）。

    Returns
    -------
    list[str]
        本次实际删除的旧版目录路径列表；无需删除返回空列表。
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
    # 越界防护：prefix 异常时 glob 可能匹配到 MD_RESULT_DIR 之外，
    # 删除前校验候选目录确在 MD_RESULT_DIR 内，阻断路径穿越删除
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
    """清理 docking 缓存目录下超过 TTL 的中间物文件。

    TTL 分级：
      - raw 子目录（原始 SDF/PDB/vina 输出，可随时重新下载生成）→
        config.DOCKING_RAW_TTL_DAYS（24 小时）
      - ligand / receptor / box（加工后的 pdbqt 等确定性产物）→
        config.DOCKING_CACHE_TTL_DAYS（7 天）

    Parameters
    ----------
    ttl_days : int | None
        显式指定时对全部目录统一使用该 TTL（便于手动全清，如传 0）；
        默认 None 按上面的分级策略清理。

    Returns
    -------
    int
        删除的文件数。
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
