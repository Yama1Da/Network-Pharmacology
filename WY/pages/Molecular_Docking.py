import os
import re
import shutil
from datetime import datetime

import pandas as pd
import streamlit as st

import config
import Utils

st.set_page_config(page_title="Molecular Docking", layout="wide")

ENRICHMENT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Static", "data", "cache", "enrichment",
)


# ---------- Hide auto navigation sidebar for full-screen display ----------
st.markdown("""
<style>
    [data-testid="stSidebar"] {
        display: none;
    }
    /* Star Management 按钮（#6347A1 底、#E0D8F2 字，与 Data Network 的 Data Management 一致） */
    .st-key-star_mgmt_btn button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }
</style>
""", unsafe_allow_html=True)

st.title("⚗️ Molecular Docking")
st.caption("AutoDock Vina v1.2.7 · ligand pubchempy/rdkit · receptor openbabel · results archived without overwrite")

col_back, col_star = st.columns([1, 1])
with col_back:
    if st.button("← Back to Main", key="back_btn"):
        st.switch_page("App.py")
with col_star:
    if st.button("Star Management", key="star_mgmt_btn"):
        st.switch_page("pages/Star_Management.py")

st.divider()

# ---------- session state（防 rerun 重跑重型计算） ----------
# 从主入口（App 按钮）进入时整体重置：不保留上次选择的文件与结果。
# 页面内 rerun 不受影响（_docking_entry 只由入口按钮设置一次）。
if st.session_state.pop("_docking_entry", False):
    st.session_state.pop("docking", None)
    for _k in ("sel_ppi_dir", "sel_enr", "gene_sel", "pdb_sel", "lig_source",
               "lig_query", "lig_smiles", "drug_name", "disease_name",
               "use_struct", "pose_sel", "conclusion", "manual_genes"):
        st.session_state.pop(_k, None)

st.session_state.setdefault("docking", {})
dock = st.session_state["docking"]

# 页面启动时清理一次过期缓存（每个浏览器会话首次进入该页时触发，后续 rerun 不重复）
if not st.session_state.get("_docking_cache_cleaned"):
    Utils.clean_docking_cache()
    st.session_state["_docking_cache_cleaned"] = True


# ==================== 辅助函数 ====================

def _extract_label(file_name: str) -> str:
    """从富集缓存文件名解析来源标识（如 gProfiler_aspirin_cancer_xxx.csv -> aspirin_cancer）。"""
    m = re.match(r"^(?:gProfiler|Enrichr)_(.+)_[0-9a-f]{32}\.csv$", file_name)
    if m:
        return m.group(1)
    m2 = re.match(r"^(?:gProfiler|Enrichr)_([0-9a-f]{32})\.csv$", file_name)
    return m2.group(1) if m2 else ""


def _parse_genes(raw: str) -> list[str]:
    """解析富集结果 genes 字段为基因符号列表（大写去重前）。

    g:Profiler 的 intersection 落盘为 Python 列表字符串（如 "['EGFR', 'ERBB2']"），
    也兼容纯空格/逗号分隔格式（"EGFR ERBB2"）。统一剥离方括号与引号再分割，
    避免 'ERBB2' 与 ERBB2 被当作两个不同靶点。
    """
    raw = (raw or "").strip()
    if not raw:
        return []
    # Python 列表字符串 → ast 解析
    if raw.startswith("[") and raw.endswith("]"):
        try:
            import ast
            items = ast.literal_eval(raw)
            if isinstance(items, list):
                return [str(x).strip() for x in items if str(x).strip()]
        except Exception:
            pass
    # 通用分隔格式：剥离残留括号/引号后按逗号、分号、空白分割
    cleaned = raw.strip("[]").replace("'", "").replace('"', "")
    return [t.strip() for t in re.split(r"[,;\s]+", cleaned) if t.strip()]


def _genes_from_df(df: pd.DataFrame) -> set[str]:
    """提取富集结果表的 genes 列（各 term 交集基因），跨行取并集（大写去重）。"""
    genes: set[str] = set()
    for raw in df.get("genes", pd.Series(dtype=str)).fillna("").astype(str):
        genes.update(t.upper() for t in _parse_genes(raw))
    return genes


def _reset_docking_downstream() -> None:
    """数据源（候选基因）变化时，清空配体/受体/盒/对接结果等全部下游状态。

    避免切换 Step 1 数据源后残留上一份文件的 Step 2-6 结果。
    """
    for _k in ("lig_pdbqt", "cid", "lig_auto_done",
               "gene", "rec_pdbqt", "pdb_id", "pdb_path",
               "box", "result", "result_dir", "starred",
               "gene_structures", "candidate_pre"):
        dock.pop(_k, None)


def _read_ppi_dataset(dir_name: str) -> tuple[set[str], set[str]]:
    """从 downloads/{dir_name}/ 读取 PPI 结果：(HUB 基因集合, 富集验证基因集合)。

    目录结构（save_ppi_results 约定）：{prefix}_data/{prefix}_nodes.csv（含 is_hub）
    与 {prefix}_enrichment.csv（含 genes 列）；文件缺失时对应集合为空。
    """
    base = os.path.join(config.DOWNLOADS_DIR, dir_name)
    prefix = dir_name[:-5] if dir_name.endswith("_data") else dir_name
    hub: set[str] = set()
    nodes_path = os.path.join(base, f"{prefix}_nodes.csv")
    if os.path.exists(nodes_path):
        nodes = pd.read_csv(nodes_path, encoding="utf-8-sig", comment="#")
        hub = set(nodes.loc[nodes["is_hub"].fillna(False), "gene"]
                  .dropna().astype(str).str.upper())
    enr: set[str] = set()
    enr_path = os.path.join(base, f"{prefix}_enrichment.csv")
    if os.path.exists(enr_path):
        enr = _genes_from_df(pd.read_csv(enr_path, encoding="utf-8-sig", comment="#"))
    return hub, enr


def _resolve_cid(query: str) -> int | None:
    """药物名（或 CID 数字）→ PubChem CID。"""
    query = query.strip()
    if query.isdigit():
        return int(query)
    import pubchempy
    try:
        comps = pubchempy.get_compounds(query, "name")
        return comps[0].cid if comps else None
    except Exception:
        return None


def _smiles_to_sdf_3d(smiles: str, out_path: str, seed: int = 0, n_confs: int = 20) -> str:
    """SMILES → RDKit 生成 3D 构象（固定 seed 枚举 + MMFF 取最低能），写 SDF 落盘。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if not cids:
        raise RuntimeError("RDKit failed to generate 3D conformers")
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


def _download_pdb(pdb_id: str) -> str:
    """按 PDB id 下载原始 PDB 到 raw 缓存目录（已存在则复用）。"""
    import requests

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


def _extract_model(out_pdbqt: str, mode: int, out_path: str) -> str:
    """从 vina 多模型输出提取第 mode 个 MODEL 块写入 out_path（不依赖 split 命名）。"""
    blocks, current = [], None
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
    if not (1 <= mode <= len(blocks)):
        raise RuntimeError(f"mode {mode} out of range ({len(blocks)} total)")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("".join(blocks[mode - 1]))
    return out_path


def _save_docking_result(prefix: str, gene: str, pdb_id: str, cid: int | None, box: dict,
                         out_pdbqt: str, scores: pd.DataFrame, conclusion: str) -> str:
    """成果归档：MDresult/{prefix}_{ts}/ 下写 summary / best pose / meta，并做版本控制与缓存清理。"""
    gene = Utils.safe_name(gene)  # gene 可能来自手动输入，进入文件名前需净化
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    result_dir = os.path.join(config.MD_RESULT_DIR, f"{prefix}_{ts}")
    os.makedirs(result_dir, exist_ok=True)

    # 打分表 + 结论
    summary = scores.copy()
    summary["conclusion"] = conclusion
    summary_path = os.path.join(result_dir, f"summary_{ts}.csv")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    # 最佳 pose（直接按 MODEL 序号从原输出提取，最可靠）
    best_mode = int(summary.iloc[0]["mode"])
    best_path = os.path.join(result_dir, f"best_{gene}_{ts}.pdbqt")
    _extract_model(out_pdbqt, best_mode, best_path)

    # 元数据（# 头，对齐富集/PPI 缓存风格）
    Utils.write_docking_metadata(os.path.join(result_dir, f"meta_{ts}.txt"), {
        "drug": prefix.split("_")[0] if "_" in prefix else prefix,
        "disease": prefix.split("_", 1)[1] if "_" in prefix else "",
        "gene": gene,
        "pdb_id": pdb_id,
        "cid": cid or "",
        "box_center": str(box["center"]),
        "box_size": str(box["size"]),
        "exhaustiveness": config.EXHAUSTIVENESS,
        "num_modes": config.NUM_MODES,
        "seed": config.DOCKING_SEED,
        "best_mode": best_mode,
        "best_affinity": float(summary.iloc[0]["affinity"]),
        "conclusion": conclusion,
    })

    # 版本控制（同前缀限 N 版，旧版直接删除）+ 中间物 TTL 清理
    Utils.prune_docking_results(result_dir)
    Utils.clean_docking_cache()
    return result_dir


def _render_3d(pdb_path: str, pose_pdbqt: str) -> None:
    """py3Dmol 内嵌展示受体（卡通）+ 最佳 pose（棍状）；未安装时提示并降级。"""
    try:
        import py3Dmol
    except ImportError:
        st.info("3D inline visualization requires py3Dmol (`pip install py3Dmol`). "
                "Meanwhile, download the pose file and open it with PyMOL / ChimeraX / VMD.")
        return
    from openbabel import openbabel as ob

    pose_pdb = os.path.join(os.path.dirname(pose_pdbqt), "pose_view.pdb")
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("pdbqt", "pdb")
    mol = ob.OBMol()
    if not (conv.ReadFile(mol, pose_pdbqt) and conv.WriteFile(mol, pose_pdb)):
        st.warning("Pose-to-PDB conversion failed; skipping 3D display")
        return
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        receptor_pdb = f.read()
    with open(pose_pdb, encoding="utf-8") as f:
        ligand_pdb = f.read()

    view = py3Dmol.view(width=900, height=520)
    view.addModel(receptor_pdb, "pdb")
    view.addModel(ligand_pdb, "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
    view.setStyle({"model": 1}, {"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()
    try:
        html = view._make_html()
    except AttributeError:
        html = view.write()
    st.components.v1.html(html, height=540)


@st.cache_data(ttl=3600)
def _cached_pdb_candidates(gene: str) -> pd.DataFrame:
    """RCSB 候选结构查询缓存 1 小时，避免每次 rerun 重复请求。"""
    return Utils.fetch_receptor_pdb_candidates(gene, max_hits=20)


# ==================== Step 1 候选基因筛选（HUB ∩ 富集验证 ∩ 结构可用性） ====================
st.subheader("1️⃣ Candidate gene filtering (HUB ∩ enrichment-validated ∩ structure-available)")

src_mode = st.radio("Data source",
                    ["PPI result directory (HUB + enrichment)", "Enrichment cache file (enrichment only)"],
                    horizontal=True, key="gene_src_mode")

# 切换数据源模式 → 重置数据源状态与全部下游，强制新模式重新读取（避免残留旧文件候选基因）
if dock.get("gene_src_mode") != src_mode:
    dock["gene_src_mode"] = src_mode
    for _k in ("ppi_dir", "enr_file", "hub_genes", "enrichment_genes",
               "label", "drug", "disease", "genes"):
        dock.pop(_k, None)
    _reset_docking_downstream()

if src_mode.startswith("PPI"):
    # 仅列出 PPI 结果目录（*_data 后缀），排除 MDresult 等非分析结果目录
    dirs = Utils.list_download_dirs()
    ppi_dirs = dirs[dirs["Directory"].str.endswith("_data", na=False)] if not dirs.empty else dirs
    if ppi_dirs.empty:
        st.warning("No PPI result directories found (downloads/*_data/). Run the Data Network analysis first, or switch to the enrichment cache.")
    else:
        sel_dir = st.selectbox("Select PPI result directory", ppi_dirs["Directory"].tolist(), key="sel_ppi_dir",
                               index=None, placeholder="Select a PPI result directory…")
        # 自动读取：目录变化即加载（无需点按钮）；选择为空时跳过，避免误读旧值
        if sel_dir and dock.get("ppi_dir") != sel_dir:
            dock["ppi_dir"] = sel_dir
            dock["hub_genes"], dock["enrichment_genes"] = _read_ppi_dataset(sel_dir)
            # 从 {药物名}_{疾病名}_data 解析药物名/疾病名，供 Step 2 自动获取配体
            prefix = sel_dir[:-5] if sel_dir.endswith("_data") else sel_dir
            dock["label"] = prefix
            parts = prefix.split("_", 1)
            dock["drug"] = parts[0] if parts[0] else ""
            dock["disease"] = parts[1] if len(parts) > 1 else ""
            _reset_docking_downstream()  # 数据源变化 → 清空 Step 2-6 下游结果
else:
    enr_files = Utils.list_cache_files("enrichment")
    if enr_files.empty:
        st.warning("No enrichment cache files found. Run the enrichment analysis first, or enter genes manually.")
    else:
        sel_file = st.selectbox("Select enrichment result file", enr_files["File name"], key="sel_enr",
                                index=None, placeholder="Select an enrichment result file…")
        # 自动读取：文件变化即加载；选择为空时跳过
        if sel_file and dock.get("enr_file") != sel_file:
            dock["enr_file"] = sel_file
            df = pd.read_csv(os.path.join(ENRICHMENT_DIR, sel_file),
                             encoding="utf-8-sig", comment="#")
            dock["enrichment_genes"] = sorted(_genes_from_df(df))
            dock["hub_genes"] = []  # 无 HUB 信息时仅按富集验证筛选
            label = _extract_label(sel_file)
            dock["label"] = label
            parts = label.split("_", 1)
            dock["drug"] = parts[0] if parts[0] else ""
            dock["disease"] = parts[1] if len(parts) > 1 else ""
            _reset_docking_downstream()  # 数据源变化 → 清空 Step 2-6 下游结果

if dock.get("enrichment_genes") is not None:
    hub_set = set(dock.get("hub_genes") or [])
    enr_set = set(dock["enrichment_genes"])
    candidate = sorted(hub_set & enr_set) if hub_set else sorted(enr_set)
    dock["candidate_pre"] = candidate

    col_h, col_e, col_c = st.columns(3)
    col_h.metric("HUB genes", len(hub_set))
    col_e.metric("Enrichment-validated genes", len(enr_set))
    col_c.metric("Intersection candidates", len(candidate))
    if candidate:
        st.write(", ".join(candidate[:60]) + (" …" if len(candidate) > 60 else ""))
    else:
        st.warning("No intersection between HUB genes and enrichment-validated genes; the candidate list is empty. Switch the PPI result directory / enrichment file, or enter genes manually.")

    # 结构可用性筛选（默认开启；RCSB 查询失败/漏检时用户可取消，保留全部候选）
    st.markdown("**Structure availability filter (RCSB human structures)**")
    use_struct = st.checkbox("Enable structure availability filter", key="use_struct",
                             disabled=not candidate, value=True)
    if use_struct:
        if st.button("Check RCSB structures gene by gene (slow)", key="check_struct"):
            with st.spinner("Querying RCSB gene by gene…"):
                try:
                    found: dict[str, int] = {}
                    progress = st.progress(0)
                    for i, g in enumerate(candidate):
                        cand = _cached_pdb_candidates(g)
                        if not cand.empty:
                            found[g] = len(cand)
                        progress.progress((i + 1) / len(candidate))
                    dock["gene_structures"] = found
                except Exception as exc:
                    # 中途失败：保留已查结果，明确提示错误原因
                    dock.pop("gene_structures", None)
                    st.error(f"Error: {exc}")
        if dock.get("gene_structures") is not None:
            struct_ok = sorted(dock["gene_structures"])
            dock["genes"] = struct_ok
            no_struct = [g for g in candidate if g not in dock["gene_structures"]]
            st.success(f"Genes with available structures ({len(struct_ok)}): "
                       + ", ".join(struct_ok[:40]) + (" …" if len(struct_ok) > 40 else ""))
            if no_struct:
                st.caption(f"Excluded (no structure found) ({len(no_struct)}): "
                           + ", ".join(no_struct[:40]) + (" …" if len(no_struct) > 40 else ""))
        else:
            dock.pop("genes", None)
            st.info("Click the button above to run the structure availability check.")
    else:
        # 未启用结构筛选 → 直接确定候选（空也赋值，避免 Step 3 误报"未完成 Step 1"）
        dock["genes"] = candidate
        if candidate:
            st.caption(f"Final candidate genes ({len(candidate)}): "
                       + ", ".join(candidate[:60]) + (" …" if len(candidate) > 60 else ""))

# 未加载任何 PPI/富集数据时，提供手动输入兜底
if dock.get("enrichment_genes") is None:
    manual = st.text_area("Or enter genes manually (comma / space separated)", key="manual_genes")
    if st.button("Use manual genes", key="load_genes_manual"):
        dock["genes"] = sorted({g.strip().upper() for g in re.split(r"[,;\s]+", manual) if g.strip()})
        dock["label"] = ""
        dock["drug"] = ""
        dock["disease"] = ""
        dock["hub_genes"], dock["enrichment_genes"] = [], []
        _reset_docking_downstream()  # 数据源变化 → 清空 Step 2-6 下游结果

if dock.get("genes"):
    shown = ", ".join(dock["genes"][:60]) + (" …" if len(dock["genes"]) > 60 else "")
    st.info(f"✅ Candidate genes confirmed ({len(dock['genes'])}): {shown}")

st.divider()

# ==================== Step 2 配体 ====================
st.subheader("2️⃣ Ligand preparation")

# 药物名/疾病名取自 Step 1 所选 {药物名}_{疾病名}_data 目录（或富集 label）
default_drug = dock.get("drug", "")
default_disease = dock.get("disease", "")

lig_source = st.radio("Ligand source", ["A. Ready-made 3D from PubChem (default)", "B. SMILES + RDKit generated 3D"],
                      horizontal=True, key="lig_source")
if lig_source.startswith("A"):
    # 自动获取：已解析出药物名且尚未准备配体时自动执行一次（先执行后渲染，便于按结果决定是否显示手动输入）
    if default_drug and not dock.get("lig_pdbqt") and not dock.get("lig_auto_done"):
        with st.spinner(f"Auto-fetching 3D structure of {default_drug} and preparing ligand…"):
            cid = _resolve_cid(default_drug)
            if cid is None:
                dock["lig_auto_done"] = True
                st.error(f"PubChem returned no results for {default_drug}. Enter a drug name or CID below to fetch manually.")
            else:
                sdf_path = Utils.fetch_ligand_sdf_3d(cid)
                if not sdf_path:
                    dock["lig_auto_done"] = True
                    st.error(f"CID {cid} has no 3D structure. Try another drug name or CID below.")
                else:
                    dock["cid"] = cid
                    dock["lig_pdbqt"] = Utils.prepare_ligand_pdbqt(sdf_path)

    # 手动输入框仅在自动获取失败（或未识别到药物名）时显示，其余情况隐藏以免用户误输入
    lig_ready = bool(dock.get("lig_pdbqt"))
    show_lig_input = not lig_ready and (not default_drug or dock.get("lig_auto_done"))
    if show_lig_input:
        if not default_drug:
            st.caption("No drug name recognized from the data source. Enter a drug name or PubChem CID manually.")
        lig_query = st.text_input("Drug name or PubChem CID (manual fetch)", value=default_drug, key="lig_query")
        if st.button("Fetch & prepare ligand", key="prep_lig"):
            dock["lig_auto_done"] = True  # 手动指定后不再自动覆盖
            with st.spinner("Fetching 3D structure from PubChem and converting to pdbqt…"):
                cid = _resolve_cid(lig_query)
                if cid is None:
                    st.error(f"PubChem returned no compound for: {lig_query}. Check the name or enter a CID directly.")
                else:
                    sdf_path = Utils.fetch_ligand_sdf_3d(cid)
                    if not sdf_path:
                        st.error(f"CID {cid} has no 3D structure. Try another drug.")
                    else:
                        dock["cid"] = cid
                        dock["lig_pdbqt"] = Utils.prepare_ligand_pdbqt(sdf_path)
else:
    smiles = st.text_area("SMILES", height=90, key="lig_smiles",
                          placeholder="e.g. aspirin: CC(=O)Oc1ccccc1C(=O)O")
    if st.button("Generate & prepare ligand", key="prep_lig_b"):
        with st.spinner("Generating 3D conformers with RDKit and converting to pdbqt…"):
            try:
                tmp_dir = config.DOCKING_RAW_CACHE_DIR
                os.makedirs(tmp_dir, exist_ok=True)
                sdf_path = _smiles_to_sdf_3d(
                    smiles,
                    os.path.join(tmp_dir, f"ligand_smiles_{datetime.now():%Y%m%d_%H%M%S}.sdf"))
                dock["cid"] = None
                dock["lig_pdbqt"] = Utils.prepare_ligand_pdbqt(sdf_path)
            except Exception as exc:
                st.error(f"Ligand preparation failed: {exc}")

if dock.get("lig_pdbqt"):
    st.success(f"✅ Ligand pdbqt: `{dock['lig_pdbqt']}`"
               + (f" (CID {dock.get('cid')})" if dock.get("cid") else " (generated from SMILES)"))

st.divider()

# ==================== Step 3 受体（基因 + PDB，决策点） ====================
st.subheader("3️⃣ Receptor selection (pick gene → auto-query candidate PDBs → pick structure)")

if dock.get("genes"):
    gene = st.selectbox("Target gene (receptor protein)", dock["genes"], key="gene_sel")
    # 基因变化 → 清空旧受体/盒/对接结果
    if dock.get("gene") != gene:
        dock["gene"] = gene
        dock.pop("rec_pdbqt", None)
        dock.pop("pdb_path", None)
        dock.pop("pdb_id", None)
        dock.pop("box", None)
        dock.pop("result", None)

    # 自动查询该基因的 RCSB 候选结构（@st.cache_data 命中时零请求）
    try:
        cand = _cached_pdb_candidates(gene)
    except Exception as exc:
        st.warning(f"Error: {exc} (retry later or pick another gene)")
        cand = None
    if cand is None:
        pass  # 已提示错误
    elif cand.empty:
        st.warning(f"No human structure found on RCSB for {gene}. Pick another gene.")
    else:
        st.dataframe(cand, hide_index=True, use_container_width=True)
        pdb_id = st.selectbox("Select PDB structure (research decision point)", cand["pdb_id"].tolist(), key="pdb_sel")
        # PDB 选择变化 → 作废旧受体与 Step 4 盒/预览（防止 Step 4 残留旧结构的预览图）
        if dock.get("pdb_sel") != pdb_id:
            dock["pdb_sel"] = pdb_id
            for _k in ("rec_pdbqt", "pdb_id", "pdb_path", "box", "result"):
                dock.pop(_k, None)
        if st.button("Download & prepare receptor", key="prep_rec"):
            with st.spinner("Downloading PDB and converting to pdbqt…"):
                try:
                    pdb_path = _download_pdb(pdb_id)
                    dock["pdb_id"] = pdb_id
                    dock["pdb_path"] = pdb_path
                    dock["rec_pdbqt"] = Utils.prepare_receptor_pdbqt(pdb_path)
                    dock.pop("result", None)
                except Exception as exc:
                    st.error(f"Receptor preparation failed: {exc}")
else:
    if dock.get("genes") is None:
        st.info("Complete Step 1 first: select a PPI result directory / enrichment file (then run the gene-by-gene RCSB structure check).")
    else:
        st.warning("The candidate gene list is empty (no HUB/enrichment intersection, or no structures passed the filter); a receptor cannot be selected. Return to Step 1 to adjust the data source.")

if dock.get("rec_pdbqt"):
    st.success(f"✅ Receptor pdbqt: `{dock['rec_pdbqt']}` (PDB: {dock.get('pdb_id')})")

st.divider()

# ==================== Step 4 对接盒辅助 ====================


def _parse_hetatm(pdb_path: str) -> list[tuple[str, str, str, tuple]]:
    """解析 PDB 中非水 HETATM（配体/辅因子）：[(resname, chain, resseq, (x,y,z)), ...]。"""
    het = []
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            resn = line[17:20].strip()
            if resn in ("HOH", "H2O", "WAT", "DOD"):
                continue
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            het.append((resn, line[21], line[22:26].strip(), xyz))
    return het


def _residues_in_box(pdb_path: str, box: dict) -> list[str]:
    """盒内蛋白残基列表（任意原子落入盒内即计入），返回 ['A123:TYR', ...]。"""
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


def _render_box_3d(pdb_path: str, box: dict) -> None:
    """py3Dmol 内嵌：受体卡通 + 共结晶配体（绿色棍状）+ 对接盒（橙色半透明线框）。

    蛋白（ATOM）与配体（HETATM）分成两个 model 渲染，确保配体棍状可靠显示；
    apo 结构（无 HETATM 配体）时只有蛋白卡通与盒框。
    """
    try:
        import py3Dmol
    except ImportError:
        st.info("3D box preview requires py3Dmol (`pip install py3Dmol`).")
        return
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM", "TER", "END"))]
    het_lines = [ln for ln in lines
                 if ln.startswith("HETATM") and ln[17:20].strip() not in ("HOH", "H2O", "WAT", "DOD")]

    c, s = box["center"], box["size"]
    view = py3Dmol.view(width=900, height=520)
    view.addModel("\n".join(atom_lines) + "\n", "pdb")  # model 0：蛋白
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
    if het_lines:
        view.addModel("\n".join(het_lines) + "\n", "pdb")  # model 1：共结晶配体
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "greenCarbon"}})
    # 对接盒：3Dmol.js 原生 box 图元（橙色半透明线框）
    view.addBox({
        "center": {"x": c[0], "y": c[1], "z": c[2]},
        "dimensions": {"w": s[0], "h": s[1], "d": s[2]},
        "color": "orange",
        "opacity": 0.6,
    })
    view.zoomTo()
    try:
        html = view._make_html()
    except AttributeError:
        html = view.write()
    st.components.v1.html(html, height=540)


# ==================== Step 4 对接盒（决策点） ====================
st.subheader("4️⃣ Docking box (manual confirmation required after inference)")

if dock.get("pdb_path"):
    if st.button("Infer docking box from ligand/cofactor", key="infer_box"):
        try:
            dock["box"] = Utils.infer_box_from_cofactor(dock["pdb_path"])
        except Exception as exc:
            st.error(f"Docking box inference failed: {exc}")
    if dock.get("box"):
        # 推断依据说明
        het = _parse_hetatm(dock["pdb_path"])
        if het:
            ligs = sorted({f"{r} (chain {c}{n})" for r, c, n, _ in het})
            st.caption(f"Inferred from: centroid of co-crystallized ligand/cofactor {', '.join(ligs[:5])}, expanded by 8 Å to cover the pocket.")
        else:
            st.caption("Inferred from: no co-crystallized ligand (apo structure), falling back to a bounding box of the whole protein. Adjust the center manually based on literature / known active sites.")

        # 微调指南
        st.info("Fine-tuning guide: ① center on the active pocket (see the co-crystal ligand or literature sites); "
                "② size should fully enclose the ligand plus surrounding residues (typical edge length 20–25 Å); "
                "③ after adjusting, check \"Residues in box\" below and the 3D box preview to confirm coverage of target residues.")

        c, s = dock["box"]["center"], dock["box"]["size"]

        # 左：3D 盒框预览；右：微调选项（并排，方便边看边调）
        col_box3d, col_tune = st.columns([2, 1])
        with col_box3d:
            st.markdown("**3D box preview (orange wireframe = docking box; green sticks = co-crystal ligand)**")
            _render_box_3d(dock["pdb_path"], dock["box"])
        with col_tune:
            st.caption("Fine-tune the docking box (values are editable):")
            cc1, cc2, cc3 = st.columns(3)
            with cc1:
                c1 = st.number_input("cx", value=float(c[0]), step=0.5, format="%.3f", key="cx_x")
            with cc2:
                c2 = st.number_input("cy", value=float(c[1]), step=0.5, format="%.3f", key="cx_y")
            with cc3:
                c3 = st.number_input("cz", value=float(c[2]), step=0.5, format="%.3f", key="cx_z")
            cs1, cs2, cs3 = st.columns(3)
            with cs1:
                s1 = st.number_input("sx", value=float(s[0]), step=1.0, format="%.3f", key="sx_x")
            with cs2:
                s2 = st.number_input("sy", value=float(s[1]), step=1.0, format="%.3f", key="sx_y")
            with cs3:
                s3 = st.number_input("sz", value=float(s[2]), step=1.0, format="%.3f", key="sx_z")
            dock["box"] = {"center": (c1, c2, c3), "size": (s1, s2, s3)}
            st.code(f"center = {dock['box']['center']}\nsize   = {dock['box']['size']}")

        # 盒内残基（辅助判断口袋覆盖；默认只显示前 30 个，避免刷屏）
        residues = _residues_in_box(dock["pdb_path"], dock["box"])
        if residues:
            st.caption(f"Residues covered by the current box ({len(residues)} total, showing first 30): "
                       + ", ".join(residues[:30])
                       + (" …" if len(residues) > 30 else ""))
            if len(residues) > 30:
                with st.expander(f"View all {len(residues)} residues"):
                    st.write(", ".join(residues))
        else:
            st.caption("No protein residues inside the current box — it may be off the protein. Adjust the center to cover the binding site.")
else:
    st.info("Complete Step 3 to prepare the receptor first.")

st.divider()

# ==================== Step 5 运行对接 ====================
st.subheader("5️⃣ Run AutoDock Vina")

if dock.get("rec_pdbqt") and dock.get("lig_pdbqt") and dock.get("box"):
    st.caption(f"Current parameters: exhaustiveness={config.EXHAUSTIVENESS} · "
               f"num_modes={config.NUM_MODES} · cpu={config.DOCKING_CPU} · "
               f"seed={config.DOCKING_SEED} · timeout={config.DOCKING_TIMEOUT}s")
    if st.button("🚀 Run docking", key="run_dock"):
        with st.spinner("Vina searching (stochastic, may take several minutes)…"):
            try:
                out_pdbqt = os.path.join(config.DOCKING_RAW_CACHE_DIR,
                                         f"vina_out_{datetime.now():%Y%m%d_%H%M%S}.pdbqt")
                proc = Utils.run_vina(dock["rec_pdbqt"], dock["lig_pdbqt"], dock["box"],
                                      out_path=out_pdbqt)
                poses = Utils.split_vina(out_pdbqt, out_dir=os.path.dirname(out_pdbqt))
                scores = Utils.parse_vina_scores(out_pdbqt)
                dock["result"] = {"proc": proc, "out": out_pdbqt, "poses": poses, "scores": scores}
            except Exception as exc:
                st.error(f"Docking failed: {exc}")

    if dock.get("result"):
        res = dock["result"]
        st.caption("Vina run log (tail)")
        st.code((res["proc"].stdout or "")[-1500:], language="text")
        st.dataframe(res["scores"], use_container_width=True, hide_index=True)
else:
    st.info("Complete ligand, receptor, and docking box preparation first.")

st.divider()

# ==================== Step 6 结果解读与归档（决策点） ====================
st.subheader("6️⃣ Result interpretation & archiving")

if dock.get("result"):
    scores = dock["result"]["scores"]
    if scores.empty:
        st.warning("No binding modes parsed; cannot archive")
    else:
        pose_sel = st.selectbox("Select pose to save (best by default)",
                                scores["mode"].tolist(), index=0, key="pose_sel")
        conclusion = st.text_area("Conclusion (scientific interpretation, written to metadata and summary)",
                                  key="conclusion", height=90,
                                  placeholder=f"e.g. forms stable hydrogen bonds with {dock.get('gene', 'the target protein')}, "
                                              f"predicted binding affinity {scores.iloc[0]['affinity']:.2f} kcal/mol …")
        drug = st.text_input("Drug (used as archive prefix)", value=default_drug, key="drug_name")
        disease = st.text_input("Disease (used as archive prefix)", value=default_disease, key="disease_name")

        col_save, col_star = st.columns([3, 1])
        with col_save:
            if st.button("💾 Save & archive result", key="save_result"):
                gene = dock.get("gene", "gene")
                prefix = f"{drug}_{disease}_{gene}".strip("_").replace(" ", "_")
                try:
                    result_dir = _save_docking_result(
                        prefix, gene, dock.get("pdb_id", ""), dock.get("cid"),
                        dock["box"], dock["result"]["out"], scores, conclusion)
                    dock["result_dir"] = result_dir
                    dock.pop("starred", None)  # 新成果 → 重置收藏状态
                    st.success("✅ Saved and archived")
                except Exception as exc:
                    st.error(f"Save & archive failed: {exc}")
        with col_star:
            # ☆/★ 收藏切换：☆ 收藏到 star/（不受自动删除影响）；★ 点击即取消收藏（删除副本）
            if dock.get("result_dir"):
                if dock.get("starred"):
                    if st.button("★ Starred", key="star_on"):
                        try:
                            star_base = os.path.basename(dock["result_dir"])
                            import glob
                            for p in glob.glob(os.path.join(config.MD_STAR_DIR, star_base) + "*"):
                                if os.path.isdir(p):
                                    shutil.rmtree(p, ignore_errors=True)
                            dock.pop("starred", None)
                            st.rerun()  # 立即更新图标为 ☆
                        except Exception as exc:
                            st.error(f"Unstar failed: {exc}")
                else:
                    if st.button("☆ Star", key="star_off"):
                        try:
                            Utils.star_docking_result(dock["result_dir"])
                            dock["starred"] = True
                            st.rerun()  # 立即更新图标为 ★
                        except Exception as exc:
                            st.error(f"Starring failed: {exc}")

        if dock.get("result_dir"):
            st.write(f"Result directory: `{dock['result_dir']}`")
            st.download_button("⬇️ Download all results (zip)",
                               data=Utils.make_dir_zip(dock["result_dir"]),
                               file_name=os.path.basename(dock["result_dir"]) + ".zip",
                               mime="application/zip", key="dl_zip")

        # 3D 内嵌可视化（py3Dmol 可选），渲染用户选择的 pose
        st.markdown("**3D structure preview**")
        if dock.get("pdb_path"):
            try:
                pose_path = os.path.join(config.DOCKING_RAW_CACHE_DIR,
                                         f"pose_{int(pose_sel)}.pdbqt")
                _extract_model(dock["result"]["out"], int(pose_sel), pose_path)
                _render_3d(dock["pdb_path"], pose_path)
            except Exception as exc:
                st.warning(f"3D preview failed: {exc}")
else:
    st.info("Run a docking to interpret and archive results here.")

st.divider()
st.caption("Cache policy: raw files (SDF/PDB/vina output) are auto-cleaned with a 24-hour TTL; "
           "processed ligand/receptor pdbqt files with a 7-day TTL; results under "
           "downloads/MDresult keep at most 3 versions (oldest auto-deleted). "
           "Starred results can be viewed and managed on the Star Management page.")
