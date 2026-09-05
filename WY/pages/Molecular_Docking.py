import os
import re
import shutil
from datetime import datetime

import pandas as pd
import streamlit as st

import config
import Utils
from i18n import t
from core.docking import (
    clamp_box,
    download_pdb,
    extract_model,
    fetch_alphafold_receptor,
    interactive_3d_html,
    parse_hetatm,
    pose_pdbqt_to_pdb,
    residues_in_box,
    resolve_cid,
    save_interactive_3d_html,
    smiles_to_sdf_3d,
)

st.set_page_config(page_title=t("dock.page_title"), layout="wide")

# ---------- 隐藏侧边栏,全屏显示页面;Star Management 按钮(#6347A1 底、#E0D8F2 字) ----------
Utils.hide_sidebar("""
    .st-key-star_mgmt_btn button {
        background-color: #6347A1;
        color: #E0D8F2;
        border: 1px solid #6347A1;
    }
""")

st.title(t("dock.title"))
st.caption(t("dock.caption"))

col_back, col_star = st.columns([1, 1])
with col_back:
    Utils.back_button()
with col_star:
    if st.button(t("dock.star_mgmt"), key="star_mgmt_btn"):
        st.switch_page("pages/Star_Management.py")

st.divider()

# 页面状态:从主页重新进入本页时全部重置,页内操作则保留(不重跑已完成步骤)。
if st.session_state.pop("_docking_entry", False):
    st.session_state.pop("docking", None)
    for _k in ("sel_ppi_dir", "sel_enr", "gene_sel", "pdb_sel", "lig_source",
               "lig_query", "lig_smiles", "drug_name", "disease_name",
               "use_struct", "pose_sel", "conclusion", "manual_genes"):
        st.session_state.pop(_k, None)

st.session_state.setdefault("docking", {})
dock = st.session_state["docking"]

# 每次进入本页时清理一次过期缓存(页内操作不重复清理)
if not st.session_state.get("_docking_cache_cleaned"):
    Utils.clean_docking_cache()
    st.session_state["_docking_cache_cleaned"] = True


# ===== 辅助函数 =====

def _extract_label(file_name: str) -> str:
    """从富集结果文件名中解析出"药物_疾病"标识,用于命名后续结果。"""
    m = re.match(r"^(?:gProfiler|Enrichr)_(.+)_[0-9a-f]{32}\.csv$", file_name)
    if m:
        return m.group(1)
    m2 = re.match(r"^(?:gProfiler|Enrichr)_([0-9a-f]{32})\.csv$", file_name)
    return m2.group(1) if m2 else ""


def _parse_genes(raw: str) -> list[str]:
    """把富集结果里的基因列表字段解析成基因符号列表,兼容列表/空格/逗号等格式,避免同一基因被重复计数。"""
    raw = (raw or "").strip()
    if not raw:
        return []
    if raw.startswith("[") and raw.endswith("]"):
        try:
            import ast
            items = ast.literal_eval(raw)
            if isinstance(items, list):
                return [str(x).strip() for x in items if str(x).strip()]
        except Exception:
            pass
    cleaned = raw.strip("[]").replace("'", "").replace('"', "")
    return [t.strip() for t in re.split(r"[,;\s]+", cleaned) if t.strip()]


def _genes_from_df(df: pd.DataFrame) -> set[str]:
    """汇总富集表中各条目涉及的全部基因,合并去重后返回。"""
    genes: set[str] = set()
    for raw in df.get("genes", pd.Series(dtype=str)).fillna("").astype(str):
        genes.update(g.upper() for g in _parse_genes(raw))
    return genes


def _reset_docking_downstream() -> None:
    """切换候选基因数据源后,清空配体、受体、对接盒和对接结果,
    避免上一份数据的旧结果混入新分析。"""
    for _k in ("lig_pdbqt", "cid", "lig_auto_done",
               "gene", "rec_pdbqt", "pdb_id", "pdb_path",
               "box", "result", "result_dir", "starred",
               "gene_structures", "candidate_pre", "batch_results"):
        dock.pop(_k, None)


def _read_ppi_dataset(dir_name: str) -> tuple[set[str], set[str]]:
    """读取所选 PPI 结果目录,返回 HUB 基因集合与富集验证基因集合;
    对应文件缺失时返回空集合。"""
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


def _save_docking_result(prefix: str, gene: str, pdb_id: str, cid: int | None, box: dict,
                         out_pdbqt: str, scores: pd.DataFrame, conclusion: str,
                         receptor_source: str = "experiment") -> str:
    """归档本次对接成果:打分表、最佳构象和参数信息存入结果文件夹,并清理旧版本。"""
    gene = Utils.safe_name(gene)  # 清理基因名中的非法字符,再用于文件命名
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    result_dir = os.path.join(config.MD_RESULT_DIR, f"{prefix}_{ts}")
    os.makedirs(result_dir, exist_ok=True)

    # 打分表和结论
    summary = scores.copy()
    summary["conclusion"] = conclusion
    summary_path = os.path.join(result_dir, f"summary_{ts}.csv")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    # 最佳构象:按序号从 Vina 原始输出中提取
    best_mode = int(summary.iloc[0]["mode"])
    best_path = os.path.join(result_dir, f"best_{gene}_{ts}.pdbqt")
    extract_model(out_pdbqt, best_mode, best_path)

    # 参数信息(药物、靶点、对接盒、打分等),便于日后复查
    Utils.write_docking_metadata(os.path.join(result_dir, f"meta_{ts}.txt"), {
        "drug": prefix.split("_")[0] if "_" in prefix else prefix,
        "disease": prefix.split("_", 1)[1] if "_" in prefix else "",
        "gene": gene,
        "pdb_id": pdb_id,
        "receptor_source": receptor_source,
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

    # 同名结果只保留最近几版,并清理过期中间文件
    Utils.prune_docking_results(result_dir)
    Utils.clean_docking_cache()
    return result_dir


def _render_3d(pdb_path: str, pose_pdbqt: str) -> None:
    """在页面内 3D 展示受体(卡通)与所选对接构象(棍状);缺少组件时仅给出提示。"""
    try:
        import py3Dmol
    except ImportError:
        st.info(t("dock.py3dmol_pose"))
        return
    pose_pdb = os.path.join(os.path.dirname(pose_pdbqt), "pose_view.pdb")
    try:
        pose_pdbqt_to_pdb(pose_pdbqt, pose_pdb)
    except Exception:
        st.warning(t("dock.pose_pdb_failed"))
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


def _save_3d_html_ui(dock: dict, pdb_path: str, pose_pdbqt: str, pose_sel) -> None:
    """「保存交互式 3D 视图」:把受体与当前构象打包成单个自包含 HTML(内嵌 3Dmol.js 与坐标),
    可下载或存入结果目录;之后双击即可在任意浏览器打开,可拖动旋转,离线可用。"""
    gene = Utils.safe_name(str(dock.get("gene") or "target"))
    pose_n = int(pose_sel)
    caption = t("dock.save3d_caption", gene=dock.get("gene") or "?", pose=pose_n,
                drug=dock.get("drug") or "-", disease=dock.get("disease") or "-")
    try:
        html_text, offline_ok = interactive_3d_html(pdb_path, pose_pdbqt, caption)
    except Exception as exc:
        st.warning(t("dock.save3d_failed", exc=exc))
        return
    file_name = f"3d_view_{gene}_pose{pose_n}.html"

    def _download(col) -> None:
        with col:
            st.download_button(t("dock.save3d_dl"), data=html_text.encode("utf-8"),
                               file_name=file_name, mime="text/html", key="save3d_dl")

    if dock.get("result_dir"):
        col_arc, col_dl = st.columns(2)
        with col_arc:
            if st.button(t("dock.save3d_btn"), key="save3d_btn"):
                try:
                    out_path, offline_ok = save_interactive_3d_html(
                        pdb_path, pose_pdbqt,
                        os.path.join(dock["result_dir"], file_name), caption)
                    key = "dock.save3d_ok" if offline_ok else "dock.save3d_ok_online"
                    st.success(t(key, path=out_path))
                except Exception as exc:
                    st.error(t("dock.save3d_failed", exc=exc))
        _download(col_dl)
    else:
        # 尚未归档时没有结果目录可存,只提供直接下载
        _download(st.container())


@st.cache_data(ttl=3600)
def _cached_pdb_candidates(gene: str, allow_alphafold: bool = False) -> pd.DataFrame:
    """查询某基因可用的受体结构列表(近期查过则复用);无实验结构时可附加 AlphaFold 预测结构。"""
    return Utils.fetch_receptor_pdb_candidates(gene, max_hits=20,
                                               allow_alphafold=allow_alphafold)


@st.cache_data(ttl=3600)
def _cached_pdb_candidates_batch(genes: tuple, allow_alphafold: bool = False) -> dict[str, pd.DataFrame]:
    """一次批量查询多个基因的受体结构(近期查过则复用);返回 {基因: 候选表}。

    搜索合并成少数几次批量请求 + GraphQL 批量归因,速度与基因数基本无关;
    查无实验结构的基因按开关附加 AlphaFold 预测结构。
    """
    return Utils.fetch_receptor_candidates_batch(list(genes), max_hits=20,
                                                 allow_alphafold=allow_alphafold)


# ===== Step 1 候选基因筛选 =====
st.subheader(t("dock.step1"))

src_mode = st.radio(t("dock.src_mode"),
                    ["ppi", "enr"],
                    format_func=lambda v: t("dock.src_ppi") if v == "ppi" else t("dock.src_enr"),
                    horizontal=True, key="gene_src_mode")

# 切换数据来源时清空旧选择和后续结果,按新来源重新读取
if dock.get("gene_src_mode") != src_mode:
    dock["gene_src_mode"] = src_mode
    for _k in ("ppi_dir", "enr_file", "hub_genes", "enrichment_genes",
               "label", "drug", "disease", "genes"):
        dock.pop(_k, None)
    _reset_docking_downstream()

if src_mode == "ppi":
    # 只列出 PPI 分析结果目录,排除其他类型的文件夹
    dirs = Utils.list_download_dirs()
    ppi_dirs = dirs[dirs["Directory"].str.endswith("_data", na=False)] if not dirs.empty else dirs
    if ppi_dirs.empty:
        st.warning(t("dock.no_ppi_dirs"))
    else:
        sel_dir = st.selectbox(t("dock.sel_ppi_dir"), ppi_dirs["Directory"].tolist(), key="sel_ppi_dir",
                               index=None, placeholder=t("dock.sel_ppi_dir_ph"))
        # 选定目录后自动读取,无需再点按钮
        if sel_dir and dock.get("ppi_dir") != sel_dir:
            dock["ppi_dir"] = sel_dir
            dock["hub_genes"], dock["enrichment_genes"] = _read_ppi_dataset(sel_dir)
            # 从目录名拆出药物名和疾病名,供 Step 2 自动获取配体
            prefix = sel_dir[:-5] if sel_dir.endswith("_data") else sel_dir
            dock["label"] = prefix
            parts = prefix.split("_", 1)
            dock["drug"] = parts[0] if parts[0] else ""
            dock["disease"] = parts[1] if len(parts) > 1 else ""
            _reset_docking_downstream()  # 数据源变了,清空后面步骤的旧结果
else:
    enr_files = Utils.list_cache_files("enrichment")
    if enr_files.empty:
        st.warning(t("dock.no_enr_files"))
    else:
        sel_file = st.selectbox(t("dock.sel_enr_file"), enr_files["File name"], key="sel_enr",
                                index=None, placeholder=t("dock.sel_enr_file_ph"))
        # 选定文件后自动读取
        if sel_file and dock.get("enr_file") != sel_file:
            dock["enr_file"] = sel_file
            df = pd.read_csv(os.path.join(config.ENRICHMENT_CACHE_DIR, sel_file),
                             encoding="utf-8-sig", comment="#")
            dock["enrichment_genes"] = sorted(_genes_from_df(df))
            dock["hub_genes"] = []  # 没有 HUB 信息,只按富集基因筛选
            label = _extract_label(sel_file)
            dock["label"] = label
            parts = label.split("_", 1)
            dock["drug"] = parts[0] if parts[0] else ""
            dock["disease"] = parts[1] if len(parts) > 1 else ""
            _reset_docking_downstream()  # 数据源变了,清空后面步骤的旧结果

if dock.get("enrichment_genes") is not None:
    hub_set = set(dock.get("hub_genes") or [])
    enr_set = set(dock["enrichment_genes"])
    candidate = sorted(hub_set & enr_set) if hub_set else sorted(enr_set)
    dock["candidate_pre"] = candidate

    col_h, col_e, col_c = st.columns(3)
    col_h.metric(t("dock.hub_metric"), len(hub_set))
    col_e.metric(t("dock.enr_metric"), len(enr_set))
    col_c.metric(t("dock.cand_metric"), len(candidate))
    if candidate:
        st.write(", ".join(candidate[:60]) + (" …" if len(candidate) > 60 else ""))
    else:
        st.warning(t("dock.no_intersection"))

    # 结构可用性筛选:默认开启;若想保留查不到结构的基因,可关闭
    st.markdown(t("dock.struct_filter_title"))
    use_struct = st.checkbox(t("dock.enable_struct_filter"), key="use_struct",
                             disabled=not candidate, value=True)
    # AF 预测结构开关(默认关):没有实验结构的基因也能用预测结构对接
    allow_af = config.ALPHAFOLD_ENABLED and st.checkbox(
        t("dock.af_toggle"), key="allow_af", value=False, help=t("dock.af_toggle_help"))
    if use_struct:
        if st.button(t("dock.check_struct_btn"), key="check_struct"):
            with st.spinner(t("dock.checking_rcsb")):
                try:
                    # 批量检索:4 路并发逐基因搜索(召回与旧版一致),
                    # 元数据合并为 GraphQL 批量请求;进度条按两阶段推进
                    progress = st.progress(0.0)
                    results = _cached_pdb_candidates_batch(
                        tuple(candidate), allow_alphafold=allow_af)
                    progress.progress(1.0)
                    dock["gene_structures"] = {g: len(df) for g, df in results.items()
                                               if not df.empty}
                    dock["batch_results"] = results
                except Exception as exc:
                    # 查询失败时提示错误原因
                    dock.pop("gene_structures", None)
                    st.error(t("dock.error_occurred", exc=exc))
        if dock.get("gene_structures") is not None:
            struct_ok = sorted(dock["gene_structures"])
            dock["genes"] = struct_ok
            no_struct = [g for g in candidate if g not in dock["gene_structures"]]
            st.success(t("dock.genes_with_struct", n=len(struct_ok))
                       + ", ".join(struct_ok[:40]) + (" …" if len(struct_ok) > 40 else ""))
            if no_struct:
                st.caption(t("dock.excluded_no_struct", n=len(no_struct))
                           + ", ".join(no_struct[:40]) + (" …" if len(no_struct) > 40 else ""))
        else:
            dock.pop("genes", None)
            st.info(t("dock.click_to_check"))
    else:
        # 不做结构筛选时,候选基因直接作为最终靶点列表
        dock["genes"] = candidate
        if candidate:
            st.caption(t("dock.final_candidates", n=len(candidate))
                       + ", ".join(candidate[:60]) + (" …" if len(candidate) > 60 else ""))

# 没有可选数据文件时,允许手动输入基因列表
if dock.get("enrichment_genes") is None:
    manual = st.text_area(t("dock.manual_genes_label"), key="manual_genes")
    if st.button(t("dock.use_manual_btn"), key="load_genes_manual"):
        dock["genes"] = sorted({g.strip().upper() for g in re.split(r"[,;\s]+", manual) if g.strip()})
        dock["label"] = ""
        dock["drug"] = ""
        dock["disease"] = ""
        dock["hub_genes"], dock["enrichment_genes"] = [], []
        _reset_docking_downstream()  # 数据源变了,清空后面步骤的旧结果

if dock.get("genes"):
    shown = ", ".join(dock["genes"][:60]) + (" …" if len(dock["genes"]) > 60 else "")
    st.info(t("dock.candidates_confirmed", n=len(dock["genes"]), lst=shown))

st.divider()

# ===== Step 2 配体 =====
st.subheader(t("dock.step2"))

# 药物名/疾病名默认取自 Step 1 所选的数据
default_drug = dock.get("drug", "")
default_disease = dock.get("disease", "")

lig_source = st.radio(t("dock.lig_source"),
                      ["pubchem", "smiles"],
                      format_func=lambda v: t("dock.lig_src_pubchem") if v == "pubchem" else t("dock.lig_src_smiles"),
                      horizontal=True, key="lig_source")
if lig_source == "pubchem":
    # 已识别到药物名时自动从 PubChem 获取配体,只自动执行一次
    if default_drug and not dock.get("lig_pdbqt") and not dock.get("lig_auto_done"):
        with st.spinner(t("dock.auto_fetching", drug=default_drug)):
            cid = resolve_cid(default_drug)
            if cid is None:
                dock["lig_auto_done"] = True
                st.error(t("dock.pubchem_not_found", drug=default_drug))
            else:
                sdf_path = Utils.fetch_ligand_sdf_3d(cid)
                if not sdf_path:
                    dock["lig_auto_done"] = True
                    st.error(t("dock.cid_no_3d", cid=cid))
                else:
                    dock["cid"] = cid
                    dock["lig_pdbqt"] = Utils.prepare_ligand_pdbqt(sdf_path)

    # 只有自动获取失败或没有药物名时,才显示手动输入框
    lig_ready = bool(dock.get("lig_pdbqt"))
    show_lig_input = not lig_ready and (not default_drug or dock.get("lig_auto_done"))
    if show_lig_input:
        if not default_drug:
            st.caption(t("dock.no_drug_recognized"))
        lig_query = st.text_input(t("dock.lig_query_label"), value=default_drug, key="lig_query")
        if st.button(t("dock.fetch_ligand_btn"), key="prep_lig"):
            dock["lig_auto_done"] = True  # 手动指定后不再自动获取
            with st.spinner(t("dock.fetching_pubchem")):
                cid = resolve_cid(lig_query)
                if cid is None:
                    st.error(t("dock.pubchem_no_compound", q=lig_query))
                else:
                    sdf_path = Utils.fetch_ligand_sdf_3d(cid)
                    if not sdf_path:
                        st.error(t("dock.cid_try_another", cid=cid))
                    else:
                        dock["cid"] = cid
                        dock["lig_pdbqt"] = Utils.prepare_ligand_pdbqt(sdf_path)
else:
    smiles = st.text_area("SMILES", height=90, key="lig_smiles",
                          placeholder=t("dock.smiles_placeholder"))
    if st.button(t("dock.gen_ligand_btn"), key="prep_lig_b"):
        with st.spinner(t("dock.gen_3d_spinner")):
            try:
                tmp_dir = config.DOCKING_RAW_CACHE_DIR
                os.makedirs(tmp_dir, exist_ok=True)
                sdf_path = smiles_to_sdf_3d(
                    smiles,
                    os.path.join(tmp_dir, f"ligand_smiles_{datetime.now():%Y%m%d_%H%M%S}.sdf"))
                dock["cid"] = None
                dock["lig_pdbqt"] = Utils.prepare_ligand_pdbqt(sdf_path)
            except Exception as exc:
                st.error(t("dock.ligand_failed", exc=exc))

if dock.get("lig_pdbqt"):
    st.success(t("dock.ligand_ready", path=dock["lig_pdbqt"])
               + (t("dock.cid_suffix", cid=dock.get("cid")) if dock.get("cid") else t("dock.smiles_suffix")))

st.divider()

# ===== Step 3 受体(基因 + PDB,决策点) =====
st.subheader(t("dock.step3"))

if dock.get("genes"):
    gene = st.selectbox(t("dock.gene_label"), dock["genes"], key="gene_sel")
    # 换靶点基因后,清空旧受体和对接结果
    if dock.get("gene") != gene:
        dock["gene"] = gene
        dock.pop("rec_pdbqt", None)
        dock.pop("pdb_path", None)
        dock.pop("pdb_id", None)
        dock.pop("box", None)
        dock.pop("result", None)
        dock.pop("receptor_source", None)

    # 自动查询该基因的候选结构:优先复用 Step 1 批量结果,查不到再单基因查
    allow_af = st.session_state.get("allow_af", False) and config.ALPHAFOLD_ENABLED
    try:
        cand = (dock.get("batch_results") or {}).get(gene)
        if cand is None:
            cand = _cached_pdb_candidates(gene, allow_alphafold=allow_af)
    except Exception as exc:
        st.warning(t("dock.rcsb_error", exc=exc))
        cand = None
    if cand is None:
        pass  # 已提示错误
    elif cand.empty:
        st.warning(t("dock.no_human_struct", gene=gene))
    else:
        st.dataframe(cand, hide_index=True, use_container_width=True)
        pdb_id = st.selectbox(t("dock.sel_pdb"), cand["pdb_id"].tolist(), key="pdb_sel")
        # 更换 PDB 结构后,作废旧受体和后续步骤的旧结果
        if dock.get("pdb_sel") != pdb_id:
            dock["pdb_sel"] = pdb_id
            for _k in ("rec_pdbqt", "pdb_id", "pdb_path", "box", "result", "receptor_source"):
                dock.pop(_k, None)
        if st.button(t("dock.dl_receptor_btn"), key="prep_rec"):
            with st.spinner(t("dock.dl_pdb_spinner")):
                try:
                    if pdb_id.startswith("AFDB:"):
                        # 通过该基因对应的 UniProt 条目获取 AlphaFold 预测模型
                        pdb_path, af_acc = fetch_alphafold_receptor(
                            gene, taxon=int(config.SPECIES_MAP["Human"]["ncbi_taxon"]))
                        dock["pdb_id"] = f"AFDB:{af_acc}"
                        dock["receptor_source"] = "alphafold"
                    else:
                        pdb_path = download_pdb(pdb_id)
                        dock["pdb_id"] = pdb_id
                        dock["receptor_source"] = "experiment"
                    dock["pdb_path"] = pdb_path
                    dock["rec_pdbqt"] = Utils.prepare_receptor_pdbqt(pdb_path)
                    dock.pop("result", None)
                except Exception as exc:
                    st.error(t("dock.receptor_failed", exc=exc))
        if dock.get("receptor_source") == "alphafold" and dock.get("rec_pdbqt"):
            # 使用预测结构的提示(只要用的是 AF 受体就持续显示)
            st.warning(t("dock.af_disclaimer",
                         acc=str(dock.get("pdb_id", "")).split(":", 1)[-1],
                         version=config.ALPHAFOLD_MODEL_VERSION))
else:
    if dock.get("genes") is None:
        st.info(t("dock.complete_step1"))
    else:
        st.warning(t("dock.no_candidates"))

if dock.get("rec_pdbqt"):
    st.success(t("dock.receptor_ready", path=dock["rec_pdbqt"], pdb=dock.get("pdb_id")))

st.divider()

# ===== Step 4 对接盒辅助函数 =====


def _render_box_3d(pdb_path: str, box: dict) -> None:
    """页面内 3D 预览:受体(卡通)、共结晶配体(绿棍)、对接盒(橙框);无共晶配体时只显示受体和盒子。"""
    try:
        import py3Dmol
    except ImportError:
        st.info(t("dock.py3dmol_box"))
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
    # 对接盒:橙色半透明线框
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


# ===== Step 4 对接盒(决策点) =====
st.subheader(t("dock.step4"))

if dock.get("pdb_path"):
    if st.button(t("dock.infer_box_btn"), key="infer_box"):
        try:
            dock["box"] = Utils.infer_box_from_cofactor(dock["pdb_path"])
            # 预测结构没有共晶配体可参考,按整蛋白范围估盒子;过大时收缩到可计算大小
            if (dock.get("receptor_source") == "alphafold"
                    and max(dock["box"]["size"]) > config.DOCKING_BOX_MAX_SIZE_AF):
                dock["box"] = clamp_box(dock["box"], config.DOCKING_BOX_MAX_SIZE_AF)
                st.warning(t("dock.af_box_clamped", m=config.DOCKING_BOX_MAX_SIZE_AF))
        except Exception as exc:
            st.error(t("dock.box_infer_failed", exc=exc))
    if dock.get("box"):
        # 说明盒子是参考哪个配体推断的
        het = parse_hetatm(dock["pdb_path"])
        if het:
            ligs = sorted({f"{r} (chain {c}{n})" for r, c, n, _ in het})
            st.caption(t("dock.inferred_ligand", ligs=", ".join(ligs[:5])))
        else:
            st.caption(t("dock.inferred_apo"))

        # 微调指南
        st.info(t("dock.finetune_guide"))

        c, s = dock["box"]["center"], dock["box"]["size"]

        # 左侧 3D 预览盒子,右侧微调参数,方便对照调整
        col_box3d, col_tune = st.columns([2, 1])
        with col_tune:
            st.caption(t("dock.finetune_caption"))
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
        with col_box3d:
            st.markdown(t("dock.box_preview_title"))
            _render_box_3d(dock["pdb_path"], dock["box"])

        # 列出盒内残基,帮助判断盒子是否覆盖结合口袋
        residues = residues_in_box(dock["pdb_path"], dock["box"])
        if residues:
            st.caption(t("dock.residues_caption", n=len(residues))
                       + ", ".join(residues[:30])
                       + (" …" if len(residues) > 30 else ""))
            if len(residues) > 30:
                with st.expander(t("dock.view_all_residues", n=len(residues))):
                    st.write(", ".join(residues))
        else:
            st.caption(t("dock.no_residues"))
else:
    st.info(t("dock.complete_step3"))

st.divider()

# ===== Step 5 运行对接 =====
st.subheader(t("dock.step5"))

if dock.get("rec_pdbqt") and dock.get("lig_pdbqt") and dock.get("box"):
    st.caption(t("dock.current_params", e=config.EXHAUSTIVENESS, m=config.NUM_MODES,
                 c=config.DOCKING_CPU, s=config.DOCKING_SEED, timeout=config.DOCKING_TIMEOUT))
    if st.button(t("dock.run_dock_btn"), key="run_dock"):
        with st.spinner(t("dock.vina_spinner")):
            try:
                out_pdbqt = os.path.join(config.DOCKING_RAW_CACHE_DIR,
                                         f"vina_out_{datetime.now():%Y%m%d_%H%M%S}.pdbqt")
                proc = Utils.run_vina(dock["rec_pdbqt"], dock["lig_pdbqt"], dock["box"],
                                      out_path=out_pdbqt)
                poses = Utils.split_vina(out_pdbqt, out_dir=os.path.dirname(out_pdbqt))
                scores = Utils.parse_vina_scores(out_pdbqt)
                dock["result"] = {"proc": proc, "out": out_pdbqt, "poses": poses, "scores": scores}
            except Exception as exc:
                st.error(t("dock.dock_failed", exc=exc))

    if dock.get("result"):
        res = dock["result"]
        st.caption(t("dock.vina_log"))
        st.code((res["proc"].stdout or "")[-1500:], language="text")
        st.dataframe(res["scores"], use_container_width=True, hide_index=True)
else:
    st.info(t("dock.complete_prep"))

st.divider()

# ===== Step 6 结果解读与归档(决策点) =====
st.subheader(t("dock.step6"))

if dock.get("result"):
    scores = dock["result"]["scores"]
    if scores.empty:
        st.warning(t("dock.no_modes"))
    else:
        pose_sel = st.selectbox(t("dock.sel_pose"),
                                scores["mode"].tolist(), index=0, key="pose_sel")
        conclusion = st.text_area(t("dock.conclusion_label"),
                                  key="conclusion", height=90,
                                  placeholder=t("dock.conclusion_ph",
                                                gene=dock.get("gene", t("dock.default_target")),
                                                aff=f"{scores.iloc[0]['affinity']:.2f}"))
        drug = st.text_input(t("dock.drug_prefix_label"), value=default_drug, key="drug_name")
        disease = st.text_input(t("dock.disease_prefix_label"), value=default_disease, key="disease_name")

        col_save, col_star = st.columns([3, 1])
        with col_save:
            if st.button(t("dock.save_btn"), key="save_result"):
                gene = dock.get("gene", "gene")
                prefix = f"{drug}_{disease}_{gene}".strip("_").replace(" ", "_")
                try:
                    result_dir = _save_docking_result(
                        prefix, gene, dock.get("pdb_id", ""), dock.get("cid"),
                        dock["box"], dock["result"]["out"], scores, conclusion,
                        receptor_source=dock.get("receptor_source", "experiment"))
                    dock["result_dir"] = result_dir
                    dock.pop("starred", None)  # 新保存的结果视为未收藏
                    st.success(t("dock.saved"))
                except Exception as exc:
                    st.error(t("dock.save_failed", exc=exc))
        with col_star:
            # 收藏开关:收藏的副本会被保留、不参与自动清理;再点一次即取消收藏
            if dock.get("result_dir"):
                if dock.get("starred"):
                    if st.button(t("dock.starred_btn"), key="star_on"):
                        try:
                            star_base = os.path.basename(dock["result_dir"])
                            import glob
                            for p in glob.glob(os.path.join(config.MD_STAR_DIR, star_base) + "*"):
                                if os.path.isdir(p):
                                    shutil.rmtree(p, ignore_errors=True)
                            dock.pop("starred", None)
                            st.rerun()  # 立即更新图标为 ☆
                        except Exception as exc:
                            st.error(t("dock.unstar_failed", exc=exc))
                else:
                    if st.button(t("dock.star_btn"), key="star_off"):
                        try:
                            Utils.star_docking_result(dock["result_dir"])
                            dock["starred"] = True
                            st.rerun()  # 立即更新图标为 ★
                        except Exception as exc:
                            st.error(t("dock.star_failed", exc=exc))

        if dock.get("result_dir"):
            st.write(t("dock.result_dir", dir=dock["result_dir"]))
            st.download_button(t("dock.dl_all_zip"),
                               data=Utils.make_dir_zip(dock["result_dir"]),
                               file_name=os.path.basename(dock["result_dir"]) + ".zip",
                               mime="application/zip", key="dl_zip")

        # 3D 预览所选构象与受体的结合情况
        st.markdown(t("dock.preview3d"))
        if dock.get("pdb_path"):
            try:
                pose_path = os.path.join(config.DOCKING_RAW_CACHE_DIR,
                                         f"pose_{int(pose_sel)}.pdbqt")
                extract_model(dock["result"]["out"], int(pose_sel), pose_path)
                _render_3d(dock["pdb_path"], pose_path)
                _save_3d_html_ui(dock, dock["pdb_path"], pose_path, pose_sel)
            except Exception as exc:
                st.warning(t("dock.preview3d_failed", exc=exc))
else:
    st.info(t("dock.run_to_interpret"))

st.divider()
st.caption(t("dock.cache_policy"))
