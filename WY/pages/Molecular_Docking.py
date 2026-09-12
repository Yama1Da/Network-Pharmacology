import logging
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
    cofactor_hetatm_pdb,
    download_pdb,
    extract_model,
    fetch_alphafold_receptor,
    interactions_near_ligand,
    interactive_3d_html,
    parse_hetatm,
    pose_pdbqt_to_pdb,
    residues_in_box,
    residues_near_ligand,
    resolve_cid,
    save_interactive_3d_html,
    smiles_to_sdf_3d,
)
from core.enrichment import parse_enrichment_filename

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

# 每次进入本页时清理一次过期缓存(页内操作不重复清理);
# 顺手把旧版完整哈希命名的富集缓存改成短哈希命名,让下方文件选择列表命名统一可读
if not st.session_state.get("_docking_cache_cleaned"):
    Utils.clean_docking_cache()
    Utils.migrate_legacy_cache_names()
    st.session_state["_docking_cache_cleaned"] = True


# ===== 辅助函数 =====

def _extract_label(file_name: str) -> str:
    """从富集结果文件名中解析出"药物_疾病"标识,用于命名后续结果。"""
    info = parse_enrichment_filename(file_name)
    return info["label"] if info else ""


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
    """切换候选基因数据源后,清空配体、受体、对接盒和对接结果,以及由这些数据派生的
    会话控件键(基因/结构/构象下拉、盒子微调数字输入、药物/疾病名输入与结论):
    残留的控件会话值不在新选项里会让 selectbox 抛异常,旧盒会静默覆盖新推断盒。"""
    for _k in ("lig_pdbqt", "cid", "lig_auto_done",
               "gene", "rec_pdbqt", "pdb_id", "pdb_path",
               "box", "result", "result_dir", "starred",
               "gene_structures", "candidate_pre", "batch_results"):
        dock.pop(_k, None)
    for _k in ("gene_sel", "pdb_sel", "pose_sel", "conclusion",
               "drug_name", "disease_name",
               "cx_x", "cx_y", "cx_z", "sx_x", "sx_y", "sx_z"):
        st.session_state.pop(_k, None)


def _fetch_ligand_into_dock(query: str, not_found_key: str, no_3d_key: str,
                            auto_mode: bool = False) -> None:
    """Step 2 的 PubChem 配体获取公共链路:药物名/CID → 3D SDF → pdbqt,成功后写入 dock。
    查无化合物/无 3D 结构两种失败提示的文案由 i18n 键区分;auto_mode 为真时失败也置
    lig_auto_done,本次会话不再自动重试。"""
    cid = resolve_cid(query)
    if cid is None:
        if auto_mode:
            dock["lig_auto_done"] = True
        st.error(t(not_found_key, drug=query, q=query))
        return
    sdf_path = Utils.fetch_ligand_sdf_3d(cid)
    if not sdf_path:
        if auto_mode:
            dock["lig_auto_done"] = True
        st.error(t(no_3d_key, cid=cid))
        return
    dock["cid"] = cid
    dock["lig_pdbqt"] = Utils.prepare_ligand_pdbqt(sdf_path)


def _read_ppi_dataset(dir_name: str) -> tuple[set[str], set[str]]:
    """读取所选 PPI 结果目录,返回 HUB 基因集合与富集验证基因集合;
    对应文件缺失时返回空集合。历史版本产出的 nodes 表可能缺 is_hub/gene 列
    (downloads 目录跨版本留存),缺列或读不了时按空集合处理,不让页面崩溃。"""
    base = os.path.join(config.DOWNLOADS_DIR, dir_name)
    prefix = dir_name[:-5] if dir_name.endswith("_data") else dir_name
    hub: set[str] = set()
    nodes_path = os.path.join(base, f"{prefix}_nodes.csv")
    if os.path.exists(nodes_path):
        try:
            nodes = pd.read_csv(nodes_path, encoding="utf-8-sig", comment="#")
            if {"is_hub", "gene"}.issubset(nodes.columns):
                hub = set(nodes.loc[nodes["is_hub"].fillna(False), "gene"]
                          .dropna().astype(str).str.upper())
            else:
                logging.warning("nodes file lacks is_hub/gene columns, hub list "
                                "skipped: %s", nodes_path)
        except Exception:
            logging.exception("Failed to read nodes file: %s", nodes_path)
    enr: set[str] = set()
    enr_path = os.path.join(base, f"{prefix}_enrichment.csv")
    if os.path.exists(enr_path):
        try:
            enr = _genes_from_df(pd.read_csv(enr_path, encoding="utf-8-sig", comment="#"))
        except Exception:
            logging.exception("Failed to read enrichment file: %s", enr_path)
    return hub, enr


def _docking_result_prefix(drug: str, disease: str, gene: str) -> str:
    """归档结果目录名前缀:三段各自过 safe_name 消毒,防止用户输入中的路径分隔符等
    特殊字符进入文件路径;空段必须先剔除再消毒(safe_name 对空串返回 "unnamed",
    直接套会凭空出现 unnamed 段),全空回退 "docking"。safe_name 只用于目录命名,
    meta 等人读字段仍用原始输入。"""
    parts = []
    for seg in (drug, disease, gene):
        seg = (seg or "").strip()
        if seg:
            parts.append(Utils.safe_name(seg))
    return "_".join(parts) or "docking"


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


# 相互作用虚线圆柱配色与 3Dmol 视图/导出 HTML 保持一致:蓝=氢键 紫=盐桥 灰=疏水
_INTERACTION_COLORS = {"hbond": "#0f62fe", "salt": "#c026d3", "hydro": "#808080"}


def _render_3d(pdb_path: str, pose_pdbqt: str, pocket: list[dict] | None = None,
               interactions: list[dict] | None = None, focus_pocket: bool = False,
               pose_pdb: str | None = None, receptor_text: str | None = None) -> None:
    """在页面内 3D 展示受体(卡通)、所选对接构象(绿棍)与结合位点残基
    (橙棍并标注"残基名-编号"),相互作用以虚线画出(蓝=氢键 紫=盐桥 灰=疏水,
    氢键与盐桥在中点标注距离),并叠加共晶配体(灰色半透明)作天然口袋参照;
    focus_pocket 为真时以配体为中心聚焦口袋并把周围蛋白卡通虚化(半透明),
    否则展示整体视图;缺少组件时仅给出提示。
    pose_pdb 传入已转换好的构象 PDB、receptor_text 传入已读取的受体 PDB 全文时直接复用,
    与「保存交互式 3D 视图」共用一次转换、一次读盘,不传时保持原有的自足行为。"""
    try:
        import py3Dmol
    except ImportError:
        st.info(t("dock.py3dmol_pose"))
        return
    if pose_pdb is None:
        pose_pdb = os.path.join(os.path.dirname(pose_pdbqt), "pose_view.pdb")
        try:
            pose_pdbqt_to_pdb(pose_pdbqt, pose_pdb)
        except Exception:
            st.warning(t("dock.pose_pdb_failed"))
            return
    if receptor_text is None:
        with open(pdb_path, encoding="utf-8", errors="replace") as f:
            receptor_text = f.read()
    with open(pose_pdb, encoding="utf-8") as f:
        ligand_pdb = f.read()
    # 受体中的共晶配体(非水 HETATM),叠加作天然口袋参照
    het_pdb = cofactor_hetatm_pdb(pdb_path, receptor_text=receptor_text)

    view = py3Dmol.view(width=900, height=520)
    view.addModel(receptor_text, "pdb")
    view.addModel(ligand_pdb, "pdb")
    # 聚焦口袋时卡通整体半透明(虚化周围蛋白),突出口袋棍状、相互作用虚线与距离标注;
    # 虚化只是显示开关,信息不删减,切回开关即还原
    cartoon: dict = {"color": "spectrum"}
    if focus_pocket:
        cartoon["opacity"] = 0.3
    view.setStyle({"model": 0}, {"cartoon": cartoon})
    view.setStyle({"model": 1}, {"stick": {"colorscheme": "greenCarbon"}})
    if het_pdb:
        view.addModel(het_pdb + "\n", "pdb")  # model 2:共晶配体参照
        view.setStyle({"model": 2}, {"stick": {"colorscheme": "greyCarbon",
                                               "opacity": 0.5, "radius": 0.15}})
    # 结合位点残基:橙色棍状高亮,按链分组选择;无链号时不限定 chain。
    # 卡通与棍状在一条语句里显式声明完整样式,不依赖多次 setStyle 的合并行为
    by_chain: dict[str, list] = {}
    for r in pocket or []:
        by_chain.setdefault(r["chain"], []).append(r["resi"])
    for chain, resis in by_chain.items():
        sel: dict = {"model": 0, "resi": resis}
        if chain:
            sel["chain"] = chain
        view.setStyle(sel, {"cartoon": cartoon,
                            "stick": {"colorscheme": "orangeCarbon", "radius": 0.2}})
    # 残基名标签,位置取该残基 CA 原子(缺 CA 时用第一个原子)
    for r in pocket or []:
        x, y, z = r["pos"]
        view.addLabel(f"{r['resn']}-{r['resi']}", {
            "position": {"x": x, "y": y, "z": z},
            "fontSize": 11, "fontColor": "black",
            "showBackground": True, "backgroundColor": "white",
            "backgroundOpacity": 0.7,
        })
    # 相互作用虚线圆柱(1px 线白底下看不清,故用带半径的虚线圆柱);
    # 氢键与盐桥在中点标注距离,疏水接触只画线不标注
    for it in interactions or []:
        view.addCylinder({"dashed": True, "radius": 0.06, "dashLength": 0.3, "gapLength": 0.22,
                          "color": _INTERACTION_COLORS.get(it["kind"], "#777777"),
                          "start": {"x": it["pos1"][0], "y": it["pos1"][1], "z": it["pos1"][2]},
                          "end": {"x": it["pos2"][0], "y": it["pos2"][1], "z": it["pos2"][2]}})
        if it["kind"] != "hydro":
            mid = [(a + b) / 2 for a, b in zip(it["pos1"], it["pos2"])]
            view.addLabel(f'{it["dist"]:.1f} Å', {
                "position": {"x": mid[0], "y": mid[1], "z": mid[2]},
                "fontSize": 10, "fontColor": "#333333",
                "showBackground": True, "backgroundColor": "white",
                "backgroundOpacity": 0.6,
            })
    if focus_pocket:
        # 以配体为中心取景,再拉远把口袋残基框进来;整体视图仍可拖动/滚轮还原
        view.zoomTo({"model": 1})
        view.zoom(0.55)
    else:
        view.zoomTo()
    try:
        html = view._make_html()
    except AttributeError:
        html = view.write()
    legend = t("dock.legend_pose" if het_pdb else "dock.legend_pose_apo")
    if interactions:
        legend += "  |  " + t("dock.legend_inter")
    st.caption(legend)
    st.components.v1.html(html, height=540)


def _interactions_table(interactions: list[dict]) -> pd.DataFrame:
    """把相互作用原子对列表转成展示表格(类型本地化,残基带链号),列顺序固定。"""
    kind_label = {"hbond": t("dock.inter_hbond"), "salt": t("dock.inter_salt"),
                  "hydro": t("dock.inter_hydro")}
    rows = [{
        t("dock.inter_col_kind"): kind_label.get(it["kind"], it["kind"]),
        t("dock.inter_col_residue"): (it["chain"] + "/" if it["chain"] else "")
                                     + f'{it["resn"]}{it["resi"]}',
        t("dock.inter_col_res_atom"): it["res_atom"],
        t("dock.inter_col_lig_atom"): it["lig_atom"],
        t("dock.inter_col_dist"): it["dist"],
    } for it in interactions]
    return pd.DataFrame(rows)


def _save_3d_html_ui(dock: dict, pdb_path: str, pose_pdbqt: str, pose_sel,
                     pocket: list[dict] | None = None,
                     interactions: list[dict] | None = None,
                     pose_pdb: str | None = None, receptor_text: str | None = None) -> None:
    """「保存交互式 3D 视图」:把受体、当前构象、结合位点残基标注、相互作用虚线
    与共晶配体参照打包成单个自包含 HTML(内嵌 3Dmol.js 与坐标),可下载或存入结果目录;
    之后双击即可在任意浏览器打开,可拖动旋转、聚焦口袋,离线可用。
    pose_pdb / receptor_text 传入 3D 预览已转换/已读取的结果时直接复用,不再重复转换与读盘。"""
    gene = Utils.safe_name(str(dock.get("gene") or "target"))
    pose_n = int(pose_sel)
    # 图例随视图一并写入,离线打开时颜色含义不靠猜;共晶配体判定与 HTML 生成共用同一份受体文本
    has_het = bool(cofactor_hetatm_pdb(pdb_path, receptor_text=receptor_text))
    caption = (t("dock.save3d_caption", gene=dock.get("gene") or "?", pose=pose_n,
                 drug=dock.get("drug") or "-", disease=dock.get("disease") or "-")
               + "  |  " + t("dock.legend_pose" if has_het else "dock.legend_pose_apo")
               + (("  |  " + t("dock.legend_inter")) if interactions else ""))
    focus_label = t("dock.focus_btn")
    try:
        html_text, offline_ok = interactive_3d_html(pdb_path, pose_pdbqt, caption, pocket=pocket,
                                                    receptor_text=receptor_text, pose_pdb=pose_pdb,
                                                    interactions=interactions,
                                                    focus_label=focus_label)
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
                        os.path.join(dock["result_dir"], file_name), caption, pocket=pocket,
                        interactions=interactions, focus_label=focus_label)
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
    查无实验结构的基因按开关附加 AlphaFold 预测结构。"""
    return Utils.fetch_receptor_candidates_batch(list(genes), max_hits=20,
                                                 allow_alphafold=allow_alphafold)


# ===== Step 1 候选基因筛选 =====
def _step1_select_source(dock: dict) -> None:
    """数据来源选择:PPI 分析结果目录 或 富集缓存文件,选中后自动读取候选基因集。"""
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
        # 只列出 PPI 分析结果目录,排除其他类型的文件夹;目录名展示为"药物 → 疾病"
        dirs = Utils.list_download_dirs()
        ppi_dirs = dirs[dirs["Directory"].str.endswith("_data", na=False)] if not dirs.empty else dirs
        if ppi_dirs.empty:
            st.warning(t("dock.no_ppi_dirs"))
        else:
            dir_disp = {d: Utils.humanize_analysis_label(d[:-5] if d.endswith("_data") else d)
                        for d in ppi_dirs["Directory"]}
            sel_dir = st.selectbox(t("dock.sel_ppi_dir"), ppi_dirs["Directory"].tolist(), key="sel_ppi_dir",
                                   index=None, placeholder=t("dock.sel_ppi_dir_ph"),
                                   format_func=lambda d: dir_disp.get(d, d))
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
        enr_files = Utils.list_cache_files("enrichment", exts=(".csv",))
        # 只保留普通富集(ORA)结果供选择:GSEA 文件没有逐条目的基因列,选中也提不出候选基因;
        # 按最后修改时间倒序,最近一次分析排最前;文件名展示为"药物 → 疾病 · 数据源"的友好名称
        if not enr_files.empty:
            is_ora = enr_files["File name"].map(
                lambda n: (parse_enrichment_filename(n) or {}).get("kind") == "ora")
            enr_files = (enr_files[is_ora]
                         .sort_values("Modified time", ascending=False))
        if enr_files.empty:
            st.warning(t("dock.no_enr_files"))
        else:
            # 同名分析可能有多份不同参数的缓存,展示名带上时间戳以便区分
            file_disp = {row["File name"]:
                         f"{Utils.enrichment_display_name(row['File name'])} · {row['Modified time']}"
                         for _, row in enr_files.iterrows()}
            sel_file = st.selectbox(t("dock.sel_enr_file"), enr_files["File name"].tolist(), key="sel_enr",
                                    index=None, placeholder=t("dock.sel_enr_file_ph"),
                                    format_func=lambda n: file_disp.get(n, n))
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


def _step1_filter_genes(dock: dict) -> None:
    """候选基因确认:HUB∩富集交集展示、RCSB 结构可用性批量筛选(可选 AlphaFold 兜底)。"""
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
                    # 被剔除的基因里有 hub 时要特别提醒:hub 是对接的核心对象,
                    # 无实验结构只能靠 AlphaFold 预测模型兜底,不提醒容易被忽略
                    hub_no_struct = sorted(hub_set & set(no_struct))
                    if hub_no_struct and not allow_af:
                        st.warning(t("dock.hub_no_struct_warn", n=len(hub_no_struct),
                                     genes=", ".join(hub_no_struct)))
            else:
                dock.pop("genes", None)
                st.info(t("dock.click_to_check"))
        else:
            # 不做结构筛选时,候选基因直接作为最终靶点列表
            dock["genes"] = candidate
            if candidate:
                st.caption(t("dock.final_candidates", n=len(candidate))
                           + ", ".join(candidate[:60]) + (" …" if len(candidate) > 60 else ""))


def _step1_manual_genes(dock: dict) -> None:
    """没有可选数据文件时的兜底入口:手动粘贴基因列表。"""
    manual = st.text_area(t("dock.manual_genes_label"), key="manual_genes")
    if st.button(t("dock.use_manual_btn"), key="load_genes_manual"):
        dock["genes"] = sorted({g.strip().upper() for g in re.split(r"[,;\s]+", manual) if g.strip()})
        dock["label"] = ""
        dock["drug"] = ""
        dock["disease"] = ""
        dock["hub_genes"], dock["enrichment_genes"] = [], []
        _reset_docking_downstream()  # 数据源变了,清空后面步骤的旧结果


st.subheader(t("dock.step1"))
_step1_select_source(dock)
_step1_filter_genes(dock)

# 没有可选数据文件时,允许手动输入基因列表
if dock.get("enrichment_genes") is None:
    _step1_manual_genes(dock)

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
            _fetch_ligand_into_dock(default_drug, "dock.pubchem_not_found",
                                    "dock.cid_no_3d", auto_mode=True)

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
                _fetch_ligand_into_dock(lig_query, "dock.pubchem_no_compound",
                                        "dock.cid_try_another")
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
    # dock["genes"] 可能被重新筛选(重跑结构批量检查等),残留会话值不在新列表里
    # 会让 selectbox 抛异常,先按新列表校验
    if st.session_state.get("gene_sel") not in dock["genes"]:
        st.session_state.pop("gene_sel", None)
    gene = st.selectbox(t("dock.gene_label"), dock["genes"], key="gene_sel")
    # 换靶点基因后,清空旧受体和对接结果,以及派生的控件会话键
    if dock.get("gene") != gene:
        dock["gene"] = gene
        for _k in ("rec_pdbqt", "pdb_path", "pdb_id", "box", "result",
                   "result_dir", "starred", "receptor_source"):
            dock.pop(_k, None)
        for _k in ("pose_sel", "conclusion",
                   "cx_x", "cx_y", "cx_z", "sx_x", "sx_y", "sx_z"):
            st.session_state.pop(_k, None)

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
        # 换基因后新候选列表不含残留的旧 PDB id,先校验再渲染,避免 selectbox 抛异常
        if st.session_state.get("pdb_sel") not in [str(x) for x in cand["pdb_id"].tolist()]:
            st.session_state.pop("pdb_sel", None)
        pdb_id = st.selectbox(t("dock.sel_pdb"), cand["pdb_id"].tolist(), key="pdb_sel")
        # 更换 PDB 结构后,作废旧受体和后续步骤的旧结果,以及派生的控件会话键
        if dock.get("pdb_sel") != pdb_id:
            dock["pdb_sel"] = pdb_id
            for _k in ("rec_pdbqt", "pdb_id", "pdb_path", "box", "result",
                       "result_dir", "starred", "receptor_source"):
                dock.pop(_k, None)
            for _k in ("pose_sel", "conclusion",
                       "cx_x", "cx_y", "cx_z", "sx_x", "sx_y", "sx_z"):
                st.session_state.pop(_k, None)
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
        text = f.read()
    atom_lines = [ln for ln in text.splitlines() if ln.startswith(("ATOM", "TER", "END"))]
    # 共晶配体(非水 HETATM)叠加作参照,过滤规则与 Step 6 的 3D 预览一致
    het_pdb = cofactor_hetatm_pdb(pdb_path, receptor_text=text)

    c, s = box["center"], box["size"]
    view = py3Dmol.view(width=900, height=520)
    view.addModel("\n".join(atom_lines) + "\n", "pdb")  # model 0：蛋白
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
    if het_pdb:
        view.addModel(het_pdb + "\n", "pdb")  # model 1：共结晶配体
        # 灰色=共晶配体(参考),与 Step 6 的叠加参照同色;绿色留给对接结果
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "greyCarbon"}})
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
    st.caption(t("dock.legend_box" if het_pdb else "dock.legend_box_apo"))
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
                # 新一轮对接:清掉上一轮残留的构象选择、结论与归档状态
                st.session_state.pop("pose_sel", None)
                st.session_state.pop("conclusion", None)
                dock.pop("result_dir", None)
                dock.pop("starred", None)
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
        # 新一轮对接的 mode 数可能少于旧会话残留的选中值,先按新打分表校验再渲染
        if st.session_state.get("pose_sel") not in scores["mode"].tolist():
            st.session_state.pop("pose_sel", None)
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
                prefix = _docking_result_prefix(drug, disease, gene)
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
                # 半径以主页齿轮的全局参数(config.POCKET_CUTOFF)为主导:齿轮改动后
                # 本页滑条自动跟随;页内 −/+/滑条仅本次会话生效,不回写配置
                # (回写会让对接页的旧值覆盖齿轮设置)。按钮必须写在滑条实例化之前
                # (控件键在其控件创建后不可再写),点击当轮即生效;步进 0.5 Å,边界夹紧
                if st.session_state.get("_pocket_cutoff_sync") != float(config.POCKET_CUTOFF):
                    st.session_state["_pocket_cutoff_sync"] = float(config.POCKET_CUTOFF)
                    st.session_state["pocket_cutoff"] = float(config.POCKET_CUTOFF)
                cbtn, cslider = st.columns([1.2, 6])
                with cbtn:
                    bminus, bplus = st.columns(2)
                    with bminus:
                        minus = st.button("−", key="cutoff_minus", use_container_width=True)
                    with bplus:
                        plus = st.button("+", key="cutoff_plus", use_container_width=True)
                if minus:
                    st.session_state.pocket_cutoff = max(
                        3.0, st.session_state.pocket_cutoff - 0.5)
                if plus:
                    st.session_state.pocket_cutoff = min(
                        6.0, st.session_state.pocket_cutoff + 0.5)
                with cslider:
                    cutoff = st.slider(t("dock.cutoff_label"), 3.0, 6.0,
                                       step=0.5, key="pocket_cutoff")
                pocket = residues_near_ligand(dock["pdb_path"], pose_path, cutoff=cutoff)
                interactions = interactions_near_ligand(dock["pdb_path"], pose_path)
                # 构象 pdbqt→pdb 只转换一次、受体 PDB 全文只读一次,
                # 3D 预览与保存 HTML 复用同一份结果(每次刷新少一次转换、少两次读盘)
                pose_pdb = os.path.join(os.path.dirname(pose_path), "pose_view.pdb")
                try:
                    pose_pdbqt_to_pdb(pose_path, pose_pdb)
                except Exception:
                    pose_pdb = None  # 转换失败:预览按原样提示,保存时由 core 内部重试
                with open(dock["pdb_path"], encoding="utf-8", errors="replace") as f:
                    receptor_text = f.read()
                focus_pocket = st.toggle(t("dock.focus_toggle"), key="focus_pocket")
                _render_3d(dock["pdb_path"], pose_path, pocket,
                           interactions=interactions, focus_pocket=focus_pocket,
                           pose_pdb=pose_pdb, receptor_text=receptor_text)
                if pocket:
                    names = [f"{r['resn']}{r['resi']}" for r in pocket]
                    st.caption(t("dock.pocket_caption", n=len(names), cutoff=f"{cutoff:g}")
                               + ", ".join(names[:40])
                               + (" …" if len(names) > 40 else ""))
                if interactions:
                    with st.expander(t("dock.inter_expander", n=len(interactions))):
                        st.dataframe(_interactions_table(interactions),
                                     use_container_width=True, hide_index=True)
                else:
                    st.caption(t("dock.inter_none"))
                _save_3d_html_ui(dock, dock["pdb_path"], pose_path, pose_sel, pocket,
                                 interactions=interactions,
                                 pose_pdb=pose_pdb, receptor_text=receptor_text)
            except Exception as exc:
                st.warning(t("dock.preview3d_failed", exc=exc))
else:
    st.info(t("dock.run_to_interpret"))

st.divider()
st.caption(t("dock.cache_policy"))
