"""分子对接(阶段 3b:Molecular_Docking 页)。"""

from fastapi import APIRouter
import os
import re
from datetime import datetime
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel
from WY import config
from WY.core import jobs, transfer
from WY.core.enrichment import parse_enrichment_filename
from WY.core.files import (list_cache_files, list_download_dirs, list_log_files,
                           make_bytes_zip, make_dirs_zip, make_ppi_zip)
from WY.core.management import (enrichment_display_name, save_api_env,
                                save_ui_settings, update_all_cache_files)
from WY.core.text_utils import safe_name
from functools import lru_cache
from WY.core import docking_session


router = APIRouter()


# 多步向导的无状态 API:会话状态存 core.docking_session(镜像原 dock 会话键),
# 每个端点对应一步;网络重的步骤(结构批量检查/受体准备/vina)经 job 提交-轮询。
# 3D 预览走产物直出:py3Dmol 在服务端渲染为 HTML,前端 iframe 内嵌。


def _embed_centered(html: str) -> str:
    """py3Dmol 生成的 HTML 注入 iframe 适配样式:去 body 边距、禁滚动条、
    视图在框内水平垂直居中(用户拖转/滚轮缩放不受影响)。"""
    style = ("<style>html,body{margin:0;padding:0;overflow:hidden;height:100%}"
             "body{display:flex;align-items:center;justify-content:center;"
             "background:#ffffff}</style>")
    if "<head>" in html:
        return html.replace("<head>", "<head>" + style, 1)
    return style + html


@lru_cache(maxsize=256)
def _pdb_candidates_cached(gene: str, allow_af: bool):
    """单基因受体候选查询的进程内缓存(原页面 @st.cache_data(ttl=1h) 的等价物)。"""
    from WY.core.docking import fetch_receptor_pdb_candidates

    return fetch_receptor_pdb_candidates(gene, max_hits=20, allow_alphafold=allow_af)


def _session(sid: str) -> dict:
    dock = docking_session.get(sid)
    if dock is None:
        raise HTTPException(status_code=404, detail="docking session expired, create a new one")
    return dock


@router.post("/api/docking/session")
def docking_create_session() -> dict:
    return {"sid": docking_session.create()}


@router.get("/api/docking/sources")
def docking_sources() -> dict:
    """Step 1 的候选来源清单:PPI 分析目录(展示名 = 药物 → 疾病)与 ORA 富集文件。

    手动基因分析(manual_* 桶及其富集缓存)不进这两个清单 —— 它们由
    「手动基因对接」专属入口接管(见 manual_dirs 与 source 的 manual 模式)。
    """
    from WY.core.docking import split_result_prefix
    from WY.core.enrichment import parse_enrichment_filename
    from WY.core.files import list_cache_files, list_download_dirs
    from WY.core.management import enrichment_display_name

    dirs = list_download_dirs()
    ppi_dirs = dirs[dirs["Directory"].str.endswith("_data", na=False)] if not dirs.empty else dirs
    ppi_out = []
    manual_out = []
    for d in ppi_dirs["Directory"].tolist():
        prefix = d[:-5] if d.endswith("_data") else d
        # 空壳目录(中断运行遗留、无任何文件)不进选择清单,避免出现 "empty →" 这类无效项
        try:
            if not any(e.is_file() for e in os.scandir(os.path.join(config.DOWNLOADS_DIR, d))):
                continue
        except OSError:
            continue
        if prefix.startswith("manual"):
            # 手动基因桶:进专属入口的清单(基因列表在 analysis_meta.json),
            # 不混入 PPI 分析下拉 —— 哈希桶名对用户无意义,且药病反解必得 "manual"
            meta = _read_manual_meta(d)
            if meta["genes"]:
                manual_out.append(meta)
            continue
        drug, disease = split_result_prefix(prefix)
        ppi_out.append({"name": d, "display": f"{drug} → {disease}" or d})

    enr = list_cache_files("enrichment", exts=(".csv",))
    if not enr.empty:
        is_ora = enr["File name"].map(
            lambda n: (parse_enrichment_filename(n) or {}).get("kind") == "ora")
        enr = enr[is_ora].sort_values("Modified time", ascending=False)
        # 手动基因模式的富集缓存(label=manual_<hash>)同样归专属入口,不进本清单
        enr = enr[~enr["File name"].map(
            lambda n: ((parse_enrichment_filename(n) or {}).get("label") or "")
            .startswith("manual"))]
    enr_out = [{"name": row["File name"],
                "display": f"{enrichment_display_name(row['File name'])} · {row['Modified time']}"}
               for _, row in enr.iterrows()]
    return {"ppi_dirs": ppi_out, "enr_files": enr_out, "manual_dirs": manual_out}


def _read_manual_meta(dir_name: str) -> dict:
    """读手动基因桶的 analysis_meta.json,返回 {name, n_genes, mtime, genes}。

    基因清单由 save_ppi_results 写入(core/ppi.py,manual 前缀桶专属);
    缺文件/无 genes(旧版桶)返回空 genes,调用方据此跳过。
    """
    import json as _json

    d = os.path.join(config.DOWNLOADS_DIR, dir_name)
    info = {"name": dir_name, "n_genes": 0, "mtime": "", "genes": []}
    try:
        info["mtime"] = datetime.fromtimestamp(
            os.path.getmtime(d)).strftime("%Y-%m-%d %H:%M")
    except OSError:
        pass
    meta_path = os.path.join(d, "analysis_meta.json")
    try:
        with open(meta_path, encoding="utf-8") as f:
            genes = (_json.load(f).get("genes") or [])
        genes = sorted({str(g).strip().upper() for g in genes if str(g).strip()})
        info["genes"] = genes
        info["n_genes"] = len(genes)
    except (OSError, ValueError):
        pass
    return info


class DockSourceBody(BaseModel):
    mode: str            # "ppi" | "enr" | "manual" | "manual-file"
    name: str = ""       # ppi 目录名 / enr 文件名
    manual_genes: str = ""
    manual_file: str = ""  # manual-file 模式:downloads 下 manual_*_data 目录名


@router.post("/api/docking/{sid}/source")
def docking_source(sid: str, body: DockSourceBody) -> dict:
    dock = _session(sid)
    import pandas as pd

    from WY.core.docking import split_result_prefix
    from WY.core.enrichment import enrichment_filename_label
    from WY.core.files import list_cache_files
    from WY.core.ppi import read_ppi_dataset
    from WY.core.text_utils import genes_from_dataframe

    hub_genes: list[str] = []
    if body.mode == "ppi":
        dock["ppi_dir"] = body.name
        dock["hub_genes"], dock["enrichment_genes"] = read_ppi_dataset(body.name)
        prefix = body.name[:-5] if body.name.endswith("_data") else body.name
        dock["label"] = prefix
        dock["drug"], dock["disease"] = split_result_prefix(prefix)
    elif body.mode == "enr":
        dock["enr_file"] = body.name
        df = pd.read_csv(os.path.join(config.ENRICHMENT_CACHE_DIR, body.name),
                         encoding="utf-8-sig", comment="#")
        dock["enrichment_genes"] = sorted(genes_from_dataframe(df))
        dock["hub_genes"] = []
        label = enrichment_filename_label(body.name)
        dock["label"] = label
        dock["drug"], dock["disease"] = split_result_prefix(label)
    elif body.mode in ("manual", "manual-file"):
        # 手动基因对接专属入口:基因列表来自输入框(manual)或历史桶元数据
        # (manual-file,读 analysis_meta.json)。无药病语境,label 只作来源标识,
        # 以 manual 开头 —— ligand 端点据此拦截"按药名查 PubChem"。
        if body.mode == "manual":
            genes = sorted({g.strip().upper() for g in re.split(r"[,;\s]+", body.manual_genes)
                            if g.strip()})
            if not genes:
                raise HTTPException(status_code=400, detail="manual_genes is required")
            dock["label"] = ""
        else:
            safe = os.path.basename(body.manual_file.strip())
            if not (safe.startswith("manual") and safe.endswith("_data")):
                raise HTTPException(status_code=400, detail="manual_file must be a manual_*_data dir")
            genes = _read_manual_meta(safe)["genes"]
            if not genes:
                raise HTTPException(status_code=400,
                                    detail=f"{safe}: analysis_meta.json 无基因清单(旧版桶)")
            dock["label"] = safe[:-5] if safe.endswith("_data") else safe
        dock["drug"] = ""
        dock["disease"] = ""
        dock["hub_genes"], dock["enrichment_genes"] = [], []
        dock["genes"] = genes
    else:
        raise HTTPException(status_code=400, detail="mode must be ppi|enr|manual|manual-file")

    docking_session.reset_downstream(dock)
    if body.mode not in ("ppi", "enr"):
        # 手动模式:全部输入基因即候选(无网络/无交集计算),跳过结构检查直接
        # 进 Step 3 也能选基因;check-structures 会把 genes 替换为"有结构"子集
        dock["candidate_pre"] = list(dock["genes"])
    else:
        hub_set = set(dock.get("hub_genes") or [])
        enr_set = set(dock.get("enrichment_genes") or [])
        dock["candidate_pre"] = sorted(hub_set & enr_set) if hub_set else sorted(enr_set)
        # genes 先按候选全集初始化:跳过结构检查直接进 Step 3 也能选基因;
        # 之后 check-structures 会把它替换为"有结构"的子集
        dock["genes"] = list(dock["candidate_pre"])
    return {
        "label": dock.get("label", ""), "drug": dock.get("drug", ""),
        "disease": dock.get("disease", ""),
        "hub_count": len(dock.get("hub_genes") or []),
        "enr_count": len(dock.get("enrichment_genes") or []),
        "candidates": dock.get("candidate_pre") or dock.get("genes") or [],
    }


class DockStructBody(BaseModel):
    use_struct: bool = True
    allow_af: bool = False


@router.post("/api/docking/{sid}/check-structures")
def docking_check_structures(sid: str, body: DockStructBody) -> dict:
    """Step 1 结构可用性批量筛选(RCSB 4 路并发 + GraphQL 批量归因):经 job 提交。"""
    dock = _session(sid)
    candidate = dock.get("candidate_pre") or []
    if not candidate:
        raise HTTPException(status_code=400, detail="no candidates: select a source first")
    if body.use_struct:
        from WY.core.docking import fetch_receptor_candidates_batch

        def task(progress) -> dict:
            progress(fraction=0.15, text="RCSB 批量检索中…")
            results = fetch_receptor_candidates_batch(list(candidate), max_hits=20,
                                                      allow_alphafold=body.allow_af)
            gene_structures = {g: len(df) for g, df in results.items() if not df.empty}
            dock["batch_results"] = results
            dock["gene_structures"] = gene_structures
            dock["genes"] = sorted(gene_structures)
            hub_set = set(dock.get("hub_genes") or [])
            no_struct = [g for g in candidate if g not in gene_structures]
            return {"genes": dock["genes"], "no_struct": no_struct,
                    "hub_no_struct": sorted(hub_set & set(no_struct)),
                    "structures": {g: gene_structures[g] for g in dock["genes"]}}
        return {"job_id": jobs.submit(task, "docking structure check")}
    dock["genes"] = list(candidate)
    return {"job_id": None, "direct": {"genes": dock["genes"], "no_struct": [],
                                       "hub_no_struct": [], "structures": {}}}


class DockLigandBody(BaseModel):
    source: str          # "pubchem" | "smiles"
    query: str = ""
    smiles: str = ""


@router.post("/api/docking/{sid}/ligand")
def docking_ligand(sid: str, body: DockLigandBody) -> dict:
    """Step 2 配体获取:PubChem 名/CID → 3D SDF → pdbqt,或 SMILES → 3D → pdbqt。
    同步端点(线程池执行);失败以 error 字段返回(与原页面分类一致),不抛 500。"""
    dock = _session(sid)
    from datetime import datetime

    from WY.core.docking import (PubChemUnavailable, fetch_ligand_sdf_3d,
                                 prepare_ligand_pdbqt, resolve_cid, smiles_to_sdf_3d)

    if body.source == "pubchem":
        query = body.query.strip()
        if not query:
            return {"ok": False, "error": "未提供药物名或 CID"}
        # 手动基因会话无药名语境:按名称盲查必然是"把来源标识当药名"的误用
        #(曾把 "manual" 直接送进 PubChem 报"查无化合物");CID 直查(纯数字)放行。
        # 判据 = 会话由手动入口建立(label 为空串或 manual_*,且未挂 PPI/富集来源)
        is_manual_session = (not (dock.get("ppi_dir") or dock.get("enr_file"))
                             and not (dock.get("label") or "").startswith(("gProfiler", "Enrichr"))
                             and bool(dock.get("genes")))
        if is_manual_session and not query.isdigit():
            return {"ok": False,
                    "error": "手动基因模式无关联药物:请输入 CID(纯数字)直接下载,或改用 SMILES"}
        try:
            cid = resolve_cid(query)
        except PubChemUnavailable as exc:
            return {"ok": False, "error": f"PubChem 暂不可用({exc}),可稍后重试"}
        if cid is None:
            return {"ok": False, "error": f"PubChem 查无化合物:{query}"}
        try:
            sdf_path = fetch_ligand_sdf_3d(cid)
        except Exception as exc:
            return {"ok": False, "error": f"配体获取失败:{exc}"}
        if not sdf_path:
            return {"ok": False, "error": f"CID {cid} 无 3D 结构,请换一个条目或用 SMILES"}
        dock["cid"] = cid
        try:
            dock["lig_pdbqt"] = prepare_ligand_pdbqt(sdf_path)
        except Exception as exc:
            return {"ok": False, "error": f"配体准备失败:{exc}"}
    elif body.source == "smiles":
        if not body.smiles.strip():
            return {"ok": False, "error": "未提供 SMILES"}
        try:
            os.makedirs(config.DOCKING_RAW_CACHE_DIR, exist_ok=True)
            sdf_path = smiles_to_sdf_3d(
                body.smiles,
                os.path.join(config.DOCKING_RAW_CACHE_DIR,
                             f"ligand_smiles_{datetime.now():%Y%m%d_%H%M%S}.sdf"))
            dock["cid"] = None
            dock["lig_pdbqt"] = prepare_ligand_pdbqt(sdf_path)
        except Exception as exc:
            return {"ok": False, "error": f"配体生成失败:{exc}"}
    else:
        raise HTTPException(status_code=400, detail="source must be pubchem|smiles")
    dock["lig_auto_done"] = True
    return {"ok": True, "cid": dock.get("cid"), "lig_pdbqt": dock.get("lig_pdbqt")}


def _candidate_rows(cand) -> list[dict]:
    """候选结构表 → JSON 行(含核酸防呆标注)。"""
    flag_notes = {"dna_only": "纯核酸结构,不含蛋白,不能对接",
                  "has_nucleic": "含核酸组分,注意确认靶点"}
    rows = []
    for _, r in cand.iterrows():
        flag = str(r.get("flag") or "")
        rows.append({"pdb_id": r["pdb_id"], "resolution": str(r.get("resolution", "")),
                     "method": str(r.get("method", "")), "title": str(r.get("title", "")),
                     "source": str(r.get("source", "")), "flag": flag,
                     "note": flag_notes.get(flag, "")})
    return rows


class DockReceptorQueryBody(BaseModel):
    gene: str
    allow_af: bool = False


@router.post("/api/docking/{sid}/receptor/candidates")
def docking_receptor_candidates(sid: str, body: DockReceptorQueryBody) -> dict:
    """Step 3:该基因的受体候选列表;优先复用 Step 1 批量结果,查不到再单基因查(带缓存)。"""
    dock = _session(sid)
    if body.gene not in (dock.get("genes") or []):
        raise HTTPException(status_code=400, detail=f"gene {body.gene} not in confirmed genes")
    try:
        cand = (dock.get("batch_results") or {}).get(body.gene)
        if cand is None:
            cand = _pdb_candidates_cached(body.gene, bool(body.allow_af and config.ALPHAFOLD_ENABLED))
    except Exception as exc:
        return {"error": f"RCSB 检索失败:{exc}", "rows": []}
    # 换基因:失效下游(与原页面一致)
    if dock.get("gene") != body.gene:
        dock["gene"] = body.gene
        for _k in ("rec_pdbqt", "pdb_path", "pdb_id", "box", "result", "result_dir",
                   "starred", "receptor_source", "pdb_sel", "pose_sel", "conclusion"):
            dock.pop(_k, None)
    return {"rows": _candidate_rows(cand) if cand is not None else [],
            "af_enabled": bool(config.ALPHAFOLD_ENABLED)}


class DockReceptorPrepareBody(BaseModel):
    gene: str
    pdb_id: str


@router.post("/api/docking/{sid}/receptor/prepare")
def docking_receptor_prepare(sid: str, body: DockReceptorPrepareBody) -> dict:
    """下载受体 PDB(或 AlphaFold 兜底)+ 转 pdbqt:网络+本地转换,经 job 提交。"""
    dock = _session(sid)
    if body.gene not in (dock.get("genes") or []):
        raise HTTPException(status_code=400, detail=f"gene {body.gene} not in confirmed genes")
    if body.pdb_id.startswith("AFDB:") and not config.ALPHAFOLD_ENABLED:
        raise HTTPException(status_code=400, detail="AlphaFold fallback is disabled")

    def task(progress) -> dict:
        from WY.core.docking import (download_pdb, fetch_alphafold_receptor,
                                     prepare_receptor_pdbqt)

        progress(fraction=0.2, text="下载受体结构…")
        if body.pdb_id.startswith("AFDB:"):
            pdb_path, af_acc = fetch_alphafold_receptor(
                body.gene, taxon=int(config.SPECIES_MAP["Human"]["ncbi_taxon"]))
            dock["pdb_id"] = f"AFDB:{af_acc}"
            dock["receptor_source"] = "alphafold"
        else:
            pdb_path = download_pdb(body.pdb_id)
            dock["pdb_id"] = body.pdb_id
            dock["receptor_source"] = "experiment"
        progress(fraction=0.7, text="受体转 pdbqt…")
        dock["pdb_path"] = pdb_path
        dock["rec_pdbqt"] = prepare_receptor_pdbqt(pdb_path)
        dock["pdb_sel"] = body.pdb_id
        dock.pop("result", None)
        return {"pdb_id": dock["pdb_id"], "receptor_source": dock["receptor_source"],
                "pdb_path": pdb_path}

    return {"job_id": jobs.submit(task, f"docking receptor {body.gene}/{body.pdb_id}")}


class DockBoxBody(BaseModel):
    infer: bool = False
    center: list[float] | None = None
    size: list[float] | None = None


@router.post("/api/docking/{sid}/box")
def docking_box(sid: str, body: DockBoxBody) -> dict:
    dock = _session(sid)
    from WY.core.docking import clamp_box, infer_box_from_cofactor, parse_hetatm

    if not dock.get("pdb_path"):
        raise HTTPException(status_code=400, detail="receptor not prepared")
    if body.infer:
        box = infer_box_from_cofactor(dock["pdb_path"])
        clamped = False
        if (dock.get("receptor_source") == "alphafold"
                and max(box["size"]) > config.DOCKING_BOX_MAX_SIZE_AF):
            box = clamp_box(box, config.DOCKING_BOX_MAX_SIZE_AF)
            clamped = True
        dock["box"] = box
    else:
        if not body.center or not body.size or len(body.center) != 3 or len(body.size) != 3:
            raise HTTPException(status_code=400, detail="center/size must be [x,y,z]")
        dock["box"] = {"center": tuple(body.center), "size": tuple(body.size)}
    het = parse_hetatm(dock["pdb_path"])
    ligs = sorted({f"{r} (chain {c}{n})" for r, c, n, _ in het}) if het else []
    from WY.core.docking import residues_in_box

    residues = residues_in_box(dock["pdb_path"], dock["box"]) if dock.get("box") else []
    return {"box": {"center": list(dock["box"]["center"]), "size": list(dock["box"]["size"])},
            "clamped": clamped if body.infer else False,
            "inferred_ligands": ligs[:5] if body.infer else None,
            "residues": residues}


@router.get("/api/docking/{sid}/box-3d")
def docking_box_3d(sid: str) -> Response:
    """Step 4 盒子 3D 预览:服务端 py3Dmol 渲染为自包含 HTML,前端 iframe 内嵌。"""
    dock = _session(sid)
    from WY.core.docking import cofactor_hetatm_pdb

    if not (dock.get("pdb_path") and dock.get("box")):
        raise HTTPException(status_code=400, detail="receptor/box not ready")
    try:
        import py3Dmol
    except ImportError:
        raise HTTPException(status_code=503, detail="py3Dmol not installed")
    pdb_path = dock["pdb_path"]
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    atom_lines = [ln for ln in text.splitlines() if ln.startswith(("ATOM", "TER", "END"))]
    het_pdb = cofactor_hetatm_pdb(pdb_path, receptor_text=text)
    c, s = dock["box"]["center"], dock["box"]["size"]
    view = py3Dmol.view(width=980, height=560)
    view.addModel("\n".join(atom_lines) + "\n", "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})
    if het_pdb:
        view.addModel(het_pdb + "\n", "pdb")
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "greyCarbon"}})
    view.addBox({"center": {"x": c[0], "y": c[1], "z": c[2]},
                 "dimensions": {"w": s[0], "h": s[1], "d": s[2]},
                 "color": "orange", "opacity": 0.6})
    view.zoomTo()
    try:
        html = view._make_html()
    except AttributeError:
        html = view.write()
    return Response(content=_embed_centered(html), media_type="text/html")


class DockRunBody(BaseModel):
    pass


@router.post("/api/docking/{sid}/run")
def docking_run(sid: str, body: DockRunBody) -> dict:
    dock = _session(sid)
    from datetime import datetime

    if not (dock.get("rec_pdbqt") and dock.get("lig_pdbqt") and dock.get("box")):
        raise HTTPException(status_code=400, detail="ligand/receptor/box not ready")

    def task(progress) -> dict:
        from WY.core.docking import docking_timeout, parse_vina_scores, run_vina, split_vina

        progress(fraction=None, text=f"vina 运行中(exhaustiveness={config.EXHAUSTIVENESS},"
                                     f" 最长 {docking_timeout()}s)…")
        out_pdbqt = os.path.join(config.DOCKING_RAW_CACHE_DIR,
                                 f"vina_out_{datetime.now():%Y%m%d_%H%M%S}.pdbqt")
        proc = run_vina(dock["rec_pdbqt"], dock["lig_pdbqt"], dock["box"], out_path=out_pdbqt)
        progress(fraction=0.9, text="解析构象与打分…")
        poses = split_vina(out_pdbqt, out_dir=os.path.dirname(out_pdbqt))
        scores = parse_vina_scores(out_pdbqt)
        # 新一轮对接:清掉上一轮残留(与原页面一致)
        for _k in ("pose_sel", "conclusion", "result_dir", "starred"):
            dock.pop(_k, None)
        dock["result"] = {"log_tail": (proc.stdout or "")[-1500:], "out": out_pdbqt,
                          "poses": poses, "scores_rows": scores.to_dict(orient="records")}
        return {"log_tail": dock["result"]["log_tail"],
                "scores": dock["result"]["scores_rows"], "poses": poses}

    return {"job_id": jobs.submit(task, "docking vina run")}


class DockSaveBody(BaseModel):
    pose_mode: int
    conclusion: str = ""
    drug: str = ""
    disease: str = ""


@router.post("/api/docking/{sid}/save")
def docking_save(sid: str, body: DockSaveBody) -> dict:
    dock = _session(sid)
    import pandas as pd

    from WY.core.docking import docking_result_prefix, save_docking_result

    if not dock.get("result"):
        raise HTTPException(status_code=400, detail="no vina result to save")
    drug = body.drug.strip() or dock.get("drug", "")
    disease = body.disease.strip() or dock.get("disease", "")
    gene = dock.get("gene", "gene")
    prefix = docking_result_prefix(drug, disease, gene)
    scores = pd.DataFrame(dock["result"]["scores_rows"])
    result_dir = save_docking_result(
        prefix, gene, dock.get("pdb_id", ""), dock.get("cid"), dock["box"],
        dock["result"]["out"], scores, body.conclusion,
        receptor_source=dock.get("receptor_source", "experiment"),
        pose_mode=int(body.pose_mode))
    dock["result_dir"] = result_dir
    dock["drug_name"], dock["disease_name"] = drug, disease
    dock.pop("starred", None)
    return {"result_dir": result_dir}


class DockStarBody(BaseModel):
    on: bool


@router.post("/api/docking/{sid}/star")
def docking_star(sid: str, body: DockStarBody) -> dict:
    dock = _session(sid)
    import glob
    import shutil as _shutil

    if not dock.get("result_dir"):
        raise HTTPException(status_code=400, detail="save the result first")
    if body.on:
        from WY.core.docking import star_docking_result

        star_docking_result(dock["result_dir"])
        dock["starred"] = True
        return {"starred": True}
    star_base = os.path.basename(dock["result_dir"])
    for p in glob.glob(os.path.join(config.MD_STAR_DIR, star_base) + "*"):
        if os.path.isdir(p):
            _shutil.rmtree(p, ignore_errors=True)
    dock.pop("starred", None)
    return {"starred": False}


@router.get("/api/docking/{sid}/view3d")
def docking_view3d(sid: str, pose: int, cutoff: float = 4.0, focus: int = 0) -> Response:
    """Step 6 3D 预览:受体卡通 + 构象绿棍 + 口袋残基橙棍标注 + 相互作用虚线 +
    共晶配体参照;服务端 py3Dmol 渲染,图例由前端展示。"""
    dock = _session(sid)
    if not dock.get("result") or not dock.get("pdb_path"):
        raise HTTPException(status_code=400, detail="no docking result")
    try:
        import py3Dmol
    except ImportError:
        raise HTTPException(status_code=503, detail="py3Dmol not installed")
    from WY.core.docking import (cofactor_hetatm_pdb, extract_model,
                                 interactions_near_ligand, pose_pdbqt_to_pdb,
                                 residues_near_ligand)

    out_stem = os.path.splitext(os.path.basename(dock["result"]["out"]))[0]
    pose_path = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{out_stem}_pose_{int(pose)}.pdbqt")
    extract_model(dock["result"]["out"], int(pose), pose_path)
    pose_pdb = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{out_stem}_pose_view.pdb")
    try:
        pose_pdbqt_to_pdb(pose_path, pose_pdb)
    except Exception:
        raise HTTPException(status_code=500, detail="pose pdbqt→pdb conversion failed")
    with open(dock["pdb_path"], encoding="utf-8", errors="replace") as f:
        receptor_text = f.read()
    with open(pose_pdb, encoding="utf-8") as f:
        ligand_pdb = f.read()
    het_pdb = cofactor_hetatm_pdb(dock["pdb_path"], receptor_text=receptor_text)
    pocket = residues_near_ligand(dock["pdb_path"], pose_path, cutoff=float(cutoff))
    interactions = interactions_near_ligand(dock["pdb_path"], pose_path)

    view = py3Dmol.view(width=1100, height=560)
    view.addModel(receptor_text, "pdb")
    view.addModel(ligand_pdb, "pdb")
    cartoon: dict = {"color": "spectrum"}
    if focus:
        cartoon["opacity"] = 0.3
    view.setStyle({"model": 0}, {"cartoon": cartoon})
    view.setStyle({"model": 1}, {"stick": {"colorscheme": "greenCarbon"}})
    if het_pdb:
        view.addModel(het_pdb + "\n", "pdb")
        view.setStyle({"model": 2}, {"stick": {"colorscheme": "greyCarbon",
                                               "opacity": 0.5, "radius": 0.15}})
    by_chain: dict[str, list] = {}
    for r in pocket or []:
        by_chain.setdefault(r["chain"], []).append(r["resi"])
    for chain, resis in by_chain.items():
        sel: dict = {"model": 0, "resi": resis}
        if chain:
            sel["chain"] = chain
        view.setStyle(sel, {"cartoon": cartoon,
                            "stick": {"colorscheme": "orangeCarbon", "radius": 0.2}})
    for r in pocket or []:
        x, y, z = r["pos"]
        view.addLabel(f"{r['resn']}-{r['resi']}", {
            "position": {"x": x, "y": y, "z": z},
            "fontSize": 11, "fontColor": "black",
            "showBackground": True, "backgroundColor": "white",
            "backgroundOpacity": 0.7,
        })
    for it in interactions or []:
        view.addCylinder({"dashed": True, "radius": 0.06, "dashLength": 0.3, "gapLength": 0.22,
                          "color": {"hbond": "#0f62fe", "salt": "#c026d3",
                                    "hydro": "#808080"}.get(it["kind"], "#777777"),
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
    if focus:
        view.zoomTo({"model": 1})
        view.zoom(0.55)
    else:
        view.zoomTo()
    try:
        html = view._make_html()
    except AttributeError:
        html = view.write()
    return Response(content=_embed_centered(html), media_type="text/html")


@router.get("/api/docking/{sid}/interactions")
def docking_interactions(sid: str, pose: int) -> dict:
    dock = _session(sid)
    from WY.core.docking import extract_model, interactions_near_ligand

    if not dock.get("result") or not dock.get("pdb_path"):
        raise HTTPException(status_code=400, detail="no docking result")
    out_stem = os.path.splitext(os.path.basename(dock["result"]["out"]))[0]
    pose_path = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{out_stem}_pose_{int(pose)}.pdbqt")
    extract_model(dock["result"]["out"], int(pose), pose_path)
    return {"rows": interactions_near_ligand(dock["pdb_path"], pose_path)}


class DockSave3dBody(BaseModel):
    pose: int
    caption: str = ""
    cutoff: float = 4.0


@router.post("/api/docking/{sid}/save3d")
def docking_save3d(sid: str, body: DockSave3dBody) -> dict:
    """保存交互式 3D 视图 HTML 到结果目录(caption 由前端传入,离线可用)。
    口袋残基与相互作用随 cutoff 一并嵌入,保证导出文件与页面预览所见一致。"""
    dock = _session(sid)
    from WY.core.docking import (interactions_near_ligand, residues_near_ligand,
                                 save_interactive_3d_html)

    if not (dock.get("result") and dock.get("pdb_path") and dock.get("result_dir")):
        raise HTTPException(status_code=400, detail="save the result first")
    out_stem = os.path.splitext(os.path.basename(dock["result"]["out"]))[0]
    pose_path = os.path.join(config.DOCKING_RAW_CACHE_DIR,
                             f"{out_stem}_pose_{int(body.pose)}.pdbqt")
    if not os.path.exists(pose_path):
        from WY.core.docking import extract_model

        extract_model(dock["result"]["out"], int(body.pose), pose_path)
    pocket = residues_near_ligand(dock["pdb_path"], pose_path, cutoff=float(body.cutoff))
    interactions = interactions_near_ligand(dock["pdb_path"], pose_path)
    file_name = f"3d_view_{safe_name(dock.get('gene', 'target'))}_pose{int(body.pose)}.html"
    out_path, offline_ok = save_interactive_3d_html(
        dock["pdb_path"], pose_path,
        os.path.join(dock["result_dir"], file_name), body.caption,
        pocket=pocket, interactions=interactions)
    rel = os.path.relpath(out_path, config.DOWNLOADS_DIR).replace(os.sep, "/")
    return {"path": out_path, "artifact_rel": rel, "offline_ok": offline_ok}


@router.get("/api/docking/{sid}/zip")
def docking_zip(sid: str) -> Response:
    dock = _session(sid)
    from WY.core.files import make_dir_zip

    if not dock.get("result_dir"):
        raise HTTPException(status_code=400, detail="save the result first")
    data = make_dir_zip(dock["result_dir"])
    return Response(content=data, media_type="application/zip",
                    headers={"Content-Disposition":
                             f'attachment; filename="{os.path.basename(dock["result_dir"])}.zip"'})
