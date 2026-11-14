"""分子对接模块(AutoDock Vina)。

流程:把受体蛋白和配体小分子转成对接所需格式,确定结合口袋范围,
运行 Vina 对接并解析打分。中间文件会缓存一段时间以便复用,对接结果只归档不缓存。

结构解析/相互作用分析与交互式 3D 视图导出分别拆在 docking_interactions 与 docking_3d。
"""

import hashlib
import logging
import os
from datetime import datetime
from typing import Iterator

import pandas as pd
import requests

from .caching import load_json_file, save_json_file, save_text_atomic
from .text_utils import load_name_map, safe_name
from WY.i18n import t

# 结构解析/相互作用(docking_interactions)与 3D 视图导出(docking_3d)已拆分为子模块;
# 此处再导入以维持 WY.core.docking 既有对外接口(api/cli/tests 均从本模块导入这些名字)。
from .docking_3d import (build_interactive_3d_html, cofactor_hetatm_pdb,
                         pose_pdbqt_to_pdb, save_interactive_3d_html)
from .docking_interactions import (_pdb_coords, _run_obabel, _strip_receptor_pdb,
                                   interactions_near_ligand,
                                   parse_hetatm, residues_in_box, residues_near_ligand)


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


def download_pdb(pdb_id: str) -> str:
    """按 PDB 编号下载受体结构文件;之前下载过则直接复用本地文件。返回本地路径。

    RCSB 对部分条目不提供旧版 PDB 格式(逐条目判断,如 8VHE),.pdb 404 时
    自动改下 mmCIF 并按 _atom_site 原样转写成 PDB(cif_atom_site_to_pdb),
    ATOM/HETATM 语义不变,后续清洗/转 pdbqt/盒子推断链路无感。
    """
    from WY import config
    from .docking_interactions import cif_atom_site_to_pdb

    pdb_id = pdb_id.strip().lower()
    os.makedirs(config.DOCKING_RAW_CACHE_DIR, exist_ok=True)
    path = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{pdb_id}.pdb")
    if os.path.exists(path):
        return path
    cif_path = os.path.join(config.DOCKING_RAW_CACHE_DIR, f"{pdb_id}.cif")
    used_cif = os.path.exists(cif_path)  # 上次兜底已下载:条目发布后内容冻结,直接复用
    if not used_cif:
        pdb_url = config.PDB_DOWNLOAD_URL.format(pdb_id=pdb_id)
        for url, dest in ((pdb_url, path), (pdb_url[:-4] + ".cif", cif_path)):
            resp = requests.get(url, timeout=config.API_TIMEOUT)
            if resp.status_code == 404:
                continue  # 该条目没有这种格式,换下一个
            resp.raise_for_status()
            # 原子写:同一 PDB 被两个会话并发处理时,半截文件不会被"存在即复用"判定永久信任
            save_text_atomic(dest, resp.text)
            used_cif = dest == cif_path
            break
        else:
            raise RuntimeError(
                f"RCSB 未提供 {pdb_id.upper()} 的结构文件(PDB 与 mmCIF 格式均不可得)")
    if used_cif and cif_atom_site_to_pdb(cif_path, path) == 0:
        raise RuntimeError(f"mmCIF 转写得到 0 个原子:{cif_path}")
    return path


# -------------------- 配体链 --------------------


class VinaInterrupted(RuntimeError):
    """vina 进程被控制台 Ctrl+C 终止(Windows rc=0xC000013A,STATUS_CONTROL_C_EXIT):
    中断信号同时送达子进程,父进程的 KeyboardInterrupt 可能不触发,须按返回码识别并按"中断"处理。"""


class PubChemUnavailable(RuntimeError):
    """PubChem 暂不可用:网络故障或服务端错误(503 重试后仍失败),与"查无化合物"区分。"""


def resolve_cid(query: str) -> int | None:
    """把药物名换成 PubChem 化合物编号(CID);输入本身就是编号时直接使用。

    查无化合物(404)返回 None;网络故障或服务端错误抛 PubChemUnavailable——
    其中 503 ServerBusy 为瞬时抖动,按退避自动重试 PUBCHEM_BUSY_RETRIES 次,
    其余不重试直接抛出。调用方据此给出与"查无化合物"不同的提示。
    """
    import time

    from WY import config
    import pubchempy

    query = query.strip()
    if query.isdigit():
        return int(query)
    attempts = config.PUBCHEM_BUSY_RETRIES + 1
    for attempt in range(attempts):
        try:
            comps = pubchempy.get_compounds(query, "name")
            return comps[0].cid if comps else None
        except pubchempy.NotFoundError:          # 404:真·查无此药
            return None
        except pubchempy.PubChemHTTPError as exc:
            # 503 服务端繁忙为瞬时抖动:退避重试;其余 HTTP 错误重试无意义
            if exc.code != 503 or attempt >= attempts - 1:
                raise PubChemUnavailable(f"HTTP {exc.code} {exc.msg}") from exc
            time.sleep(config.PUBCHEM_BUSY_BACKOFF * (2 ** attempt))
        except Exception as exc:                 # requests 网络层异常(断网/DNS/超时)
            raise PubChemUnavailable(str(exc) or type(exc).__name__) from exc


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
    # MMFF 力场的第二个位置参数(分子参数表)必需:个别分子无 MMFF 参数时跳过
    # 能量最小化,直接用第一个构象,不让配体准备因此失败
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    for cid in cids:
        ff = (AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid)
              if mmff_props is not None else None)
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
    from WY import config
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
    from WY import config
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
    rcsb_entry_info {
      resolution_combined
      polymer_entity_count_protein
      polymer_entity_count_DNA
      polymer_entity_count_RNA
    }
  }
}
"""


def composition_flag(n_protein, n_dna, n_rna) -> str:
    """按聚合物组分给受体检索结果加防呆标注(2026-10-15 AGT 三重撞名事故的防再犯措施)。

    基因符号全文检索会混入撞名条目,其中最危险的是不含任何蛋白质的纯核酸结构
    (如事故中的 2MMQ)。返回值:"dna_only"(无蛋白,不能作蛋白对接受体)、
    "has_nucleic"(蛋白+核酸复合物,受体清洗不剥离核酸原子,警示)、""(正常或组分未知)。
    某项计数缺失(旧条目无该字段)时按保守策略处理:不因未知而阻断,已知含核酸即警示。
    """

    def _i(v):
        return int(v) if isinstance(v, (int, float)) else 0

    has_nucl = _i(n_dna) + _i(n_rna) > 0
    if not has_nucl:
        return ""
    if n_protein == 0:
        return "dna_only"
    return "has_nucleic"


def _rcsb_search_entry_ids(gene_symbol: str, max_hits: int) -> list[str]:
    """在 RCSB 中用基因符号精确检索结构,返回 entry 编号列表(最多 max_hits 条)。
    与旧的逐基因检索完全同一查询,保证召回一致;失败向上抛错。
    """
    from WY import config
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
        cache_ttl_hours=config.RCSB_SEARCH_CACHE_TTL_HOURS,
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
    from WY import config
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
            info = e.get("rcsb_entry_info") or {}
            res = info.get("resolution_combined")
            res = res[0] if isinstance(res, list) and res else res
            out[eid] = {
                "resolution": res,
                "method": ", ".join(m.get("method", "") for m in e.get("exptl") or []),
                "title": (e.get("struct") or {}).get("title", "") or "",
                # 聚合物组分计数:供 composition_flag 做防呆标注(缺失字段按未知处理)
                "n_protein": info.get("polymer_entity_count_protein"),
                "n_dna": info.get("polymer_entity_count_DNA"),
                "n_rna": info.get("polymer_entity_count_RNA"),
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
    from WY import config

    uniq = sorted({g.strip().upper() for g in genes if g and g.strip()})
    empty = pd.DataFrame(columns=["pdb_id", "resolution", "method", "title", "source", "flag"])
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

    # 3) 按基因组装结果表(与旧版字段完全一致,新增 flag 防呆标注列;个别 entry 元数据缺失按缺信息展示)
    for g, ids in gene_ids.items():
        rows = [{"pdb_id": eid,
                 "resolution": meta.get(eid, {}).get("resolution"),
                 "method": meta.get(eid, {}).get("method", ""),
                 "title": meta.get(eid, {}).get("title", ""),
                 "source": "experiment",
                 "flag": composition_flag(meta.get(eid, {}).get("n_protein"),
                                          meta.get(eid, {}).get("n_dna"),
                                          meta.get(eid, {}).get("n_rna"))}
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
                        "flag": "",
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
    from WY import config

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

    # -xr:刚性受体模式(不加扭转树),加氢与 AutoDock 原子类型由 pdbqt 写出器完成
    _run_obabel(source, out_path, "-xr")
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


def docking_timeout() -> int:
    """vina 有效超时(秒):随搜索深度线性放宽——EXHAUSTIVENESS 每放大,墙钟预算同比增加,
    单条采样链的预算不变。EXHAUSTIVENESS<=8 时即基线 DOCKING_TIMEOUT(默认行为不变,
    16→2 倍,32→4 倍);钳制在基线的 8 倍,防手改设置文件把超时写成天文数字。"""
    from WY import config

    base = max(int(config.DOCKING_TIMEOUT), 1)
    exh = max(int(config.EXHAUSTIVENESS), 1)
    return min(max(base, base * exh // 8), base * 8)


def run_vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    box: dict,
    out_path: str | None = None,
) -> "subprocess.CompletedProcess":
    """调用 AutoDock Vina,在指定口袋范围内对受体和配体执行分子对接;超时或运行失败时报错,运行日志(含打分)保留在返回结果中。"""
    import subprocess
    from WY import config

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
    timeout = docking_timeout()
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=timeout,
            encoding="utf-8", errors="replace",
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"vina 对接超时（>{timeout} 秒，搜索深度 {config.EXHAUSTIVENESS}）；"
            "可在设置页调低搜索深度后重试")
    if proc.returncode == 3221225786:  # 0xC000013A:Ctrl+C 经控制台进程组杀死了 vina
        raise VinaInterrupted(f"vina 被控制台 Ctrl+C 终止(rc={proc.returncode})")
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
    from WY import config

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
                # 排除原始多模型输出自身:glob 的 f"{stem}*" 会把它一并收进来,
                # 且字典序还排在各 pose 文件之前,混进 pose 列表会误导调用方
                out_name = os.path.basename(out_pdbqt)
                files = sorted(
                    f for f in glob.glob(os.path.join(out_dir, f"{stem}*"))
                    if os.path.basename(f) != out_name)
                if files:
                    return files
        except (OSError, subprocess.SubprocessError):
            # OSError:可执行文件缺失/无法执行;SubprocessError:含 TimeoutExpired
            # (vina_split 卡死超时)。两者都回落到内置拆分,不让超时直接炸到页面
            logging.warning("vina_split failed (%s), fallback to python split", split_exe)

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
    from WY import config

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
    from WY import config

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
    from WY import config

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


def docking_result_prefix(drug: str, disease: str, gene: str) -> str:
    """归档结果目录名前缀:三段各自过 safe_name 消毒,防止用户输入中的路径分隔符等
    特殊字符进入文件路径;空段必须先剔除再消毒(safe_name 对空串返回 "unnamed",
    直接套会凭空出现 unnamed 段),全空回退 "docking"。safe_name 只用于目录命名,
    meta 等人读字段仍用原始输入。"""
    parts = []
    for seg in (drug, disease, gene):
        seg = (seg or "").strip()
        if seg:
            parts.append(safe_name(seg))
    return "_".join(parts) or "docking"


def split_result_prefix(prefix: str) -> tuple[str, str]:
    """把结果目录名/富集缓存名前缀拆回 (药物名, 疾病名)。

    前缀是 f"{safe(drug)}_{safe(disease)}" 拼接,而药物名自身可能含空格
    (如 "Salicylic acid" → "Salicylic_acid"),简单 split("_", 1) 会把药物名截断,
    导致后续用半截名字去 PubChem 取错配体、归档前缀张冠李戴。
    这里借 name_map(safe 名 → 原始名,各数据源抓取时都有记录)贪心匹配:
    从最长的药物名前缀试起,药物与疾病两段都能还原才采用;识别不出时退回旧拆法。
    """
    name_map = load_name_map()
    parts = prefix.split("_")
    for i in range(len(parts) - 1, 0, -1):
        drug = name_map.get("_".join(parts[:i]))
        disease = name_map.get("_".join(parts[i:]))
        if drug and disease:
            return str(drug), str(disease)
    parts = prefix.split("_", 1)
    return (parts[0], parts[1] if len(parts) > 1 else "")


def save_docking_result(prefix: str, gene: str, pdb_id: str, cid: int | None, box: dict,
                        out_pdbqt: str, scores: pd.DataFrame, conclusion: str,
                        receptor_source: str = "experiment",
                        pose_mode: int | None = None) -> str:
    """归档本次对接成果:打分表、所选构象和参数信息存入结果文件夹,并清理旧版本。
    pose_mode 提供时归档用户选定的 pose,未提供(默认)归档打分最佳的构象。"""
    from WY import config

    gene = safe_name(gene)  # 清理基因名中的非法字符,再用于文件命名
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    result_dir = os.path.join(config.MD_RESULT_DIR, f"{prefix}_{ts}")
    os.makedirs(result_dir, exist_ok=True)

    # 打分表和结论
    summary = scores.copy()
    summary["conclusion"] = conclusion
    summary_path = os.path.join(result_dir, f"summary_{ts}.csv")
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")

    # 归档构象:选定了 pose 就按选择归档,否则按亲和力排序取最佳;
    # 选择值不在打分表里时回退最佳,保证归档必有构象
    best_mode = int(summary.iloc[0]["mode"])
    save_mode = best_mode
    if pose_mode is not None and int(pose_mode) in set(summary["mode"].tolist()):
        save_mode = int(pose_mode)
    best_path = os.path.join(result_dir, f"best_{gene}_{ts}.pdbqt")
    extract_model(out_pdbqt, save_mode, best_path)

    # 参数信息(药物、靶点、对接盒、打分等),便于日后复查
    write_docking_metadata(os.path.join(result_dir, f"meta_{ts}.txt"), {
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
        "best_mode": save_mode,
        "best_affinity": float(
            summary.loc[summary["mode"] == save_mode, "affinity"].iloc[0]),
        "conclusion": conclusion,
    })

    # 同名结果只保留最近几版,并清理过期中间文件
    prune_docking_results(result_dir)
    clean_docking_cache()
    return result_dir


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
    from WY import config
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
    from WY import config
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
    from WY import config

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

    from WY import config
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
