"""网络药理学流水线 CLI:core/ 之上的薄封装。

与网页共享缓存与设置(ui_settings.json / 各级请求缓存),不改 core 一行代码;
处理函数一律返回退出码(0 成功 / 1 失败)。
"""

import argparse
import hashlib
import logging
import os
import re
import sys
import time
from datetime import datetime

# 导入统一走 WY 包绝对导入,可见性由可编辑安装(pip install -e .)提供,
# 全仓不再有 sys.path 补丁(selfcheck_architecture 守护)。入口:npap / python -m WY.cli。
from WY import i18n

i18n.set_lang_provider(lambda: "en")  # 纪律 1:core 内 t() 输出钉死英文,先于一切 core 调用

import pandas as pd

from WY.core import interactome
from WY.core.caching import load_json_file, save_json_file
from WY.core.enrichment import (TRUSTED_REASONS, enrichment_retry_hint_minutes,
                             run_enrichment)
from WY.core.logging_setup import setup_logging
from WY.core.network_algorithms import run_proximity_cached, run_rwr_cached
from WY.core.ppi import (build_ppi_graph, compute_mcc_scores, compute_ppi_metrics,
                      fetch_string_ppi, mark_hub_genes, pick_top_mcc_genes)
from WY.core.targets import (get_disease_targets_merged, get_drug_targets_merged,
                          get_targets_concurrent, intersect_target_symbols)
from WY.core.analysis import expand_targets_with_rwr, pair_prefix
from WY.core.text_utils import safe_name
from WY import config

# _slug 结果目录名的最大长度:Windows 路径 260 字符上限下给深层文件名留足余量
_SLUG_MAX_LEN = 120
# analyze 与网页的 RWR 扩展/邻近度/RWR 缓存 label 保持同一格式,跑过的对才能跨端命中缓存


# ======== 通用小工具 ========

def _slug(*parts: str) -> str:
    """把药物/疾病等多段名字拼成合法目录名:各段先 safe_name,再以双下划线连接;
    超长截断并追加 6 位 MD5 短哈希防撞(analyze/batch 的每对结果子目录名)。"""
    name = "__".join(safe_name(p) for p in parts)
    if len(name) > _SLUG_MAX_LEN:
        digest = hashlib.md5(name.encode("utf-8")).hexdigest()[:6]
        name = f"{name[:_SLUG_MAX_LEN - 7]}~{digest}"
    return name


def _resolve_out_dir(out_arg: str | None) -> str | None:
    """--out 落点三级解析(先命中先用):命令行参数 > API.env 的 CLI_RESULTS_DIR > 兜底 ./results。
    第二级取 config.CLI_RESULTS_DIR(config 导入时经 load_dotenv 读到的本机路径,与 VINA_DIR 同类,不入库)。
    --out 指向已存在的文件、或目录创建失败(权限/盘符不存在)时警告并返回 None,由调用方按失败退出。"""
    if out_arg:
        path = os.path.abspath(out_arg)
        if os.path.exists(path) and not os.path.isdir(path):
            _warn(f"--out must be a directory, but a file already exists at: {path}")
            return None
    else:
        path = os.path.abspath(config.CLI_RESULTS_DIR if config.CLI_RESULTS_DIR else "results")
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        _warn(f"cannot create output directory '{path}': {e}")
        return None
    return path


def _save_csv(df, path: str) -> str | None:
    """结果表统一 utf-8-sig 落盘(Excel 直接打开不乱码);空表不落盘(与 core 约定一致),返回实际路径或 None。
    目标被其他程序占用(如 Excel 开着该文件)时警告并返回 None,不让已完成的整次分析以崩溃收场。"""
    if df is None or df.empty:
        return None
    try:
        df.to_csv(path, index=False, encoding="utf-8-sig")
    except PermissionError:
        _warn(f"cannot write '{path}': file is locked by another program "
              "(close it, e.g. Excel, and rerun this command to refresh the file)")
        return None
    return path


def _save_json(data: dict, path: str) -> str | None:
    """JSON 结果落盘(复用 core 的原子写,目录不存在自动创建),返回路径;被其他程序占用时警告并返回 None。"""
    try:
        save_json_file(path, data, indent=2, sort_keys=True)
    except PermissionError:
        _warn(f"cannot write '{path}': file is locked by another program "
              "(close it and rerun this command to refresh the file)")
        return None
    return path


def _read_genes(value: str) -> list[str]:
    """解析基因清单:存在的路径按文件读(.csv 取 gene/approvedSymbol 列,其余按一行一个基因读,
    # 注释忽略,兼容单行逗号分隔——方案端到端示例直接喂 common_targets.csv);否则按
    逗号/分号/空白分隔的字面基因列表解析。结果去重(保序)。"""
    if os.path.isdir(value):
        _warn(f"expected a gene-list file or a literal gene list, but this is a directory: {value}")
        return []
    if os.path.isfile(value):
        if value.lower().endswith(".csv"):
            df = pd.read_csv(value, encoding="utf-8-sig", dtype=str)
            if df.empty:
                return []
            col = next((c for c in ("gene", "approvedSymbol") if c in df.columns),
                       df.columns[0])
            items = df[col].dropna().astype(str).tolist()
            text = "\n".join(items)
        else:
            with open(value, "r", encoding="utf-8-sig") as f:
                text = "\n".join(line.partition("#")[0] for line in f)
    else:
        text = value
    genes = [g for g in (tok.strip() for tok in re.split(r"[\s,;]+", text)) if g]
    return list(dict.fromkeys(genes))


def _parse_sources(value: str) -> tuple[str, ...] | None:
    """解析 --sources 逗号分隔列表(大小写不敏感,去重保序);不在 config.ENRICHMENT_SOURCES
    白名单内的值报错返回 None,避免拼错的数据源静默走空结果。"""
    toks = [t.strip().upper() for t in value.split(",") if t.strip()]
    bad = [t for t in toks if t not in config.ENRICHMENT_SOURCES]
    if bad:
        _warn(f"unknown enrichment source(s): {', '.join(bad)}; "
              f"valid: {', '.join(config.ENRICHMENT_SOURCES)}")
        return None
    return tuple(dict.fromkeys(toks))


def _pvalue(value: str) -> float:
    """--pvalue 的类型校验:必须是 (0, 1] 内的数值,越界按用法错误退出(argparse 退出码 2)。
    在进任何网络请求前拦下,避免非法值导致富集接口失败后误触发与网页共享的富集熔断器。"""
    try:
        v = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid float value: '{value}'")
    if not 0 < v <= 1:
        raise argparse.ArgumentTypeError(f"pvalue must be in (0, 1]; got {v}")
    return v


def _warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


# ======== 共享分析步骤(独立子命令与 run_pair 共用,纪律 6 的测试接缝) ========

def _targets_step(drug: str, disease: str | None, out_dir: str,
                  progress_cb=None) -> tuple:
    """获取药物/疾病靶点(四数据源,降噪整合),明细表落盘;
    返回 (药物表, 疾病表或 None, 失败源列表, 实际落盘路径列表)。"""
    written: list[str] = []
    if disease:
        drug_df, disease_df, failed = get_targets_concurrent(drug, disease,
                                                             progress_cb=progress_cb)
    else:
        drug_df, disease_df, failed = get_drug_targets_merged(drug), None, []
    path = _save_csv(drug_df, os.path.join(out_dir, "drug_targets.csv"))
    if path:
        written.append(path)
    if disease:
        path = _save_csv(disease_df, os.path.join(out_dir, "disease_targets.csv"))
        if path:
            written.append(path)
    return drug_df, disease_df, failed, written


def _ppi_step(genes: list[str], out_dir: str, label: str | None) -> dict:
    """STRING PPI 网络 → 拓扑指标 → MCC 核心靶点,边表/指标表/核心靶点表落盘。
    返回 dict(ok, detail, edges, metrics, core):ok=核心靶点可用;边表/指标表非空也会写出,
    供失败时排查(基因不被识别、无达标互作等)。"""
    truncated = len(genes) > config.STRING_QUERY_LIMIT
    if truncated:
        genes = genes[:config.STRING_QUERY_LIMIT]
    edges_df = fetch_string_ppi(
        genes, species=config.SPECIES_MAP["Human"]["ncbi_taxon"],
        score=config.PPI_CONFIDENCE, network_type="full", label=label)
    if edges_df is None or edges_df.empty:
        return {"ok": False, "truncated": truncated, "written": [], "detail":
                "STRING returned no interactions (genes unrecognized, or nothing "
                "above the confidence threshold)", "edges": None, "metrics": None,
                "core": []}
    written: list[str] = []
    p = _save_csv(edges_df, os.path.join(out_dir, "ppi_edges.csv"))
    if p:
        written.append(p)

    graph = build_ppi_graph(edges_df, genes)
    metrics_df = compute_ppi_metrics(graph)
    try:
        mcc_map = compute_mcc_scores(graph)
    except Exception:
        logging.exception("MCC computation failed; core targets skipped")
        mcc_map = {}
    if mcc_map and "mcc" not in metrics_df.columns:
        # hub_method 为 degree 时 mark_hub_genes 不带 mcc 列,提前补上保证指标表完整
        metrics_df["mcc"] = metrics_df["gene"].map(mcc_map).fillna(0.0)
    metrics_df = mark_hub_genes(
        metrics_df, top_n=config.HUB_TOP_N, method=config.HUB_METHOD,
        mcc_map=mcc_map, small_net_factor=config.HUB_SMALL_NET_FACTOR)
    p = _save_csv(metrics_df, os.path.join(out_dir, "ppi_metrics.csv"))
    if p:
        written.append(p)

    core_genes = pick_top_mcc_genes(metrics_df, mcc_map, config.MCC_TOP_N)
    if core_genes:
        # 核心靶点按 MCC 排名顺序输出完整指标行
        rank = {g: i for i, g in enumerate(core_genes)}
        core_df = metrics_df[metrics_df["gene"].isin(rank)].copy()
        core_df = core_df.sort_values("gene", key=lambda s: s.map(rank))
        p = _save_csv(core_df, os.path.join(out_dir, "core_targets.csv"))
        if p:
            written.append(p)

    detail = f"{len(edges_df)} edges, {len(metrics_df)} nodes " \
             f"({int((metrics_df['degree'] > 0).sum())} with interactions)"
    if truncated:
        detail += f"; input truncated to {config.STRING_QUERY_LIMIT}"
    if not core_genes:
        detail += "; core targets unavailable (MCC failed or no interactions)"
    ok = bool(core_genes) and len(written) == 3
    if ok is False and core_genes:
        detail += "; some result files could NOT be written (file locked)"
    return {"ok": ok, "truncated": truncated, "detail": detail,
            "edges": edges_df, "metrics": metrics_df, "core": core_genes,
            "written": written}


def _enrichment_step(genes: list[str], out_dir: str, label: str | None,
                     organism: str | None = None, pvalue: float | None = None,
                     sources: tuple | None = None) -> dict:
    """GO/KEGG 富集(g:Profiler 主 + Enrichr 自动降级),结果表落盘。
    状态语义必须区分:no_significant=两源都成功但真无显著项(正常,不算失败);
    both_sources_down / gprofiler_down_no_fallback=服务不可用,结果不可信需冷却后重跑。"""
    status: dict = {}
    df, provider = run_enrichment(
        genes, organism=organism or config.SPECIES_MAP["Human"]["gprofiler"],
        pvalue=pvalue if pvalue is not None else config.ENRICHMENT_PVALUE,
        sources=sources, label=label, status=status)
    if df is not None and not df.empty:
        path = _save_csv(df, os.path.join(out_dir, "enrichment.csv"))
        if path:
            return {"ok": True, "status": "ok", "reason": "success", "path": path,
                    "detail": f"{len(df)} significant terms via {provider}"}
        return {"ok": False, "status": "failed", "reason": "write_locked", "path": None,
                "detail": f"{len(df)} significant terms computed via {provider}, but "
                          "enrichment.csv could NOT be written (file locked)"}
    reason = status.get("reason", "")
    if reason in TRUSTED_REASONS:
        return {"ok": True, "status": "ok", "reason": reason,
                "detail": "no significant terms (both sources succeeded)"}
    hint = ""
    minutes = enrichment_retry_hint_minutes()
    if minutes:
        hint = f"; circuit cooling, retry in ~{minutes} min"
    return {"ok": False, "status": "failed", "reason": reason or "empty",
            "detail": f"enrichment unavailable ({reason or 'empty result'}){hint}; "
                      "result NOT usable for literature comparison, rerun later"}


def _proximity_step(drug_syms: list[str], disease_syms: list[str],
                    out_dir: str, label: str | None) -> dict:
    """网络邻近度 Z-score(依赖背景互作组),结果落盘 proximity.json。
    背景互作组未构建时跳过并说明,绝不自动触发数 GB 下载(构建走网页或脚本)。"""
    if not interactome.is_built():
        return {"ok": False, "status": "skipped",
                "detail": "background interactome not built; build it in the web UI "
                          "or via scripts/build_interactome.py, then rerun"}
    result = run_proximity_cached(drug_syms, disease_syms, label=label)
    if not result.get("ok"):
        return {"ok": False, "status": "failed",
                "detail": f"proximity failed: {result.get('error') or 'unknown error'}"}
    path = _save_json(result, os.path.join(out_dir, "proximity.json"))
    cached = " (cache hit)" if result.get("cached") else ""
    z, p = result.get("z"), result.get("p_empirical")
    z_txt = f"{z:.3f}" if z is not None else "n/a"
    p_txt = f"{p:.4g}" if p is not None else "n/a"
    if not path:
        return {"ok": False, "status": "failed", "result": result,
                "detail": f"Z={z_txt}, p={p_txt} computed{cached}, but proximity.json "
                          "could NOT be written (file locked)"}
    return {"ok": True, "status": "ok", "result": result, "path": path,
            "detail": f"Z={z_txt}, p={p_txt}{cached}"}


def _save_rwr_csv(result: dict, path: str) -> str | None:
    """把 RWR 结果写成 gene,score,is_seed 长表(得分降序即优先级;种子全量,非种子已由 core 截断)。"""
    if not result.get("ok"):
        return None
    rows = [{"gene": r["gene"], "score": r["score"], "is_seed": 1}
            for r in result.get("seeds") or []]
    rows += [{"gene": r["gene"], "score": r["score"], "is_seed": 0}
             for r in result.get("top_nonseed") or []]
    rows.sort(key=lambda r: -float(r["score"]))
    return _save_csv(pd.DataFrame(rows), path)


def _rwr_both_sides_step(drug_syms: list[str], disease_syms: list[str],
                         out_dir: str, label: str) -> dict:
    """双側 RWR(网页同款):疾病靶点为种子(药物靶点标注)+药物靶点为种子(疾病靶点标注),
    各落一张排序表;与网页缓存同 label,跑过的对秒回。"""
    if not interactome.is_built():
        return {"ok": False, "status": "skipped",
                "detail": "background interactome not built; build it in the web UI "
                          "or via scripts/build_interactome.py, then rerun"}
    res_d = run_rwr_cached(list(disease_syms), label=f"{label}_disease_seeds",
                           annotate=set(drug_syms))
    res_g = run_rwr_cached(list(drug_syms), label=f"{label}_drug_seeds",
                           annotate=set(disease_syms))
    p_d = _save_rwr_csv(res_d, os.path.join(out_dir, "rwr_disease_seeds.csv"))
    p_g = _save_rwr_csv(res_g, os.path.join(out_dir, "rwr_drug_seeds.csv"))
    parts = []
    for side, res, path in (("disease-seeded", res_d, p_d), ("drug-seeded", res_g, p_g)):
        if not res.get("ok"):
            parts.append(f"{side}: failed ({res.get('error') or 'unknown error'})")
        elif not path:
            parts.append(f"{side}: computed but table NOT written (file locked)")
        else:
            parts.append(f"{side}: {len(res.get('top_nonseed') or [])} ranked"
                         + (" (cache hit)" if res.get("cached") else ""))
    ok = bool(res_d.get("ok")) and bool(res_g.get("ok")) and bool(p_d) and bool(p_g)
    return {"ok": ok, "status": "ok" if ok else ("failed" if not (p_d or p_g) else "degraded"),
            "detail": "; ".join(parts)}


# ======== analyze 单对全流程(纪律 6:batch 复用此函数) ========

def _jsonable_steps(steps: dict) -> dict:
    """挑出各步结果里可 JSON 序列化的字段(summary.json 用):步骤字典里夹带的
    DataFrame 等内部大对象(run_pair 步骤间传数据用)不进摘要。"""
    out = {}
    for name, res in steps.items():
        out[name] = {k: v for k, v in res.items() if not hasattr(v, "to_csv")}
    return out


def run_pair(drug: str, disease: str, out_dir: str) -> dict:
    """单对(药物,疾病)全流程:靶点 → 交集 →(RWR 扩展)→ PPI → 富集 → 邻近度 → RWR,
    全部产出落盘 out_dir,参数快照与各步状态/耗时写入 summary.json 并返回。
    某一步失败不阻断后续步骤(相互独立的步骤照跑),状态如实记录,供调用方判定退出码。"""
    os.makedirs(out_dir, exist_ok=True)
    steps: dict[str, dict] = {}
    t0 = time.monotonic()
    prefix_full = pair_prefix(drug, disease)

    def _run(name: str, fn) -> dict:
        t = time.monotonic()
        res = fn()
        res["seconds"] = round(time.monotonic() - t, 1)
        if not res.get("status"):
            res["status"] = "ok" if res.get("ok") else "failed"
        steps[name] = res
        print(f"[analyze] {name:<11s} {res['status'] or 'ok':<8s} {res['seconds']:>6.1f}s  "
              f"{res.get('detail', '')}")
        return res

    # ---- 靶点 + 交集 ----
    def _targets_fn() -> dict:
        drug_df, disease_df, failed, _written = _targets_step(drug, disease, out_dir)
        detail_bits = []
        if drug_df is not None and not drug_df.empty:
            detail_bits.append(f"drug {len(drug_df)}")
        if disease_df is not None and not disease_df.empty:
            detail_bits.append(f"disease {len(disease_df)}")
        if failed:
            detail_bits.append(f"failed sources: {', '.join(failed)}")
        common: list[str] = []
        if drug_df is not None and not drug_df.empty \
                and disease_df is not None and not disease_df.empty:
            common = intersect_target_symbols(drug_df, disease_df,
                                              config.STRING_QUERY_LIMIT)
            _save_csv(pd.DataFrame({"approvedSymbol": common}),
                      os.path.join(out_dir, "common_targets.csv"))
            detail_bits.append(f"common {len(common)}")
        ok = bool(common) and not failed and len(_written) == 2
        if len(_written) < 2:
            detail_bits.append("some target tables NOT written (file locked)")
        return {"ok": ok, "status": "ok" if ok else "failed", "failed_sources": list(failed),
                "detail": ", ".join(detail_bits) or "no targets",
                "drug": drug_df, "disease": disease_df, "common": common}

    t_res = _run("targets", _targets_fn)
    common = t_res["common"]

    # ---- RWR 扩展(与网页 Build Report 共用 core.analysis.expand_targets_with_rwr:
    #      交集靶点在背景网络上扩展查询集;扩展失败按原始交集继续,不让分析挂掉) ----
    genes = list(common)
    prefix = prefix_full
    expansion: dict | None = None
    if int(config.RWR_EXPAND_TOP_N) > 0 and common and interactome.is_built():
        genes, prefix, rwr_exp = expand_targets_with_rwr(
            common, int(config.RWR_EXPAND_TOP_N), prefix_full)
        added = [str(r["gene"]) for r in rwr_exp["added"]]
        if rwr_exp["ok"]:
            expansion = {"added": len(added), "genes": added, "cached": rwr_exp["cached"]}
            cache_txt = " (cache hit)" if rwr_exp["cached"] else ""
            print(f"[analyze] expansion  ok        query set +{len(added)} genes "
                  f"via RWR (total {len(genes)}){cache_txt}")
        else:
            # 扩展失败按原始交集继续(与网页一致),失败原因记进 summary 供追溯
            # (旧实现在此分支引用仅成功分支定义的 added,直接 NameError,已随收敛修复)
            expansion = {"added": 0, "error": rwr_exp["error"]}
            print(f"[analyze] expansion  skipped   RWR expansion failed "
                  f"({rwr_exp['error']}); continuing with raw intersection")
            steps["targets"]["detail"] += f"; RWR expansion failed ({rwr_exp['error']})"

    # ---- PPI / 富集 / 邻近度 / RWR(相互独立,某步失败不影响其余) ----
    if common:
        p_res = _run("ppi", lambda: _ppi_step(genes, out_dir, label=prefix))
        _run("enrichment", lambda: _enrichment_step(genes, out_dir, label=prefix))
        drug_syms = sorted({str(g).strip().upper()
                            for g in t_res["drug"]["approvedSymbol"].dropna()})
        disease_syms = sorted({str(g).strip().upper()
                               for g in t_res["disease"]["approvedSymbol"].dropna()})
        _run("proximity", lambda: _proximity_step(drug_syms, disease_syms, out_dir,
                                                  label=prefix))
        _run("rwr", lambda: _rwr_both_sides_step(drug_syms, disease_syms, out_dir,
                                                 label=prefix))
    else:
        print("[analyze] empty intersection: ppi/enrichment/proximity/rwr skipped")
        for name in ("ppi", "enrichment", "proximity", "rwr"):
            steps[name] = {"ok": False, "status": "skipped", "seconds": 0.0,
                           "detail": "skipped: empty target intersection"}

    summary = {
        "drug": drug, "disease": disease,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "out_dir": os.path.abspath(out_dir),
        "failed_sources": list(steps.get("targets", {}).get("failed_sources") or []),
        "rwr_expansion": expansion,
        "params": {
            "organism": config.SPECIES_MAP["Human"]["gprofiler"],
            "enrichment_pvalue": config.ENRICHMENT_PVALUE,
            "enrichment_sources": list(config.ENRICHMENT_SOURCES),
            "string_version": config.STRING_VERSION,
            "ppi_confidence": config.PPI_CONFIDENCE,
            "target_filter_version": config.TARGET_FILTER_VERSION,
            "hub_method": config.HUB_METHOD, "hub_top_n": config.HUB_TOP_N,
            "mcc_top_n": config.MCC_TOP_N,
            "rwr_alpha": config.RWR_ALPHA,
            "rwr_expand_top_n": config.RWR_EXPAND_TOP_N,
            "interactome_built": interactome.is_built(),
        },
        "steps": _jsonable_steps(steps),
        "elapsed_seconds": round(time.monotonic() - t0, 1),
        "ok": all(s.get("status") == "ok" for s in steps.values()),
    }
    _save_json(summary, os.path.join(out_dir, "summary.json"))
    return summary


# ======== 子命令:targets ========

def _setup_targets(p: argparse.ArgumentParser) -> None:
    p.add_argument("--drug", required=True, help="Drug name, e.g. aspirin")
    p.add_argument("--disease", default=None,
                   help="Disease name in English, e.g. 'colorectal cancer' (omit to fetch drug targets only)")
    p.add_argument("--out", default=None,
                   help="Output directory (default: CLI_RESULTS_DIR from API.env, else ./results)")


def cmd_targets(args) -> int:
    """并发请求 ChEMBL/DGIdb/Open Targets/DisGeNET,降噪整合去重后落盘靶点明细表。
    失败数据源按空表处理并打警告,不中断(DisGeNET 报空先检查 API.env 密钥)。"""
    out_dir = _resolve_out_dir(args.out)
    if out_dir is None:
        return 1
    print(f"[targets] output directory: {out_dir}")

    t0 = time.monotonic()
    shown: list[str] = []

    def _progress(done: list[str], total: int) -> None:
        # core 每 0.5 秒回调一次,只在有新完成源时打印,避免刷屏
        if len(done) != len(shown):
            shown[:] = done
            print(f"[targets] sources done {len(done)}/{total}: {', '.join(done)}")

    if not args.disease:
        print("[targets] no --disease given: fetching drug targets only")
    drug_df, disease_df, failed, written = _targets_step(args.drug, args.disease, out_dir,
                                                         progress_cb=_progress)

    for src in failed:
        hint = ""
        if src == "disgenet" and not config.DisGeNET_API_KEY:
            hint = " (DisGeNET_API_KEY not set in API.env)"
        _warn(f"data source '{src}' failed and was treated as empty{hint}")

    elapsed = time.monotonic() - t0
    empty_sides = []
    if drug_df is None or drug_df.empty:
        empty_sides.append(f"drug '{args.drug}'")
    else:
        print(f"[targets] drug '{args.drug}': {len(drug_df)} targets")
    if args.disease:
        if disease_df is None or disease_df.empty:
            empty_sides.append(f"disease '{args.disease}'")
        else:
            print(f"[targets] disease '{args.disease}': {len(disease_df)} targets")
    for side in empty_sides:
        _warn(f"no targets found for {side}")

    print(f"[targets] done in {elapsed:.1f}s: {len(written)} file(s) written"
          + ("".join(f"\n  {p}" for p in written)))
    # 失败源意味着明细表不完整(文献对照不可用);明确请求的一侧查无靶点按失败计(防拼写错误静默通过);
    # 该写出的文件没写成(如被 Excel 锁定)同样按失败计,避免旧文件被当成新结果
    expected = 1 if not args.disease else 2
    return 1 if (failed or empty_sides or len(written) < expected) else 0


# ======== 子命令:ppi ========

def _setup_ppi(p: argparse.ArgumentParser) -> None:
    p.add_argument("--genes", required=True,
                   help="Gene list: path to a .txt (one gene per line) or .csv (gene/approvedSymbol column), or a comma-separated list")
    p.add_argument("--out", default=None,
                   help="Output directory (default: CLI_RESULTS_DIR from API.env, else ./results)")


def cmd_ppi(args) -> int:
    """STRING PPI 网络 → 网络指标 → MCC 排序选核心靶点,落盘边表/指标表/核心靶点表。
    基因数超 STRING_QUERY_LIMIT 时截断并警告;与网页手动基因模式共享同一份边表缓存。"""
    out_dir = _resolve_out_dir(args.out)
    if out_dir is None:
        return 1
    genes = _read_genes(args.genes)
    if not genes:
        _warn(f"no genes parsed from --genes '{args.genes}'")
        return 1
    print(f"[ppi] {len(genes)} genes, STRING confidence >= {config.PPI_CONFIDENCE}, "
          f"output directory: {out_dir}")

    t0 = time.monotonic()
    # label="manual" 与网页手动基因模式一致:同参数共享 StringPPI_manual_* 边表缓存
    res = _ppi_step(genes, out_dir, label="manual")
    if res["truncated"]:
        _warn(f"{len(genes)} genes exceed STRING query limit "
              f"{config.STRING_QUERY_LIMIT}; truncated")
    if res["edges"] is not None and len(res["metrics"]):
        print(f"[ppi] hub rule: {res['metrics']['hub_rule'].iloc[0]} "
              f"(method={config.HUB_METHOD}, top {config.HUB_TOP_N})")
    if res["core"]:
        print(f"[ppi] core targets (top {len(res['core'])} by MCC): "
              f"{', '.join(res['core'])}")
    else:
        _warn(res["detail"])
    written = res["written"]
    print(f"[ppi] done in {time.monotonic() - t0:.1f}s: {len(written)} file(s) written"
          + ("".join(f"\n  {p}" for p in written)))
    return 0 if res["ok"] else 1


# ======== 子命令:enrich ========

def _setup_enrich(p: argparse.ArgumentParser) -> None:
    p.add_argument("--genes", required=True,
                   help="Gene list: path to a .txt (one gene per line) or .csv (gene/approvedSymbol column), or a comma-separated list")
    p.add_argument("--pvalue", type=_pvalue, default=config.ENRICHMENT_PVALUE,
                   help=f"Significance threshold in (0, 1] (default: {config.ENRICHMENT_PVALUE})")
    p.add_argument("--sources", default=",".join(config.ENRICHMENT_SOURCES),
                   help=f"Comma-separated term sources (default: {','.join(config.ENRICHMENT_SOURCES)})")
    p.add_argument("--organism", default=config.SPECIES_MAP["Human"]["gprofiler"],
                   help=f"Organism code for g:Profiler (default: {config.SPECIES_MAP['Human']['gprofiler']})")
    p.add_argument("--out", default=None,
                   help="Output directory (default: CLI_RESULTS_DIR from API.env, else ./results)")


def cmd_enrich(args) -> int:
    """GO/KEGG 富集:先 g:Profiler,失败自动切 Enrichr,出富集条目表 + 终端状态行。
    状态语义必须区分:no_significant=正常;both_sources_down=服务挂,结果不可信需重跑。"""
    sources = _parse_sources(args.sources)
    if sources is None:
        return 1
    genes = _read_genes(args.genes)
    if not genes:
        _warn(f"no genes parsed from --genes '{args.genes}'")
        return 1
    out_dir = _resolve_out_dir(args.out)
    if out_dir is None:
        return 1
    print(f"[enrich] {len(genes)} genes, p<={args.pvalue}, sources={','.join(sources)}, "
          f"organism={args.organism}, output directory: {out_dir}")

    t0 = time.monotonic()
    res = _enrichment_step(genes, out_dir, label=None, organism=args.organism,
                           pvalue=args.pvalue, sources=sources)
    tag = "ok" if res["ok"] else "FAILED"
    print(f"[enrich] {res['detail']}")
    written = [res["path"]] if res.get("path") else []
    print(f"[enrich] done in {time.monotonic() - t0:.1f}s [{tag}]: "
          f"{len(written)} file(s) written" + ("".join(f"\n  {p}" for p in written)))
    # both_sources_down 不是普通失败:结果不可用于文献对照,必须按非零退出暴露给 CI/脚本
    return 0 if res["ok"] else 1


# ======== 子命令:proximity ========

def _setup_proximity(p: argparse.ArgumentParser) -> None:
    p.add_argument("--drug", required=True, help="Drug name, e.g. aspirin")
    p.add_argument("--disease", required=True,
                   help="Disease name in English, e.g. 'colorectal cancer'")
    p.add_argument("--out", default=None,
                   help="Output directory (default: CLI_RESULTS_DIR from API.env, else ./results)")


def cmd_proximity(args) -> int:
    """药物靶点模块与疾病模块在背景互作组上的邻近度显著性,落盘 proximity.json。
    背景互作组未构建时警告并跳过,CLI 绝不自动触发数 GB 下载。"""
    out_dir = _resolve_out_dir(args.out)
    if out_dir is None:
        return 1
    print(f"[proximity] drug '{args.drug}' vs disease '{args.disease}', "
          f"output directory: {out_dir}")
    if not interactome.is_built():
        _warn("background interactome not built; build it in the web UI or via "
              "scripts/build_interactome.py, then rerun")
        return 1
    t0 = time.monotonic()
    drug_df = get_drug_targets_merged(args.drug)
    disease_df = get_disease_targets_merged(args.disease)
    if drug_df is None or drug_df.empty or disease_df is None or disease_df.empty:
        _warn("no targets on at least one side; cannot compute proximity")
        return 1
    drug_syms = sorted({str(g).strip().upper() for g in drug_df["approvedSymbol"].dropna()})
    disease_syms = sorted({str(g).strip().upper() for g in disease_df["approvedSymbol"].dropna()})
    res = _proximity_step(drug_syms, disease_syms, out_dir,
                          label=pair_prefix(args.drug, args.disease))
    print(f"[proximity] {res['detail']}")
    if not res["ok"]:
        _warn(res["detail"])
        return 1
    written = [res["path"]] if res.get("path") else []
    print(f"[proximity] done in {time.monotonic() - t0:.1f}s: "
          f"{len(written)} file(s) written" + ("".join(f"\n  {p}" for p in written)))
    return 0


# ======== 子命令:rwr ========

def _setup_rwr(p: argparse.ArgumentParser) -> None:
    p.add_argument("--seeds", required=True,
                   help="Seed genes (usually core targets): path to a .txt/.csv, or a comma-separated list")
    p.add_argument("--out", default=None,
                   help="Output directory (default: CLI_RESULTS_DIR from API.env, else ./results)")


def cmd_rwr(args) -> int:
    """以种子基因在背景网络上做随机游走重启,输出全网络基因排序(gene,score,is_seed,得分降序即优先级)。
    同样依赖背景互作组,未构建时行为与 proximity 相同:警告并跳过。"""
    out_dir = _resolve_out_dir(args.out)
    if out_dir is None:
        return 1
    seeds = _read_genes(args.seeds)
    if not seeds:
        _warn(f"no seeds parsed from --seeds '{args.seeds}'")
        return 1
    print(f"[rwr] {len(seeds)} seed genes, alpha={config.RWR_ALPHA}, "
          f"output directory: {out_dir}")
    if not interactome.is_built():
        _warn("background interactome not built; build it in the web UI or via "
              "scripts/build_interactome.py, then rerun")
        return 1

    t0 = time.monotonic()
    result = run_rwr_cached(seeds, label=None)  # 裸查询:缓存文件名用基因集短哈希,互不覆盖
    if not result.get("ok"):
        _warn(f"RWR failed: {result.get('error') or 'unknown error'}")
        return 1
    path = _save_rwr_csv(result, os.path.join(out_dir, "rwr.csv"))
    cov = result.get("coverage") or {}
    top = [r["gene"] for r in (result.get("top_nonseed") or [])[:10]]
    print(f"[rwr] converged={result.get('converged')} in {result.get('n_iter')} iters, "
          f"seed coverage {cov.get('seed_in_lcc', '?')}/{cov.get('seed_total', '?')} in LCC")
    if top:
        print(f"[rwr] top non-seed genes: {', '.join(top)}")
    print(f"[rwr] done in {time.monotonic() - t0:.1f}s: "
          f"{1 if path else 0} file(s) written" + (f"\n  {path}" if path else ""))
    return 0 if path else 1


# ======== 子命令:analyze ========

def _setup_analyze(p: argparse.ArgumentParser) -> None:
    p.add_argument("--drug", required=True, help="Drug name, e.g. aspirin")
    p.add_argument("--disease", required=True,
                   help="Disease name in English, e.g. 'colorectal cancer'")
    p.add_argument("--out", default=None,
                   help="Output directory (default: CLI_RESULTS_DIR from API.env, else ./results); the pair gets a subdirectory named after it")


def cmd_analyze(args) -> int:
    """单对全流程:靶点 → 交集 → PPI → 富集 → 邻近度 → RWR,落一个完整结果目录 + summary.json。
    某一步失败不阻断后续步骤;summary.json 如实记录每步状态、耗时与当次生效的全部参数。"""
    out_dir = _resolve_out_dir(args.out)
    if out_dir is None:
        return 1
    pair_dir = os.path.join(out_dir, _slug(args.drug, args.disease))
    os.makedirs(pair_dir, exist_ok=True)
    print(f"[analyze] drug '{args.drug}' vs disease '{args.disease}'")
    print(f"[analyze] pair directory: {pair_dir}")

    t0 = time.monotonic()
    summary = run_pair(args.drug, args.disease, pair_dir)
    written = sorted(f for f in os.listdir(pair_dir) if f != "summary.json")
    failed = [n for n, s in summary["steps"].items() if s.get("status") != "ok"]
    print(f"[analyze] done in {time.monotonic() - t0:.1f}s: "
          f"{len(written)} file(s) + summary.json"
          + ("".join(f"\n  {os.path.join(pair_dir, f)}" for f in written) if written else ""))
    if failed:
        _warn(f"steps not ok: {', '.join(failed)} (see summary.json)")
    else:
        print(f"[analyze] all {len(summary['steps'])} steps ok")
    return 0 if summary["ok"] else 1


# ======== 子命令:batch ========

def _setup_batch(p: argparse.ArgumentParser) -> None:
    p.add_argument("--pairs", required=True,
                   help="Pairs CSV with 'drug' and 'disease' columns (case-insensitive, any position; extra columns like row numbers are ignored; utf-8-sig, Excel-friendly)")
    p.add_argument("--out", default=None,
                   help="Output directory (default: CLI_RESULTS_DIR from API.env, else ./results); one subdirectory per pair")
    p.add_argument("--resume", action="store_true",
                   help="Skip pairs whose subdirectory already holds a successful summary.json (checkpoint resume; failed pairs are retried)")


def _read_pairs(path: str) -> tuple[list[tuple[str, str]], list[tuple[int, str]]]:
    """读 pairs.csv:表头须含 drug、disease 两列(大小写不敏感、列位置不限;Excel 常见的
    No. 序号列等多余列自动忽略);utf-8-sig 读(兼容 Excel 导出)。
    行 = 一对,一对一配对;同药多病写多行,不做笛卡尔积。
    返回 (合法对列表, 非法行列表[(行号, 原因)]):缺列抛 ValueError;drug/disease 有空缺的行记非法。"""
    df = pd.read_csv(path, encoding="utf-8-sig", dtype=str)
    by_lower = {str(c).strip().lower(): c for c in df.columns}
    missing = [name for name in ("drug", "disease") if name not in by_lower]
    if missing:
        raise ValueError(f"pairs CSV must have 'drug' and 'disease' columns; "
                         f"got: {', '.join(str(c) for c in df.columns)}")
    drug_col, disease_col = by_lower["drug"], by_lower["disease"]
    pairs: list[tuple[str, str]] = []
    bad: list[tuple[int, str]] = []
    for i, row in df.iterrows():
        drug_raw, disease_raw = row[drug_col], row[disease_col]
        drug = str(drug_raw).strip() if pd.notna(drug_raw) else ""
        disease = str(disease_raw).strip() if pd.notna(disease_raw) else ""
        if drug and disease:
            pairs.append((drug, disease))
        else:
            bad.append((i + 2, "empty drug or disease"))  # +2:文件首行为表头,行号从 1 计
    return pairs, bad


_BATCH_SUMMARY_COLUMNS = ("drug", "disease", "status", "resumed",
                          "failed_steps", "detail", "seconds", "directory")


def cmd_batch(args) -> int:
    """批量跑批(文献对照主力):逐行读药病对循环调用 run_pair,每对一个子目录,
    最后出 batch_summary.csv(每对一行:状态/耗时/失败原因)。单对失败不中断;
    --resume 跳过已成功的对(断点续跑);结束打印 ok=N failed=M,failed>0 退出码 1。"""
    out_dir = _resolve_out_dir(args.out)
    if out_dir is None:
        return 1
    try:
        pairs, bad_rows = _read_pairs(args.pairs)
    except (OSError, ValueError) as e:
        _warn(f"failed to read pairs CSV '{args.pairs}': {e}")
        return 1
    if not pairs and not bad_rows:
        _warn("pairs CSV has no data rows")
        return 1
    print(f"[batch] {len(pairs)} pair(s)"
          + (f", {len(bad_rows)} invalid row(s)" if bad_rows else "")
          + f" from {args.pairs}\n[batch] output directory: {out_dir}"
          + (" (resume on)" if args.resume else ""))

    t0 = time.monotonic()
    rows: list[dict] = []
    for n, (drug, disease) in enumerate(pairs, 1):
        slug = _slug(drug, disease)
        pair_dir = os.path.join(out_dir, slug)
        print(f"\n[batch] ({n}/{len(pairs)}) {drug} vs {disease}")
        record = {"drug": drug, "disease": disease, "status": "failed", "resumed": 0,
                  "failed_steps": "", "detail": "", "seconds": "", "directory": slug}
        rows.append(record)
        summary_path = os.path.join(pair_dir, "summary.json")
        prev = load_json_file(summary_path, None) if args.resume else None
        if prev is not None and prev.get("ok"):
            # 断点续跑:已成功(且落盘 summary.json)的对直接跳过,状态沿用落盘记录
            print("[batch] skipped: successful summary.json found")
            record.update(status="ok", resumed=1, detail="resumed from existing summary",
                          seconds=prev.get("elapsed_seconds", ""))
            continue
        if args.resume and os.path.exists(summary_path):
            print("[batch] previous run not ok: retrying")
        try:
            summary = run_pair(drug, disease, pair_dir)
            failed_steps = [name for name, st in summary["steps"].items()
                            if st.get("status") != "ok"]
            detail = ("; ".join(f"{name}: {summary['steps'][name].get('detail', '')}"
                                for name in failed_steps)
                      if failed_steps else "all steps ok")
            record.update(status="ok" if summary["ok"] else "failed",
                          failed_steps=", ".join(failed_steps), detail=detail,
                          seconds=summary["elapsed_seconds"])
        except KeyboardInterrupt:
            raise  # 主动中断:整批停止;已完成的对已各自落盘,--resume 可续跑
        except Exception as e:
            logging.exception("pair '%s' vs '%s' crashed", drug, disease)
            record.update(failed_steps="crash", detail=f"{type(e).__name__}: {e}")

    # 结构非法的行(空 drug/disease)也按失败入汇总,保证"每行一条记录"对得上账
    for lineno, reason in bad_rows:
        rows.append({"drug": "", "disease": "", "status": "failed", "resumed": 0,
                     "failed_steps": "invalid_row", "detail": f"line {lineno}: {reason}",
                     "seconds": "", "directory": ""})

    summary_csv = os.path.join(out_dir, "batch_summary.csv")
    try:
        pd.DataFrame(rows, columns=_BATCH_SUMMARY_COLUMNS).to_csv(
            summary_csv, index=False, encoding="utf-8-sig")
    except PermissionError:
        _warn(f"cannot write '{summary_csv}': file is locked by another program "
              "(close it, e.g. Excel, and rerun this command to refresh the summary)")
        return 1

    ok_n = sum(1 for r in rows if r["status"] == "ok")
    failed_n = len(rows) - ok_n
    print(f"\n[batch] done in {time.monotonic() - t0:.1f}s: ok={ok_n} failed={failed_n}")
    for r in rows:
        if r["status"] != "ok":
            print(f"[batch] FAILED: {r['drug'] or '(row)'} | {r['disease'] or '(row)'}"
                  f" -> {r['detail']}"[:220])
    print(f"[batch] summary: {summary_csv}")
    return 1 if failed_n else 0


def _setup_serve(p) -> None:
    p.add_argument("--host", default="127.0.0.1", help="bind address (default 127.0.0.1)")
    p.add_argument("--port", type=int, default=8000, help="port (default 8000)")
    p.add_argument("--no-browser", action="store_true",
                   help="do not auto-open the browser after startup")


def cmd_serve(args) -> int:
    """启动 Web 服务:FastAPI(/api)伺服已构建的 TS 前端(WY/web/dist)与 API,
    可选自动打开浏览器。前端改代码后需在 WY/web 下 npm run build 重建 dist。"""
    import threading
    import webbrowser

    import uvicorn
    from WY.api import app

    url = f"http://{args.host}:{args.port}"
    # 端口预检:占用时给友好指引而不是 winerror 10048 裸异常
    # (双击 bat 场景最常见的错误 = 上一个实例还在运行,或别的程序占了端口)。
    import socket

    probe = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        probe.bind((args.host, args.port))
    except OSError:
        print(f"[serve] 端口 {args.port} 已被占用 —— 可能已有一个 NPAP 在运行,"
              f"浏览器打开 {url} 即是;若为其他程序,请关闭它或加 --port 换端口。")
        return 1
    finally:
        probe.close()
    print(f"[serve] {url}  (Ctrl+C 停止)")
    if not args.no_browser:
        threading.Timer(1.5, webbrowser.open, args=[url]).start()
    uvicorn.run(app, host=args.host, port=args.port, log_level="warning")
    return 0


# ======== 命令注册表:名字 → (参数装配, 处理函数) ========

HANDLER_HELP = {
    "serve": "Launch the local web server (FastAPI API + built TS frontend) and open the browser",
    "targets": "Fetch drug/disease targets from 4 databases (ChEMBL/DGIdb/OpenTargets/DisGeNET)",
    "ppi": "Build STRING PPI network, compute topology metrics and MCC core targets",
    "enrich": "GO/KEGG enrichment (g:Profiler with automatic Enrichr fallback)",
    "proximity": "Network proximity Z-score between drug and disease target modules",
    "rwr": "Random-walk-with-restart ranking of the whole network from seed genes",
    "analyze": "Full pipeline for one drug-disease pair (all steps + summary.json)",
    "batch": "Run the full pipeline for many drug-disease pairs from a CSV (batch summary + resume)",
}

COMMANDS = {
    "serve": (_setup_serve, cmd_serve),
    "targets": (_setup_targets, cmd_targets),
    "ppi": (_setup_ppi, cmd_ppi),
    "enrich": (_setup_enrich, cmd_enrich),
    "proximity": (_setup_proximity, cmd_proximity),
    "rwr": (_setup_rwr, cmd_rwr),
    "analyze": (_setup_analyze, cmd_analyze),
    "batch": (_setup_batch, cmd_batch),
}


def main(argv=None) -> int:
    setup_logging()  # core 日志同时进终端(stderr)与按天日志文件,进度可见
    parser = argparse.ArgumentParser(
        prog="npap",
        description="Network pharmacology pipeline CLI "
                    "(shares cache & settings with the web UI)")
    sub = parser.add_subparsers(dest="command", required=True)
    for name, (setup, handler) in COMMANDS.items():
        p = sub.add_parser(name, help=HANDLER_HELP[name])
        setup(p)
        p.set_defaults(func=handler)
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        return 130
    except Exception:
        logging.exception("command '%s' crashed", args.command)
        return 1


if __name__ == "__main__":
    sys.exit(main())
