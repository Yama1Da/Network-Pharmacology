"""主页页脚统计:本机数据规模实数(全部由缓存数出来,无装饰值)。"""

from fastapi import APIRouter
import csv
import os
import threading
import time
from WY import config


router = APIRouter()


# 视为"数据库"的来源文件前缀(其余前缀如 GSEA 是方法名,不计)
_DB_FILE_PREFIXES = ("ChEMBL", "DGIdb", "OpenTargets", "DisGeNET", "gProfiler", "Enrichr")
_stats_cache: dict = {"at": 0.0, "data": None}
_stats_lock = threading.Lock()


def _csv_rows(path: str) -> list[list[str]]:
    """读 CSV 数据行(跳过 # 注释行与表头行),BOM 兼容。"""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = [r for r in csv.reader(fh) if r and not r[0].startswith("#")]
    return rows[1:] if rows else []


def _gene_union(dirs: list[str]) -> set[str]:
    """各来源/整合表首列均为 approvedSymbol,并集即本机缓存靶点。"""
    out: set[str] = set()
    for d in dirs:
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if f.endswith(".csv"):
                out |= {r[0] for r in _csv_rows(os.path.join(d, f)) if r[0]}
    return out


def _compound_set() -> set[str]:
    """化合物 = 药物缓存文件名的查询词元 ∪ DGIdb 行内 drug_name 化合物条目。"""
    out: set[str] = set()
    if not os.path.isdir(config.DRUG_CACHE_DIR):
        return out
    for f in os.listdir(config.DRUG_CACHE_DIR):
        if not f.endswith(".csv"):
            continue
        prefix, _, term = f[:-4].partition("_")
        if prefix in _DB_FILE_PREFIXES and term:
            out.add(term)
        if prefix == "DGIdb":
            out |= {r[1] for r in _csv_rows(os.path.join(config.DRUG_CACHE_DIR, f))
                    if len(r) > 1 and r[1]}
    return out


def _pathway_union() -> set[str]:
    """通路 = 富集缓存各行去重:gProfiler/Enrichr 取 term_id,GSEA 只有 term_name。"""
    out: set[str] = set()
    if not os.path.isdir(config.ENRICHMENT_CACHE_DIR):
        return out
    for f in os.listdir(config.ENRICHMENT_CACHE_DIR):
        if not f.endswith(".csv"):
            continue
        with open(os.path.join(config.ENRICHMENT_CACHE_DIR, f),
                  newline="", encoding="utf-8-sig") as fh:
            rows = [r for r in csv.reader(fh) if r and not r[0].startswith("#")]
        if len(rows) < 2:
            continue
        header = rows[0]
        for key in ("term_id", "term_name"):
            if key in header:
                i = header.index(key)
                out |= {r[i] for r in rows[1:] if len(r) > i and r[i]}
                break
    return out


def _dir_has(path: str, suffix: str) -> bool:
    return os.path.isdir(path) and any(f.endswith(suffix) for f in os.listdir(path))


def _compute_stats() -> dict:
    names: set[str] = set()
    for d in (config.DRUG_CACHE_DIR, config.DISEASE_CACHE_DIR, config.ENRICHMENT_CACHE_DIR):
        if os.path.isdir(d):
            names |= {f.split("_", 1)[0] for f in os.listdir(d)
                      if f.endswith(".csv")} & set(_DB_FILE_PREFIXES)
    if _dir_has(config.INTERACTOME_CACHE_DIR, ".parquet"):
        names.add("STRING")
    if os.path.isdir(os.path.join(config.CACHE_DIR, config.GEO_SERIES_CACHE_SUBDIR)) \
            and os.listdir(os.path.join(config.CACHE_DIR, config.GEO_SERIES_CACHE_SUBDIR)):
        names.add("GEO")
    if _dir_has(config.GSEA_GMT_DIR, ".gmt"):
        names.add("GMT")  # 基因集库(KEGG/MSigDB GMT)
    return {
        "databases": len(names),
        "compounds": len(_compound_set()),
        "targets": len(_gene_union([config.DRUG_CACHE_DIR, config.DISEASE_CACHE_DIR,
                                    config.MERGED_CACHE_DIR])),
        "pathways": len(_pathway_union()),
    }


@router.get("/api/stats")
def home_stats() -> dict:
    """主页页脚统计:数据库(有本机缓存数据的外部库)/化合物/靶点/通路,60s 进程内缓存。

    过期重算走 single-flight(2026-11-10 审计):并发请求只放一个线程进锁重算,
    其余等锁后直接拿到新值——避免过期瞬间多线程重复全目录扫描(重算要遍历
    drug/disease/merged/enrichment 四个目录的 CSV)。双检锁与 http_client.get_session
    同款模式;锁外的快路径读 dict 字段无锁(GIL 下字段赋值原子,最坏读到刚过期的旧值)。"""
    if _stats_cache["data"] is None or time.monotonic() - _stats_cache["at"] > 60:
        with _stats_lock:
            if _stats_cache["data"] is None or time.monotonic() - _stats_cache["at"] > 60:
                _stats_cache["data"] = _compute_stats()
                _stats_cache["at"] = time.monotonic()
    return _stats_cache["data"]
