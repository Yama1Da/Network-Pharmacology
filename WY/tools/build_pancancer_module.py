#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""构建泛癌对照基因集（二期 1.9）：Open Targets 癌症靶点 ∪ CancerMine 癌症驱动基因。

背景（清单 §1.9 / 方案 §3.3 B5 定稿）：DrugComb 药对评测里，细胞系无法映射到具体
癌种的组合**不剔除**，改挂"泛癌对照模块"分层报告；该对照必须是一个"泛癌"基因集，
而不是全基因组——全基因组会把疾病模块信号稀释成噪声，对照就失去意义。

数据源（均为开放许可，且都不含 COSMIC CGC / MSigDB）：
  1. Open Targets Platform GraphQL（CC BY 4.0）：cancer 的关联靶点。
     疾病编号用 MONDO_0004992——OT 平台已迁 MONDO 主键，老 EFO 编号（如
     EFO_0000311）查询返回 HTTP 200 + disease:null 而非报错（2026-11-18 实证），
     网络重试层拦不住，会静默拿到空臂；
  2. CancerMine（CC0）：文献挖掘的 driver / oncogene / tumor_suppressor。
原定 NCG 因官网停服（2026-11-18 实测维护中；最新为 2022 年 NCGHD 7.0、无 8.0）
改用 CancerMine。

**CancerMine 取数策略（2026-11-18 定案：不依赖 Zenodo 可达性）**
Zenodo 在国内 DNS 被污染，在线下载经常失败；因此主路径是**手动下载 + 本地快照**：
    1) 在能访问 Zenodo 的环境（VPN）下载 cancermine_collated.tsv；
    2) 跑一次带 --cancermine-file（指向下载的 tsv）与 --cancermine-update 的构建
       —— 脚本过滤/降噪后写 WY/resources/cancermine_cancer_genes.csv（仓库内
       快照，带来源与日期）；之后一律读快照，--force 重建也不再联网，只有
       --refresh 才会重新下载。
顺序：--cancermine-file ＞ 快照（--refresh 时跳过）＞ raw 缓存（--offline）＞ 在线。
各参数语义见 argparse --help。

产物（版本化，落 config.CACHE_DIR/pancancer/）：
    pancancer_genes_v<YYYYMMDD>.parquet   列：gene / sources
    meta.json                             来源、计数、疾病编号、快照版本、拉取时间
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
from datetime import date as _date
from datetime import datetime, timezone

import pandas as pd

from WY import config
from WY.core.http_client import request_with_retry

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
# OT 已迁 MONDO 主键（2026-11-18 实证）：老 EFO_0000311 查询返回 200+disease:null，
# 与四库检索 find_disease_id 优先 MONDO 的取向一致，勿改回 EFO 编号。
OT_DISEASE_ID = "MONDO_0004992"  # cancer
ZENODO_CONCEPT_DOI = "10.5281/zenodo.1156241"
ZENODO_RECORDS_API = "https://zenodo.org/api/records"
CANCERMINE_FALLBACK_URL = (
    "https://zenodo.org/records/1156241/files/cancermine_collated.tsv?download=1")
CANCERMINE_FILE = "cancermine_collated.tsv"
CANCERMINE_SNAPSHOT = "cancermine_cancer_genes.csv"
OT_PAGE_SIZE = 500
ROLES = {"driver", "oncogene", "tumor_suppressor", "tumour_suppressor"}

DISEASE_TARGETS_QUERY = """
query ($diseaseId: String!, $index: Int!, $size: Int!) {
  disease(efoId: $diseaseId) {
    name
    associatedTargets(page: {index: $index, size: $size}) {
      count
      rows { target { approvedSymbol } }
    }
  }
}
"""


class PancancerError(RuntimeError):
    """构建失败（缺数据/上游不可用/验收不过）。"""


def data_dir() -> str:
    return os.path.join(config.CACHE_DIR, "pancancer")


def raw_dir() -> str:
    return os.path.join(data_dir(), "raw")


def resources_dir() -> str:
    """仓库源数据目录 WY/resources（与运行时缓存 Static 区分）。"""
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "resources")


def default_snapshot_path() -> str:
    return os.path.join(resources_dir(), CANCERMINE_SNAPSHOT)


# ---------------------------------------------------------------------------
# Open Targets
# ---------------------------------------------------------------------------
def _ot_post(query: str, variables: dict) -> dict:
    resp = request_with_retry(
        OT_GRAPHQL, method="post", json={"query": query, "variables": variables},
        max_retries=config.LIGHT_API_MAX_RETRIES, timeout=config.API_TIMEOUT,
        label="OT pancancer targets", throttle_key="opentargets",
    )
    if resp is None:
        raise PancancerError("Open Targets 请求失败（重试耗尽）")
    payload = resp.json()
    if payload.get("errors"):
        raise PancancerError(f"Open Targets 返回错误：{payload['errors'][:1]}")
    return payload.get("data") or {}


def fetch_ot_targets(disease_id: str = OT_DISEASE_ID) -> tuple[list[str], str]:
    """分页拉取该疾病关联靶点，返回 (approvedSymbol 列表, 疾病名)。"""
    symbols: list[str] = []
    disease_name = ""
    index = 0
    total = None
    while True:
        data = _ot_post(DISEASE_TARGETS_QUERY,
                        {"diseaseId": disease_id, "index": index, "size": OT_PAGE_SIZE})
        disease = data.get("disease") or {}
        disease_name = str(disease.get("name") or disease_name)
        block = disease.get("associatedTargets") or {}
        if total is None:
            total = int(block.get("count") or 0)
        rows = block.get("rows") or []
        for row in rows:
            symbol = str(((row.get("target") or {}).get("approvedSymbol")) or "").strip()
            if symbol:
                symbols.append(symbol)
        if not rows or len(symbols) >= (total or 0):
            break
        index += 1
        if index > 200:  # 防御：分页异常时不死循环
            break
    return sorted(set(symbols)), disease_name


# ---------------------------------------------------------------------------
# CancerMine：TSV 解析 / 快照读写 / 取数
# ---------------------------------------------------------------------------
def resolve_cancermine_url() -> str:
    """按概念 DOI 解析 CancerMine 最新版本的 collated TSV 直链（需能访问 Zenodo）。"""
    resp = request_with_retry(
        ZENODO_RECORDS_API, method="get",
        params={"q": f'conceptdoi:"{ZENODO_CONCEPT_DOI}"',
                "sort": "mostrecent", "size": 1},
        max_retries=config.LIGHT_API_MAX_RETRIES, timeout=config.API_TIMEOUT,
        label="Zenodo CancerMine record",
    )
    if resp is None:
        return CANCERMINE_FALLBACK_URL
    try:
        hits = (resp.json() or {}).get("hits", {}).get("hits") or []
    except ValueError:
        return CANCERMINE_FALLBACK_URL
    for hit in hits:
        for item in hit.get("files") or []:
            if str(item.get("key")) == CANCERMINE_FILE:
                link = ((item.get("links") or {}).get("self")
                        or (item.get("links") or {}).get("download"))
                if link:
                    return str(link)
    return CANCERMINE_FALLBACK_URL


def download_cancermine(url: str | None = None) -> tuple[str, str]:
    """下载 CancerMine TSV 原文，返回 (文本, 实际 URL)。"""
    target = url or resolve_cancermine_url()
    resp = request_with_retry(
        target, method="get", max_retries=config.LIGHT_API_MAX_RETRIES,
        timeout=max(config.API_TIMEOUT, 60), label="CancerMine TSV",
    )
    if resp is None:
        raise PancancerError(
            f"CancerMine 下载失败：{target}\n"
            "Zenodo 在国内常被 DNS 污染。请在可访问 Zenodo 的环境（VPN）手动下载 "
            f"{CANCERMINE_FILE}，然后运行：\n"
            "  python -m WY.tools.build_pancancer_module "
            "--cancermine-file <下载的 tsv 路径> --cancermine-update")
    resp.encoding = resp.encoding or "utf-8"
    return resp.text, target


def parse_cancermine_rows(text: str, *, min_citations: int = 2) -> list[dict]:
    """解析 collated TSV：保留 driver/oncogene/TSG，按引用数过滤，返回行。"""
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not reader.fieldnames:
        raise PancancerError("CancerMine TSV 无表头")
    gene_col = next((c for c in ("gene_normalized", "gene", "Gene")
                     if c in reader.fieldnames), None)
    role_col = next((c for c in ("role", "Role") if c in reader.fieldnames), None)
    cite_col = next((c for c in ("citation_count", "citationCount")
                     if c in reader.fieldnames), None)
    if not gene_col or not role_col:
        raise PancancerError(f"CancerMine TSV 缺列（现有：{reader.fieldnames}）")
    rows: dict[str, dict] = {}
    for row in reader:
        raw_roles = str(row.get(role_col) or "").replace(",", ";")
        roles = {part.strip().lower() for part in raw_roles.split(";") if part.strip()}
        if not (roles & ROLES):
            continue
        citations = 0
        if cite_col:
            try:
                citations = int(float(row.get(cite_col) or 0))
            except ValueError:
                citations = 0
            if citations < int(min_citations):
                continue
        symbol = str(row.get(gene_col) or "").strip()
        if not symbol:
            continue
        item = rows.setdefault(symbol, {"gene": symbol, "roles": set(), "citations": 0})
        item["roles"] |= roles & ROLES
        item["citations"] = max(item["citations"], citations)
    return [{"gene": item["gene"],
             "roles": ";".join(sorted(item["roles"])),
             "citations": item["citations"]}
            for item in sorted(rows.values(), key=lambda x: x["gene"])]


def parse_cancermine(text: str, *, min_citations: int = 2) -> list[str]:
    """兼容接口：只返回基因符号列表。"""
    return [row["gene"] for row in parse_cancermine_rows(text, min_citations=min_citations)]


def write_snapshot(path: str, rows: list[dict], *, source: str,
                   min_citations: int) -> None:
    """把过滤后的 CancerMine 基因表写成仓库快照（`#` 注释头 + CSV）。"""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        f.write("# cancermine_cancer_genes.csv —— CancerMine 泛癌基因快照（CC0）\n")
        f.write(f"# snapshot_date: {_date.today().isoformat()}\n")
        f.write(f"# source: {source}\n")
        f.write(f"# filter: role in driver/oncogene/tumor_suppressor; "
                f"citation_count >= {min_citations}\n")
        f.write("# 由 WY/tools/build_pancancer_module.py 生成，勿手改；"
                "更新请重跑 --cancermine-update\n")
        writer = csv.DictWriter(f, fieldnames=["gene", "roles", "citations"])
        writer.writeheader()
        for row in rows:
            writer.writerow({"gene": row["gene"], "roles": row["roles"],
                             "citations": row["citations"]})
    os.replace(tmp, path)


def read_snapshot(path: str) -> list[dict]:
    """读快照；文件缺失返回空列表，格式坏抛 PancancerError。"""
    if not os.path.isfile(path):
        return []
    with open(path, encoding="utf-8") as f:
        lines = [line for line in f if not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    if not reader.fieldnames or "gene" not in reader.fieldnames:
        raise PancancerError(f"快照格式损坏（缺 gene 列）：{path}")
    rows = []
    for row in reader:
        gene = str(row.get("gene") or "").strip()
        if gene:
            rows.append({"gene": gene, "roles": str(row.get("roles") or ""),
                         "citations": int(float(row.get("citations") or 0))})
    return rows


def _cancermine_rows(*, force: bool, refresh: bool, offline: bool,
                     cancermine_file: str | None, cancermine_url: str | None,
                     min_citations: int, snapshot_path: str,
                     update_snapshot: bool) -> tuple[list[dict], str]:
    """按优先级取 CancerMine 行：本地文件 ＞ 快照 ＞ raw 缓存 ＞ 在线下载。"""
    source = ""
    rows: list[dict] = []

    if cancermine_file:
        if not os.path.isfile(cancermine_file):
            raise PancancerError(
                f"--cancermine-file 不存在：{cancermine_file}"
                f"（按当前目录解析为 {os.path.abspath(cancermine_file)}；"
                f"仓库源数据目录在 {resources_dir()}，请传实际路径或先 cd 到文件所在目录）")
        with open(cancermine_file, encoding="utf-8") as f:
            text = f.read()
        rows = parse_cancermine_rows(text, min_citations=min_citations)
        source = f"file:{os.path.basename(cancermine_file)}"
        if update_snapshot:
            write_snapshot(snapshot_path, rows, source=source,
                           min_citations=min_citations)
            source += f"|snapshot:{os.path.basename(snapshot_path)}"
        return rows, source

    snapshot_rows = read_snapshot(snapshot_path)
    if snapshot_rows and not refresh:
        return snapshot_rows, f"snapshot:{os.path.basename(snapshot_path)}"

    cm_cache = os.path.join(raw_dir(), CANCERMINE_FILE)
    if os.path.isfile(cm_cache) and not refresh and not force:
        with open(cm_cache, encoding="utf-8") as f:
            text = f.read()
        rows = parse_cancermine_rows(text, min_citations=min_citations)
        source = "cache"
        if update_snapshot:
            write_snapshot(snapshot_path, rows, source=source,
                           min_citations=min_citations)
        return rows, source

    if offline:
        raise PancancerError(
            "离线且没有 CancerMine 快照/缓存。请在 VPN 环境手动下载 "
            f"{CANCERMINE_FILE} 后运行 --cancermine-file <tsv> --cancermine-update")

    try:
        text, url = download_cancermine(cancermine_url)
    except PancancerError as exc:
        if snapshot_rows:
            raise PancancerError(f"{exc}\n（检测到旧快照，可去掉 --refresh 直接用）")
        raise
    rows = parse_cancermine_rows(text, min_citations=min_citations)
    source = url
    os.makedirs(raw_dir(), exist_ok=True)
    with open(cm_cache, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    if update_snapshot:
        write_snapshot(snapshot_path, rows, source=source, min_citations=min_citations)
        source += f"|snapshot:{os.path.basename(snapshot_path)}"
    return rows, source


# ---------------------------------------------------------------------------
# 构建
# ---------------------------------------------------------------------------
def _write_json(path: str, payload: dict) -> None:
    tmp = path + ".tmp"
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


def build(*, force: bool = False, offline: bool = False,
          disease_id: str = OT_DISEASE_ID, min_citations: int = 2,
          cancermine_url: str | None = None, cancermine_file: str | None = None,
          snapshot_path: str | None = None, update_snapshot: bool = False,
          refresh: bool = False, out: str | None = None) -> dict:
    snapshot = snapshot_path or default_snapshot_path()

    # ---- Open Targets ----
    ot_cache = os.path.join(raw_dir(), f"ot_{disease_id}.json")
    os.makedirs(raw_dir(), exist_ok=True)
    ot_symbols, ot_name = [], ""
    cache_hit = False
    if os.path.isfile(ot_cache) and not (force or refresh):
        with open(ot_cache, encoding="utf-8") as f:
            cached = json.load(f)
        ot_symbols = list(cached.get("symbols") or [])
        ot_name = str(cached.get("disease_name") or "")
        cache_hit = bool(ot_symbols)  # 空缓存视作未命中（历史坏缓存），走重抓
    if not cache_hit:
        if offline:
            raise PancancerError(f"离线且无有效 OT 缓存：{ot_cache}")
        ot_symbols, ot_name = fetch_ot_targets(disease_id)
        if not ot_symbols:
            # 空结果当场抛错、不写缓存：OT 迁 MONDO 后老编号返回 200+null（无
            # errors），重试层看是成功——2026-11-18 实证空缓存落盘会让对照静默缺臂
            raise PancancerError(
                f"Open Targets 疾病 {disease_id} 一条靶点都没拿到（HTTP 层成功但"
                "疾病查空，多为疾病编号失效/上游接口改动）——不写缓存，核对编号后重试")
        _write_json(ot_cache, {"symbols": ot_symbols, "disease_name": ot_name,
                               "fetched_at": datetime.now(timezone.utc).isoformat()})

    # ---- CancerMine ----
    cm_rows, cm_source = _cancermine_rows(
        force=force, refresh=refresh, offline=offline,
        cancermine_file=cancermine_file, cancermine_url=cancermine_url,
        min_citations=min_citations, snapshot_path=snapshot,
        update_snapshot=update_snapshot)
    cm_symbols = [row["gene"] for row in cm_rows]

    # ---- 并集 ----
    def norm(values: list[str]) -> list[str]:
        return sorted({value.strip().upper() for value in values
                       if value and value.strip()})

    ot_set, cm_set = set(norm(ot_symbols)), set(norm(cm_symbols))
    union = sorted(ot_set | cm_set)
    if not union:
        raise PancancerError("并集为空：两个来源都没拿到基因")
    # 验收门 1（2026-11-18 增）：两源各至少命中一条——单臂缺失时并集数量可能
    # 依然可观（CancerMine 单臂 ~4000），不拦住则对照悄悄降级（OT 空缓存事故实证）
    empty_arms = [name for name, subset in
                  (("Open Targets", ot_set), ("CancerMine", cm_set)) if not subset]
    if empty_arms:
        raise PancancerError(f"验收失败：{'/'.join(empty_arms)} 一条基因都没拿到，泛癌对照缺臂")
    # 验收门 2：关键癌基因必须在集内（防两源口径漂移）
    missing = [g for g in ("BRAF", "TP53", "PIK3CA") if g not in set(union)]
    if missing:
        raise PancancerError(f"验收失败：关键基因缺失 {missing}")

    rows = [{"gene": gene,
             "sources": ";".join(s for s, subset in
                                 (("OT", ot_set), ("CANCERMINE", cm_set))
                                 if gene in subset)}
            for gene in union]
    df = pd.DataFrame(rows, columns=["gene", "sources"])

    out_dir = out or data_dir()
    os.makedirs(out_dir, exist_ok=True)
    version = _date.today().strftime("%Y%m%d")
    parquet_path = os.path.join(out_dir, f"pancancer_genes_v{version}.parquet")
    tmp = parquet_path + ".tmp"
    df.to_parquet(tmp, index=False)
    os.replace(tmp, parquet_path)

    meta = {
        "module": "pancancer_genes",
        "version": version,
        "built_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "ot_disease_id": disease_id,
        "ot_disease": ot_name,
        "ot_targets": len(ot_set),
        "cancermine_source": cm_source,
        "cancermine_min_citations": min_citations,
        "cancermine_genes": len(cm_set),
        "union_genes": len(union),
        "overlap": len(ot_set & cm_set),
        "markers_ok": True,
        "parquet": os.path.basename(parquet_path),
        "note": "NCG 官网 2026-11 停服，改用 CancerMine（CC0）；"
                "CancerMine 走本地快照，不依赖 Zenodo 可达性。",
    }
    _write_json(os.path.join(out_dir, "meta.json"), meta)
    return meta


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="构建泛癌对照基因集（二期 1.9）")
    parser.add_argument("--force", action="store_true",
                        help="忽略 OT raw 缓存/旧产物重建（仍优先本地快照）")
    parser.add_argument("--refresh", action="store_true",
                        help="VPN 下重新下载 CancerMine 并更新快照")
    parser.add_argument("--offline", action="store_true", help="只用本地数据")
    parser.add_argument("--ot-disease-id", default=OT_DISEASE_ID,
                        help=f"Open Targets 疾病编号（MONDO 主键，默认 {OT_DISEASE_ID} cancer；"
                             "老 EFO 编号已查空勿用）")
    parser.add_argument("--min-citations", type=int, default=2,
                        help="CancerMine 引用数下限（降噪，默认 2）")
    parser.add_argument("--cancermine-file", default=None,
                        help="手动下载的 cancermine_collated.tsv（优先于快照/网络）")
    parser.add_argument("--cancermine-update", action="store_true",
                        help="把当前 CancerMine 数据写入 WY/resources 快照")
    parser.add_argument("--cancermine-url", default=None,
                        help="覆盖 CancerMine TSV 直链（需能访问 Zenodo）")
    parser.add_argument("--snapshot-path", default=None,
                        help="快照路径（默认 WY/resources/cancermine_cancer_genes.csv）")
    parser.add_argument("--out", default=None, help="输出目录（默认 CACHE_DIR/pancancer）")
    args = parser.parse_args(argv)
    try:
        meta = build(force=args.force, offline=args.offline, disease_id=args.ot_disease_id,
                     min_citations=args.min_citations,
                     cancermine_url=args.cancermine_url,
                     cancermine_file=args.cancermine_file,
                     snapshot_path=args.snapshot_path,
                     update_snapshot=args.cancermine_update,
                     refresh=args.refresh, out=args.out)
    except PancancerError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 2
    print("[done] %s/%s" % (meta["version"], meta["parquet"]))
    print("[done] OT=%d CancerMine=%d overlap=%d union=%d source=%s"
          % (meta["ot_targets"], meta["cancermine_genes"],
             meta["overlap"], meta["union_genes"], meta["cancermine_source"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
