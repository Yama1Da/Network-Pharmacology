#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""DailyMed 说明书机制靶点精选工具（二期 1.10；独立单文件）。

三段式流水线：
    prefill   三级预填（OT MOA → FDA 药理类别 → DailyMed MoA 原文），生成审查队列
    review    终端逐条审查（y 确认 / t 待定 / x 否决），每键立即写盘、可断点续传
    validate  CSV 校验（列/枚举/日期/来源 section/基因符号），可离线降级

独立工具约束（清单 §1.10）：单文件、只依赖标准库 + requests，**import 区不出现
WY**——内部使用与对外分享是同一份代码。缓存/限速自带，不依赖平台 http_client。

数据来源与许可：
  - Open Targets Platform GraphQL（CC BY 4.0，无需 key）：药物 MOA 及靶点基因符号；
  - openFDA drug/ndc（公开域，无需 key，限速 240 次/分、1000 次/天）：
    FDA SPL 派生的药理类别（EPC / MoA / CS 三类原文）；
  - DailyMed SPL（公有领域）：说明书 XML 的 12.1 Mechanism of Action
    （LOINC 43679-0）原文。
  ATC 已按 2026-11-18 定案砍掉：openFDA 与 SPL XML 均无该字段，WHO 亦不开放
  再分发；类别佐证用 EPC/MoA/CS 已足够（详见清单 §1.10）。

离线语义（2026-11-18 定案）：有缓存用缓存；无缓存则**警告并跳过该药/该校验**，
不静默失败、也不擅自联网重试到卡死。

循环性说明（清单 §1.10 陷阱 6）：OT MOA 与 R 值计算输入（ChEMBL 活性 + PK 暴露）
无交集，预填候选不构成循环论证；真值仍以人工核对说明书原文为准。
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import xml.etree.ElementTree as ET
from datetime import date as _date
from hashlib import sha1
from typing import Any, Iterable

import requests

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------
OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
OPENFDA_NDC = "https://api.fda.gov/drug/ndc.json"
DAILYMED_SPLS = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json"
DAILYMED_SPL = "https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{setid}.xml"

MOA_SECTION_CODE = "43679-0"          # LOINC：12.1 Mechanism of Action
PENDING = "待定"                       # 人工无法判定时的唯一合法占位
COLUMNS = ("drug", "label_mechanism_targets", "ot_moa_prefill",
           "spl_pharm_class", "label_section", "decision", "curator",
           "date", "notes")
DECISIONS = ("", "y", "t", "x")
USER_AGENT = "npap-dailymed-curate/1.0 (+https://github.com/; curation helper)"
THROTTLE_SECONDS = 0.34               # openFDA 240/min 与 DailyMed 礼貌间隔
RETRIES = 3
TIMEOUT = 30

CSV_HEADER_COMMENT = (
    "# dailymed_primary_targets.csv —— FDA 说明书机制靶点人工精选（二期 1.10）\n"
    "# 版本 1 | 2026-11-18 | 列定义与填表规则见 WY/resources/README.md\n"
    "# 许可：DailyMed 公有领域 / Open Targets CC BY 4.0 / openFDA 公开域\n"
)


class CurateError(RuntimeError):
    """可读的业务错误（CLI 打印后以非零码退出）。"""


# ---------------------------------------------------------------------------
# 轻量缓存 + 限速（独立工具自带，不 import WY）
# ---------------------------------------------------------------------------
def default_cache_dir() -> str:
    env = os.environ.get("NPAP_CURATE_CACHE", "").strip()
    if env:
        return env
    return os.path.join(os.path.expanduser("~"), ".cache", "npap", "dailymed_curate")


def _cache_file(cache_dir: str, kind: str, key: str, ext: str) -> str:
    digest = sha1(key.encode("utf-8")).hexdigest()[:16]
    sub = os.path.join(cache_dir, kind)
    os.makedirs(sub, exist_ok=True)
    return os.path.join(sub, f"{digest}.{ext}")


def _read_cache(path: str) -> Any:
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def _write_cache(path: str, value: Any) -> None:
    tmp = path + ".tmp"
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(value, f, ensure_ascii=False)
        os.replace(tmp, path)
    except OSError:
        pass  # 缓存失败不影响主流程


_THROTTLE: dict[str, float] = {}


def _throttle(url: str) -> None:
    host = urllib.parse.urlsplit(url).netloc
    last = _THROTTLE.get(host, 0.0)
    wait = THROTTLE_SECONDS - (time.monotonic() - last)
    if wait > 0:
        time.sleep(wait)
    _THROTTLE[host] = time.monotonic()


def _request(method: str, url: str, **kwargs) -> requests.Response | None:
    """带重试与礼貌间隔的请求；重试耗尽返回 None（由调用方决定降级）。"""
    kwargs.setdefault("timeout", TIMEOUT)
    kwargs.setdefault("headers", {"User-Agent": USER_AGENT})
    last_exc: Exception | None = None
    for attempt in range(RETRIES):
        _throttle(url)
        try:
            resp = requests.request(method, url, **kwargs)
            if resp.status_code in (429, 500, 502, 503, 504):
                last_exc = RuntimeError(f"HTTP {resp.status_code}")
                time.sleep(2.0 * (attempt + 1))
                continue
            return resp
        except requests.RequestException as exc:  # 网络层失败
            last_exc = exc
            time.sleep(2.0 * (attempt + 1))
    _warn(f"请求失败（重试 {RETRIES} 次）：{url} —— {last_exc}")
    return None


# ---------------------------------------------------------------------------
# 三个数据源
# ---------------------------------------------------------------------------
def ot_post(query: str, variables: dict, *, cache_dir: str,
            offline: bool = False) -> dict | None:
    """Open Targets GraphQL（只接受 POST）。离线/失败返回 None。"""
    key = "OT:" + json.dumps({"q": query, "v": variables}, sort_keys=True)
    path = _cache_file(cache_dir, "ot", key, "json")
    cached = _read_cache(path)
    if cached is not None:
        return cached
    if offline:
        _warn("离线且无缓存，跳过 Open Targets")
        return None
    resp = _request("POST", OT_GRAPHQL,
                    json={"query": query, "variables": variables})
    if resp is None or resp.status_code != 200:
        return None
    try:
        payload = resp.json()
    except ValueError:
        return None
    if payload.get("errors"):
        _warn(f"Open Targets 返回错误：{payload['errors'][:1]}")
        return None
    data = payload.get("data") or {}
    _write_cache(path, data)
    return data


OT_SEARCH_DRUG = """
query ($q: String!) {
  search(queryString: $q, entityNames: ["drug"]) {
    hits { id name entity }
  }
}
"""

OT_DRUG_MOA = """
query ($id: String!) {
  drug(chemblId: $id) {
    name
    mechanismsOfAction {
      rows { mechanismOfAction actionType targets { approvedSymbol } }
    }
  }
}
"""

OT_TARGET_LOOKUP = """
query ($q: String!) {
  search(queryString: $q, entityNames: ["target"]) {
    hits { id name entity }
  }
}
"""


def ot_drug_moa(drug: str, *, cache_dir: str, offline: bool = False) -> dict:
    """返回 {"chembl_id","moa":[...],"targets":[...]}；无结果时字段为空。"""
    out = {"chembl_id": "", "moa": [], "targets": []}
    data = ot_post(OT_SEARCH_DRUG, {"q": drug}, cache_dir=cache_dir, offline=offline)
    if not data:
        return out
    hits = ((data.get("search") or {}).get("hits") or [])
    drug_hits = [h for h in hits if (h.get("entity") or "").lower() == "drug"]
    if not drug_hits:
        return out
    chembl_id = str(drug_hits[0].get("id") or "")
    out["chembl_id"] = chembl_id
    if not chembl_id:
        return out
    data = ot_post(OT_DRUG_MOA, {"id": chembl_id}, cache_dir=cache_dir, offline=offline)
    rows = (((data or {}).get("drug") or {}).get("mechanismsOfAction") or {}).get("rows") or []
    targets: set[str] = set()
    moa: list[str] = []
    for row in rows:
        mech = str(row.get("mechanismOfAction") or "").strip()
        if mech and mech not in moa:
            moa.append(mech)
        for tgt in row.get("targets") or []:
            sym = str(tgt.get("approvedSymbol") or "").strip()
            if sym:
                targets.add(sym)
    out["moa"] = moa
    out["targets"] = sorted(targets)
    return out


def openfda_class(drug: str, *, cache_dir: str, offline: bool = False) -> list[str]:
    """openFDA NDC 的 FDA SPL 药理类别（EPC/MoA/CS 原文），按去重顺序返回。"""
    key = "OPENFDA:" + drug.lower()
    path = _cache_file(cache_dir, "openfda", key, "json")
    cached = _read_cache(path)
    if cached is not None:
        return list(cached)
    if offline:
        _warn("离线且无缓存，跳过 openFDA")
        return []
    params = {"search": f'openfda.generic_name:"{drug}"', "limit": 5}
    resp = _request("GET", OPENFDA_NDC, params=params)
    if resp is None or resp.status_code != 200:
        return []
    try:
        results = (resp.json() or {}).get("results") or []
    except ValueError:
        return []
    classes: list[str] = []
    for record in results:
        openfda = record.get("openfda") or {}
        for field in ("pharm_class_epc", "pharm_class_moa", "pharm_class_cs"):
            for value in openfda.get(field) or []:
                text = str(value).strip()
                if text and text not in classes:
                    classes.append(text)
    _write_cache(path, classes)
    return classes


def dailymed_moa(drug: str, *, cache_dir: str, offline: bool = False) -> dict:
    """取该药说明书 XML 的 12.1 机制章节原文；返回 {"setid","text"}。"""
    key = "DAILYMED:" + drug.lower()
    path = _cache_file(cache_dir, "dailymed", key, "json")
    cached = _read_cache(path)
    if cached is not None:
        return dict(cached)
    if offline:
        _warn("离线且无缓存，跳过 DailyMed")
        return {"setid": "", "text": ""}
    resp = _request("GET", DAILYMED_SPLS, params={"drug_name": drug, "pagesize": 5})
    if resp is None or resp.status_code != 200:
        return {"setid": "", "text": ""}
    try:
        data = (resp.json() or {}).get("data") or []
    except ValueError:
        return {"setid": "", "text": ""}
    if not data:
        return {"setid": "", "text": ""}
    # 优先标题包含药名的条目（多药复方标签噪声大），否则取第一条
    chosen = data[0]
    for item in data:
        if drug.lower() in str(item.get("title") or "").lower():
            chosen = item
            break
    setid = str(chosen.get("setid") or "")
    if not setid:
        return {"setid": "", "text": ""}
    resp = _request("GET", DAILYMED_SPL.format(setid=setid))
    text = ""
    if resp is not None and resp.status_code == 200:
        text = _extract_moa_section(resp.text)
    payload = {"setid": setid, "text": text}
    _write_cache(path, payload)
    return payload


def _extract_moa_section(xml_text: str) -> str:
    """从 SPL XML 抽取 LOINC 43679-0（12.1 Mechanism of Action）正文。"""
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return ""
    for section in root.iter():
        if not str(section.tag).endswith("section"):
            continue
        codes = [child for child in section if str(child.tag).endswith("code")]
        if not any(str(c.get("code")) == MOA_SECTION_CODE for c in codes):
            continue
        chunks = [chunk.strip() for chunk in section.itertext() if chunk.strip()]
        # 去掉标题与通用章节名，保留正文
        body = [c for c in chunks if not re.fullmatch(r"(12(\.1)?[ .]?)?Mechanism of Action", c, re.I)]
        return " ".join(body).strip()
    return ""


# ---------------------------------------------------------------------------
# CSV 读写
# ---------------------------------------------------------------------------
def read_rows(path: str) -> list[dict]:
    """读 CSV，跳过 `#` 注释头；列名去空白。"""
    if not os.path.isfile(path):
        raise CurateError(f"文件不存在：{path}")
    rows: list[dict] = []
    with open(path, encoding="utf-8", newline="") as f:
        lines = [line for line in f if not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise CurateError(f"缺少表头：{path}")
    for row in reader:
        rows.append({str(k).strip(): (v or "").strip() for k, v in row.items() if k is not None})
    return rows


def write_rows(path: str, rows: Iterable[dict]) -> None:
    """原子写 CSV（tmp + os.replace），保留固定注释头。"""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        f.write(CSV_HEADER_COMMENT)
        writer = csv.DictWriter(f, fieldnames=list(COLUMNS))
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in COLUMNS})
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# prefill
# ---------------------------------------------------------------------------
def _read_drug_names(args) -> list[str]:
    names: list[str] = []
    if args.drugs:
        with open(args.drugs, encoding="utf-8") as f:
            for line in f:
                text = line.strip()
                if text and not text.startswith("#"):
                    names.append(text)
    else:
        names = [row.get("drug", "") for row in read_rows(args.csv)]
    seen: set[str] = set()
    unique: list[str] = []
    for name in names:
        key = name.lower()
        if name and key not in seen:
            seen.add(key)
            unique.append(name)
    return unique


def prefill(args) -> int:
    rows = {row.get("drug", ""): row for row in read_rows(args.csv)} if os.path.isfile(args.csv) else {}
    names = _read_drug_names(args)
    if args.limit:
        names = names[: args.limit]
    if not names:
        raise CurateError("没有可预填的药物名（--drugs 或现有 CSV 的 drug 列）")

    built: list[dict] = []
    for index, drug in enumerate(names, 1):
        print(f"[prefill {index}/{len(names)}] {drug} …", file=sys.stderr)
        ot = ot_drug_moa(drug, cache_dir=args.cache_dir, offline=args.offline)
        classes = openfda_class(drug, cache_dir=args.cache_dir, offline=args.offline)
        dm = dailymed_moa(drug, cache_dir=args.cache_dir, offline=args.offline)

        row = dict(rows.get(drug) or {col: "" for col in COLUMNS})
        row["drug"] = drug
        row["ot_moa_prefill"] = ";".join(ot["targets"])
        row["spl_pharm_class"] = ";".join(classes)
        row["label_section"] = f"{dm['setid']}#{MOA_SECTION_CODE}" if dm.get("setid") else ""
        if dm.get("text"):
            snippet = dm["text"] if len(dm["text"]) <= 400 else dm["text"][:400] + "…"
            row["notes"] = _merge_note(row.get("notes", ""), f"MoA原文：{snippet}")
        if ot.get("moa"):
            row["notes"] = _merge_note(row.get("notes", ""), "OT机制：" + "；".join(ot["moa"]))
        if not row.get("curator"):
            row["curator"] = args.curator
        if not row.get("date"):
            row["date"] = _date.today().isoformat()
        built.append(row)

    # 一致性排序：有基因级答案且有原文证据的排最前，矛盾/缺失的垫后
    def rank(row: dict) -> tuple[int, str]:
        has_genes = bool(row.get("ot_moa_prefill"))
        has_text = "MoA原文" in (row.get("notes") or "")
        return (0 if (has_genes and has_text) else 1 if has_genes else 2), row["drug"].lower()

    merged = {row["drug"]: row for row in built}
    ordered = sorted(merged.values(), key=rank)
    write_rows(args.csv, ordered)
    with_genes = sum(1 for row in ordered if row.get("ot_moa_prefill"))
    print(f"[done] {len(ordered)} 行写入 {args.csv}（有 OT 基因级候选 {with_genes} 行）")
    return 0


def _merge_note(existing: str, extra: str) -> str:
    parts = [p for p in (existing or "").split(" | ") if p]
    if extra not in parts:
        parts.append(extra)
    return " | ".join(parts)


# ---------------------------------------------------------------------------
# review
# ---------------------------------------------------------------------------
def review(args) -> int:
    rows = read_rows(args.csv)
    pending = [row for row in rows if not row.get("decision")]
    if not pending:
        print("[review] 没有待审查行（decision 均已有值）")
        return 0
    print(f"[review] 待审查 {len(pending)} 行；y=确认 t=待定 x=否决 q=退出（每键立即写盘）")
    _print_cheatsheet()
    for index, row in enumerate(pending, 1):
        _render_card(index, len(pending), row)
        while True:
            answer = input("决策 [y/t/x/q]: ").strip().lower()
            if answer in ("y", "t", "x", "q"):
                break
        if answer == "q":
            print("[review] 已退出，进度已保存")
            return 0
        if answer == "y":
            genes = row.get("label_mechanism_targets") or row.get("ot_moa_prefill", "")
            if not genes:
                note = input("无 OT 候选，请手填靶点基因（分号分隔，留空则按待定）: ").strip()
                genes = note or PENDING
            row["decision"] = "y"
            row["label_mechanism_targets"] = genes
        elif answer == "t":
            row["decision"] = "t"
            row["label_mechanism_targets"] = PENDING
        else:
            note = input("否决原因（可留空）: ").strip()
            row["decision"] = "x"
            row["label_mechanism_targets"] = ""
            if note:
                row["notes"] = _merge_note(row.get("notes", ""), f"否决：{note}")
        row["curator"] = args.curator or row.get("curator", "")
        row["date"] = _date.today().isoformat()
        write_rows(args.csv, rows)   # 立即写盘（断点续传）
    print(f"[review] 完成，{len(pending)} 行已处理")
    return 0


def _render_card(index: int, total: int, row: dict) -> None:
    line = "=" * 72
    print(f"\n{line}\n[{index}/{total}] {row.get('drug', '')}\n{line}")
    print(f"① OT MOA 候选（基因级）：{row.get('ot_moa_prefill') or '（无）'}")
    print(f"② FDA 药理类别（EPC/MoA/CS）：{row.get('spl_pharm_class') or '（无）'}")
    print(f"③ DailyMed 来源：{row.get('label_section') or '（无）'}")
    notes = row.get("notes") or ""
    if notes:
        print(f"证据摘录：{notes}")


def _print_cheatsheet() -> None:
    print("常见药理类别 → 靶点基因速查（仅供提示，判定以说明书原文为准）：")
    print("  双胍类 → 无单一基因靶点（AMPK 通路，通常填待定）")
    print("  ACE 抑制剂 → ACE；ARB → AGTR1；他汀 → HMGCR；β 受体阻滞剂 → ADRB1/ADRB2")
    print("  质子泵抑制剂 → ATP4A/ATP4B；COX 抑制剂 → PTGS1/PTGS2")
    print("  DPP-4 抑制剂 → DPP4；SGLT2 抑制剂 → SLC5A2；磺脲类 → ABCC8/KCNJ11")


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------
_GENE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9\-]{0,19}$")


def validate(args) -> int:
    rows = read_rows(args.csv)
    errors: list[str] = []
    warnings: list[str] = []
    seen: set[str] = set()
    for idx, row in enumerate(rows, start=2):  # 头部注释行不计入，仅作提示
        drug = row.get("drug", "")
        if not drug:
            errors.append(f"第 {idx} 行：drug 为空")
            continue
        if drug.lower() in seen:
            errors.append(f"第 {idx} 行：药物重复 {drug}")
        seen.add(drug.lower())

        decision = row.get("decision", "")
        if decision not in DECISIONS:
            errors.append(f"第 {idx} 行：decision={decision!r} 不在 y/t/x/空 之内")
        if decision in ("y", "x", "t"):
            if not row.get("curator"):
                errors.append(f"第 {idx} 行：已判定但 curator 为空")
            if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", row.get("date", "")):
                errors.append(f"第 {idx} 行：date 需为 YYYY-MM-DD")
        section = row.get("label_section", "")
        if decision == "y":
            targets = row.get("label_mechanism_targets", "")
            if not targets or targets == PENDING:
                errors.append(f"第 {idx} 行：decision=y 但靶点为空或待定")
            else:
                for gene in [g.strip() for g in targets.split(";") if g.strip()]:
                    if not _GENE_RE.match(gene):
                        errors.append(f"第 {idx} 行：基因符号可疑 {gene!r}")
            if not re.fullmatch(r"[0-9a-fA-F\-]{8,}#43679-0", section):
                warnings.append(f"第 {idx} 行：label_section 未按 <setid>#43679-0 记录")
        if decision == "t" and row.get("label_mechanism_targets") != PENDING:
            errors.append(f"第 {idx} 行：decision=t 时靶点应填“待定”")
        if decision == "":
            warnings.append(f"第 {idx} 行：尚未审查（decision 为空）")

    if args.online and not args.offline:
        errors.extend(_check_genes_online(rows, args))

    for msg in warnings:
        _warn(msg)
    for msg in errors:
        print(f"[error] {msg}", file=sys.stderr)
    if errors:
        print(f"[validate] {len(errors)} 个错误、{len(warnings)} 个警告")
        return 1
    print(f"[validate] 通过（{len(rows)} 行，{len(warnings)} 个警告）")
    return 0


def _check_genes_online(rows: list[dict], args) -> list[str]:
    errors: list[str] = []
    for idx, row in enumerate(rows, start=2):
        if row.get("decision") != "y":
            continue
        for gene in [g.strip() for g in row.get("label_mechanism_targets", "").split(";") if g.strip()]:
            data = ot_post(OT_TARGET_LOOKUP, {"q": gene}, cache_dir=args.cache_dir)
            hits = ((data or {}).get("search") or {}).get("hits") or []
            if not any(str(h.get("id", "")).startswith("ENSG") for h in hits):
                errors.append(f"第 {idx} 行：OT 查不到靶点 {gene}")
    return errors


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _warn(msg: str) -> None:
    print(f"[warn] {msg}", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--csv", default=os.path.join("WY", "resources",
                                                      "dailymed_primary_targets.csv"),
                        help="精选 CSV 路径（默认 WY/resources/dailymed_primary_targets.csv）")
    common.add_argument("--cache-dir", default=default_cache_dir(), help="本地缓存目录")
    common.add_argument("--curator", default="", help="填表人（写入 curator 列）")

    parser = argparse.ArgumentParser(
        description="DailyMed 说明书机制靶点精选工具（prefill / review / validate）")
    sub = parser.add_subparsers(dest="command", required=True)

    p_pre = sub.add_parser("prefill", parents=[common], help="三级预填并生成审查队列")
    p_pre.add_argument("--drugs", help="药物名清单文件（每行一个；省略则读现有 CSV 的 drug 列）")
    p_pre.add_argument("--limit", type=int, default=0, help="只预填前 N 个（冒烟用）")
    p_pre.add_argument("--offline", action="store_true", help="只用缓存，不联网")

    p_rev = sub.add_parser("review", parents=[common], help="终端逐条审查（y/t/x，立即写盘）")
    p_rev.add_argument("--offline", action="store_true", help="保留参数（review 不联网）")

    p_val = sub.add_parser("validate", parents=[common], help="校验 CSV")
    p_val.add_argument("--offline", action="store_true", help="跳过在线基因校验")
    p_val.add_argument("--online", action="store_true", help="在线校验基因符号存在性（需网络）")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "prefill":
            return prefill(args)
        if args.command == "review":
            return review(args)
        return validate(args)
    except CurateError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
