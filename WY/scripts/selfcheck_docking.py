#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""M6 自检：AlphaFold 兜底链路的离线不变式（纯本地，无网络请求）。

运行方式（在 WY 目录下）：
    python scripts/selfcheck_docking.py

覆盖点：
  1. clamp_box：未超限原样返回 / 超限以质心为中心收拢
  2. _strip_receptor_pdb 兼容 AF 风格 PDB（REMARK 头 + 仅 ATOM + altloc 去重）
  3. EBI 条目基因名解析（name.value + synonyms[].value）
  4. fetch_uniprot_accession：负缓存（查无不再重复请求）/ 网络失败不缓存且上抛
  5. fetch_alphafold_pdb：文件已存在时跳过下载（不发任何请求）
  6. fetch_receptor_pdb_candidates：RCSB 空结果 + AF 开启 → 追加 AFDB 行（source=alphafold）
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import requests

import config
import core.docking as dock
from core.http_client import throttle


def check(name, fn):
    try:
        fn()
        print(f"PASS  {name}")
        return True
    except AssertionError as e:
        print(f"FAIL  {name}: {e}")
        return False
    except Exception as e:
        print(f"ERROR {name}: {type(e).__name__}: {e}")
        return False


def t_clamp_box():
    box = {"center": (1.0, 2.0, 3.0), "size": (30.0, 35.0, 60.0)}
    under = dock.clamp_box(box, 40.0)
    assert under["size"] == (30.0, 35.0, 40.0), under  # 只有超限的边被收拢
    assert under["center"] == (1.0, 2.0, 3.0), "center must be preserved"
    small = dock.clamp_box({"center": (0, 0, 0), "size": (20, 20, 20)}, 40.0)
    assert small["size"] == (20.0, 20.0, 20.0), "under-limit box must be unchanged"


def t_strip_af_style_pdb():
    """AF 模型为 REMARK 头 + 仅 ATOM（无 HETATM），清洗后应保留 ATOM 并去 altloc。"""
    src = os.path.join(tempfile.gettempdir(), "selfcheck_af_like.pdb")
    out = src + ".clean.pdb"
    with open(src, "w", encoding="utf-8") as f:
        f.write("REMARK   1 ALPHAFOLD MONOMER V2.0 PREDICTION\n")
        f.write("TITLE    ALPHAFOLD MONOMER PREDICTION FOR TEST\n")
        f.write("ATOM      1  N   MET A   1      63.063 -13.528  -5.589  1.00 42.72           N\n")
        f.write("ATOM      2  CA  MET A   1      62.454 -12.831  -4.439  1.00 42.72           C\n")
        f.write("ATOM      3  CA  MET A   1      62.454 -12.831  -4.439  1.00 42.72 B\n")  # altloc B → 去重
        f.write("HETATM    4  O   HOH A 999      10.000  10.000  10.000  1.00 20.00           O\n")
        f.write("END\n")
    try:
        dock._strip_receptor_pdb(src, out)
        lines = open(out, encoding="utf-8").read().splitlines()
        assert all(ln.startswith(("ATOM", "TER", "END")) for ln in lines), lines
        assert sum(1 for ln in lines if ln.startswith("ATOM")) == 2, "altloc B must be dropped"
    finally:
        for p in (src, out):
            if os.path.exists(p):
                os.remove(p)


def t_ebi_gene_names():
    entry = {"gene": [{"name": {"value": "ANO9"},
                       "synonyms": [{"value": "PIG5"}, {"value": "TP53I5"}]}]}
    names = dock._gene_names_from_ebi_entry(entry)
    assert names == {"ANO9", "PIG5", "TP53I5"}, names
    assert dock._gene_names_from_ebi_entry({}) == set()


def t_uniprot_negative_cache():
    """查无（API 正常返回但无精确匹配）→ 负缓存；第二次调用不再发起请求。"""
    orig_cache = config.UNIPROT_LOOKUP_CACHE
    orig_ebi = dock._lookup_uniprot_ebi
    calls = {"n": 0}

    def fake_ebi(gene, taxon):
        calls["n"] += 1
        return ""  # 确定性"查无"

    tmp = os.path.join(tempfile.gettempdir(), "selfcheck_uniprot_map.json")
    if os.path.exists(tmp):
        os.remove(tmp)
    config.UNIPROT_LOOKUP_CACHE = tmp
    dock._lookup_uniprot_ebi = fake_ebi
    try:
        assert dock.fetch_uniprot_accession("NOPE_X", 9606) == ""
        assert dock.fetch_uniprot_accession("NOPE_X", 9606) == ""
        assert calls["n"] == 1, f"negative cache failed: {calls['n']} calls"
    finally:
        dock._lookup_uniprot_ebi = orig_ebi
        config.UNIPROT_LOOKUP_CACHE = orig_cache
        if os.path.exists(tmp):
            os.remove(tmp)


def t_uniprot_network_failure_not_cached():
    """全部源网络失败 → 抛 RuntimeError 且不写缓存（抖动不得固化为"查无"）。"""
    orig_cache = config.UNIPROT_LOOKUP_CACHE
    orig = (dock._lookup_uniprot_ebi, dock._lookup_uniprot_rest)

    def boom(gene, taxon):
        raise ConnectionError("simulated network outage")

    tmp = os.path.join(tempfile.gettempdir(), "selfcheck_uniprot_map2.json")
    if os.path.exists(tmp):
        os.remove(tmp)
    config.UNIPROT_LOOKUP_CACHE = tmp
    dock._lookup_uniprot_ebi = boom
    dock._lookup_uniprot_rest = boom
    try:
        try:
            dock.fetch_uniprot_accession("TP53", 9606)
            raise AssertionError("all sources down must raise RuntimeError")
        except RuntimeError:
            pass
        assert not os.path.exists(tmp) or not _json_nonempty(tmp), "failure must not be cached"
    finally:
        dock._lookup_uniprot_ebi, dock._lookup_uniprot_rest = orig
        config.UNIPROT_LOOKUP_CACHE = orig_cache
        if os.path.exists(tmp):
            os.remove(tmp)


def _json_nonempty(path):
    import json
    try:
        with open(path, encoding="utf-8") as f:
            return bool(json.load(f))
    except Exception:
        return False


def t_af_pdb_skip_existing():
    """AF 模型文件已存在时不发任何请求（直接复用）。"""
    tmpdir = tempfile.mkdtemp()
    orig_dir, orig_get = config.DOCKING_RECEPTOR_CACHE_DIR, requests.get
    config.DOCKING_RECEPTOR_CACHE_DIR = tmpdir
    pre = os.path.join(tmpdir, f"AF_P99999_{config.ALPHAFOLD_MODEL_VERSION}.pdb")
    with open(pre, "w", encoding="utf-8") as f:
        f.write("REMARK x\n" + "ATOM " * 300)
    requests.get = lambda *a, **k: (_ for _ in ()).throw(AssertionError("network call made"))
    try:
        path = dock.fetch_alphafold_pdb("P99999")
        assert path == pre
    finally:
        requests.get = orig_get
        config.DOCKING_RECEPTOR_CACHE_DIR = orig_dir
        os.remove(pre)
        os.rmdir(tmpdir)


def t_candidates_af_fallback_row():
    """RCSB 空结果 + AF 开启 → 追加 AFDB 行（source=alphafold）；关闭时不追加。"""
    orig_acc = dock.fetch_uniprot_accession
    orig_post = requests.post
    dock.fetch_uniprot_accession = lambda gene, taxon=9606: "P99999"

    class _FakeResp:
        def raise_for_status(self):
            return None

        def json(self):
            return {"result_set": []}  # RCSB 确定无实验结构

    requests.post = lambda *a, **k: _FakeResp()
    try:
        df = dock.fetch_receptor_pdb_candidates("FAKEGENE", allow_alphafold=True)
        assert len(df) == 1 and df.iloc[0]["pdb_id"] == "AFDB:P99999", df
        assert df.iloc[0]["source"] == "alphafold"
        df2 = dock.fetch_receptor_pdb_candidates("FAKEGENE", allow_alphafold=False)
        assert df2.empty and "source" in df2.columns
    finally:
        dock.fetch_uniprot_accession = orig_acc
        requests.post = orig_post


def t_throttle_spacing():
    """主机级节流：同键串行调用间隔 ≥ 最小间隔；未知键不限速；并发线程正确排队。"""
    import threading
    import time as _t

    import core.http_client as hc

    orig = config.HOST_MIN_REQUEST_INTERVAL
    config.HOST_MIN_REQUEST_INTERVAL = dict(orig or {})
    config.HOST_MIN_REQUEST_INTERVAL["selftest"] = 0.3
    with hc._throttle_lock:
        hc._next_allowed_at.pop("selftest", None)
    try:
        # 串行 3 次 → 至少跨越 2 个完整间隔
        t0 = _t.monotonic()
        for _ in range(3):
            throttle("selftest")
        el = _t.monotonic() - t0
        assert el >= 0.6, f"3 calls in {el:.3f}s, expected >= 0.6 (2 intervals)"

        # 未知键 = 不限速
        t0 = _t.monotonic()
        for _ in range(3):
            throttle("unknown_host_x")
        assert _t.monotonic() - t0 < 0.05

        # 并发预约：4 线程同时抢 → 总耗时 ≥ 3 个间隔（不会同时睡醒连发）
        barrier = threading.Barrier(4)
        t0 = _t.monotonic()
        threads = [threading.Thread(target=lambda: (barrier.wait(), throttle("selftest")))
                   for _ in range(4)]
        [t.start() for t in threads]
        [t.join() for t in threads]
        el = _t.monotonic() - t0
        assert el >= 0.9, f"4 concurrent calls in {el:.3f}s, expected >= 0.9 (3 intervals)"
    finally:
        config.HOST_MIN_REQUEST_INTERVAL = orig
        with hc._throttle_lock:
            hc._next_allowed_at.pop("selftest", None)


def main():
    tests = [
        ("clamp_box under/over limit", t_clamp_box),
        ("throttle spacing & concurrency", t_throttle_spacing),
        ("strip AF-style PDB (altloc/HETATM)", t_strip_af_style_pdb),
        ("EBI entry gene-name parsing", t_ebi_gene_names),
        ("uniprot negative cache", t_uniprot_negative_cache),
        ("uniprot network failure not cached", t_uniprot_network_failure_not_cached),
        ("alphafold pdb skip-if-exists", t_af_pdb_skip_existing),
        ("candidates AF fallback row", t_candidates_af_fallback_row),
    ]
    print(f"==== docking AlphaFold selfcheck: {len(tests)} tests ====")
    ok = all([check(name, fn) for name, fn in tests])
    print("==== ALL PASS ====" if ok else "==== FAILED ====")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
