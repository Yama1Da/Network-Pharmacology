#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""分子对接模块自检脚本:检查 AlphaFold 结构兜底相关功能是否正常。

纯本地运行,不发网络请求。运行方式(在 WY 目录下):
    python scripts/selfcheck_docking.py
"""

import json
import os
import sys
import tempfile
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from selfcheck_common import _tmp_uniprot_cache, patch_attr, run

import requests

import config
import core.docking as dock
import core.http_client as hc
from core.http_client import throttle


def t_clamp_box():
    box = {"center": (1.0, 2.0, 3.0), "size": (30.0, 35.0, 60.0)}
    under = dock.clamp_box(box, 40.0)
    assert under["size"] == (30.0, 35.0, 40.0), under  # 只有超限的边被收拢
    assert under["center"] == (1.0, 2.0, 3.0), "center must be preserved"
    small = dock.clamp_box({"center": (0, 0, 0), "size": (20, 20, 20)}, 40.0)
    assert small["size"] == (20.0, 20.0, 20.0), "under-limit box must be unchanged"


def t_strip_af_style_pdb():
    """AlphaFold 模型的 PDB 格式特殊,清洗后应只保留原子记录并去掉重复原子。"""
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
    """数据库里查无此基因时记入缓存,第二次查询不再重复发请求。"""
    calls = {"n": 0}

    def fake_ebi(gene, taxon):
        calls["n"] += 1
        return ""  # 确定性"查无"

    with _tmp_uniprot_cache(config, "selfcheck_uniprot_map.json"), \
            patch_attr(dock, "_lookup_uniprot_ebi", fake_ebi):
        assert dock.fetch_uniprot_accession("NOPE_X", 9606) == ""
        assert dock.fetch_uniprot_accession("NOPE_X", 9606) == ""
        assert calls["n"] == 1, f"negative cache failed: {calls['n']} calls"


def t_uniprot_network_failure_not_cached():
    """网络全部失败时应当报错且不写缓存,避免把临时故障当成"查无此基因"。"""
    def boom(gene, taxon):
        raise ConnectionError("simulated network outage")

    with _tmp_uniprot_cache(config, "selfcheck_uniprot_map2.json") as tmp, \
            patch_attr(dock, "_lookup_uniprot_ebi", boom), \
            patch_attr(dock, "_lookup_uniprot_rest", boom):
        try:
            dock.fetch_uniprot_accession("TP53", 9606)
            raise AssertionError("all sources down must raise RuntimeError")
        except RuntimeError:
            pass
        assert not os.path.exists(tmp) or not _json_nonempty(tmp), "failure must not be cached"


def _json_nonempty(path):
    try:
        with open(path, encoding="utf-8") as f:
            return bool(json.load(f))
    except Exception:
        return False


def t_af_pdb_skip_existing():
    """AlphaFold 模型文件已存在时直接复用,不再发网络请求。"""
    tmpdir = tempfile.mkdtemp()
    pre = os.path.join(tmpdir, f"AF_P99999_{config.ALPHAFOLD_MODEL_VERSION}.pdb")
    with open(pre, "w", encoding="utf-8") as f:
        f.write("REMARK x\n" + "ATOM " * 300)
    try:
        with patch_attr(config, "DOCKING_RECEPTOR_CACHE_DIR", tmpdir), \
                patch_attr(requests, "get",
                           lambda *a, **k: (_ for _ in ()).throw(AssertionError("network call made"))):
            path = dock.fetch_alphafold_pdb("P99999")
            assert path == pre
    finally:
        os.remove(pre)
        os.rmdir(tmpdir)


def t_candidates_af_fallback_row():
    """实验结构查不到且允许 AlphaFold 时,追加预测结构作候选;关闭时不追加。"""
    # 单基因入口是批量版的薄包装,检索统一走 _rcsb_search_entry_ids(带 HTTP 缓存会话),
    # 打桩层级在检索函数而不是 requests.post
    with patch_attr(dock, "fetch_uniprot_accession", lambda gene, taxon=9606: "P99999"), \
            patch_attr(dock, "_rcsb_search_entry_ids", lambda gene, max_hits: []):
        df = dock.fetch_receptor_pdb_candidates("FAKEGENE", allow_alphafold=True)
        assert len(df) == 1 and df.iloc[0]["pdb_id"] == "AFDB:P99999", df
        assert df.iloc[0]["source"] == "alphafold"
        df2 = dock.fetch_receptor_pdb_candidates("FAKEGENE", allow_alphafold=False)
        assert df2.empty and "source" in df2.columns


def t_candidates_search_failure_friendly_error():
    """检索服务不可达时,单基因入口抛出普通用户可读的中文错误(而不是静默按无结构处理)。"""
    def _boom(gene, max_hits):
        raise requests.exceptions.ConnectionError("simulated outage")

    with patch_attr(dock, "_rcsb_search_entry_ids", _boom):
        try:
            dock.fetch_receptor_pdb_candidates("TP53")
            raise AssertionError("expected RuntimeError")
        except RuntimeError as e:
            msg = str(e)
            assert "RCSB 搜索失败（TP53）" in msg, msg
            assert "无法连接 RCSB 服务器" in msg, msg
    # 批量模式保持原语义:单基因失败不中断整批,按无结构处理
    with patch_attr(dock, "_rcsb_search_entry_ids", _boom):
        results = dock.fetch_receptor_candidates_batch(["TP53", "EGFR"])
        assert results["TP53"].empty and results["EGFR"].empty


def t_throttle_spacing():
    """对同一网站的请求有最小时间间隔;未配置的请求不限速;并发调用也正确排队。"""
    orig = config.HOST_MIN_REQUEST_INTERVAL
    config.HOST_MIN_REQUEST_INTERVAL = dict(orig or {})
    config.HOST_MIN_REQUEST_INTERVAL["selftest"] = 0.3
    with hc._throttle_lock:
        hc._next_allowed_at.pop("selftest", None)
    try:
        # 串行 3 次 → 至少跨越 2 个完整间隔
        t0 = time.monotonic()
        for _ in range(3):
            throttle("selftest")
        el = time.monotonic() - t0
        assert el >= 0.6, f"3 calls in {el:.3f}s, expected >= 0.6 (2 intervals)"

        # 未知键 = 不限速
        t0 = time.monotonic()
        for _ in range(3):
            throttle("unknown_host_x")
        assert time.monotonic() - t0 < 0.05

        # 4 个线程同时请求 → 总耗时至少 3 个间隔(不会一起放行)
        barrier = threading.Barrier(4)
        t0 = time.monotonic()
        threads = [threading.Thread(target=lambda: (barrier.wait(), throttle("selftest")))
                   for _ in range(4)]
        [t.start() for t in threads]
        [t.join() for t in threads]
        el = time.monotonic() - t0
        assert el >= 0.9, f"4 concurrent calls in {el:.3f}s, expected >= 0.9 (3 intervals)"
    finally:
        config.HOST_MIN_REQUEST_INTERVAL = orig
        with hc._throttle_lock:
            hc._next_allowed_at.pop("selftest", None)


def t_interactive_3d_html():
    """交互式 3D HTML:坐标须以合法 JSON 内嵌,离线/在线两种模式标记正确。"""
    receptor = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
    ligand = "ATOM      1  N   LIG B   1       1.000   1.000   1.000  1.00  0.00           N\n"
    # 离线模式:传入库源码,HTML 直接内嵌且标记为离线可用
    html, offline = dock.build_interactive_3d_html(
        receptor, ligand, caption="gene · pose 1", lib_js="var $3Dmol = {};")
    assert offline and "$3Dmol.createViewer" in html
    assert 'var receptor = "ATOM' in html, "receptor must be embedded as a JSON string"
    # 内嵌的坐标文本必须是合法 JSON(转义正确,浏览器端才能还原)
    embedded = html.split("var receptor = ", 1)[1].split(";\n", 1)[0]
    assert json.loads(embedded) == receptor, "embedded receptor must round-trip"
    # 在线模式:不传库源码,应改用 CDN 引用并明确标记非离线
    html2, offline2 = dock.build_interactive_3d_html(receptor, ligand)
    assert not offline2 and 'src="https://3dmol.org' in html2
    # 默认导出:无口袋残基、无共晶配体参照,浏览器端不渲染对应元素
    assert "var pocket = []" in html and "var het = null" in html
    # 传入口袋残基与共晶配体:两类数据须以合法 JSON 嵌入(浏览器端据此还原样式与标注)
    pocket = [{"chain": "A", "resi": "1", "resn": "ALA", "pos": (0.0, 0.0, 0.0)}]
    html3, _ = dock.build_interactive_3d_html(
        receptor, ligand, pocket=pocket,
        het_pdb="HETATM    1  O   LIG B   2       2.000   1.000   1.000  1.00  0.00           O\n")
    assert 'var pocket = [{"chain": "A"' in html3, "pocket must be embedded as JSON"
    assert 'var het = "HETATM' in html3, "co-crystal ligand text must be embedded as JSON"


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
        ("candidates search failure friendly error", t_candidates_search_failure_friendly_error),
        ("interactive 3D html embed", t_interactive_3d_html),
    ]
    run(tests, "docking AlphaFold selfcheck")


if __name__ == "__main__":
    main()
