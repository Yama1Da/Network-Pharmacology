"""core.interactome 的测试:蛋白信息列挑选、边表转换、分段下载、互作组加载与新鲜度。

迁移自 selfcheck_network.py(44 项中的 8 项);下载一律打桩(requests.get 假响应),
真实互作组未构建时 pytest.skip(CI 机器上同样跳过,跳过不是失败)。
"""
import os

import pandas as pd
import pytest
import requests

from WY import config
from WY.core import interactome as itx


def test_interactome_info_cols():
    """_info_column_names 挑出 STRING 蛋白信息表的必备列;列名不符时必须报 ValueError。"""
    ok = pd.DataFrame({"protein_external_id": ["9606.ENSP1"],
                       "preferredName": ["TP53"], "other": [1]})
    assert itx._info_column_names(ok) == ("protein_external_id", "preferredName")
    bad = pd.DataFrame({"foo": [1], "bar": [2]})
    try:
        itx._info_column_names(bad)
        raise AssertionError("unexpected columns must raise ValueError")
    except ValueError:
        pass


def test_interactome_links_to_edges():
    """links → edges:自环去除、同对蛋白取最高分、置信度阈值过滤、端点按字典序排列。"""
    links = pd.DataFrame({
        "protein1": ["9606.E1", "9606.E2", "9606.E3", "9606.E4",
                     "9606.E5", "9606.E6", "9606.E7"],
        "protein2": ["9606.E2", "9606.E1", "9606.E3", "9606.E5",
                     "9606.E6", "9606.E7", "9606.E8"],
        # 前两行是同一对蛋白,分数 700 和 800
        "combined_score": [700, 800, 699, 950,
                           900, 850, 400],
    })
    id2sym = {"9606.E1": "TP53", "9606.E2": "EGFR", "9606.E3": "EGFR",
              "9606.E4": "KRAS", "9606.E5": "BRAF", "9606.E6": "BRCA1",
              "9606.E7": "BRCA1", "9606.E8": "ABL1"}
    edges = itx._links_to_edges(links, id2sym, confidence=700)
    records = {(r.protein_a, r.protein_b): int(r.score) for r in edges.itertuples()}
    # 无向边的两个端点统一按基因名字典序排列
    # 1) 基因塌缩产生的自环(EGFR-EGFR)必须去掉
    assert ("EGFR", "EGFR") not in records, "self-loop must be dropped"
    # 2) 同一对蛋白出现两次时保留最高分 800
    assert records.get(("EGFR", "TP53")) == 800, records
    # 3) E7-E8 分数 400 低于阈值,被过滤
    assert ("ABL1", "BRCA1") not in records and ("BRCA1", "ABL1") not in records
    # 4) E4-E5 (KRAS-BRAF, 950) 与 E5-E6 (BRAF-BRCA1, 900) 存在
    assert records.get(("BRAF", "KRAS")) == 950
    assert records.get(("BRAF", "BRCA1")) == 900
    assert len(edges) == 3, f"expect 3 edges, got {len(edges)}"
    # 5) 置信度换算一致:0.7 * 1000 == 700
    assert config.INTERACTOME_CONFIDENCE == int(config.PPI_CONFIDENCE * 1000)


class _FakeResp:
    """_download_gz 打桩用的最小响应对象:支持 with 上下文与流式 iter_content。"""

    def __init__(self, status, headers=None, body=b""):
        self.status_code = status
        self.headers = headers or {}
        self._body = body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def close(self):
        pass

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"http {self.status_code}")

    def iter_content(self, chunk_size=1 << 20):
        for i in range(0, len(self._body), chunk_size):
            yield self._body[i:i + chunk_size]


def _fake_string_get(payload, fail_ranges=False, delay=0.0):
    """伪装 STRING 下载源:带 Range 头按区间切片回 206,Range 探测回总大小,
    不带 Range 回 200 全量(模拟不支持分段的源);fail_ranges=True 时区间请求一律失败;
    delay>0 时每段响应前等待(把下载时长拖过 _download_gz_ranged 的 0.5 秒进度轮询窗口,
    使进度事件触发成为确定性——真实 GB 级文件必然跨过,内存切片的假下载否则瞬时完成)。"""
    import time

    def fake_get(url, headers=None, stream=True, timeout=None, **kw):
        rng = (headers or {}).get("Range", "")
        if rng == "bytes=0-0":
            return _FakeResp(206, {"Content-Range": f"bytes 0-0/{len(payload)}"},
                             payload[:1])
        if rng.startswith("bytes="):
            if fail_ranges:
                raise IOError("simulated segment failure")
            if delay:
                time.sleep(delay)
            s, e = rng[6:].split("-")
            return _FakeResp(206, {}, payload[int(s):int(e) + 1])
        return _FakeResp(200, {}, payload)

    return fake_get


def test_interactome_ranged_download(monkeypatch, tmp_path):
    """Range 可用时按段并发下载并拼接出与原文件一致的产物,无 .part 残留。"""
    payload = bytes(range(256)) * 64  # 16 KiB,把分段下限调小以便离线测试
    events = []
    # delay=0.3:分段下载跨越 0.5s 进度轮询窗口,进度事件至少触发一次(CI 曾因
    # 假下载瞬时完成而零事件,属时序侥幸翻车)
    monkeypatch.setattr(requests, "get", _fake_string_get(payload, delay=0.3))
    monkeypatch.setattr(itx, "_RANGE_MIN_SIZE", 4096)
    monkeypatch.setattr(itx, "_RANGE_SEG_MIN_SIZE", 4096)
    dest = os.path.join(str(tmp_path), "links.txt.gz")
    out = itx._download_gz("/x.gz", dest,
                           progress_cb=lambda *a: events.append(a),
                           progress_label="download_links")
    assert out == dest
    with open(dest, "rb") as f:
        assert f.read() == payload
    leftover = [n for n in os.listdir(str(tmp_path)) if ".part" in n]
    assert not leftover, leftover
    assert any(ev[0] == "download_links" for ev in events), events


def test_interactome_fallback_notice(monkeypatch, tmp_path):
    """大文件的分段在重试轮次耗尽后仍失败:回退单流,文件完整,且只发一次 fallback 提醒。"""
    payload = bytes(range(256)) * 64  # 16 KiB
    events = []
    monkeypatch.setattr(requests, "get", _fake_string_get(payload, fail_ranges=True))
    monkeypatch.setattr(itx, "_RANGE_MIN_SIZE", 4096)
    monkeypatch.setattr(itx, "_RANGE_SEG_MIN_SIZE", 4096)
    monkeypatch.setattr(config, "API_MAX_RETRIES", 1)
    monkeypatch.setattr(config, "API_RETRY_DELAY", 0)
    dest = os.path.join(str(tmp_path), "links.txt.gz")
    itx._download_gz("/x.gz", dest, progress_cb=lambda *a: events.append(a))
    with open(dest, "rb") as f:
        assert f.read() == payload
    assert sum(1 for ev in events if ev[0] == "fallback") == 1, events


def test_interactome_small_file_no_notice(monkeypatch, tmp_path):
    """小文件(低于分段阈值)走单流属正常路径,不应触发 fallback 误报。"""
    payload = b"tiny" * 512  # 2 KiB

    def no_range_get(url, headers=None, stream=True, timeout=None, **kw):
        return _FakeResp(200, {}, payload)

    events = []
    monkeypatch.setattr(requests, "get", no_range_get)
    monkeypatch.setattr(itx, "_RANGE_MIN_SIZE", 4096)
    monkeypatch.setattr(itx, "_RANGE_SEG_MIN_SIZE", 4096)
    dest = os.path.join(str(tmp_path), "info.txt.gz")
    itx._download_gz("/x.gz", dest, progress_cb=lambda *a: events.append(a))
    with open(dest, "rb") as f:
        assert f.read() == payload
    assert not any(ev[0] == "fallback" for ev in events), events


def test_interactome_segment_resume(monkeypatch, tmp_path):
    """源每次只给出区间的一小段就提前收流:靠段内续传补完,全程不应触发 fallback。"""
    payload = bytes(range(256)) * 64  # 16 KiB
    calls = {"range": 0}

    def stingy_get(url, headers=None, stream=True, timeout=None, **kw):
        rng = (headers or {}).get("Range", "")
        if rng == "bytes=0-0":
            return _FakeResp(206, {"Content-Range": f"bytes 0-0/{len(payload)}"},
                             payload[:1])
        if rng.startswith("bytes="):
            s, e = (int(x) for x in rng[6:].split("-"))
            calls["range"] += 1
            cut = min(s + 512, e + 1)  # 每次只给 512 字节就提前收流
            return _FakeResp(206, {}, payload[s:cut])
        return _FakeResp(200, {}, payload)

    events = []
    monkeypatch.setattr(requests, "get", stingy_get)
    monkeypatch.setattr(itx, "_RANGE_MIN_SIZE", 4096)
    monkeypatch.setattr(itx, "_RANGE_SEG_MIN_SIZE", 4096)
    monkeypatch.setattr(config, "API_MAX_RETRIES", 64)
    monkeypatch.setattr(config, "API_RETRY_DELAY", 0)
    dest = os.path.join(str(tmp_path), "links.txt.gz")
    itx._download_gz("/x.gz", dest, progress_cb=lambda *a: events.append(a))
    with open(dest, "rb") as f:
        assert f.read() == payload
    assert not any(ev[0] == "fallback" for ev in events), events
    assert calls["range"] > 4, calls  # 每段都要多次续传才能补完


def test_interactome_load_or_skip():
    """互作组未构建时跳过;已构建时正常加载并检查基因覆盖情况。"""
    if not itx.is_built("Human"):
        pytest.skip("Human 背景互作组未构建")
    graph = itx.load_interactome("Human")
    assert graph.number_of_nodes() >= 10000, f"human interactome too small: {graph.number_of_nodes()}"
    inside, outside = itx.get_network_coverage(["tp53 ", "EGFR", "FAKE_GENE_X", "tp53"])
    assert inside == ["TP53", "EGFR"], inside
    assert outside == ["FAKE_GENE_X"], outside
    deg = itx.load_degree_table("Human")
    assert set(deg.columns) == {"gene", "degree"} and len(deg) == graph.number_of_nodes()


def test_interactome_staleness():
    """检查互作组数据的"过期"判断:未构建时 age 为 None;已构建时 age 为非负整数,
    且 is_stale 的阈值比较是严格大于。"""
    if not itx.is_built("Human"):
        pytest.skip("Human 背景互作组未构建")
    age = itx.interactome_age_days("Human")
    assert isinstance(age, int) and age >= 0, age
    # 判断规则:age 严格大于 days 才算过期
    assert itx.is_stale("Human", days=-1) is True, "negative threshold must be stale"
    assert itx.is_stale("Human", days=age) is False, "boundary: age > days must be strict"
    assert itx.is_stale("Human", days=10_000) is False, "huge threshold must not be stale"
    # 默认阈值来自 config,且为正数
    assert config.INTERACTOME_STALE_DAYS > 0
