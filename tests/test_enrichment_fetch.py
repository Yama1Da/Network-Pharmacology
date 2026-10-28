"""core.enrichment 抓取链补测(假请求,零网络):缓存命中/键漂移/空结果/归一化空/
skip 守卫/g:Profiler 与 Enrichr 入口。

深水区冲线批 II。契约:request 以零参闭包被 _call_with_request_guard 调用。
"""
import pandas as pd
import pytest

import WY.core.enrichment as en


def _normalized(terms):
    return pd.DataFrame([{"source": "KEGG", "term_id": t, "term_name": f"N{t}",
                          "p_value": 0.01, "adjusted_p_value": 0.01,
                          "intersection_size": 2, "query_size": 3,
                          "genes": "A;B"} for t in terms])


def _run_fetch(genes, *, request, normalize=None, label="lbl", status=None,
               skip_if=None):
    return en._fetch_enrichment_cached(
        genes, "hsapiens", 0.05, label, status,
        source_name="gProfiler", cache_key_extra=("GO:BP",), input_ok=bool(genes),
        request=request, normalize=normalize or (lambda res: _normalized(["hsa04151"])),
        skip_if=skip_if)


def test_fetch_cached_hit():
    """第二次同参数查询:命中同名覆盖式缓存(头部键一致),状态 cache,零请求。"""
    calls = {"n": 0}

    def request():
        calls["n"] += 1
        return [{"x": 1}]

    st1: dict = {}
    df1 = _run_fetch(["TP53", "EGFR"], request=request, status=st1)
    assert st1["reason"] == "success" and calls["n"] == 1
    st2: dict = {}
    df2 = _run_fetch(["TP53", "EGFR"], request=request, status=st2)
    assert st2["reason"] == "cache" and calls["n"] == 1  # 零新请求
    pd.testing.assert_frame_equal(df1, df2)


def test_fetch_key_mismatch_recomputes():
    """同标识不同基因集:头部键不一致 → 视为未命中,重算并覆盖缓存。"""
    _run_fetch(["TP53"], request=lambda: [{"x": 1}])
    calls = {"n": 0}

    def request():
        calls["n"] += 1
        return [{"x": 1}]

    st: dict = {}
    df = _run_fetch(["KRAS"], request=request, status=st)
    assert st["reason"] == "success" and calls["n"] == 1


def test_fetch_empty_input():
    """空基因输入:请求前直接拦下,状态 empty_input。"""
    st: dict = {}
    df = _run_fetch([], status=st,
                    request=lambda: (_ for _ in ()).throw(AssertionError("called")))
    assert df.empty and st["reason"] == "empty_input"


def test_fetch_request_failure():
    """请求保护返回 None(重试耗尽):状态 request_failed,不写缓存。"""
    st: dict = {}
    df = _run_fetch(["TP53"], request=lambda: None, status=st)
    assert df.empty and st["reason"] == "request_failed"


def test_fetch_no_significant_not_cached():
    """请求成功但零条目:状态 no_significant(可信空结果),不写缓存(下次仍真请求)。"""
    calls = {"n": 0}

    def request():
        calls["n"] += 1
        return []

    st: dict = {}
    df = _run_fetch(["TP53"], request=request, status=st)
    assert df.empty and st["reason"] == "no_significant"
    st2: dict = {}
    _run_fetch(["TP53"], request=request, status=st2)
    assert calls["n"] == 2  # 空结果不留缓存


def test_fetch_normalize_empty_not_cached():
    """归一化后为空(异常响应):状态 normalize_empty,不把异常空表写进缓存。"""
    st: dict = {}
    df = _run_fetch(["TP53"], request=lambda: [{"x": 1}], status=st,
                    normalize=lambda res: pd.DataFrame())
    assert df.empty and st["reason"] == "normalize_empty"


def test_fetch_skip_guard():
    """skip_if 命中(如熔断开启):不发请求,状态透传跳过原因。"""
    st: dict = {}
    df = _run_fetch(["TP53"], status=st,
                    request=lambda: (_ for _ in ()).throw(AssertionError("called")),
                    skip_if=lambda: "circuit_open")
    assert df.empty and st["reason"] == "circuit_open"


def test_gprofiler_and_enrichr_entries(monkeypatch):
    """两个数据源入口:g:Profiler 正常解析落缓存(假请求返回 DataFrame,与真实
    _profile_gprofiler 同构);Enrichr 空结果记 no_significant。"""
    monkeypatch.setattr(en, "_profile_gprofiler", lambda *a, **k: pd.DataFrame([
        {"source": "KEGG", "native": "hsa04151", "name": "PI3K-Akt",
         "p_value": 0.01, "intersection_size": 2, "query_size": 3,
         "intersections": "TP53;EGFR"}]))
    st: dict = {}
    df = en.fetch_enrichment_gprofiler(["TP53"], status=st)
    assert st["reason"] == "success" and list(df["term_id"]) == ["hsa04151"]

    monkeypatch.setattr(en, "_query_enrichr", lambda *a, **k: pd.DataFrame())
    st2: dict = {}
    df2 = en.fetch_enrichment_enrichr(["TP53"], libs=["KEGG_2021_Human"], status=st2)
    assert df2.empty and st2["reason"] == "no_significant"
