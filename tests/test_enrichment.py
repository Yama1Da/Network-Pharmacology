"""core.enrichment 与富集编排(na.enrich_modules)的测试:模块富集、熔断器、GSEA、通路图。

迁移自 selfcheck_network.py(44 项中的 14 项)。富集服务一律打桩(假函数,不联网);
熔断标记文件随 config.ENRICHMENT_CACHE_DIR 被 conftest 沙箱隔离到每个测试的
临时目录——自检脚本里的手工沙箱(_circuit_sandbox)由此整个省去。
"""
import glob
import os
import time

import pandas as pd

import WY.core.network_algorithms as na
from WY import config
from WY.core import enrichment as en
from WY.core.caching import save_json_file
from WY.core.enrichment import (
    TRUSTED_REASONS,
    _gsea_cache_file,
    _normalize_enrichment_gsea,
    build_ranked_list,
    fetch_enrichment_gsea_prerank,
    render_drug_target_pathway_plot,
)


def _fake_fetch(reason, df):
    """假的富集服务函数工厂(g:Profiler/Enrichr 共用,调用方均按关键字传参):
    reason=None 表示成功,否则按给定 reason 记录状态。"""
    def f(genes, libs=None, organism="hsapiens", sources=("GO:BP",), pvalue=0.05,
          label=None, status=None):
        if status is not None:
            if reason is None:
                status.update(ok=True, reason="success")
            else:
                status.update(ok=reason in TRUSTED_REASONS, reason=reason)
        return df.copy()
    return f


def _module5(genes=("A", "B", "C", "D", "E"), module_id=1):
    """测试用的 5 基因模块 dict(多个富集测试共用同一形状的数据)。"""
    return {"module_id": module_id, "genes": list(genes),
            "size": len(genes), "avg_degree": 2.0}


def _purge(pattern):
    """按 glob 模式删除旧缓存文件,文件不存在时为无害空操作。"""
    for old in glob.glob(pattern):
        try:
            os.remove(old)
        except OSError:
            pass


def _expire_circuit(circuit, stage: int) -> None:
    """写入一条"已过冷却期"的失败标记,让下一次请求作为探测放行。"""
    save_json_file(circuit._file(),
                   {"failed_at": time.time() - circuit.cooldown(stage) - 1,
                    "stage": stage})


# ---------- 模块富集编排(两个服务均打桩) ----------


def test_enrich_fail_both(monkeypatch):
    """两个富集服务都失败 → 结果表为空、失败记录在案、且不算"确无显著项"。"""
    monkeypatch.setattr(na, "fetch_enrichment_gprofiler",
                        _fake_fetch("request_failed", pd.DataFrame()))
    monkeypatch.setattr(na, "fetch_enrichment_enrichr",
                        _fake_fetch("request_failed", pd.DataFrame()))
    modules = [_module5()]
    long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
    assert long_df.empty, "failed module must produce empty table"
    assert len(failures) == 1 and failures[0]["module_id"] == 1, "failure must be recorded"
    assert not no_sig, "failed module must not be counted as no-significant"


def test_enrich_fail_gp_ok_enr(monkeypatch):
    """一个服务被熔断跳过、另一个确实无显著项 → 组合结果不可信,按失败处理。"""
    monkeypatch.setattr(na, "fetch_enrichment_gprofiler",
                        _fake_fetch("circuit_open", pd.DataFrame()))
    monkeypatch.setattr(na, "fetch_enrichment_enrichr",
                        _fake_fetch("no_significant", pd.DataFrame()))
    modules = [_module5()]
    long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
    assert len(failures) == 1, "circuit_open + empty enrichr must be a failure"
    assert not no_sig


def test_enrich_genuine_no_sig(monkeypatch):
    """两条链路都成功但确实无显著项 → 不算失败,记入 no_sig 列表。"""
    monkeypatch.setattr(na, "fetch_enrichment_gprofiler",
                        _fake_fetch("no_significant", pd.DataFrame()))
    monkeypatch.setattr(na, "fetch_enrichment_enrichr",
                        _fake_fetch("no_significant", pd.DataFrame()))
    modules = [_module5()]
    long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
    assert long_df.empty and not failures and no_sig == [1]


def test_enrich_success(monkeypatch):
    """模块 1 富集成功、有结果;模块 2 两个服务都失败,应记入 failures。"""
    ok_df = pd.DataFrame([
        {"source": "KEGG", "term_id": "hsa04151", "term_name": "PI3K-Akt",
         "p_value": 0.001, "adjusted_p_value": 0.001,
         "intersection_size": 3, "query_size": 5, "genes": "A;B;C"},
    ])

    def gp(genes, organism="hsapiens", sources=("GO:BP",), pvalue=0.05, label=None, status=None):
        if "A" in genes:
            if status is not None:
                status.update(ok=True, reason="success")
            return ok_df.copy()
        if status is not None:
            status.update(ok=False, reason="request_failed")
        return pd.DataFrame()

    def enr(genes, libs, organism="human", pvalue=0.05, label=None, status=None):
        if status is not None:
            status.update(ok=False, reason="request_failed")
        return pd.DataFrame()

    monkeypatch.setattr(na, "fetch_enrichment_gprofiler", gp)
    monkeypatch.setattr(na, "fetch_enrichment_enrichr", enr)
    modules = [
        _module5(),
        _module5(genes=("X", "Y", "Z", "W", "V"), module_id=2),
    ]
    long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens")
    assert not no_sig
    assert not long_df.empty and (long_df["module_id"] == 1).all()
    assert "module_size" in long_df.columns
    assert [f["module_id"] for f in failures] == [2], f"module 2 should fail, got {failures}"


def test_enrich_min_size_filter(monkeypatch):
    """规模不足的模块不参与富集(根本不应发起请求)。"""
    calls = {"n": 0}

    def counting_gp(genes, **kw):
        calls["n"] += 1
        return pd.DataFrame()

    monkeypatch.setattr(na, "fetch_enrichment_gprofiler", counting_gp)
    monkeypatch.setattr(na, "fetch_enrichment_enrichr", lambda genes, **kw: pd.DataFrame())
    modules = [{"module_id": 1, "genes": ["A", "B"], "size": 2, "avg_degree": 1.0}]
    long_df, failures, no_sig = na.enrich_modules(modules, organism="hsapiens",
                                                  min_module_size=5)
    assert calls["n"] == 0, "module below min size must not be enriched"
    assert long_df.empty and not failures and not no_sig


# ---------- g:Profiler / Enrichr 熔断器(标记文件由 conftest 沙箱隔离) ----------


def test_circuit_ladder_advance():
    """阶梯语义:常规失败保持第 1 档;冷却期内的失败(并发探测撞车)不推进;
    到期后的探测失败推进一档;冷却时长按阶梯加长且最后一档封顶。"""
    en._gp_circuit.trip(advance=False)
    en._gp_circuit.trip(advance=False)  # 常规失败不推进档位
    assert en._gp_circuit.read()["stage"] == 1
    en._gp_circuit.trip(advance=True)   # 冷却期内:探测撞车去抖,不推进
    assert en._gp_circuit.read()["stage"] == 1
    # 熔断冷却期内,重试提示应给出大致分钟数(向上取整,不足 1 分钟按 1 分钟)
    assert en.enrichment_retry_hint_minutes() == 1
    _expire_circuit(en._gp_circuit, 1)
    en._gp_circuit.trip(advance=True)   # 确实到期后的失败才推进
    assert en._gp_circuit.read()["stage"] == 2
    assert en._gp_circuit.state() == "open"
    cooldowns = [en._gp_circuit.cooldown(s) for s in (1, 2, 3, 4, 9)]
    assert cooldowns[0] < cooldowns[1] < cooldowns[2] <= cooldowns[3] == cooldowns[4], cooldowns


def test_circuit_probe_empty_result_resets(monkeypatch):
    """冷却到期后的探测拿到空结果(服务可用、恰好无显著项)也必须完全复位,不再保持熔断。"""
    en._gp_circuit.trip(advance=False)
    _expire_circuit(en._gp_circuit, 1)
    monkeypatch.setattr(en, "_profile_gprofiler", lambda *a, **k: pd.DataFrame())
    st: dict = {}
    df = en.fetch_enrichment_gprofiler(["TP53", "TNF"], status=st)
    assert df.empty and st.get("reason") == "no_significant", st
    assert en._gp_circuit.state() == "closed", \
        "empty probe result must reset the breaker"
    assert en.enrichment_retry_hint_minutes() is None  # 复位后无冷却信息


def test_circuit_probe_fail_advances(monkeypatch):
    """冷却到期后的探测失败:status 记 probe_failed,冷却阶梯推进一档;
    冷却期内后续请求直接跳过(circuit_open),不再发起。"""
    en._gp_circuit.trip(advance=False)
    _expire_circuit(en._gp_circuit, 1)

    def _down(*a, **k):
        raise ConnectionError("g:Profiler still down")

    monkeypatch.setattr(en, "_profile_gprofiler", _down)
    st: dict = {}
    df = en.fetch_enrichment_gprofiler(["TP53", "TNF"], status=st)
    assert df.empty and st.get("reason") == "probe_failed", st
    assert en._gp_circuit.read()["stage"] == 2
    assert en._gp_circuit.state() == "open"
    st2: dict = {}
    df2 = en.fetch_enrichment_gprofiler(["TP53", "TNF"], status=st2)
    assert df2.empty and st2.get("reason") == "circuit_open", st2


def test_enrichr_breaker_single_attempt(monkeypatch):
    """Enrichr 兜底源熔断:单次尝试失败即跳闸(请求内不重试);冷却期内跳过;探测成功复位。"""
    calls = {"n": 0}

    def _down(*a, **k):
        calls["n"] += 1
        raise ConnectionError("Enrichr down")

    monkeypatch.setattr(en, "_query_enrichr", _down)
    st: dict = {}
    df = en.fetch_enrichment_enrichr(["TP53", "TNF"],
                                     libs=["KEGG_2021_Human"], status=st)
    assert df.empty and st.get("reason") == "request_failed", st
    assert calls["n"] == 1, "fallback source must not retry in-request"
    assert en._enrichr_circuit.read()["stage"] == 1
    st2: dict = {}
    en.fetch_enrichment_enrichr(["TP53", "TNF"],
                                libs=["KEGG_2021_Human"], status=st2)
    assert st2.get("reason") == "circuit_open", st2
    assert calls["n"] == 1, "breaker open must skip the request"
    _expire_circuit(en._enrichr_circuit, 1)
    monkeypatch.setattr(en, "_query_enrichr", lambda *a, **k: pd.DataFrame())
    st3: dict = {}
    df3 = en.fetch_enrichment_enrichr(["TP53", "TNF"],
                                      libs=["KEGG_2021_Human"], status=st3)
    assert df3.empty and st3.get("reason") == "no_significant", st3
    assert en._enrichr_circuit.state() == "closed"


# ---------- GSEA prerank(真实 gseapy 离线运行) ----------


def test_gsea_build_ranked():
    """排序列表构建:基因符号统一大写、同基因取最高分、缺分基因剔除、按分数降序。"""
    df = pd.DataFrame({
        "approvedSymbol": ["tp53", "TP53", "EGFR", "BRCA1", "KRAS", "MT-CO1"],
        "score_pchembl": [5.0, 9.0, 7.0, None, 6.0, float("nan")],
    })
    rnk, stats = build_ranked_list(df, "approvedSymbol", "score_pchembl")
    # TP53 大小写重复取最高分;缺分的基因剔除 → 保留 3 个
    assert stats["n_input"] == 6 and stats["n_kept"] == 3, stats
    assert abs(stats["drop_ratio"] - 0.5) < 1e-9
    assert list(rnk["gene"]) == ["TP53", "EGFR", "KRAS"], rnk
    assert list(rnk["score"]) == [9.0, 7.0, 6.0], "must be sorted desc"


def test_gsea_normalize():
    """gseapy res2d 结果归一化:按 |NES| 降序、保留来源列;缺必需列时返回空表而非报错。"""
    # 按真实 gseapy 输出的列名构造测试数据
    res2d = pd.DataFrame([
        {"Name": "prerank", "Term": "SET_A", "ES": 0.6, "NES": 2.1,
         "NOM p-val": 0.01, "FDR q-val": 0.05, "FWER p-val": 0.02,
         "Tag %": "10/15", "Gene %": "33%", "Lead_genes": "G1;G2", "source": "LIB_X"},
        {"Name": "prerank", "Term": "SET_B", "ES": -0.5, "NES": -1.8,
         "NOM p-val": 0.08, "FDR q-val": 0.2, "FWER p-val": 0.1,
         "Tag %": "8/15", "Gene %": "27%", "Lead_genes": "G9;G8", "source": "LIB_X"},
    ])
    out = _normalize_enrichment_gsea(res2d)
    assert list(out["term_name"]) == ["SET_A", "SET_B"], "must be sorted by |NES| desc"
    assert out["source"].iloc[0] == "LIB_X"
    assert {"term_name", "es", "nes", "nominal_p_value", "fdr_q_value",
            "tag_pct", "gene_pct", "lead_genes", "source"} <= set(out.columns)
    # 缺少必需列时应返回空表而不是报错
    assert _normalize_enrichment_gsea(pd.DataFrame({"foo": [1]})).empty


def _tiny_gmt(directory: str) -> str:
    """40 基因、两组基因集的微型 GMT:HIGH_SET 取高分段 15 个,LOW_SET 取低分段 15 个。"""
    genes = [f"G{i:02d}" for i in range(40)]
    path = os.path.join(directory, "selfcheck_gsea_tiny.gmt")
    with open(path, "w", encoding="utf-8") as f:
        f.write("HIGH_SET\tna\t" + "\t".join(genes[:15]) + "\n")
        f.write("LOW_SET\tna\t" + "\t".join(genes[25:]) + "\n")
    return path


def _rnk40() -> pd.DataFrame:
    """40 基因、分数从 2.0 线性降到 -1.9 的排序列表。"""
    return pd.DataFrame({
        "gene": [f"G{i:02d}" for i in range(40)],
        "score": [2.0 - 0.1 * i for i in range(40)],
    })


def test_gsea_prerank_local(tmp_path):
    """用真实 gseapy 离线跑 GSEA:排位靠前的基因集方向应正确;第二次调用应命中缓存。"""
    libs = {"SELFTEST_LIB": _tiny_gmt(str(tmp_path))}
    label = "selfcheck_gsea"
    # 先删除旧测试缓存,避免直接命中
    _purge(os.path.join(config.ENRICHMENT_CACHE_DIR, f"GSEA_{label}_*.csv"))

    st1: dict = {}
    df1 = fetch_enrichment_gsea_prerank(_rnk40(), libs=libs, label=label, status=st1)
    assert st1.get("reason") == "success", st1
    assert not df1.empty and {"source", "term_name", "nes", "fdr_q_value"} <= set(df1.columns)
    high = df1[df1["term_name"] == "HIGH_SET"].iloc[0]
    low = df1[df1["term_name"] == "LOW_SET"].iloc[0]
    assert high["nes"] > 0 and low["nes"] < 0, "enrichment directions are wrong"

    st2: dict = {}
    df2 = fetch_enrichment_gsea_prerank(_rnk40(), libs=libs, label=label, status=st2)
    assert st2.get("reason") == "cache", st2
    # 经缓存读写后内容应不变(允许微小格式差异)
    pd.testing.assert_frame_equal(
        df1.reset_index(drop=True), df2.reset_index(drop=True),
        check_dtype=False, check_exact=False, rtol=1e-6,
    )

    try:
        cf = _gsea_cache_file(_rnk40(), libs, label)
        if os.path.exists(cf):
            os.remove(cf)
    except OSError:
        pass
    assert config.GSEA_MIN_RANKED_GENES <= 40, "test rnk must pass the gate"


def test_gsea_insufficient_no_cache(tmp_path):
    """排位基因数低于门槛:直接拒绝(insufficient_input),且失败不得留下缓存文件。"""
    st1: dict = {}
    df = fetch_enrichment_gsea_prerank(_rnk40().head(10), libs={"L": _tiny_gmt(str(tmp_path))},
                                       label="selfcheck_gsea_short", status=st1)
    assert df.empty and st1.get("reason") == "insufficient_input", st1
    left = glob.glob(os.path.join(config.ENRICHMENT_CACHE_DIR,
                                  "GSEA_selfcheck_gsea_short_*.csv"))
    assert not left, f"failure must not leave cache files: {left}"


# ---------- 药物-靶点-通路网络图 ----------


def test_drug_target_pathway_plot():
    """药物-靶点-通路网络图:KEGG 行参与构图、GO 行忽略;无 KEGG/无命中返回 None;
    g:Profiler 列表样式 / 分号 / 逗号三种基因格式都兼容;
    无药名与 drug_name=None 均按手动模式渲染。"""
    import matplotlib.pyplot as plt

    targets = ["TP53", "VEGFA", "TNF", "IL6", "AKT1", "STAT3", "CASP3", "MYC",
               "EGFR", "PTGS2"]
    enr = pd.DataFrame([
        {"source": "KEGG", "term_id": "hsa04151", "term_name": "PI3K-Akt signaling pathway",
         "p_value": 1e-4, "adjusted_p_value": 1e-4, "intersection_size": 4,
         "query_size": 10, "genes": "['TP53', 'AKT1', 'MYC', 'EGFR']"},  # g:Profiler 格式
        {"source": "KEGG", "term_id": "hsa04066", "term_name": "HIF-1 signaling pathway",
         "p_value": 2e-4, "adjusted_p_value": 2e-4, "intersection_size": 3,
         "query_size": 10, "genes": "VEGFA;STAT3;TP53"},                 # Enrichr 格式
        {"source": "KEGG", "term_id": "hsa04917", "term_name": "Insulin signaling pathway",
         "p_value": 3e-4, "adjusted_p_value": 3e-4, "intersection_size": 2,
         "query_size": 10, "genes": "IL6,PTGS2"},                        # 逗号分隔
        {"source": "GO:BP", "term_id": "GO:0000001", "term_name": "apoptotic process",
         "p_value": 1e-5, "adjusted_p_value": 1e-5, "intersection_size": 2,
         "query_size": 10, "genes": "TP53;TNF"},
    ])
    fig = render_drug_target_pathway_plot("metformin", targets, enr)
    assert fig is not None
    plt.close(fig)
    assert render_drug_target_pathway_plot(
        "metformin", targets, enr[enr["source"] == "GO:BP"]) is None
    bad = enr.iloc[[0]].copy()
    bad["genes"] = "ZZZ1;ZZZ2"  # 命中基因与输入靶点无交集
    assert render_drug_target_pathway_plot("metformin", targets, bad) is None
    for no_drug in ("", None):  # 手动模式(含 drug_name 缺省容错)照常出图
        fig2 = render_drug_target_pathway_plot(no_drug, targets, enr)
        assert fig2 is not None
        plt.close(fig2)
