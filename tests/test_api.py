"""WY.api 的冒烟测试(FastAPI TestClient,进程内,不占端口):50 个端点逐个「路由通」。

冒烟口径(不是端到端,业务语义的唯一定义在 core,由各自模块测试覆盖):
- 读端点:断言 200 + 响应形状;
- 写端点:断言校验防线(400/404/穿越拦截);
- job 端点:代表性走完「提交 → 轮询 → done/error」管线,任务体内的 core 重活
  一律打桩(打在任务函数实际解析的命名空间:api 模块命名空间或 core 源模块);
- 安全红线单列:路径穿越(../)与越权名必须在 400/404 被堵住。

不触碰真实落盘:conftest 沙箱接管全部缓存/下载/日志目录。唯一例外是
save_api_env 会写 <repo>/WY/API.env(config.BASE_DIR 不在沙箱清单),该端点测试打桩替身。
"""
import os
import time
from types import SimpleNamespace

import pandas as pd
import pytest

from WY import config


@pytest.fixture()
def client():
    from fastapi.testclient import TestClient

    from WY.api import app

    return TestClient(app)


def _wait_job(client, job_id, timeout=15.0):
    """轮询 job 到终态(线程在 daemon 池里跑,轮询带超时防挂死)。"""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        r = client.get(f"/api/jobs/{job_id}")
        assert r.status_code == 200
        body = r.json()
        if body["status"] in ("done", "error"):
            return body
        time.sleep(0.02)
    raise AssertionError(f"job {job_id} 未在 {timeout}s 内结束")


# ======== 六类集合管理页:清单 / 打包 / 删除 ========

COLLECTION_KINDS = ("drug_sources", "disease_sources", "ppi", "enrichment",
                    "stars", "downloads")


@pytest.mark.parametrize("kind", COLLECTION_KINDS)
def test_collections_all_kinds(client, kind):
    r = client.get(f"/api/collections/{kind}")
    assert r.status_code == 200
    body = r.json()
    assert body["kind"] == kind and body["mode"] in ("files", "dirs")
    assert isinstance(body["items"], list)


def test_collections_unknown_kind_404(client):
    assert client.get("/api/collections/nope").status_code == 404


def test_collection_delete_rejects_unknown_and_traversal(client):
    for bad in (["../evil"], ["ghost.csv"], [".", ".."]):
        r = client.post("/api/collections/ppi/delete", json={"names": bad})
        assert r.status_code == 400, f"{bad} 应被越权校验拦下"


def test_collection_zip_rejects_empty(client):
    assert client.post("/api/collections/ppi/zip",
                       json={"names": []}).status_code == 400


def test_collection_zip_delete_roundtrip(client, monkeypatch):
    """空沙箱里造一个 PPI 缓存文件 → 清单可见 → 打包出 zip → 删除生效。

    清单走 CACHE_DIR/ppi、路径解析走 PPI_CACHE_DIR(生产里两者是同一目录,
    沙箱逐属性独立指临时目录后才会分叉),测试里把后者对齐前者保持一致。"""
    ppi_dir = os.path.join(config.CACHE_DIR, "ppi")
    monkeypatch.setattr(config, "PPI_CACHE_DIR", ppi_dir)
    os.makedirs(ppi_dir, exist_ok=True)
    p = os.path.join(ppi_dir, "aspirin_cancer.csv")
    with open(p, "w", encoding="utf-8") as f:
        f.write("gene,score\nTP53,1\n")

    items = client.get("/api/collections/ppi").json()["items"]
    assert [it["name"] for it in items] == ["aspirin_cancer.csv"]

    z = client.post("/api/collections/ppi/zip", json={"names": ["aspirin_cancer.csv"]})
    assert z.status_code == 200 and z.content[:2] == b"PK"

    d = client.post("/api/collections/ppi/delete", json={"names": ["aspirin_cancer.csv"]})
    assert d.json() == {"deleted": 1, "failed": []}
    assert not os.path.exists(p)


# ======== 来源刷新 / 背景互作组 / 任务轮询 ========

def test_sources_refresh_unknown_kind_404(client):
    assert client.post("/api/sources/geo/refresh",
                       json={"filenames": None}).status_code == 404


def test_sources_refresh_job_flow(client, monkeypatch):
    monkeypatch.setattr("WY.api.update_all_cache_files",
                        lambda *a, **k: pd.DataFrame())
    r = client.post("/api/sources/drug/refresh", json={"filenames": None})
    assert r.status_code == 200
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done"


def test_interactome_status_payload(client):
    """状态端点按真实状态返回对应形态:默认沙箱必为未构建;设
    NPAP_TEST_REAL_INTERACTOME=1 且本机已构建时为已构建形态(字段齐全)。"""
    body = client.get("/api/interactome").json()
    if not body["built"]:
        assert body["version"] == config.STRING_VERSION
    else:
        for key in ("nodes", "edges", "version", "confidence",
                    "built_at", "age", "stale"):
            assert key in body, f"已构建形态缺字段 {key}"
        assert body["nodes"] >= 10000 and body["edges"] >= 0


def test_interactome_build_job_flow(client, monkeypatch):
    monkeypatch.setattr("WY.api.build_interactome",
                        lambda progress_cb=None: {"n_nodes": 2, "n_edges": 1})
    monkeypatch.setattr("WY.api.clean_old_version_caches", lambda name: 0)
    r = client.post("/api/interactome/build")
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done" and job["result"]["nodes"] == 2


def test_jobs_unknown_id_404(client):
    assert client.get("/api/jobs/no-such-job").status_code == 404


# ======== 日志页:清单 / 内容尾部 / 级别切换(含穿越拦截) ========

def test_logs_files_empty(client):
    assert client.get("/api/logs/files").json() == {"files": []}


def test_logs_content_roundtrip_traversal_and_truncation(client):
    os.makedirs(config.LOG_DIR, exist_ok=True)
    log = os.path.join(config.LOG_DIR, "app.log")
    with open(log, "w", encoding="utf-8", newline="") as f:
        f.write("hello log\n")

    body = client.get("/api/logs/content", params={"name": "app.log"}).json()
    assert body["content"] == "hello log\n" and body["truncated"] is False

    with open(log, "w", encoding="utf-8", newline="") as f:
        f.write("x" * 3000)
    body = client.get("/api/logs/content",
                      params={"name": "app.log", "max_bytes": 1000}).json()
    assert body["truncated"] is True and len(body["content"]) <= 1000

    assert client.get("/api/logs/content",
                      params={"name": "../secrets"}).status_code == 400
    assert client.get("/api/logs/content",
                      params={"name": "ghost.log"}).status_code == 404


def test_logs_level_roundtrip_and_bad_level(client):
    original = client.get("/api/logs/level").json()["level"]
    try:
        assert client.post("/api/logs/level",
                           json={"level": "debug"}).json()["level"] == "DEBUG"
        assert client.get("/api/logs/level").json()["level"] == "DEBUG"
        assert client.post("/api/logs/level",
                           json={"level": "WARN"}).status_code == 400
    finally:
        client.post("/api/logs/level", json={"level": original})


# ======== PPI 分析:校验防线 / job 全流程(分析器打桩)/ 产物伺服 ========

def test_ppi_analyze_validation(client):
    assert client.post("/api/ppi/analyze",
                       json={"mode": "wat"}).status_code == 400
    assert client.post("/api/ppi/analyze",
                       json={"mode": "names", "drug": "a"}).status_code == 400


def _fake_ppi_result():
    d = os.path.join(config.DOWNLOADS_DIR, "fake_data")
    os.makedirs(d, exist_ok=True)
    empty = pd.DataFrame()
    return {"prefix": "fake", "saved": {"directory": d}, "params": (),
            "genes": [], "metrics_df": empty, "edges_df": empty, "hub_df": empty,
            "hub_rule": "degree", "hub_method_eff": "degree",
            "string_fig": None, "cyto_fig": None, "expansion_df": None,
            "gsea": None, "rwr": None, "modules": None,
            "enrichment_source": "", "enrichment_reason": "", "enrichment_df": None}


def test_ppi_analyze_job_flow(client, monkeypatch):
    monkeypatch.setattr("WY.api.run_ppi_analysis", lambda *a, **k: _fake_ppi_result())
    r = client.post("/api/ppi/analyze",
                    json={"mode": "names", "drug": "aspirin", "disease": "cancer"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done"
    body = job["result"]
    assert body["prefix"] == "fake" and body["nodes"] == 0 and "artifacts" in body


def test_artifacts_traversal_and_missing_404(client):
    assert client.get("/api/artifacts", params={"rel": "../API.env"}).status_code == 404
    assert client.get("/api/artifacts", params={"rel": "ghost.png"}).status_code == 404


def test_artifacts_serves_downloads_file(client):
    d = os.path.join(config.DOWNLOADS_DIR, "demo_data")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "plot.png"), "wb") as f:
        f.write(b"\x89PNG fake")
    r = client.get("/api/artifacts", params={"rel": "demo_data/plot.png"})
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"


def test_ppi_rerun_modules_no_dir_404(client):
    assert client.post("/api/ppi/rerun-modules",
                       json={"directory_rel": "../outside", "prefix": "x"}).status_code == 404


def test_ppi_rerun_modules_job_error_terminal(client):
    """目录存在但缺产物 CSV:任务体内抛错 → job 走到 error 终态(错误管道的冒烟)。"""
    d = os.path.join(config.DOWNLOADS_DIR, "broken_data")
    os.makedirs(d, exist_ok=True)
    r = client.post("/api/ppi/rerun-modules", json={"directory_rel": "broken_data",
                                                    "prefix": "broken"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "error" and "FileNotFoundError" in job["error"]


def test_ppi_rerun_gsea_validation_and_flow(client, monkeypatch):
    assert client.post("/api/ppi/rerun-gsea",
                       json={"side": "wat", "drug": "a", "disease": "b",
                             "prefix": "p"}).status_code == 400
    monkeypatch.setattr("WY.core.analysis.gsea_prerank_side",
                        lambda *a, **k: (None, {"nes": 0.0}, "skipped"))
    r = client.post("/api/ppi/rerun-gsea",
                    json={"side": "drug", "drug": "aspirin", "disease": "cancer",
                          "prefix": "p"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done" and job["result"]["reason"] == "skipped"


# ======== 主页统计与设置 ========

def test_stats_shape(client):
    body = client.get("/api/stats").json()
    assert set(body) == {"databases", "compounds", "targets", "pathways"}
    assert all(isinstance(v, int) for v in body.values())


def test_settings_snapshot_shape(client):
    body = client.get("/api/settings").json()
    for key in ("keys", "vina_dir", "string_version", "known_versions",
                "chembl_blacklist", "pocket_cutoff"):
        assert key in body, f"设置快照缺 {key}"
    assert set(body["keys"]) == {"disgenet", "ncbi"}  # 密钥只回布尔不回值


def test_settings_string_version_same_and_bad_format(client):
    r = client.post("/api/settings/string-version",
                    json={"version": config.STRING_VERSION})
    assert r.status_code == 200 and r.json()["changed"] is False
    assert client.post("/api/settings/string-version",
                       json={"version": "bogus!!"}).status_code == 400


def test_settings_ui_validation_and_noop(client):
    current = client.get("/api/settings").json()["pocket_cutoff"]
    assert client.post("/api/settings/ui", json={"changes": {}}).status_code == 400
    assert client.post("/api/settings/ui",
                       json={"changes": {"NOPE": 1}}).status_code == 400
    assert client.post("/api/settings/ui",
                       json={"changes": {"POCKET_CUTOFF": 99}}).status_code == 400
    r = client.post("/api/settings/ui", json={"changes": {"POCKET_CUTOFF": current}})
    assert r.status_code == 200 and r.json()["changed"] == []
    assert r.json()["snapshot"]["pocket_cutoff"] == current


def test_settings_ui_blacklist_normalizes_and_restores(client):
    original = list(config.CHEMBL_TARGET_BLACKLIST)
    try:
        r = client.post("/api/settings/ui",
                        json={"changes": {"CHEMBL_TARGET_BLACKLIST": [" tp53 "]}})
        assert r.status_code == 200
        assert r.json()["changed"] == ["CHEMBL_TARGET_BLACKLIST"]
        assert "TP53" in client.get("/api/settings").json()["chembl_blacklist"]
        assert client.post("/api/settings/ui",
                           json={"changes": {"CHEMBL_TARGET_BLACKLIST": [1]}}).status_code == 400
    finally:
        client.post("/api/settings/ui",
                    json={"changes": {"CHEMBL_TARGET_BLACKLIST": original}})


def test_settings_api_env_stubbed(client, monkeypatch):
    """save_api_env 会写真实 WY/API.env(BASE_DIR 不在沙箱),冒烟打桩验证 HTTP 形态。"""
    monkeypatch.setattr("WY.api.save_api_env", lambda changes: ["NCBI_API_KEY"])
    r = client.post("/api/settings/api-env",
                    json={"changes": {"NCBI_API_KEY": "x"}})
    assert r.status_code == 200 and r.json() == {"applied": ["NCBI_API_KEY"]}


# ======== 主页四库检索:job + 服务端 venn ========

def test_targets_retrieve_validation(client):
    assert client.post("/api/targets/retrieve",
                       json={"drug": " ", "disease": "cancer"}).status_code == 400


def test_targets_retrieve_and_venn(client, monkeypatch):
    # 假表带 source_count 列(2026-10 合并 master 后,端点对交集表做
    # drop(columns=["source_count"]) 无条件消费,假数据必须满足新契约)
    drug_df = pd.DataFrame({"approvedSymbol": ["TP53", "BRCA1"],
                            "score": [0.9, 0.8], "source_count": [2, 1]})
    disease_df = pd.DataFrame({"approvedSymbol": ["TP53", "EGFR"],
                               "score": [0.7, 0.6], "source_count": [3, 1]})
    monkeypatch.setattr("WY.core.targets.get_targets_concurrent",
                        lambda *a, **k: (drug_df, disease_df, []))
    r = client.post("/api/targets/retrieve",
                    json={"drug": "aspirin", "disease": "cancer"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done"
    body = job["result"]
    assert (body["drug_count"], body["disease_count"], body["common_count"]) == (2, 2, 1)
    assert body["common_genes"] == ["TP53"] and body["failed"] == []
    # 行载荷契约:三张行表齐全,交集行已去掉 source_count 冗余列
    assert len(body["drug_rows"]) == 2 and len(body["disease_rows"]) == 2
    assert [row["approvedSymbol"] for row in body["common_rows"]] == ["TP53"]
    assert "source_count" not in body["common_rows"][0]

    venn = client.get("/api/targets/venn", params={"job_id": r.json()["job_id"]})
    assert venn.status_code == 200
    assert venn.headers["content-type"] == "image/png"
    assert venn.content[:4] == b"\x89PNG"
    # 暗色版按需渲染并缓存:同 job 二次请求不再画图,但同样出 PNG
    dark = client.get("/api/targets/venn",
                      params={"job_id": r.json()["job_id"], "theme": "dark"})
    assert dark.status_code == 200 and dark.content[:4] == b"\x89PNG"


def test_targets_venn_unknown_job_404(client):
    assert client.get("/api/targets/venn",
                      params={"job_id": "ghost"}).status_code == 404


# ======== 分子对接向导:会话 + 各步防线 ========

@pytest.fixture()
def sid(client):
    return client.post("/api/docking/session").json()["sid"]


def test_docking_session_create_and_expired_404(client):
    assert client.post("/api/docking/session").json()["sid"]
    r = client.post("/api/docking/deadbeef/source",
                    json={"mode": "manual", "manual_genes": "TP53"})
    assert r.status_code == 404


def test_docking_sources_empty(client):
    body = client.get("/api/docking/sources").json()
    assert body == {"ppi_dirs": [], "enr_files": []}


def test_docking_source_manual_and_bad_mode(client, sid):
    assert client.post(f"/api/docking/{sid}/source",
                       json={"mode": "wat"}).status_code == 400
    body = client.post(f"/api/docking/{sid}/source",
                       json={"mode": "manual",
                             "manual_genes": "tp53, brca1;;EGFR"}).json()
    assert body["candidates"] == ["BRCA1", "EGFR", "TP53"]  # 大写归一 + 分隔符切分
    assert body["label"] == "" and body["hub_count"] == 0


def test_docking_check_structures_guards(client, sid):
    fresh = client.post("/api/docking/session").json()["sid"]
    assert client.post(f"/api/docking/{fresh}/check-structures",
                       json={"use_struct": True}).status_code == 400
    # manual 模式不产生候选集(candidate_pre 为空):结构检查必须 400 拒绝而非空跑
    client.post(f"/api/docking/{sid}/source",
                json={"mode": "manual", "manual_genes": "TP53"})
    assert client.post(f"/api/docking/{sid}/check-structures",
                       json={"use_struct": False}).status_code == 400


def test_docking_source_enr_and_check_structures_direct(client, monkeypatch):
    """enr 来源:富集文件 → 候选基因;跳过结构检查(use_struct=False)走 direct 分支。"""
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    with open(os.path.join(config.ENRICHMENT_CACHE_DIR,
                           "gProfiler_aspirin_cancer.csv"), "w",
              encoding="utf-8") as f:
        f.write("term_id,genes\nhsa04115,TP53;BRCA1\n")
    sid = client.post("/api/docking/session").json()["sid"]
    src = client.post(f"/api/docking/{sid}/source",
                      json={"mode": "enr", "name": "gProfiler_aspirin_cancer.csv"})
    assert src.status_code == 200
    assert src.json()["candidates"] == ["BRCA1", "TP53"]
    body = client.post(f"/api/docking/{sid}/check-structures",
                       json={"use_struct": False}).json()
    assert body["job_id"] is None
    assert body["direct"]["genes"] == ["BRCA1", "TP53"]
    assert body["direct"]["structures"] == {}


def test_docking_ligand_guards(client, sid):
    body = client.post(f"/api/docking/{sid}/ligand",
                       json={"source": "pubchem", "query": " "}).json()
    assert body["ok"] is False
    body = client.post(f"/api/docking/{sid}/ligand",
                       json={"source": "smiles", "smiles": " "}).json()
    assert body["ok"] is False
    assert client.post(f"/api/docking/{sid}/ligand",
                       json={"source": "wat"}).status_code == 400
    # SMILES 路径真实执行(RDKit 本地;obabel 缺失时按契约返回 ok=False 而非 500)
    body = client.post(f"/api/docking/{sid}/ligand",
                       json={"source": "smiles", "smiles": "CCO"}).json()
    assert "ok" in body


def test_docking_receptor_candidates_guard_and_flow(client, sid, monkeypatch):
    assert client.post(f"/api/docking/{sid}/receptor/candidates",
                       json={"gene": "EGFR"}).status_code == 400
    client.post(f"/api/docking/{sid}/source",
                json={"mode": "manual", "manual_genes": "TP53"})
    import WY.core.docking as docking

    monkeypatch.setattr(docking, "fetch_receptor_pdb_candidates",
                        lambda *a, **k: pd.DataFrame())
    _pdb_cached = __import__("WY.api", fromlist=["_pdb_candidates_cached"])
    _pdb_cached._pdb_candidates_cached.cache_clear()
    body = client.post(f"/api/docking/{sid}/receptor/candidates",
                       json={"gene": "TP53"}).json()
    assert body["rows"] == [] and isinstance(body["af_enabled"], bool)


def test_docking_downstream_guards_without_receptor(client, sid):
    """受体未就绪时,向导后续每一步都必须 400 拒绝(而非 500)。"""
    client.post(f"/api/docking/{sid}/source",
                json={"mode": "manual", "manual_genes": "TP53"})
    assert client.post(f"/api/docking/{sid}/box",
                       json={"infer": False, "center": [0, 0, 0],
                             "size": [1, 1, 1]}).status_code == 400
    assert client.post(f"/api/docking/{sid}/box",
                       json={"infer": False}).status_code == 400
    assert client.post(f"/api/docking/{sid}/receptor/prepare",
                       json={"gene": "EGFR", "pdb_id": "1CRN"}).status_code == 400
    assert client.get(f"/api/docking/{sid}/box-3d").status_code == 400
    assert client.post(f"/api/docking/{sid}/run", json={}).status_code == 400
    assert client.post(f"/api/docking/{sid}/save",
                       json={"pose_mode": 1}).status_code == 400
    assert client.post(f"/api/docking/{sid}/star",
                       json={"on": True}).status_code == 400
    assert client.get(f"/api/docking/{sid}/view3d",
                      params={"pose": 1}).status_code == 400
    assert client.get(f"/api/docking/{sid}/interactions",
                      params={"pose": 1}).status_code == 400
    assert client.post(f"/api/docking/{sid}/save3d",
                       json={"pose": 1}).status_code == 400
    assert client.get(f"/api/docking/{sid}/zip").status_code == 400


# ======== GEO 表达谱 + ML 筛选页 ========

def test_ml_diseases_empty(client):
    assert client.get("/api/ml/diseases").json() == {"diseases": [], "drugs": []}


def test_ml_search_validation_and_stub(client, monkeypatch):
    assert client.post("/api/ml/search", json={"disease": " "}).status_code == 400
    monkeypatch.setattr("WY.core.geo.search_geo_datasets",
                        lambda *a, **k: pd.DataFrame())
    body = client.post("/api/ml/search", json={"disease": "cancer"}).json()
    assert body == {"results": [], "disease": "cancer"}


def test_ml_probe_validation_and_metadata_path(client, monkeypatch):
    assert client.post("/api/ml/probe", json={"gse": " "}).status_code == 400
    monkeypatch.setattr("WY.core.geo.probe_series_matrix",
                        lambda gse: {"size_bytes": 10})
    body = client.post("/api/ml/probe", json={"gse": "gsex"}).json()
    assert body["gse"] == "GSEX" and body["err"] == ""
    assert body["assess"]["matrix_empty"] is True
    assert isinstance(body["assess"]["counts_available"], bool)


def test_ml_overview_job_flow(client, monkeypatch):
    fake_ds = SimpleNamespace(issues=[], table=pd.DataFrame(),
                              samples=pd.DataFrame({"title": ["s1"]}))
    monkeypatch.setattr("WY.core.geo.fetch_geo_dataset", lambda gse, **k: fake_ds)
    monkeypatch.setattr("WY.core.geo.group_fields_overview",
                        lambda ds: pd.DataFrame())
    monkeypatch.setattr("WY.core.expression.assess_expression_dataset",
                        lambda ds: {"usable": True})
    r = client.post("/api/ml/overview", json={"gse": "gsex"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done"
    assert job["result"]["overview"] == [] and job["result"]["issues"] == []


def test_ml_ranked_empty(client):
    assert client.post("/api/ml/ranked", json={}).json() == {"files": []}


def test_ml_run_validation(client):
    assert client.post("/api/ml/run", json={"mode": "wat"}).status_code == 400
    assert client.post("/api/ml/run", json={"mode": "gse", "gse": "G"}).status_code == 400
    assert client.post("/api/ml/run", json={"mode": "local"}).status_code == 400
    assert client.post("/api/ml/run",
                       json={"mode": "local", "expr": "a.csv", "labels": "b.csv",
                             "ranked_labels": ["ghost"]}).status_code == 400


def test_ml_run_local_job_flow(client, monkeypatch):
    # 生产中 GEO 结果目录在 downloads 之下(_rel_downloads 以它为根);沙箱里两者
    # 是兄弟目录,测试里对齐父子关系,顺带验证产物相对路径映射。
    geo_root = os.path.join(config.DOWNLOADS_DIR, "GEO_results")
    monkeypatch.setattr(config, "GEO_RESULT_DIR", geo_root)
    out_dir = os.path.join(geo_root, "GSEX_20260101_000000")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        f.write("{}")

    def fake_chain(**kwargs):
        return {"gse": None, "tag": "GSEX", "out_dir": out_dir,
                "summary": {"elapsed_total_s": 1.0, "geo_cache_fallback": False,
                            "geo_cache_date": ""},
                "meta": {},
                "funnel": pd.DataFrame([{"stage": "s", "n": 1, "note": "",
                                         "dropped_genes": ""}]),
                # consensus 以 gene 为索引(API 侧 rename_axis+reset_index 还原成列)
                "consensus": pd.DataFrame({"consensus_score": [1.0]},
                                          index=pd.Index(["TP53"], name="gene")),
                "shortlist": pd.DataFrame([{"gene": "TP53", "in_shortlist": True}]),
                "deg": pd.DataFrame([{"gene": "TP53"}]), "roc_png": ""}

    monkeypatch.setattr("WY.api.run_screen_chain", fake_chain)
    r = client.post("/api/ml/run",
                    json={"mode": "local", "expr": "expr.csv", "labels": "labels.csv"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done"
    body = job["result"]
    assert body["tag"] == "GSEX" and body["consensus_total"] == 1
    assert body["shortlist_count"] == 1
    assert "summary.json" in body["artifacts"] and body["artifacts"]["summary.json"]


def test_ml_history_and_runs_empty(client):
    assert client.get("/api/ml/history").json() == {"runs": []}
    assert client.get("/api/ml/runs").json() == {"runs": []}


def test_ml_runs_delete_rejects_outside_root(client):
    r = client.post("/api/ml/runs/delete", json={"paths": ["C:/Windows/temp"]})
    assert r.status_code == 200 and r.json()["deleted"] == []


def test_ml_runs_zip_empty_ok(client):
    r = client.post("/api/ml/runs/zip", json={"paths": []})
    assert r.status_code == 200 and r.headers["content-type"] == "application/zip"


def test_ml_runs_cleanup_empty(client):
    assert client.post("/api/ml/runs/cleanup",
                       json={"keep": 3}).json() == {"deleted": []}


def test_ml_force_refresh_flow(client, monkeypatch):
    assert client.post("/api/ml/force-refresh", json={"gse": " "}).status_code == 400
    monkeypatch.setattr("WY.core.geo.fetch_geo_dataset", lambda gse, **k: None)
    r = client.post("/api/ml/force-refresh", json={"gse": "gsex"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done" and job["result"] == {"ok": True}


# ======== 前端静态托管(SPA fallback) ========

def test_spa_root_and_fallback(client):
    r = client.get("/")
    assert r.status_code == 200
    text = r.text
    if text.lstrip().startswith("<"):
        assert "html" in text[:200].lower()
        deep = client.get("/data/network")
        assert deep.status_code == 200 and deep.text.lstrip().startswith("<")
    else:
        # dist 未构建的兜底分支
        assert "前端未构建" in text
