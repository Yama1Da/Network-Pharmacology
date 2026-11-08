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
    monkeypatch.setattr("WY.api.sources.update_all_cache_files",
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
    monkeypatch.setattr("WY.api.interactome.build_interactome",
                        lambda progress_cb=None: {"n_nodes": 2, "n_edges": 1})
    monkeypatch.setattr("WY.api.interactome.clean_old_version_caches", lambda name: 0)
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
    monkeypatch.setattr("WY.api.ppi.run_ppi_analysis", lambda *a, **k: _fake_ppi_result())
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
                "chembl_blacklist", "pocket_cutoff", "ui_lang"):
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


def test_settings_ui_lang_normalizes_and_restores(client):
    """UI_LANG:合法值归一小写并进快照与落盘,非法值 400;测后恢复原值。"""
    original = config.UI_LANG
    try:
        r = client.post("/api/settings/ui", json={"changes": {"UI_LANG": " EN "}})
        assert r.status_code == 200
        assert r.json()["changed"] == ["UI_LANG"]
        assert client.get("/api/settings").json()["ui_lang"] == "en"
        assert client.post("/api/settings/ui",
                           json={"changes": {"UI_LANG": "fr"}}).status_code == 400
    finally:
        client.post("/api/settings/ui", json={"changes": {"UI_LANG": original}})


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
    monkeypatch.setattr("WY.api.settings.save_api_env", lambda changes: ["NCBI_API_KEY"])
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
    assert body == {"ppi_dirs": [], "enr_files": [], "manual_dirs": []}


def test_docking_source_manual_and_bad_mode(client, sid):
    assert client.post(f"/api/docking/{sid}/source",
                       json={"mode": "wat"}).status_code == 400
    body = client.post(f"/api/docking/{sid}/source",
                       json={"mode": "manual",
                             "manual_genes": "tp53, brca1;;EGFR"}).json()
    assert body["candidates"] == ["BRCA1", "EGFR", "TP53"]  # 大写归一 + 分隔符切分
    assert body["label"] == "" and body["hub_count"] == 0


def test_docking_manual_dirs_listing_and_source_filter(client, sid):
    """manual_* 桶进专属清单(基因从 analysis_meta.json 读出),PPI/富集清单不再收录。"""
    import json as _json

    bucket = os.path.join(config.DOWNLOADS_DIR, "manual_abc12345_data")
    os.makedirs(bucket, exist_ok=True)
    with open(os.path.join(bucket, "analysis_meta.json"), "w", encoding="utf-8") as f:
        _json.dump({"genes": ["TP53", "EGFR", "APOE"]}, f)
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    with open(os.path.join(config.ENRICHMENT_CACHE_DIR,
                           "gProfiler_manual_abc12345.csv"), "w", encoding="utf-8") as f:
        f.write("term_id,genes\nhsa04115,TP53\n")

    body = client.get("/api/docking/sources").json()
    assert [m["name"] for m in body["manual_dirs"]] == ["manual_abc12345_data"]
    m = body["manual_dirs"][0]
    assert m["n_genes"] == 3 and m["genes"] == ["APOE", "EGFR", "TP53"]
    assert all(not d["name"].startswith("manual") for d in body["ppi_dirs"])
    assert all(not f["name"].startswith("gProfiler_manual") for f in body["enr_files"])

    # manual-file 模式:从桶元数据回填基因(与 manual 输入同构落 dock)
    r = client.post(f"/api/docking/{sid}/source",
                    json={"mode": "manual-file", "manual_file": "manual_abc12345_data"})
    assert r.status_code == 200
    assert r.json()["candidates"] == ["APOE", "EGFR", "TP53"]
    # 非法路径拒绝
    assert client.post(f"/api/docking/{sid}/source",
                       json={"mode": "manual-file", "manual_file": "../evil"}).status_code == 400
    assert client.post(f"/api/docking/{sid}/source",
                       json={"mode": "manual-file", "manual_file": "ppi_x_data"}).status_code == 400


def test_docking_check_structures_guards(client, sid):
    """未选任何来源时结构检查 400;手动模式已修通:候选即输入基因,直检走 direct 分支。"""
    fresh = client.post("/api/docking/session").json()["sid"]
    assert client.post(f"/api/docking/{fresh}/check-structures",
                       json={"use_struct": True}).status_code == 400
    client.post(f"/api/docking/{sid}/source",
                json={"mode": "manual", "manual_genes": "TP53"})
    body = client.post(f"/api/docking/{sid}/check-structures",
                       json={"use_struct": False}).json()
    assert body["job_id"] is None and body["direct"]["genes"] == ["TP53"]


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


def test_docking_ligand_manual_name_query_blocked(client, sid):
    """手动会话按药名盲查被拦(曾把 'manual' 直送 PubChem);CID 直查放行。"""
    client.post(f"/api/docking/{sid}/source",
                json={"mode": "manual", "manual_genes": "TP53"})
    body = client.post(f"/api/docking/{sid}/ligand",
                       json={"source": "pubchem", "query": "manual"}).json()
    assert body["ok"] is False and "CID" in body["error"]
    # 纯数字 CID 经 resolve_cid 直通下载层,不触名称解析
    from WY.core.docking import resolve_cid
    assert resolve_cid("2244") == 2244  # 纯数字直通,不触网


def test_docking_source_name_traversal_rejected(client, sid):
    """回归测试(P1-C):ppi/enr 来源名不得携带路径分隔符或逃逸缓存根。

    缺陷:mode=ppi 的 name 直传 read_ppi_dataset(os.path.join 无包含校验),
    mode=enr 的 name 直拼 ENRICHMENT_CACHE_DIR 读 CSV——`../` 或盘符名可读
    缓存根之外的文件内容(以 hub/enrichment 基因行集形式回显);对照
    manual-file 分支早有 basename+形态校验,三模式口径不一。
    复现要点:`name="../WY/API.env"` 与 `name="C:/Windows/win.ini"` 在修复前
    均 400 之外的状态(ppi 以空集合 200 成功建会话,enr 读不到则 500);
    修复后统一 400。合法纯名保绿(数据集不存在按空集合处理)。
    """
    for mode, name in (("ppi", "../WY/API.env"), ("ppi", "C:/Windows/win.ini"),
                       ("enr", "../WY/API.env"), ("enr", "C:/Windows/win.ini"),
                       ("ppi", "sub/dir_data"), ("enr", "a/b.csv")):
        r = client.post(f"/api/docking/{sid}/source", json={"mode": mode, "name": name})
        assert r.status_code == 400, f"{mode}/{name} 必须 400,实际 {r.status_code}"
    # 保绿:合法纯名(数据集不存在,read_ppi_dataset 按空集合处理)
    assert client.post(f"/api/docking/{sid}/source",
                       json={"mode": "ppi", "name": "aspirin_cancer_data"}
                       ).status_code == 200


def test_docking_manual_guard_survives_source_switch(client, sid, monkeypatch):
    """回归测试(P1-D):先选 ppi 再切 manual,防呆判据不得被残留键绕过。

    缺陷:ppi/enr 分支写 dock["ppi_dir"]/["enr_file"],切 manual 时不清除,
    reset_downstream 的失效清单也不含它们;防呆判据从「未挂 ppi/enr」反推,
    残留使 is_manual_session 恒为 False——手动基因清单(如基因名被当药名)
    可直送 PubChem 盲查(注释记载的事故形态)。
    复现要点:mode=ppi(数据集不存在,read_ppi_dataset 按空集合处理)→ 切
    manual → 非数字药名查询。修复前走 resolve_cid 盲查(桩替记录调用并返回
    None,响应为「查无化合物」);修复后判据只看显式 source_mode,盲查不被
    调用,返回手动模式拦截文案;manual→ppi 反向切换后药名查询恢复放行。
    """
    import WY.core.docking as docking

    calls: list[str] = []
    monkeypatch.setattr(docking, "resolve_cid",
                        lambda query: (calls.append(query), None)[1])

    # 先选 ppi 来源,再切 manual:残留 ppi_dir 曾使防呆判据失效
    assert client.post(f"/api/docking/{sid}/source",
                       json={"mode": "ppi", "name": "aspirin_cancer_data"}
                       ).status_code == 200
    assert client.post(f"/api/docking/{sid}/source",
                       json={"mode": "manual", "manual_genes": "TP53"}
                       ).status_code == 200
    body = client.post(f"/api/docking/{sid}/ligand",
                       json={"source": "pubchem", "query": "aspirin"}).json()
    assert calls == [], f"不得发生 PubChem 盲查,实际调用了 {calls}"
    assert body["ok"] is False and "CID" in body["error"], body["error"]

    # 反向切换 manual → ppi:显式来源后,药名查询应恢复放行(走盲查路径)
    assert client.post(f"/api/docking/{sid}/source",
                       json={"mode": "ppi", "name": "aspirin_cancer_data"}
                       ).status_code == 200
    body = client.post(f"/api/docking/{sid}/ligand",
                       json={"source": "pubchem", "query": "aspirin"}).json()
    assert calls == ["aspirin"], f"ppi 来源的药名查询应放行解析: calls={calls}"
    assert "查无化合物" in body["error"], body["error"]  # 桩返回 None 的固定文案


def test_docking_receptor_prepare_stale_discard(client, sid, monkeypatch):
    """回归测试(P0-B 后端兜底):受体准备任务写回前核对会话基因,被换即放弃写入。

    缺陷:prepare 任务闭包直接写共享会话 dock["rec_pdbqt"]/["pdb_id"],无代际
    核对——提交与写回之间基因被换(另一标签页/竞态),旧基因的受体会落进
    新基因的会话,canRun 满足后 vina 对错误靶点运行。
    复现要点:manual 会话种子 gene=TP53 → 提交 prepare(桩 download 2s)→
    立即把会话 gene 改为 EGFR → 任务终态 result.stale=True 且会话无受体写入。
    """
    import WY.core.docking as docking
    from WY.core import docking_session

    client.post(f"/api/docking/{sid}/source",
                json={"mode": "manual", "manual_genes": "TP53,EGFR"})
    dock = docking_session.get(sid)
    dock["gene"] = "TP53"
    dock["genes"] = ["TP53", "EGFR"]

    def slow_download(pdb_id):
        time.sleep(2)
        return f"C:/fake/{pdb_id}.pdb"

    monkeypatch.setattr(docking, "download_pdb", slow_download)
    monkeypatch.setattr(docking, "prepare_receptor_pdbqt", lambda p: "fake.pdbqt")

    job = client.post(f"/api/docking/{sid}/receptor/prepare",
                      json={"gene": "TP53", "pdb_id": "1XQ8"}).json()["job_id"]
    dock["gene"] = "EGFR"  # 提交后立刻换基因(模拟竞态)
    body = _wait_job(client, job)
    assert body["status"] == "done" and body["result"].get("stale") is True, body
    dock = docking_session.get(sid)
    assert "rec_pdbqt" not in dock and "pdb_id" not in dock, "stale 受体不得写入会话"


def test_docking_ligand_stale_discard(client, sid, monkeypatch):
    """回归测试(P0-B 后端兜底):配体获取写回前核对来源上下文,被切即丢弃。

    缺陷:配体端点直接写 dock["cid"]/["lig_pdbqt"],来源在获取期间被切换时,
    旧来源的配体落进新来源会话,canRun 满足后 vina 对错误配体运行并打 ✓。
    复现要点:manual 会话 → smiles 获取(桩 2s)→ 1 秒时后台线程把
    source_mode 翻为 ppi → 响应 ok=False/stale=True 且会话无配体写入。
    """
    import threading

    import WY.core.docking as docking
    from WY.core import docking_session

    client.post(f"/api/docking/{sid}/source",
                json={"mode": "manual", "manual_genes": "TP53"})
    assert docking_session.get(sid).get("source_mode") == "manual"

    def slow_smiles(smiles, path):
        time.sleep(2)
        return path

    monkeypatch.setattr(docking, "smiles_to_sdf_3d", slow_smiles)
    monkeypatch.setattr(docking, "prepare_ligand_pdbqt", lambda p: "fake.pdbqt")

    def flip():
        d = docking_session.get(sid)
        if d is not None:
            d["source_mode"] = "ppi"

    threading.Timer(1.0, flip).start()
    body = client.post(f"/api/docking/{sid}/ligand",
                       json={"source": "smiles", "smiles": "CCO"}).json()
    assert body.get("ok") is False and body.get("stale") is True, body
    dock = docking_session.get(sid)
    assert "lig_pdbqt" not in dock and "cid" not in dock, "stale 配体不得写入会话"


def test_docking_receptor_candidates_guard_and_flow(client, sid, monkeypatch):
    assert client.post(f"/api/docking/{sid}/receptor/candidates",
                       json={"gene": "EGFR"}).status_code == 400
    client.post(f"/api/docking/{sid}/source",
                json={"mode": "manual", "manual_genes": "TP53"})
    import WY.core.docking as docking

    monkeypatch.setattr(docking, "fetch_receptor_pdb_candidates",
                        lambda *a, **k: pd.DataFrame())
    _pdb_cached = __import__("WY.api.docking", fromlist=["_pdb_candidates_cached"])
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


def test_ml_probe_batch(client, monkeypatch):
    """批量体检:单检同款判定归纳为四态徽章;空 rows/超限 → 400;
    series matrix 404 = 确定性事实(无 counts → 红),不留灰态。"""
    def _fake_probe(gse):
        if gse == "GSE404":
            raise RuntimeError("404 Client Error: Not Found for url: .../GSE404_series_matrix.txt.gz")
        return {"size_bytes": 10}           # 其余走 <20KB 空表分支
    monkeypatch.setattr("WY.core.geo.probe_series_matrix", _fake_probe)
    monkeypatch.setattr("WY.core.geo.has_counts_supplement",
                        lambda gse: gse == "GSEOK")
    r = client.post("/api/ml/probe-batch", json={"rows": []})
    assert r.status_code == 400
    r = client.post("/api/ml/probe-batch",
                    json={"rows": [{"gse": f"G{i}", "platform": "", "type": ""}
                                   for i in range(51)]})
    assert r.status_code == 400   # 超检索上限(50)才拒;50 行合法
    body = client.post("/api/ml/probe-batch", json={"rows": [
        {"gse": "gseok", "platform": "", "type": ""},
        {"gse": "gseno", "platform": "", "type": ""},
        {"gse": "gse404", "platform": "", "type": ""},
    ]}).json()
    by = {x["gse"]: x for x in body["results"]}
    assert by["GSEOK"]["status"] == "ok" and by["GSEOK"]["assess"]["counts_available"] is True
    assert by["GSENO"]["status"] == "unusable"
    # 404:assess 走空矩阵分支(counts 无 → 红),err 保留供详情展示
    assert by["GSE404"]["status"] == "unusable"
    assert by["GSE404"]["assess"]["matrix_empty"] is True
    assert by["GSE404"]["assess"]["counts_available"] is False
    assert "404" in by["GSE404"]["err"]
    assert by["GSEOK"]["counts_confirmed"] is None   # 未落盘,文件名级判定不实锤


def test_ml_overview_job_flow(client, monkeypatch):
    samples = pd.DataFrame({
        "title": [f"S{i}" for i in range(9)],
        "characteristics_ch1_1": (["disease state: obese T2DM and dyslipidemia"] * 6
                                  + ["disease state: control"] * 3),
        "characteristics_ch1_2": ["control"] * 5 + ["control group"] * 4,   # 子串重叠对
        "gender": ["male"] * 9,           # deny 名单 → 不进选择器
    })
    fake_ds = SimpleNamespace(issues=[], table=pd.DataFrame(), samples=samples)
    monkeypatch.setattr("WY.core.geo.fetch_geo_dataset", lambda gse, **k: fake_ds)
    monkeypatch.setattr("WY.core.geo.group_fields_overview",
                        lambda ds: pd.DataFrame())
    monkeypatch.setattr("WY.core.expression.assess_expression_dataset",
                    lambda ds: {"usable": True})
    r = client.post("/api/ml/overview", json={"gse": "gsex"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done"
    assert job["result"]["overview"] == [] and job["result"]["issues"] == []
    # 分组选择器 payload:字段/取值/整行命中数/重叠对/方向建议
    gp = job["result"]["group_picker"]
    by_name = {f["name"]: f for f in gp["fields"]}
    assert set(by_name) == {"characteristics_ch1_1", "characteristics_ch1_2"}  # gender 被 deny
    f1 = by_name["characteristics_ch1_1"]
    assert f1["label"] == "disease state"
    assert f1["suggest_control"] == "disease state: control"
    assert f1["suggest_case"] == "disease state: obese T2DM and dyslipidemia"
    n_map = {v["value"]: v["n"] for v in f1["values"]}
    assert n_map["disease state: control"] == 3
    assert n_map["disease state: obese T2DM and dyslipidemia"] == 6
    # 子串重叠:整行包含语义下 "control" 也命中 "control group" 行(9+4=9?——5+4=9)
    f2 = by_name["characteristics_ch1_2"]
    assert {v["value"]: v["n"] for v in f2["values"]} == {"control": 9, "control group": 4}
    assert {"a": "control", "b": "control group", "n": 4} in f2["overlaps"]


def test_ml_overview_counts_picker(client, monkeypatch):
    """counts 模式选择器:按字段取值分组 counts 列,建议=组内列名公共前缀;
    非样本列(gene_Dbxref)不入组;取值无对照词时不给方向建议。"""
    samples = pd.DataFrame({
        "title": [f"D5_20_{i}" for i in range(1, 4)] + [f"D5_CTL_{i}" for i in range(1, 4)],
        "characteristics_ch1_3": (["treatment: HGF increased"] * 3
                                  + ["treatment: medium only"] * 3),
        "gender": ["male"] * 6,
    })
    fake_ds = SimpleNamespace(issues=["series matrix empty"], table=pd.DataFrame(),
                              samples=samples)
    monkeypatch.setattr("WY.core.geo.fetch_geo_dataset", lambda gse, **k: fake_ds)
    monkeypatch.setattr("WY.core.geo.group_fields_overview", lambda ds: pd.DataFrame())
    monkeypatch.setattr("WY.core.expression.assess_expression_dataset",
                        lambda ds: {"usable": True, "matrix_empty": True})
    cnt = pd.DataFrame({c: [1] for c in
                        ["D5_20_1", "D5_20_2", "D5_20_3",
                         "D5_CTL_1", "D5_CTL_2", "D5_CTL_3", "gene_Dbxref"]})
    monkeypatch.setattr("WY.core.geo.fetch_gse_count_matrix",
                        lambda gse, prefer=None: {"counts": cnt})
    r = client.post("/api/ml/overview", json={"gse": "gsex"})
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done"
    cp = job["result"]["counts_picker"]
    f = cp["fields"][0]
    assert f["label"] == "treatment"
    assert f["suggest_case"] is None and f["suggest_control"] is None
    groups = {g["value"]: g for g in f["groups"]}
    assert groups["treatment: HGF increased"]["suggest"] == "D5_20_"
    assert groups["treatment: HGF increased"]["columns"] == ["D5_20_1", "D5_20_2", "D5_20_3"]
    assert groups["treatment: medium only"]["columns"] == ["D5_CTL_1", "D5_CTL_2", "D5_CTL_3"]
    assert all("gene_Dbxref" not in g["columns"] for g in f["groups"])


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

    monkeypatch.setattr("WY.api.ml.run_screen_chain", fake_chain)
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


def test_ml_runs_zip_rejects_outside_root(client, tmp_path):
    """回归测试(P1-B):zip 端点不得打包登记结果根之外的目录。

    缺陷:zip_ml_runs 只查 os.path.isdir,POST /api/ml/runs/zip 传任意存在的
    目录即整树打包返回(任意读)。对照 delete_ml_runs 早有「登记根包含 +
    summary.json 存在」双重校验,zip 漏配。
    复现要点:在登记结果根之外造一个含已知内容的目录,POST 其绝对路径;
    修复前 200 且 zip 内含该文件(任意读),修复后 400 且 detail 带被拒路径。
    顺带保绿:登记根下合法 run 仍 200,空列表仍 200(空目录包)。
    """
    import io as _io
    import zipfile as _zipfile

    outside = tmp_path / "outside_secret"
    outside.mkdir()
    (outside / "secret.txt").write_text("TOP-SECRET", encoding="utf-8")

    r = client.post("/api/ml/runs/zip", json={"paths": [str(outside)]})
    assert r.status_code == 400, f"登记根外目录必须 400,实际 {r.status_code}"
    assert "outside_secret" in (r.json().get("detail") or "")

    # 保绿:登记根下的合法 run 仍 200;根外目录不会出现在任何包里
    run_dir = os.path.join(config.GEO_RESULT_DIR, "GSEX_20260101_000000")
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, "summary.json"), "w", encoding="utf-8") as f:
        f.write("{}")
    r = client.post("/api/ml/runs/zip", json={"paths": [run_dir]})
    assert r.status_code == 200
    with _zipfile.ZipFile(_io.BytesIO(r.content)) as z:
        assert "secret.txt" not in "".join(z.namelist())

    # 保绿:空列表仍 200
    assert client.post("/api/ml/runs/zip", json={"paths": []}).status_code == 200


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


def test_spa_wildcard_blocks_drive_traversal(client, monkeypatch, tmp_path):
    """回归测试(P1-A):SPA 通配路由必须挡住 Windows 盘符绝对路径穿越。

    缺陷:旧实现仅以 `if ".." in full_path` 防穿越;Windows 下
    os.path.join(_DIST, "C:/...") 会整段丢弃 dist 前缀,
    GET /C:/... 可读任意盘任意文件(2026-11-08 实测确认)。
    复现要点:沙箱 dist 外放一个已知内容的文件,请求其盘符绝对路径;
    修复前返回 200 且内容命中(任意读),修复后 404。
    修复:与 /api/artifacts 同款 realpath 包含校验;正常资产与 SPA 回退
    (index.html + no-cache)响应形态不变。
    CI 教训:盘符断言仅 Windows 可执行——POSIX 下 as_posix() 以 / 开头,
    f"/{path}" 拼出的双斜杠 URL 会被客户端规范化剥成相对名(2026-11 CI
    ubuntu runner 实测:请求落回 SPA 回退 200),故按 os.name 守卫。
    """
    import WY.api as api_mod

    dist = tmp_path / "dist"
    (dist / "assets").mkdir(parents=True)
    (dist / "index.html").write_text("<html>spa-ok</html>", encoding="utf-8")
    (dist / "app.js").write_text("console.log(1)", encoding="utf-8")
    monkeypatch.setattr(api_mod, "_DIST", str(dist))

    # 盘符绝对路径(Windows 特有:join 丢前缀)→ 修复前 200 且命中内容,修复后 404
    if os.name == "nt":
        secret = tmp_path / "secret.txt"
        secret.write_text("TOP-SECRET-OUTSIDE-DIST", encoding="utf-8")
        r = client.get(f"/{secret.as_posix()}")
        assert r.status_code == 404, f"盘符绝对路径必须 404,实际 {r.status_code}"
        assert "TOP-SECRET" not in r.text
    # URL 编码穿越 → 404(编码段不被客户端规范化,可到达服务端;裸 ../ 会被
    # httpx 规范化,TestClient 层测不到,不设断言)
    assert client.get("/..%2f..%2f/etc/passwd").status_code == 404
    # dist 内正常文件仍 200
    r = client.get("/app.js")
    assert r.status_code == 200 and "console.log" in r.text
    # SPA 回退仍 200 + index.html + no-cache(R3:响应形态零变更)
    r = client.get("/some/deep/route")
    assert r.status_code == 200 and "spa-ok" in r.text
    assert r.headers.get("cache-control") == "no-cache"


# ======== 数据迁移:导出 job / 下载端点 / 两段式导入 ========

def _mk_analysis_dir(name: str = "aspirin_cancer_data") -> None:
    d = os.path.join(config.DOWNLOADS_DIR, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "nodes.csv"), "w", encoding="utf-8") as f:
        f.write("gene\nTP53\n")


def test_transfer_export_job_and_download(client):
    """导出:提交 job → 打包 done → 下载端点取回 zip;空范围/未知范围 400。"""
    assert client.post("/api/transfer/export",
                       json={"scopes": ["analysis"]}).status_code == 400  # 空库
    assert client.post("/api/transfer/export",
                       json={"scopes": ["md"]}).status_code == 400  # 未知范围
    assert client.post("/api/transfer/export",
                       json={"scopes": []}).status_code == 400

    _mk_analysis_dir()
    r = client.post("/api/transfer/export", json={"scopes": ["analysis", "ui_settings"]})
    assert r.status_code == 200
    job = _wait_job(client, r.json()["job_id"])
    assert job["status"] == "done"
    result = job["result"]
    assert result["counts"]["analysis"] == 1
    filename = result["filename"]

    z = client.get(f"/api/transfer/export/{filename}")
    assert z.status_code == 200 and z.content[:2] == b"PK"
    # 下载端点防线:固定命名形态之外一律 400,未知文件 404
    # (含 ../ 的多段路径在 httpx 侧已被归一,到不了本端点,测的是单段坏名)
    assert client.get("/api/transfer/export/evil.zip").status_code == 400
    assert client.get("/api/transfer/export/npap_export_20260101_000000.zip").status_code == 404
    assert client.get("/api/transfer/export/npap_export_20260101_0000_extra.zip").status_code == 400


def test_transfer_import_two_stage_flow(client):
    """导入两段式:上传有效包 → 预览报告 → commit 全 skip → 会话自动销毁。"""
    from WY.core import transfer as tr

    _mk_analysis_dir()
    export = client.post("/api/transfer/export", json={"scopes": ["analysis"]})
    job = _wait_job(client, export.json()["job_id"])
    zip_path = os.path.join(config.TRANSFER_DIR, job["result"]["filename"])
    with open(zip_path, "rb") as f:
        zip_bytes = f.read()

    up = client.post("/api/transfer/import", content=zip_bytes,
                     headers={"Content-Type": "application/zip"})
    assert up.status_code == 200
    body = up.json()
    sid, report = body["sid"], body["report"]
    assert report["items"][0]["name"] == "aspirin_cancer_data"
    assert report["items"][0]["conflict"] is True  # 本机原件还在

    commit = client.post(f"/api/transfer/import/{sid}/commit", json={
        "decisions": [{"kind": "analysis", "name": "aspirin_cancer_data",
                       "action": "skip"}]})
    assert commit.status_code == 200
    assert commit.json()["report"]["skipped"] == [
        {"kind": "analysis", "name": "aspirin_cancer_data"}]
    # 无失败项 → 会话已销毁;再 commit 404
    assert client.post(f"/api/transfer/import/{sid}/commit", json={
        "decisions": []}).status_code == 404
    assert client.delete(f"/api/transfer/import/{sid}").status_code == 404


def test_transfer_import_rejects_bad_bundles(client, tmp_path):
    """导入防线:非 zip / 缺 manifest 的包 400;空上传 400;未知会话 404。"""
    assert client.post("/api/transfer/import", content=b"",
                       headers={"Content-Type": "application/zip"}).status_code == 400
    assert client.post("/api/transfer/import", content=b"not a zip",
                       headers={"Content-Type": "application/zip"}).status_code == 400

    import zipfile

    nomanifest = tmp_path / "nomanifest.zip"
    with zipfile.ZipFile(nomanifest, "w") as zf:
        zf.writestr("evil.txt", "x")
    r = client.post("/api/transfer/import", content=nomanifest.read_bytes(),
                    headers={"Content-Type": "application/zip"})
    assert r.status_code == 400 and "manifest" in r.json()["detail"]

    assert client.post("/api/transfer/import/nosid/commit",
                       json={"decisions": []}).status_code == 404
    assert client.delete("/api/transfer/import/nosid").status_code == 404


def test_transfer_import_commit_failure_keeps_session(client):
    """commit 有失败项时保留会话,用户可调整决策重试。"""
    from WY.core import transfer as tr

    _mk_analysis_dir()
    export = client.post("/api/transfer/export", json={"scopes": ["analysis"]})
    job = _wait_job(client, export.json()["job_id"])
    zip_path = os.path.join(config.TRANSFER_DIR, job["result"]["filename"])
    with open(zip_path, "rb") as f:
        zip_bytes = f.read()
    sid = client.post("/api/transfer/import", content=zip_bytes,
                      headers={"Content-Type": "application/zip"}).json()["sid"]

    # rename 撞现有名 → failed,会话保留
    r = client.post(f"/api/transfer/import/{sid}/commit", json={
        "decisions": [{"kind": "analysis", "name": "aspirin_cancer_data",
                       "action": "rename", "new_name": "aspirin_cancer_data"}]})
    assert r.status_code == 200 and r.json()["report"]["failed"]
    assert tr.get_import_session(sid) is not None

    # 改用 overwrite → 成功,会话销毁
    r2 = client.post(f"/api/transfer/import/{sid}/commit", json={
        "decisions": [{"kind": "analysis", "name": "aspirin_cancer_data",
                       "action": "overwrite"}]})
    assert r2.json()["report"]["applied"] == [
        {"kind": "analysis", "name": "aspirin_cancer_data"}]
    assert tr.get_import_session(sid) is None


def test_settings_ui_stats_sources_and_gsea_pair(client):
    """新参数面板:快照含全部新字段;SOURCES 取值表校验(空集/非法 400);
    数值范围双保险;GSEA 规模窗成对校验(min<max,含单改一边与当前值比对)。"""
    body = client.get("/api/settings").json()
    for key in ("rwr_alpha", "gsea_permutation_num", "proximity_n_perm",
                "rwr_null_perms", "gsea_set_min_size", "gsea_set_max_size",
                "gsea_fdr_screen", "enrichment_sources", "exhaustiveness",
                "docking_cpu"):
        assert key in body, f"设置快照缺 {key}"

    originals = {k: getattr(config, k) for k in
                 ("RWR_ALPHA", "GSEA_PERMUTATION_NUM", "ENRICHMENT_SOURCES",
                  "GSEA_SET_MIN_SIZE", "GSEA_SET_MAX_SIZE", "EXHAUSTIVENESS")}
    try:
        r = client.post("/api/settings/ui", json={"changes": {
            "ENRICHMENT_SOURCES": ["kegg", "GO:CC"], "RWR_ALPHA": 0.85,
            "EXHAUSTIVENESS": 16}})
        assert r.status_code == 200
        snap = r.json()["snapshot"]
        assert snap["enrichment_sources"] == ["KEGG", "GO:CC"]  # 归一:大写、去重、保序
        assert snap["rwr_alpha"] == 0.85 and snap["exhaustiveness"] == 16
        assert config.ENRICHMENT_SOURCES == ("KEGG", "GO:CC")  # 运行时为元组,下游迭代安全

        for bad in ({"ENRICHMENT_SOURCES": []},                    # 空集
                    {"ENRICHMENT_SOURCES": ["bogus"]},             # 全非法
                    {"ENRICHMENT_SOURCES": "KEGG"},                # 非列表
                    {"RWR_ALPHA": 2},                              # 范围外(α 退化区)
                    {"EXHAUSTIVENESS": 99},                        # 范围外
                    {"DOCKING_CPU": -1},                           # 范围外
                    {"GSEA_SET_MIN_SIZE": 1000}):                  # 单改下限撞当前上限
            assert client.post("/api/settings/ui",
                               json={"changes": bad}).status_code == 400, bad
        # 成对修改仍倒挂 → 400;合法配对 → 200
        assert client.post("/api/settings/ui", json={"changes": {
            "GSEA_SET_MIN_SIZE": 600, "GSEA_SET_MAX_SIZE": 100}}).status_code == 400
        assert client.post("/api/settings/ui", json={"changes": {
            "GSEA_SET_MIN_SIZE": 20, "GSEA_SET_MAX_SIZE": 600}}).status_code == 200
    finally:
        client.post("/api/settings/ui", json={"changes": {
            "RWR_ALPHA": originals["RWR_ALPHA"],
            "GSEA_PERMUTATION_NUM": originals["GSEA_PERMUTATION_NUM"],
            "ENRICHMENT_SOURCES": list(originals["ENRICHMENT_SOURCES"]),
            "GSEA_SET_MIN_SIZE": originals["GSEA_SET_MIN_SIZE"],
            "GSEA_SET_MAX_SIZE": originals["GSEA_SET_MAX_SIZE"],
            "EXHAUSTIVENESS": originals["EXHAUSTIVENESS"]}})


def test_settings_snapshot_docking_params(client):
    """快照携带对接参数行所需的实时值:num_modes 常量与随深度放宽的有效超时。"""
    body = client.get("/api/settings").json()
    assert body["num_modes"] == config.NUM_MODES
    expected = max(config.DOCKING_TIMEOUT,
                   config.DOCKING_TIMEOUT * config.EXHAUSTIVENESS // 8)
    assert body["docking_timeout"] == expected
