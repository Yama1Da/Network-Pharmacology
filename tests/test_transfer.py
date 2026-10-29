"""core.transfer 的测试:导出打包(范围/密钥排除)、两段式导入(校验/冲突/决策落盘)、
导入会话与暂存清扫。目录路径全部走 conftest 沙箱(TRANSFER_DIR 已入 _SANDBOX_ATTRS)。

恶意包(路径穿越/缺 manifest/炸弹)用 zipfile 手工构造;正常包优先走
build_export_bundle 产出,保证"自己导出的包自己能导回"这一核心承诺被测到。
"""
import json
import os
import time
import zipfile

import pytest

from WY import config
from WY.core import transfer
from WY.core.transfer import TransferError


def _mkfile(path: str, content: str = "x") -> str:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    return path


def _mk_analysis(name: str = "aspirin_cancer_data") -> str:
    return _mkfile(os.path.join(config.DOWNLOADS_DIR, name, "nodes.csv"), "a,b\n")


def _mk_star(name: str = "manual_EGFR_20261025_130808") -> str:
    return _mkfile(os.path.join(config.MD_STAR_DIR, name, "best.pdbqt"))


def _mk_geo(name: str = "GSE123_20261023_225313", summary: bool = True) -> str:
    if summary:
        _mkfile(os.path.join(config.GEO_RESULT_DIR, name, "summary.json"), "{}")
    return _mkfile(os.path.join(config.GEO_RESULT_DIR, name, "deg_results.csv"))


def _mk_ui_settings(raw: dict | None = None) -> dict:
    raw = raw if raw is not None else {"UI_THEME": "dark", "EVIL_KEY": "x"}
    import os as _os
    from WY.core.caching import save_json_file

    save_json_file(config.UI_SETTINGS_FILE, raw)
    return raw


def _valid_manifest(items: list[dict], ui: bool = False) -> dict:
    m = transfer._manifest()
    m["items"] = items
    m["ui_settings"] = ui
    return m


def _mkzip(path: str, entries: dict[str, bytes | str], manifest: dict) -> str:
    """手工构造导入包:manifest + 任意条目(测试恶意包用)。"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        for arc, content in entries.items():
            zf.writestr(arc, content)
        zf.writestr(transfer.MANIFEST_NAME,
                    json.dumps(manifest, ensure_ascii=False))
    return path


# ======== 导出 ========

def test_export_bundle_roundtrip_content():
    """导出:四类范围各造一份数据,包内条目/manifest/counts 全对;MDresult 不进包。"""
    _mkfile(os.path.join(config.MD_RESULT_DIR, "run_20260101_000000", "out.pdbqt"))
    _mk_analysis()
    _mk_star()
    _mk_geo()
    _mk_ui_settings({"UI_THEME": "dark"})

    out = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    result = transfer.build_export_bundle(transfer.ALL_SCOPES, out)

    zf = zipfile.ZipFile(out)
    names = set(zf.namelist())
    assert "items/analysis/aspirin_cancer_data/nodes.csv" in names
    assert "items/star/manual_EGFR_20261025_130808/best.pdbqt" in names
    assert "items/geo/GSE123_20261023_225313/summary.json" in names
    assert "items/geo/GSE123_20261023_225313/deg_results.csv" in names
    assert "ui_settings.json" in names
    # MDresult 对接成果目录不迁移;密钥文件绝不进包
    assert not any("MDresult" in n for n in names)
    assert not any(n.lower().endswith(".env") for n in names)
    manifest = json.loads(zf.read(transfer.MANIFEST_NAME))
    assert manifest["manifest_version"] == transfer.MANIFEST_VERSION
    assert manifest["ui_settings"] is True
    assert {(it["kind"], it["name"]) for it in manifest["items"]} == {
        ("analysis", "aspirin_cancer_data"),
        ("star", "manual_EGFR_20261025_130808"),
        ("geo", "GSE123_20261023_225313"),
    }
    assert result["counts"] == {"analysis": 1, "star": 1, "geo": 1}
    assert result["ui_settings"] is True


def test_export_analysis_excludes_product_dirs_and_env_files():
    """analysis 范围:排除 MDresult/star/GEOresult 目录、隐藏项与 API.env/.env 文件。"""
    _mkfile(os.path.join(config.DOWNLOADS_DIR, "MDresult", "r", "x.pdbqt"))
    _mkfile(os.path.join(config.DOWNLOADS_DIR, "star", "s", "x.pdbqt"))
    _mkfile(os.path.join(config.DOWNLOADS_DIR, "GEOresult", "g", "summary.json"))
    _mkfile(os.path.join(config.DOWNLOADS_DIR, ".hidden_data", "x.csv"))
    _mkfile(os.path.join(config.DOWNLOADS_DIR, "aspirin_cancer_data", "API.env"), "k=1")
    _mkfile(os.path.join(config.DOWNLOADS_DIR, "aspirin_cancer_data", "nodes.csv"))
    out = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    transfer.build_export_bundle(("analysis",), out)
    names = set(zipfile.ZipFile(out).namelist())
    assert names == {"items/analysis/aspirin_cancer_data/nodes.csv",
                     transfer.MANIFEST_NAME}


def test_export_scope_validation_and_empty():
    """导出:未知范围/空范围直接拒;范围内没有任何结果时抛错(api 层转 400)。"""
    with pytest.raises(TransferError):
        transfer.build_export_bundle(("md",), os.path.join(config.TRANSFER_DIR, "x.zip"))
    with pytest.raises(TransferError):
        transfer.build_export_bundle((), os.path.join(config.TRANSFER_DIR, "x.zip"))
    with pytest.raises(TransferError):
        transfer.build_export_bundle(("analysis",), os.path.join(config.TRANSFER_DIR, "x.zip"))


def test_export_without_ui_settings_file():
    """ui_settings 文件不存在时:包里没有该文件,manifest 如实标 false。"""
    _mk_analysis()
    out = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    result = transfer.build_export_bundle(transfer.ALL_SCOPES, out)
    names = set(zipfile.ZipFile(out).namelist())
    assert "ui_settings.json" not in names
    assert result["ui_settings"] is False


# ======== 导入·分析(正常包走"自己导出自己导回") ========

def test_import_roundtrip_report_and_commit(tmp_path):
    """导出→导入全流程:预览报告条目/冲突/设置项正确;skip 不落盘、rename 落新名。"""
    import shutil

    _mk_analysis()
    _mk_geo()
    _mk_ui_settings({"UI_THEME": "dark"})
    zip_path = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    transfer.build_export_bundle(transfer.ALL_SCOPES, zip_path)
    # 模拟新机器:删掉本机原件,导入包是唯一副本
    shutil.rmtree(os.path.join(config.DOWNLOADS_DIR, "aspirin_cancer_data"))
    shutil.rmtree(os.path.join(config.GEO_RESULT_DIR, "GSE123_20261023_225313"))
    os.remove(config.UI_SETTINGS_FILE)

    staging = os.path.join(config.TRANSFER_DIR, "import_test")
    report = transfer.analyze_bundle(zip_path, staging)
    assert [(it["kind"], it["name"], it["conflict"]) for it in report["items"]] == [
        ("analysis", "aspirin_cancer_data", False),
        ("geo", "GSE123_20261023_225313", False),
    ]
    assert report["ui_settings"]["present"] is True
    assert report["ui_settings"]["keys"] == 1  # EVIL_KEY 被白名单过滤
    assert os.path.isfile(os.path.join(staging, "items", "geo",
                                       "GSE123_20261023_225313", "summary.json"))

    # skip:目标不受影响;rename:落新名
    report1 = transfer.commit_import(staging, [
        {"kind": "analysis", "name": "aspirin_cancer_data", "action": "skip"},
        {"kind": "geo", "name": "GSE123_20261023_225313", "action": "rename",
         "new_name": "GSE123_20261029_101010"},
    ])
    assert report1["skipped"] == [{"kind": "analysis", "name": "aspirin_cancer_data"}]
    assert report1["failed"] == []
    assert os.path.isdir(os.path.join(config.GEO_RESULT_DIR,
                                      "GSE123_20261029_101010"))
    assert not os.path.isdir(os.path.join(config.GEO_RESULT_DIR,
                                          "GSE123_20261023_225313"))


def test_import_commit_conflict_overwrite_and_rename_collision():
    """冲突决策:overwrite 替换旧目录;rename 撞现有名进 failed,旧数据不受影响。"""
    _mk_analysis()
    zip_path = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    transfer.build_export_bundle(("analysis",), zip_path)
    staging = os.path.join(config.TRANSFER_DIR, "import_test")
    report = transfer.analyze_bundle(zip_path, staging)
    assert report["items"][0]["conflict"] is True

    # rename 撞现有名 → failed,磁盘不动
    bad = transfer.commit_import(staging, [
        {"kind": "analysis", "name": "aspirin_cancer_data", "action": "rename",
         "new_name": "aspirin_cancer_data"}])
    assert bad["failed"] and "exists" in bad["failed"][0]["reason"]

    ok = transfer.commit_import(staging, [
        {"kind": "analysis", "name": "aspirin_cancer_data", "action": "overwrite"}])
    assert ok["failed"] == []
    nodes = os.path.join(config.DOWNLOADS_DIR, "aspirin_cancer_data", "nodes.csv")
    assert os.path.isfile(nodes)


def test_import_commit_ui_settings_applies_whitelist_and_runtime(monkeypatch):
    """ui_settings 导入:白名单过滤后整体替换,运行时同步生效(EVIL_KEY 绝不落地)。"""
    _mk_ui_settings({"UI_THEME": "dark", "EVIL_KEY": "x"})
    zip_path = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    transfer.build_export_bundle(("ui_settings",), zip_path)
    staging = os.path.join(config.TRANSFER_DIR, "import_test")
    transfer.analyze_bundle(zip_path, staging)

    monkeypatch.setattr(config, "UI_THEME", "light")  # 导入后应变 dark,测后自动还原
    result = transfer.commit_import(staging, [
        {"kind": "ui_settings", "name": "ui_settings.json", "action": "overwrite"}])
    assert result["failed"] == []
    assert result["applied"][0]["keys"] == ["UI_THEME"]
    assert config.UI_THEME == "dark"
    with open(config.UI_SETTINGS_FILE, encoding="utf-8") as f:
        assert json.load(f) == {"UI_THEME": "dark"}


def test_import_commit_invalid_geo_rejected():
    """geo 条目缺 summary.json 或目录名无时间戳:预览标 invalid,commit 进 failed 不落盘。"""
    import shutil

    _mk_geo("GSE123_20261023_225313", summary=False)
    _mk_geo("no_timestamp_geo")
    zip_path = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    transfer.build_export_bundle(("geo",), zip_path)
    # 模拟新机器:清掉原件,只留导入包
    shutil.rmtree(config.GEO_RESULT_DIR)
    staging = os.path.join(config.TRANSFER_DIR, "import_test")
    report = transfer.analyze_bundle(zip_path, staging)
    assert all(it.get("invalid") for it in report["items"])

    result = transfer.commit_import(staging, [
        {"kind": "geo", "name": "GSE123_20261023_225313", "action": "overwrite"}])
    assert result["failed"] and "summary.json" in result["failed"][0]["reason"]
    assert not os.path.isdir(os.path.join(config.GEO_RESULT_DIR,
                                          "GSE123_20261023_225313"))


def test_import_schema_warnings():
    """manifest 的算法/数据版本与当前不同:返回 warnings 提示,不阻断导入。"""
    _mk_analysis()
    zip_path = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    transfer.build_export_bundle(("analysis",), zip_path)
    staging = os.path.join(config.TRANSFER_DIR, "import_test")
    report = transfer.analyze_bundle(zip_path, staging)
    assert report["warnings"] == []  # 本机导出的包 schema 一致

    tampered = os.path.join(config.TRANSFER_DIR, "tampered.zip")
    with zipfile.ZipFile(zip_path) as src, \
            zipfile.ZipFile(tampered, "w", zipfile.ZIP_DEFLATED) as dst:
        for info in src.infolist():
            data = src.read(info.filename)
            if info.filename == transfer.MANIFEST_NAME:
                m = json.loads(data)
                m["schema"] = {"target_filter_version": config.TARGET_FILTER_VERSION + 5,
                               "ml_screening_version": config.ML_SCREENING_VERSION + 1,
                               "string_version": "v11.0"}
                data = json.dumps(m).encode("utf-8")
            dst.writestr(info, data)
    report2 = transfer.analyze_bundle(tampered,
                                      os.path.join(config.TRANSFER_DIR, "import_t2"))
    assert len(report2["warnings"]) == 3


# ======== 导入·安全防线(手工构造恶意包) ========

def test_import_rejects_path_traversal_and_absolute_entries(tmp_path):
    """路径安检:../ 穿越(含深嵌套)、绝对路径、盘符条目一律整包拒收。

    反斜杠条目单列:Python zipfile 在读入时就把 \\ 归一成 /(Windows),归一结果
    里的 .. 仍会被组件检查拦下,Linux 上则由 _safe_target 的反斜杠检查直接拒收。"""
    manifest = _valid_manifest([{"kind": "analysis", "name": "evil_data"}])
    for bad_entry in ("../evil.txt", "items/analysis/evil_data/../../../evil.txt",
                      "/etc/passwd", "C:/evil.txt"):
        zip_path = _mkzip(str(tmp_path / "bad.zip"), {bad_entry: "x"}, manifest)
        with pytest.raises(TransferError):
            transfer.analyze_bundle(zip_path, str(tmp_path / "staging"))


def test_import_rejects_missing_or_bad_manifest(tmp_path):
    """缺 manifest、manifest 非 JSON、缺版本号、版本超前、未知 kind、伪造条目名:全拒。"""
    with zipfile.ZipFile(str(tmp_path / "nomanifest.zip"), "w") as zf:
        zf.writestr("items/analysis/foo/x.csv", "x")
    with pytest.raises(TransferError):
        transfer.analyze_bundle(str(tmp_path / "nomanifest.zip"), str(tmp_path / "s1"))

    _mkzip(str(tmp_path / "badjson.zip"), {}, {"manifest_version": 1, "items": []})
    os.replace(str(tmp_path / "badjson.zip"), str(tmp_path / "badjson.zip.bak"))
    with zipfile.ZipFile(str(tmp_path / "badjson2.zip"), "w") as zf:
        zf.writestr(transfer.MANIFEST_NAME, "{not json")
    with pytest.raises(TransferError):
        transfer.analyze_bundle(str(tmp_path / "badjson2.zip"), str(tmp_path / "s2"))

    cases = [
        {"items": []},  # 缺 manifest_version
        {"manifest_version": transfer.MANIFEST_VERSION + 1, "items": []},  # 版本超前
        {"manifest_version": 1, "items": [{"kind": "md", "name": "x"}]},  # 未知 kind
        {"manifest_version": 1, "items": [{"kind": "analysis", "name": "../x"}]},  # 伪造名
        {"manifest_version": 1, "items": [{"kind": "analysis", "name": "a"},
                                          {"kind": "analysis", "name": "a"}]},  # 重复
    ]
    for i, m in enumerate(cases):
        zip_path = _mkzip(str(tmp_path / f"case{i}.zip"), {}, m)
        with pytest.raises(TransferError):
            transfer.analyze_bundle(zip_path, str(tmp_path / f"s{i}"))


def test_import_rejects_unexpected_entries(tmp_path):
    """包内条目必须对得上 manifest:清单外的条目(夹带私货)整包拒收。"""
    manifest = _valid_manifest([{"kind": "analysis", "name": "known_data"}])
    zip_path = _mkzip(str(tmp_path / "smuggled.zip"),
                      {"items/analysis/known_data/a.csv": "x",
                       "items/analysis/ghost_data/a.csv": "x"}, manifest)
    with pytest.raises(TransferError):
        transfer.analyze_bundle(zip_path, str(tmp_path / "staging"))
    # manifest 未声明 ui_settings 却夹带该文件:拒收
    zip2 = _mkzip(str(tmp_path / "smuggled_ui.zip"),
                  {"ui_settings.json": "{}"}, _valid_manifest([]))
    with pytest.raises(TransferError):
        transfer.analyze_bundle(zip2, str(tmp_path / "staging2"))


def test_import_rejects_zip_bomb(tmp_path, monkeypatch):
    """zip 炸弹防护:声明未压缩总量超上限拒收;单条目谎报大小(实际解压更大)也拒收。"""
    monkeypatch.setattr(config, "TRANSFER_MAX_EXTRACT_BYTES", 4)
    zip_path = _mkzip(str(tmp_path / "big.zip"),
                      {"items/analysis/a_data/x.csv": "12345678"},
                      _valid_manifest([{"kind": "analysis", "name": "a_data"}]))
    with pytest.raises(TransferError):
        transfer.analyze_bundle(zip_path, str(tmp_path / "s1"))

    monkeypatch.setattr(config, "TRANSFER_MAX_EXTRACT_BYTES", 1 << 30)
    lie = os.path.join(str(tmp_path), "lie.zip")
    with zipfile.ZipFile(lie, "w", zipfile.ZIP_DEFLATED) as zf:
        info = zipfile.ZipInfo("items/analysis/a_data/x.csv")
        zf.writestr(info, "12345678", compress_type=zipfile.ZIP_DEFLATED)
        # 直接改写中央目录声明大小,模拟"小声明大内容"
        zf.infolist()[0].file_size = 1
    with pytest.raises(TransferError):
        transfer.analyze_bundle(lie, str(tmp_path / "s2"))


def test_import_not_a_zip(tmp_path):
    """非 zip 文件:BadZipFile 转 TransferError(api 层转 400)。"""
    fake = str(tmp_path / "fake.zip")
    _mkfile(fake, "definitely not a zip")
    with pytest.raises(TransferError):
        transfer.analyze_bundle(fake, str(tmp_path / "s1"))


def test_import_ignores_os_junk_entries(tmp_path):
    """打包工具/操作系统杂项(__MACOSX、.DS_Store、._*):静默跳过不拒收。"""
    manifest = _valid_manifest([{"kind": "analysis", "name": "a_data"}])
    zip_path = _mkzip(str(tmp_path / "junk.zip"),
                      {"items/analysis/a_data/x.csv": "x",
                       "__MACOSX/items/._x.csv": "junk",
                       ".DS_Store": "junk",
                       "items/analysis/a_data/._meta": "junk"}, manifest)
    report = transfer.analyze_bundle(zip_path, str(tmp_path / "staging"))
    assert report["items"][0]["files"] == 1


# ======== 会话与暂存清扫 ========

def test_import_session_lifecycle():
    """会话:创建返回报告、get 续命、drop 连暂存目录一起删。"""
    _mk_analysis()
    zip_path = os.path.join(config.TRANSFER_DIR, transfer.export_filename())
    transfer.build_export_bundle(("analysis",), zip_path)

    sid, report = transfer.create_import_session(zip_path)
    assert len(sid) == 12
    assert report["items"][0]["name"] == "aspirin_cancer_data"
    session = transfer.get_import_session(sid)
    assert session is not None
    staging = session["staging"]
    assert os.path.isdir(staging)
    assert not os.path.exists(zip_path)  # 上传的 zip 解压后即删

    transfer.drop_import_session(sid)
    assert transfer.get_import_session(sid) is None
    assert not os.path.isdir(staging)


def test_cleanup_transfer_staging_sweeps_orphans_and_exports():
    """清扫:过期会话目录/半截上传/过期导出包删除;在册会话的暂存不动。"""
    os.makedirs(config.TRANSFER_DIR, exist_ok=True)
    old = time.time() - config.TRANSFER_SESSION_TTL - 100

    orphan_dir = os.path.join(config.TRANSFER_DIR, "import_orphan")
    os.makedirs(orphan_dir)
    stale_upload = _mkfile(os.path.join(config.TRANSFER_DIR, "upload_dead.zip"))
    stale_export = _mkfile(os.path.join(config.TRANSFER_DIR,
                                        "npap_export_20260101_000000.zip"))
    keep_export = _mkfile(os.path.join(config.TRANSFER_DIR,
                                       "npap_export_29991231_235959.zip"))
    for p in (orphan_dir, stale_upload, stale_export):
        os.utime(p, (old, old))

    # 一个在册会话:目录 mtime 很老,但会话刚摸过 → 必须保留
    sid, _ = transfer.create_import_session(
        _mkzip(os.path.join(config.TRANSFER_DIR, "upload_live.zip"),
               {}, _valid_manifest([])))
    live_dir = transfer.get_import_session(sid)["staging"]
    os.utime(live_dir, (old, old))

    transfer.cleanup_transfer_staging()
    assert not os.path.isdir(orphan_dir)
    assert not os.path.exists(stale_upload)
    assert not os.path.exists(stale_export)
    assert os.path.exists(keep_export)  # 未过期导出包保留
    assert os.path.isdir(live_dir)      # 在册会话暂存保留
