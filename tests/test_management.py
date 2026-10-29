"""core.management 的测试:UI 设置落盘白名单、API.env 设置链路;env_discovery 工具链探测。

迁移自 selfcheck_docking.py(18 项中的 3 项)。ui_settings.json 路径由 conftest
沙箱接管;os.environ 操纵改用 monkeypatch.setenv/delenv(测试结束自动还原)。
"""
import json
import os

from WY import config
from WY.core.management import save_ui_settings
from WY.core.text_utils import safe_name


def test_ui_settings_roundtrip():
    """口袋半径作为全局分析参数:在齿轮白名单内,经 save_ui_settings 落盘到
    ui_settings.json,且多次写入按键合并(对接页与齿轮面板互不覆盖)。"""
    assert "POCKET_CUTOFF" in config.UI_SETTINGS_KEYS, config.UI_SETTINGS_KEYS
    save_ui_settings({"POCKET_CUTOFF": 4.5})
    save_ui_settings({"DGIDB_MIN_SCORE": 0.1})
    with open(config.UI_SETTINGS_FILE, encoding="utf-8") as f:
        saved = json.load(f)
    assert saved == {"POCKET_CUTOFF": 4.5, "DGIDB_MIN_SCORE": 0.1}, saved
    # 白名单外的键(如黑名单以外的任意字符串值)不落盘
    save_ui_settings({"NOT_IN_WHITELIST": 1.0})
    with open(config.UI_SETTINGS_FILE, encoding="utf-8") as f:
        assert "NOT_IN_WHITELIST" not in json.load(f)
    # 黑名单(列表型参数):保存侧归一化——去空格、统一大写、去重、丢空项
    save_ui_settings({"CHEMBL_TARGET_BLACKLIST": ["alb", "spike ", "", "ALB"]})
    with open(config.UI_SETTINGS_FILE, encoding="utf-8") as f:
        assert json.load(f)["CHEMBL_TARGET_BLACKLIST"] == ["ALB", "SPIKE"]


def test_api_env_roundtrip(monkeypatch, tmp_path):
    """API.env 设置链路(save_api_env):白名单过滤、值消毒(换行/井号拒收)、set_key 落盘
    保留注释与其他行、运行时回填含派生值联动(NCBI key 生效 → E-UTILS 节流间隔切 0.10)。"""
    from WY.core.management import save_api_env

    base = str(tmp_path)
    env_path = os.path.join(base, "API.env")
    with open(env_path, "w", encoding="utf-8") as f:
        f.write("# 注释行必须保留\nCLI_RESULTS_DIR=D:/old\n")
    orig = {k: getattr(config, k) for k in config.API_ENV_KEYS}
    orig_interval = config.HOST_MIN_REQUEST_INTERVAL["NCBI_EUTILS"]
    try:
        monkeypatch.setattr(config, "BASE_DIR", base)
        applied = save_api_env({
            "DisGeNET_API_KEY": "  abc123  ",   # 消毒:去首尾空白后生效
            "NCBI_API_KEY": "bad\nvalue",       # 消毒:含换行拒收
            "VINA_DIR": "",                     # 空值跳过(留空=不修改)
            "NOT_IN_WHITELIST": "x",            # 白名单外丢弃
            "OBABEL_DIR": '"D:/obabel dir"',    # 消毒:剥掉成对首尾引号(复制路径常带引号)
            "CLI_RESULTS_DIR": "D:/new",        # 已有行被覆盖
        })
        assert applied == ["DisGeNET_API_KEY", "OBABEL_DIR", "CLI_RESULTS_DIR"], applied
        assert config.DisGeNET_API_KEY == "abc123" and config.CLI_RESULTS_DIR == "D:/new"
        assert config.OBABEL_DIR == "D:/obabel dir", config.OBABEL_DIR
        with open(env_path, encoding="utf-8") as f:
            text = f.read()
        assert "# 注释行必须保留" in text
        assert "abc123" in text and "D:/new" in text and "D:/old" not in text
        assert "NOT_IN_WHITELIST" not in text and "bad" not in text
        # 引号已剥离:文件里不应出现带引号的值
        assert '"D:/obabel dir"' not in text
        # 运行时回填联动派生值:NCBI key 生效 → E-UTILS 节流间隔切 0.10
        save_api_env({"NCBI_API_KEY": "ncbi_key_456"})
        assert config.NCBI_API_KEY == "ncbi_key_456"
        assert config.HOST_MIN_REQUEST_INTERVAL["NCBI_EUTILS"] == 0.10
    finally:
        # save_api_env 会原地修改运行时 config(API_ENV_KEYS 各值与节流间隔),
        # monkeypatch 管不到,按自检脚本同款手法显式恢复
        config.apply_api_env(orig)
        config.HOST_MIN_REQUEST_INTERVAL["NCBI_EUTILS"] = orig_interval


def test_env_discovery(monkeypatch, tmp_path):
    """工具链探测(env_discovery,自 config 拆出):find_exe 的前缀匹配/exclude/非文件与
    空目录安全空串,find_vina_dir 的 环境变量>PATH 回退与"可执行文件全路径"识别、
    都找不到时返回空串。打桩 PATH 只含测试目录,不受本机真实安装干扰。"""
    from WY import env_discovery

    base = str(tmp_path)
    vina_dir = os.path.join(base, "vina_bin")
    os.makedirs(vina_dir)
    for name in ("vina_1.2.7_win.exe", "vina_split_1.2.7_win.exe",
                 "obabel.exe", "readme.txt"):
        open(os.path.join(vina_dir, name), "wb").close()
    # find_exe:前缀匹配 + exclude 排除 vina_split + 多匹配取字典序最后一个
    assert env_discovery.find_exe(vina_dir, "vina", exclude="vina_split") == \
        os.path.join(vina_dir, "vina_1.2.7_win.exe")
    assert env_discovery.find_exe(vina_dir, "vina_split") == \
        os.path.join(vina_dir, "vina_split_1.2.7_win.exe")
    assert env_discovery.find_exe(vina_dir, "obabel") == \
        os.path.join(vina_dir, "obabel.exe")
    # 目录本身(其名以 vina 开头但是目录)、空参、不存在的目录:安全返回空串
    assert env_discovery.find_exe(base, "vina") == ""
    assert env_discovery.find_exe("", "vina") == ""
    assert env_discovery.find_exe(os.path.join(base, "nope"), "vina") == ""
    # find_vina_dir:环境变量优先;填可执行文件全路径时返回其所在目录
    monkeypatch.setenv("VINA_DIR", vina_dir)
    assert env_discovery.find_vina_dir() == vina_dir
    monkeypatch.setenv("VINA_DIR", os.path.join(vina_dir, "vina_1.2.7_win.exe"))
    assert env_discovery.find_vina_dir() == vina_dir
    # PATH 回退:打桩期间 PATH 只含该目录,防本机真实安装干扰
    monkeypatch.delenv("VINA_DIR", raising=False)
    monkeypatch.setenv("PATH", vina_dir)
    assert env_discovery.find_vina_dir() == vina_dir
    # 都找不到:返回空串(运行对接时报错提示设置 VINA_DIR)
    monkeypatch.setenv("PATH", base)
    assert env_discovery.find_vina_dir() == ""
    # find_obabel_dir 同规则:环境变量指向目录即命中
    monkeypatch.setenv("OBABEL_DIR", vina_dir)
    assert env_discovery.find_obabel_dir() == vina_dir


# ---------- 阶段三第 3 批扩充:管理编排与设置链路 ----------


def _fake_spec(src_dir, merged_dir, fetch_results=None, merged_df=None):
    """管理页共用的 spec 字典假件:单来源前缀,抓取/整合行为可注入。"""
    calls = {"invalidated": [], "merged": []}

    def fake_fetch(name, cache_hours=None):
        calls["last_fetch"] = name
        if fetch_results is not None:
            return fetch_results(name)
        return __import__("pandas").DataFrame({"approvedSymbol": ["TP53"]}), "ok"

    def fake_fetcher(name, cache_hours=None):
        df, state = fake_fetch(name, cache_hours)
        if not df.empty:  # 与真实抓取器一致:成功才落盘
            import pandas as pd
            df.to_csv(os.path.join(src_dir, f"ChEMBL_{safe_name(name)}.csv"),
                      index=False)
        return df, state

    def fake_merged(name):
        calls["merged"].append(name)
        import pandas as pd
        return merged_df if merged_df is not None else pd.DataFrame(
            {"approvedSymbol": ["TP53", "EGFR"]})

    def fake_invalidate(name):
        calls["invalidated"].append(name)

    spec = {
        "cache_dir": src_dir,
        "merged_dir": merged_dir,
        "merged_prefix": "merged_drug_v1_",
        "source_prefixes": ("ChEMBL_",),
        "fetch_source": {"ChEMBL_": fake_fetcher},
        "fetch_merged": fake_merged,
        "invalidate_merged": fake_invalidate,
    }
    return spec, calls


def test_purge_merged_cache(tmp_path, monkeypatch):
    """清空整合缓存:只删 *.csv,返回删除数;目录为空/不存在返回 0。"""
    from WY.core.management import purge_merged_cache
    monkeypatch.setattr(config, "MERGED_CACHE_DIR", str(tmp_path / "merged"))
    os.makedirs(config.MERGED_CACHE_DIR)
    for name in ("a.csv", "b.csv"):
        open(os.path.join(config.MERGED_CACHE_DIR, name), "w").close()
    open(os.path.join(config.MERGED_CACHE_DIR, "keep.txt"), "w").close()
    assert purge_merged_cache() == 2
    assert os.path.exists(os.path.join(config.MERGED_CACHE_DIR, "keep.txt"))
    assert purge_merged_cache() == 0
    monkeypatch.setattr(config, "MERGED_CACHE_DIR", str(tmp_path / "nope"))
    assert purge_merged_cache() == 0


def test_humanize_analysis_label():
    """分析标识人读化:首个下划线换箭头,其余换回空格。"""
    from WY.core.management import humanize_analysis_label
    assert humanize_analysis_label("aspirin_colorectal_cancer") == \
        "aspirin → colorectal cancer"
    assert humanize_analysis_label("metformin") == "metformin"
    assert humanize_analysis_label("") == ""


def test_enrichment_display_name_formats():
    """富集文件名展示:ORA 带哈希/无哈希、GSEA 带侧别、无法识别原样返回。"""
    from WY.core.management import enrichment_display_name
    out = enrichment_display_name("gProfiler_aspirin_cancer_12345678.csv")
    assert out == "aspirin → cancer · gProfiler"
    out2 = enrichment_display_name("gProfiler_aspirin_cancer_module3.csv")
    assert out2.startswith("aspirin → cancer · gProfiler · ") and "3" in out2
    out3 = enrichment_display_name("GSEA_drug_metformin_type_2_diabetes.csv")
    assert out3.startswith("metformin → type 2 diabetes · GSEA · ")
    assert enrichment_display_name("random_name.csv") == "random_name.csv"


def test_refresh_source_ok_replaces_backup_cleared(tmp_path):
    """来源刷新成功:新文件落盘、旧备份清除、整合缓存被失效、状态 OK。"""
    from WY.core.management import refresh_cache_file
    src, merged = str(tmp_path / "src"), str(tmp_path / "merged")
    os.makedirs(src)
    old = os.path.join(src, "ChEMBL_metformin.csv")
    with open(old, "w", encoding="utf-8") as f:
        f.write("old")
    spec, calls = _fake_spec(src, merged)
    result = refresh_cache_file("ChEMBL_metformin.csv", spec, clear_http=False)
    assert result["Status"] == "OK" and result["Message"] == "1 targets"
    assert not os.path.exists(old + ".bak"), "成功后旧备份应清除"
    assert calls["invalidated"] == ["metformin"]


def test_refresh_source_failure_restores_old(tmp_path):
    """来源刷新失败:旧文件从 .bak 恢复,数据不丢,状态 FAILED。"""
    from WY.core.management import refresh_cache_file
    src, merged = str(tmp_path / "src"), str(tmp_path / "merged")
    os.makedirs(src)
    old = os.path.join(src, "ChEMBL_metformin.csv")
    with open(old, "w", encoding="utf-8") as f:
        f.write("precious old data")
    spec, _ = _fake_spec(src, merged, fetch_results=lambda n: (_ for _ in ()).throw(
        ConnectionError("simulated outage")))
    result = refresh_cache_file("ChEMBL_metformin.csv", spec, clear_http=False)
    assert result["Status"] == "FAILED"
    with open(old, encoding="utf-8") as f:
        assert f.read() == "precious old data"
    assert not os.path.exists(old + ".bak")


def test_refresh_source_empty_restores_and_reports(tmp_path):
    """抓取成功但 0 靶点(不落盘):旧文件恢复占位,状态 EMPTY(不算失败)。"""
    import pandas as pd
    from WY.core.management import refresh_cache_file
    src, merged = str(tmp_path / "src"), str(tmp_path / "merged")
    os.makedirs(src)
    old = os.path.join(src, "ChEMBL_metformin.csv")
    with open(old, "w", encoding="utf-8") as f:
        f.write("old")
    spec, _ = _fake_spec(src, merged,
                         fetch_results=lambda n: (pd.DataFrame(), "empty"))
    result = refresh_cache_file("ChEMBL_metformin.csv", spec, clear_http=False)
    assert result["Status"] == "EMPTY"
    with open(old, encoding="utf-8") as f:
        assert f.read() == "old"


def test_refresh_merged_rebuilds_after_sources(tmp_path):
    """整合文件刷新:先更新两个来源再重建;来源失败则整体 FAILED 不重建。"""
    from WY.core.management import refresh_cache_file
    src, merged = str(tmp_path / "src"), str(tmp_path / "merged")
    os.makedirs(src)
    os.makedirs(merged)
    spec, calls = _fake_spec(src, merged)
    result = refresh_cache_file("merged_drug_v1_metformin.csv", spec, clear_http=False)
    assert result["Status"] == "OK"
    assert calls["merged"] == ["metformin"] and calls["invalidated"] == ["metformin"]

    # 来源失败 → 整合流程直接 FAILED,不调 fetch_merged
    spec2, calls2 = _fake_spec(src, merged, fetch_results=lambda n: (
        __import__("pandas").DataFrame(), "failed"))
    result2 = refresh_cache_file("merged_drug_v1_metformin.csv", spec2, clear_http=False)
    assert result2["Status"] == "FAILED"
    assert calls2["merged"] == []


def test_update_all_cache_files_serial(tmp_path):
    """批量更新:来源文件并发清点、逐个刷新、merged 串行重建,进度回调走完全程。"""
    from WY.core.management import update_all_cache_files
    src, merged = str(tmp_path / "src"), str(tmp_path / "merged")
    os.makedirs(src)
    for n in ("ChEMBL_metformin.csv", "ChEMBL_aspirin.csv", "junk.txt"):
        open(os.path.join(src, n), "w").close()
    spec, calls = _fake_spec(src, merged)
    ticks = []
    out = update_all_cache_files(spec, max_workers=1,
                                 on_progress=lambda f, t: ticks.append(f))
    assert list(out["File"]) == [
        "ChEMBL_aspirin.csv", "ChEMBL_metformin.csv", "merged_drug_v1_aspirin.csv",
        "merged_drug_v1_metformin.csv"]
    assert set(out["Status"]) == {"OK"}
    assert len(calls["merged"]) == 2
    assert ticks and abs(ticks[-1] - 1.0) < 1e-9  # 进度走到 100%


def test_activate_string_version_updates_derived(tmp_path, monkeypatch):
    """切换 STRING 版本:版本常量、版本化子域名即时联动,并持久化到设置文件。"""
    import WY.core.interactome as itx_mod
    from WY.core.management import activate_string_version
    cleared = []
    monkeypatch.setattr(itx_mod, "clear_runtime_caches", lambda: cleared.append(1))
    orig_version, orig_url = config.STRING_VERSION, config.STRING_BASE_URL
    try:
        out = activate_string_version("v11.5")
        assert out == "v11.5"
        assert config.STRING_VERSION == "v11.5"
        assert config.STRING_BASE_URL == "https://version-11-5.string-db.org/api/json"
        assert cleared == [1]
        with open(config.UI_SETTINGS_FILE, encoding="utf-8") as f:
            import json
            assert json.load(f)["STRING_VERSION"] == "v11.5"
    finally:
        config.STRING_VERSION = orig_version
        config.STRING_BASE_URL = orig_url
