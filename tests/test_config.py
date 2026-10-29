"""WY.config 的测试:STRING 版本号校验归一、UI 设置覆盖加载(白名单与坏值防线)。

阶段三第 3 批新写。UI_SETTINGS_FILE 由 conftest 沙箱接管;覆盖函数会原地改写
config 模块常量,测试自行快照并恢复。
"""
import json
import os

from WY import config


def test_normalized_string_version():
    """版本号校验归一:漏 v 补 v、统一小写、允许字母后缀;非法格式返回 None。"""
    assert config.normalized_string_version("v12.0") == "v12.0"
    assert config.normalized_string_version("12.0") == "v12.0"      # 漏 v 自动补
    assert config.normalized_string_version("V10.0A") == "v10.0a"   # 大写归小写
    assert config.normalized_string_version("v10") == "v10"
    assert config.normalized_string_version(12.0) == "v12.0"        # 数字输入
    assert config.normalized_string_version("v12.0.1") is None      # 多段非法
    assert config.normalized_string_version("abc") is None
    assert config.normalized_string_version("") is None
    assert config.normalized_string_version("v") is None


def test_apply_ui_settings_overrides_whitelist_and_guards():
    """设置覆盖加载:白名单数值/合法版本/合法日志级生效;白名单外键与非法版本拒收;
    黑名单列表按字符串列表套用。函数原地改写 config 常量,测试快照恢复。"""
    saved = {k: getattr(config, k) for k in
             ("OT_TOP_N", "STRING_VERSION", "LOG_LEVEL", "CHEMBL_TARGET_BLACKLIST")}
    try:
        with open(config.UI_SETTINGS_FILE, "w", encoding="utf-8") as f:
            json.dump({
                "OT_TOP_N": 555,
                "STRING_VERSION": "11.5",            # 经 normalized_string_version 归一
                "LOG_LEVEL": "debug",                # 字符串枚举统一大写
                "CHEMBL_TARGET_BLACKLIST": ["alb", "spike ", ""],
                "BOGUS_KEY": 42,                     # 白名单外 → 忽略
            }, f)
        config.apply_ui_settings_overrides()
        assert config.OT_TOP_N == 555
        assert config.STRING_VERSION == "v11.5"
        assert config.LOG_LEVEL == "DEBUG"
        # 加载侧原样套用(列表即真值);归一化在保存侧 save_ui_settings 完成
        assert config.CHEMBL_TARGET_BLACKLIST == ("alb", "spike ", "")
    finally:
        for k, v in saved.items():
            setattr(config, k, v)


def test_apply_ui_settings_overrides_rejects_bad_values():
    """坏值防线:非法版本号与超界日志级一律不套用(按内置默认运行)。"""
    saved = {k: getattr(config, k) for k in ("STRING_VERSION", "LOG_LEVEL")}
    try:
        with open(config.UI_SETTINGS_FILE, "w", encoding="utf-8") as f:
            json.dump({"STRING_VERSION": "12.0.1", "LOG_LEVEL": "LOUD"}, f)
        config.apply_ui_settings_overrides()
        assert config.STRING_VERSION == saved["STRING_VERSION"]
        assert config.LOG_LEVEL == saved["LOG_LEVEL"]
    finally:
        for k, v in saved.items():
            setattr(config, k, v)


def test_apply_ui_settings_overrides_missing_file_is_noop(tmp_path, monkeypatch):
    """设置文件不存在:函数直接返回,不改任何常量。"""
    monkeypatch.setattr(config, "UI_SETTINGS_FILE", str(tmp_path / "nope.json"))
    before = config.OT_TOP_N
    config.apply_ui_settings_overrides()
    assert config.OT_TOP_N == before


def test_merged_prefixes_are_versioned():
    """merged 缓存前缀带过滤规则版本号(降噪规则变更后旧缓存自动失效的关键)。"""
    assert config.MERGED_DRUG_PREFIX.endswith("_") and "drug" in config.MERGED_DRUG_PREFIX
    assert config.MERGED_DISEASE_PREFIX.endswith("_") and "disease" in config.MERGED_DISEASE_PREFIX
    assert config.MERGED_DRUG_PREFIX != config.MERGED_DISEASE_PREFIX
