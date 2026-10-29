"""core.caching 的测试:JSON/文本/CSV 原子写读回环、TTL 新鲜度、损坏容错、过期清理。

阶段三第 1 批新写(原自检脚本未覆盖此模块)。
"""
import os
import time

import pandas as pd

from WY import config
from WY.core.caching import (
    clean_expired_cache_files,
    is_cache_fresh,
    load_json_file,
    open_csv_cache,
    open_named_csv_cache,
    read_cache_key,
    save_csv_cache,
    save_json_file,
    save_text_atomic,
)


def test_json_roundtrip(tmp_path):
    """存进去什么,读出来就该是什么(含中文与嵌套结构)。"""
    path = str(tmp_path / "x.json")
    save_json_file(path, {"a": 1, "b": ["汉字", 2.5]})
    assert load_json_file(path) == {"a": 1, "b": ["汉字", 2.5]}


def test_load_missing_returns_default(tmp_path):
    """读不存在的文件返回默认值而不是崩溃——这是全项目容错的起点。"""
    assert load_json_file(str(tmp_path / "nope.json"), default={}) == {}


def test_load_corrupt_returns_default(tmp_path):
    """损坏的 JSON 按不存在处理:记日志、返回 default,绝不抛异常炸掉调用方。"""
    path = str(tmp_path / "broken.json")
    with open(path, "w", encoding="utf-8") as f:
        f.write("{not valid json")
    assert load_json_file(path, default={}) == {}


def test_save_creates_missing_dirs(tmp_path):
    """目录不存在时自动创建(多级嵌套)。"""
    path = str(tmp_path / "lvl1" / "lvl2" / "x.json")
    save_json_file(path, {"ok": True})
    assert load_json_file(path) == {"ok": True}


def test_atomic_write_cleans_tmp_on_failure(tmp_path):
    """写入中途失败(值不可序列化):.tmp 必须被清掉,不能留半截文件冒充有效缓存。"""
    path = str(tmp_path / "x.json")
    try:
        save_json_file(path, {"bad": {1, 2}})  # set 不可 JSON 序列化 → dump 抛 TypeError
        raise AssertionError("non-serializable value must raise")
    except TypeError:
        pass
    leftovers = [n for n in os.listdir(str(tmp_path)) if n.endswith(".tmp")]
    assert not leftovers, leftovers
    assert not os.path.exists(path)


def test_save_text_atomic_roundtrip(tmp_path):
    """文本原子写:内容一致,且不留 .tmp。"""
    path = str(tmp_path / "rec.pdb")
    save_text_atomic(path, "ATOM      1  N   MET\n")
    assert open(path, encoding="utf-8").read() == "ATOM      1  N   MET\n"
    assert not [n for n in os.listdir(str(tmp_path)) if n.endswith(".tmp")]


def test_cache_freshness(tmp_path):
    """刚写的文件是新鲜的;把修改时间拨回 2 小时前,超过 1 小时 TTL 就不新鲜;缺文件直接不新鲜。"""
    path = str(tmp_path / "c.csv")
    pd.DataFrame({"g": ["TP53"]}).to_csv(path, index=False)
    assert is_cache_fresh(path, ttl_seconds=3600)
    old = time.time() - 2 * 3600
    os.utime(path, (old, old))
    assert not is_cache_fresh(path, ttl_seconds=3600)
    assert not is_cache_fresh(str(tmp_path / "nope.csv"), ttl_seconds=3600)


def test_open_csv_cache_hit_vs_fresh(tmp_path):
    """fresh=False:文件存在即命中;fresh=True:超过全局 TTL 视为未命中(返回 None 由调用方重新请求)。"""
    path = str(tmp_path / "drug.csv")
    pd.DataFrame({"gene": ["TP53"]}).to_csv(path, index=False)
    assert open_csv_cache(path) is not None            # 存在即命中
    # 回拨到超过全局 TTL(config.CACHE_TTL_SECONDS)之前 → 过期
    old = time.time() - (config.CACHE_TTL_SECONDS + 3600)
    os.utime(path, (old, old))
    assert open_csv_cache(path, fresh=True, label="test") is None  # 过期 → 未命中
    assert open_csv_cache(path) is not None            # fresh=False 仍按存在命中
    assert open_csv_cache(str(tmp_path / "nope.csv")) is None      # 缺文件 → 未命中


def test_save_csv_cache_skips_empty(tmp_path):
    """默认空表不落盘(避免把失败/空结果当成有效缓存),返回 False;非空表正常写出。"""
    path = str(tmp_path / "e.csv")
    assert save_csv_cache(pd.DataFrame(), path) is False
    assert not os.path.exists(path)
    assert save_csv_cache(pd.DataFrame({"g": ["TP53"]}), path) is True
    assert os.path.exists(path)


def test_csv_cache_header_lines_roundtrip(tmp_path):
    """缓存键行(# 前缀)与数据行互不干扰:read_cache_key 读到键,pandas 跳过注释行读到数据。"""
    path = str(tmp_path / "k.csv")
    save_csv_cache(pd.DataFrame({"gene": ["TP53"]}), path,
                   header_lines=["drug=metformin", "disease=T2D"])
    assert read_cache_key(path) == "drug=metformin"
    back = pd.read_csv(path, encoding="utf-8-sig", comment="#")
    assert list(back["gene"]) == ["TP53"]
    # 旧格式(无键行)→ None;缺文件 → None
    plain = str(tmp_path / "plain.csv")
    pd.DataFrame({"g": [1]}).to_csv(plain, index=False)
    assert read_cache_key(plain) is None
    assert read_cache_key(str(tmp_path / "nope.csv")) is None


def test_read_cache_csv_damaged_returns_none(tmp_path):
    """内容残缺的缓存(非法字节流):删除后按未命中(None)返回,由调用方重新抓取覆盖。"""
    path = str(tmp_path / "bad.csv")
    with open(path, "wb") as f:
        f.write(b"\xff\xfe\xff,b\n\x00\x01,2\n")  # 非法 UTF-8 → UnicodeDecodeError
    assert open_csv_cache(path) is None
    assert not os.path.exists(path), "损坏缓存必须被删除"


def test_open_named_csv_cache(tmp_path, monkeypatch):
    """按规范名拼缓存路径样板:路径含 safe_name、remember 回调先于读缓存被调、未命中返回 None。"""
    from WY.core.text_utils import safe_name

    seen = []
    csv_path, df = open_named_csv_cache(
        str(tmp_path), "drug_", "Type 2 Diabetes",
        label="test", remember=lambda original: seen.append(original))
    assert os.path.basename(csv_path) == f"drug_{safe_name('Type 2 Diabetes')}.csv"
    assert seen == ["Type 2 Diabetes"]
    assert df is None  # 未命中


def test_clean_expired_cache_files(tmp_path):
    """过期清理:超龄缓存删除、新文件与其他扩展名保留、目录不存在返回 0。"""
    old_csv = tmp_path / "old.csv"
    new_csv = tmp_path / "new.csv"
    old_json = tmp_path / "old.json"
    keep_txt = tmp_path / "old.txt"  # 非缓存扩展名,不参与清理
    for p in (old_csv, new_csv, old_json, keep_txt):
        p.write_text("x", encoding="utf-8")
    old = time.time() - 10 * 86400
    for p in (old_csv, old_json, keep_txt):
        os.utime(p, (old, old))
    removed = clean_expired_cache_files(str(tmp_path), ttl_seconds=7 * 86400)
    assert removed == 2, removed  # old.csv + old.json
    assert not old_csv.exists() and not old_json.exists()
    assert new_csv.exists() and keep_txt.exists()
    assert clean_expired_cache_files(str(tmp_path / "nope")) == 0
    # 缺省 TTL 取全局配置(为正数)
    assert config.CACHE_TTL_SECONDS > 0
