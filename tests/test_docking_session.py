"""core.docking_session 的测试:会话创建/读取/删除、下游失效、TTL 清扫。

阶段三第 3 批新写(0% 覆盖小模块)。进程内会话表经 monkeypatch 隔离。
"""
import time

from WY.core import docking_session as ds


def test_create_get_drop(monkeypatch):
    """会话生命周期:创建返回可取到的会话(仅 _touched 簿记),删除后取不到。"""
    monkeypatch.setattr(ds, "_SESSIONS", {})
    sid = ds.create()
    dock = ds.get(sid)
    assert dock is not None and set(dock) == {"_touched"}
    ds.drop(sid)
    assert ds.get(sid) is None
    assert ds.get("never") is None  # 未知 sid 安全返回 None


def test_second_create_does_not_reap_previous_session(monkeypatch):
    """回归(2026-10 API 冒烟发现):create 从未写 _touched,sweep 把所有既有
    会话判成过期——表现为「新建一个会话,之前所有会话立即 404」。"""
    monkeypatch.setattr(ds, "_SESSIONS", {})
    s1 = ds.create()
    s2 = ds.create()
    assert ds.get(s1) is not None, "新建第二个会话不得清掉第一个"
    assert ds.get(s2) is not None


def test_get_touches_session_extending_ttl(monkeypatch):
    """get 摸一次续命:TTL 按最后访问计,刚访问过的会话不被清扫。"""
    monkeypatch.setattr(ds, "_SESSIONS", {})
    sid = ds.create()
    ds._SESSIONS[sid]["_touched"] = time.time() - 3 * 3600   # 3 小时前创建
    ds.get(sid)                                              # 访问 → 续命
    ds.create()                                              # 触发清扫
    assert sid in ds._SESSIONS


def test_reset_downstream_clears_only_downstream_keys(monkeypatch):
    """下游失效:配体/受体/盒子/结果等旧步骤数据清空,会话本身与其他键保留。"""
    monkeypatch.setattr(ds, "_SESSIONS", {})
    dock = {"lig_pdbqt": "old", "gene": "TP53", "pdb_id": "2WXW", "result": "old",
            "drug_name": "metformin", "keep_me": "still here"}
    ds.reset_downstream(dock)
    assert "lig_pdbqt" not in dock and "gene" not in dock
    assert "pdb_id" not in dock and "result" not in dock and "drug_name" not in dock
    assert dock["keep_me"] == "still here"


def test_sweep_expires_stale_sessions(monkeypatch):
    """TTL 清扫:过期会话在新建会话时被清扫,未过期的保留。"""
    fresh: dict = {}
    monkeypatch.setattr(ds, "_SESSIONS", fresh)
    now = time.time()
    fresh["stale"] = {"_touched": now - 5 * 3600}   # 超过 4 小时 TTL
    fresh["alive"] = {"_touched": now - 60}
    ds.create()                                      # 新建触发清扫
    assert "stale" not in fresh
    assert "alive" in fresh
    assert len(fresh) == 2                           # alive + 刚创建的新会话
