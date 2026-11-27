# -*- coding: utf-8 -*-
"""reach_compare 单元测试(§3.2 招牌图,2026-11-27):门槛混合规则、两版排序、
demoted 名单、渲染落盘冒烟。纯函数 + Figure API,无网络/无缓存依赖。"""

import os

import pytest

from WY.core.reach_compare import (COMPARE_SCHEMA, build_reach_compare,
                                   render_reach_compare_png)


def _rows(n_unknown=160, judged=(("LMNA", "reachable"), ("ADRA1A", "reachable"),
                                ("DPP4", "partial"))):
    """metformin 形态:少量判档 + 大量 unknown。判档基因分数 0.09 = 真实案例
    形态——一期排在 top-10 之外(UNK009=0.1 压线),加权后爬进中游并把
    UNK008/UNK009 挤出前 10(demoted 名单非空的最低配置)。"""
    scores = {}
    tiers = {}
    for i in range(n_unknown):
        g = f"UNK{i:03d}"
        scores[g] = 1.0 / (i + 1)   # 等权降序
        tiers[g] = "unknown"
    for g, t in judged:
        scores[g] = 0.09
        tiers[g] = t
    return scores, tiers


def test_gate_hybrid_rule_blocks_metformin_case():
    """门槛 max(5, ceil(5%×N)):metformin 形态 3/163 < max(5, 9) → 拦(方案 §3.2 ②);
    sorafenib 形态 145/562 > 29 → 放行;全 unknown 0/8 → 拦(纯百分比会得 0,混合规则兜住)。"""
    scores, tiers = _rows()
    res = build_reach_compare(scores, tiers)
    assert res["schema"] == COMPARE_SCHEMA and res["ok"]
    assert res["gate"] is False and res["judged"] == 3

    # sorafenib 形态:562 靶点 145 判档
    big_scores = {f"G{i}": 1.0 / (i + 1) for i in range(562)}
    big_tiers = {f"G{i}": ("reachable" if i < 145 else "unknown") for i in range(562)}
    res2 = build_reach_compare(big_scores, big_tiers)
    assert res2["gate"] is True and res2["judged"] == 145 and res2["total"] == 562

    # 小样本:8 靶点 0 判档 → ceil(5%×8)=1 但绝对下限 5 → 0 < 5 拦
    small = {f"G{i}": 1.0 - i * 0.1 for i in range(8)}
    small_t = {g: "unknown" for g in small}
    assert build_reach_compare(small, small_t)["gate"] is False
    # 小样本全判档:8 ≥ max(5, 1) → 放行
    assert build_reach_compare(small, {g: "reachable" for g in small})["gate"] is True


def test_rank_math_and_demoted_list():
    """两版排序:unknown 全体同乘 0.5 保序;判档基因 ×1.0 越过 2 倍分数距内
    unknown;一期 Top-10 的 unknown 被挤出 → demoted 名单(metformin OCT1 案例形态)。"""
    scores, tiers = _rows()
    res = build_reach_compare(scores, tiers)
    rows = res["rows"]
    # 排名严格置换:两列 rank 各自 1..N 无重复
    assert sorted(rows["rank_before"]) == list(range(1, len(rows) + 1))
    assert sorted(rows["rank_after"]) == list(range(1, len(rows) + 1))
    # 判档可达基因(w=1.0):0.09×1.0 越过全部加权分 < 0.09 的 unknown(UNK005
    # 起),从一期 11/12 名爬到 6/7 名;DPP4(×0.6=0.144)爬 1 位
    lmna = rows[rows["gene"] == "LMNA"].iloc[0]
    assert lmna["rank_after"] < lmna["rank_before"] and lmna["rank_before"] >= 11
    assert lmna["w"] == 1.0
    dpp4 = rows[rows["gene"] == "DPP4"].iloc[0]
    assert dpp4["w"] == 0.6 and dpp4["rank_after"] < dpp4["rank_before"]
    # demoted:一期 Top-10 全 unknown(判档 0.09 < UNK009 0.1 压线在外),
    # 加权后 UNK008/UNK009 被挤出 → 名单非空且全 unknown、去向 >10(metformin 案例)
    assert res["demoted"], "expected UNK008/UNK009 demoted"
    assert all(r["tier"] == "unknown" for r in res["demoted"])
    assert all(r["rank_after"] > 10 for r in res["demoted"])
    assert {r["gene"] for r in res["demoted"]} == {"UNK008", "UNK009"}


def test_render_smoke(tmp_path):
    """渲染冒烟:过门槛数据 → PNG 落盘非空;图内英文无 CJK(venn 同款纪律)。"""
    big_scores = {f"G{i}": 1.0 / (i + 1) for i in range(120)}
    big_tiers = {f"G{i}": ("reachable" if i % 4 == 0 else "unknown")
                 for i in range(120)}
    res = build_reach_compare(big_scores, big_tiers)
    assert res["gate"] is True
    path = os.path.join(str(tmp_path), "reach_compare.png")
    out = render_reach_compare_png(res["rows"], "sorafenib",
                                   res["judged"], res["total"], path)
    assert out == path and os.path.getsize(path) > 10_000


def test_missing_tier_defaults_unknown():
    """tiers 缺记录的基因归 unknown(不静默丢、不虚构档位)。"""
    scores = {"A": 2.0, "B": 1.0}
    res = build_reach_compare(scores, {"A": "reachable"})
    b = res["rows"][res["rows"]["gene"] == "B"].iloc[0]
    assert b["tier"] == "unknown" and b["w"] == 0.5
