"""暴露约束前后靶点排序对比(方案 §3.2 招牌图,W7;2026-11-27 负责人核可三条款)。

前 = 一期等权 RWR 排序(药物靶点为种子,种子得分降序);后 = 同一得分 × w(R)
分档常数权重(可达 1.0 / 部分 0.6 / 未知 0.5 / 不可达 0.2——w_unknown=0.5 为
1.7 校准前临时口径,网格必测点见方案 §3.2/清单 §1.7)。
渲染门槛(混合规则,§3.2 ②):判档数 ≥ max(5, ceil(5%×N)) 才出图;低于门槛
只回结论行摘要——低覆盖药乘法权重有"少数记录碾压"病理(metformin 3/168
案例:OCT1/SLC22A1 被挤出前三,§3.2 ③)。
图内文字英文(matplotlib 默认字体无中文,与 venn 渲染同款纪律);Figure API
非 pyplot(线程安全,run_ppi_analysis 调用线程渲染)。"""

import logging
import math
import os

logger = logging.getLogger(__name__)

COMPARE_SCHEMA = "reach_compare.v1"

# w(R) 分档常数权重(B3 提案定稿初值;1.7 网格校准前临时口径)
PLAN_WEIGHTS = {"reachable": 1.0, "partial": 0.6, "unknown": 0.5, "unreachable": 0.2}
# 渲染门槛混合规则(§3.2 ②):绝对下限 + 百分比取大
GATE_FLOOR = 5
GATE_PCT = 0.05

# 图内配色(v1 demo 同款;rgba 由 matplotlib 十六进制直用)
TIER_COLORS = {"reachable": "#2e7d32", "partial": "#ef6c00",
               "unknown": "#9e9e9e", "unreachable": "#c62828"}


def _weights() -> dict:
    """w(R) 常数权重:config.REACH_TIER_WEIGHTS 已合入(提案 B),缺省按方案
    定稿值兜底(reachability.reach_weights 同款纪律,此处不重复告警)。"""
    try:
        from WY import config
        if hasattr(config, "REACH_TIER_WEIGHTS"):
            return dict(getattr(config, "REACH_TIER_WEIGHTS"))
    except Exception:
        pass
    return dict(PLAN_WEIGHTS)


def build_reach_compare(seed_scores: dict, tiers: dict, *,
                        weights: dict | None = None,
                        gate_floor: int = GATE_FLOOR,
                        gate_pct: float = GATE_PCT) -> dict:
    """等权 vs 加权两版排序(纯函数,渲染与 IO 解耦):
    seed_scores = {gene: rwr_score}(药物靶点 RWR 种子自身得分);
    tiers = {gene: tier}(reachability.v1.1 逐靶点档,无记录归 unknown)。
    返回 schema/ok/gate/judged/total/rows(按加权后名次排序的 DataFrame)/
    demoted(一期 Top-10 中被挤出且 unknown 的基因行)。"""
    import pandas as pd

    w = dict(weights) if weights else _weights()
    rows = []
    for gene, score in seed_scores.items():
        tier = str(tiers.get(gene) or "unknown")
        wt = float(w.get(tier, w.get("unknown", 0.5)))
        rows.append({"gene": str(gene), "tier": tier, "w": wt,
                     "score": float(score), "score_w": float(score) * wt})
    total = len(rows)
    judged = sum(1 for r in rows if r["tier"] != "unknown")
    gate = total > 0 and judged >= max(int(gate_floor),
                                       math.ceil(float(gate_pct) * total))
    for r in rows:
        r["rank_before"] = None
        r["rank_after"] = None
    order_b = sorted(rows, key=lambda r: r["score"], reverse=True)
    for i, r in enumerate(order_b, 1):
        r["rank_before"] = i
    order_a = sorted(rows, key=lambda r: r["score_w"], reverse=True)
    for i, r in enumerate(order_a, 1):
        r["rank_after"] = i
    rows.sort(key=lambda r: r["rank_after"])
    df = pd.DataFrame(rows, columns=["gene", "tier", "w", "score", "score_w",
                                     "rank_before", "rank_after"])
    # demoted:一期 Top-10 中被加权重排挤出前 10、且自身无证据(unknown)——
    # 结论栏的"谁让位了"名单(招牌图叙事核心,metformin 案例即 OCT1)
    demoted = [r for r in rows if r["tier"] == "unknown"
               and r["rank_before"] is not None and r["rank_before"] <= 10
               and (r["rank_after"] or 0) > 10]
    return {"schema": COMPARE_SCHEMA, "ok": True, "gate": gate,
            "judged": judged, "total": total, "rows": df, "demoted": demoted}


def render_reach_compare_png(rows, drug: str, judged: int, total: int,
                             path: str) -> str:
    """v2 斜率图(负责人认可的形态):加权重排后 Top-15 全标注 + demoted 名单 +
    How to read 结论栏;Figure API 线程安全,图内英文。"""
    from matplotlib.figure import Figure

    top = rows.head(15).to_dict("records")
    demoted = [r for r in rows.to_dict("records")
               if r["tier"] == "unknown" and r["rank_before"] <= 10
               and r["rank_after"] > 10][:4]
    fig = Figure(figsize=(11.8, 7.2), dpi=150, facecolor="white")
    ax = fig.add_axes([0.06, 0.11, 0.60, 0.75])
    ax.set_title(f"{drug}: target ranking before vs after exposure weighting",
                 fontsize=13, pad=14)
    ax.invert_yaxis()
    max_b = max((r["rank_before"] for r in top), default=1)
    for i, r in enumerate(top, 1):
        c = TIER_COLORS.get(r["tier"], "#9e9e9e")
        ax.plot([r["rank_before"], r["rank_after"]], [i, i], c=c, lw=2.6,
                alpha=0.9, marker="o", ms=7, zorder=3)
        ax.text(-0.012 * max_b, i, f"{i:>2}. {r['gene']}", ha="right",
                va="center", fontsize=10.5, color=c, fontweight="bold")
        shift = r["rank_before"] - r["rank_after"]
        tag = f"up {shift}" if shift > 0 else ("same" if shift == 0 else f"down {-shift}")
        ax.text(r["rank_before"], i - 0.28,
                f"{r['rank_before']} → {r['rank_after']}  ({tag})",
                va="center", fontsize=8.2, color="#666666")
    ax.set_xlim(-0.16 * max_b, max_b * 1.12)
    ax.set_ylim(len(top) + 0.9, 0.4)
    ax.set_xlabel("rank before → rank after (1 = best; phase-1 equal-weight "
                  "RWR vs × w(R))", fontsize=10)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)

    n_ev = sum(1 for r in top if r["tier"] != "unknown")
    takeaway = (
        f"After weighting, {n_ev}/{len(top)} top targets\n"
        "have exposure-backed activity evidence.\n\n"
        "Phase-1 equal-weight ranking puts hubs first regardless\n"
        f"of evidence — demoted out of phase-1 top-10\n"
        " (all no activity record):\n"
        + "".join(f"  {r['gene']}  #{r['rank_before']} → #{r['rank_after']}\n"
                  for r in demoted) +
        f"\nw(R): reachable 1.0 / partial 0.6 / unknown 0.5 /\n"
        f"unreachable 0.2. Judged {judged}/{total}.")
    fig.text(0.685, 0.82, "How to read", fontsize=12.5, fontweight="bold")
    fig.text(0.685, 0.76, takeaway, fontsize=9.3, color="#333333", va="top",
             linespacing=1.55)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, facecolor="white")
    logger.debug("前后对比图落盘: %s", path)
    return path


def save_reach_compare(directory: str, drug: str, seed_scores: dict,
                       tiers: dict) -> dict:
    """编排:纯函数建表 → 门槛 → 渲染落盘 → 摘要。返回
    {"ok", "gate", "judged", "total", "reach_compare_png"(过门槛才有,绝对路径),
    "reason"(不过门槛/失败令牌)}。调用方自行 try/except,本函数只 log 不抛。"""
    try:
        res = build_reach_compare(seed_scores, tiers)
        judged, total = res["judged"], res["total"]
        if not res["gate"]:
            logger.debug("前后对比图低于渲染门槛(判档 %d/%d),只回结论行", judged, total)
            return {"schema": COMPARE_SCHEMA, "ok": True, "gate": False,
                    "judged": judged, "total": total, "reason": "below_gate"}
        path = os.path.join(directory, "reach_compare.png")
        render_reach_compare_png(res["rows"], drug, judged, total, path)
        return {"schema": COMPARE_SCHEMA, "ok": True, "gate": True,
                "judged": judged, "total": total, "reach_compare_png": path}
    except Exception as e:  # 对比图失败绝不推翻主报告(与其他装饰性产物同纪律)
        logging.exception("reach compare step failed — skipped gracefully")
        return {"schema": COMPARE_SCHEMA, "ok": False, "gate": False,
                "judged": 0, "total": 0, "reason": f"exception: {type(e).__name__}"}
