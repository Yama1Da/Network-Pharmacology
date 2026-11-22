"""疾病模块多算法一致性扩展(二期方案 §4;文件清单 §1.4,队友实现)。

定位(方案 §4):产品化层空集,主张限定在"多方法一致性排序 + 置换检验
显著性 + 端到端平台化"。RWR 零改动复用一期(清单 §1.4 陷阱:禁止复制
实现),DIAMOnD 严格按原文自实现(度校正超几何),近邻度扩张为最朴素
对照。纯函数库:ctx=(nodes, idx, Mt, degree_map) 由调用方经
network_algorithms._lcc_context(species) 注入;留一法验证在
benchmark/disease_module_benchmark.py(1.8,队友)。

陷阱钉死(清单 §1.4,实现处再各钉一次):
  1) DIAMOnD 每步重算"候选节点与当前模块的连接数"超几何显著性——
     先正确后优化(逐步全候选重算,O(N·候选·度),真实图分钟级;优化走
     增量更新,改前先立基准测试防静默错)
  2) scipy.stats.hypergeom 参数顺序:M=全图节点数, n=模块内节点数,
     N=候选节点度数——填错位不报错但结果错(tests 有手算对照用例)
  3) 种子集先与互作组节点取交集,交集外基因计数入报告(不静默丢)
  4) 三方法同一个 N、同一个输入,否则共识无可比性
  5) 共识平均秩统一"秩从 1 开始"
"""

import logging
import os

from WY.core.reachability import _plan_or_config, _warn_once

logger = logging.getLogger(__name__)

# 模块输出 schema 锚点:列变更必须升版本并在清单文档留纪要(冻结纪律,
# 同 reachability.REACHABILITY_SCHEMA 模式)
CONSENSUS_SCHEMA = "consensus.v1"


def _module_size(n=None) -> int:
    """统一扩张规模 N(§4.1 定稿默认 100,可调,进参数快照;三方法共用——
    陷阱#4)。显式传参优先(敏感性/小图测试),None 时读 config,未合入按
    方案定稿值兜底。"""
    if n is not None:
        return int(n)
    return int(_plan_or_config("MODULE_EXPAND_N", 100, "模块扩张规模 N"))


def _consensus_k(k=None) -> int:
    """"被 ≥k 种方法捞入"标记的 k(§4.2 默认 2,可调);显式传参优先。"""
    if k is not None:
        return int(k)
    return int(_plan_or_config("CONSENSUS_MIN_METHODS", 2, "一致性标记 k"))


def _graph_from_ctx(ctx) -> "nx.Graph":
    """ctx → nx.Graph(仅 LCC 节点与边;STRING 边在本管线为二值——PPI_CONFIDENCE
    已在上游阈值化,近邻度扩张的\u201c加权度\u201d因此退化为连接数,§4.1 朴素对照语义)。"""
    import networkx as nx

    nodes, _idx, Mt, _degree_map = ctx
    g = nx.Graph()
    g.add_nodes_from(nodes)
    coo = Mt.tocoo()
    for i, j in zip(coo.row, coo.col):
        if i < j:  # 对称矩阵只取上三角,防重复加边
            g.add_edge(nodes[i], nodes[j])
    return g


def _seeds_intersect(ctx, seeds) -> tuple[list[str], list[str]]:
    """种子集与 LCC 取交集(陷阱#3):交集内基因序返回,交集外**计数入报告**
    不静默丢(§1.4-3;synergy coverage 同约定)。"""
    _nodes, idx, _Mt, _degree_map = ctx
    seen: set[str] = set()
    inside: list[str] = []
    outside: list[str] = []
    for s in seeds:
        g = str(s).strip().upper()
        if not g or g in seen:
            continue
        seen.add(g)
        (inside if g in idx else outside).append(g)
    return sorted(inside), sorted(outside)


def _diamond_pvalue(k_out: int, n_mod: int, k: int, M: int) -> float:
    """DIAMOnD 超几何尾部 p(陷阱#2 参数顺序手算对照的唯一定义点):
    P(X >= k_out | M=全图节点数, n=模块内节点数, N=候选节点度数 k)。
    候选的 k 条连线在\u201c随机\u201d落位下落入当前模块 s 条的概率模型。"""
    from scipy.stats import hypergeom

    return float(hypergeom.sf(k_out - 1, M, n_mod, k))


def diamond_expand(ctx, seeds, n=None) -> dict:
    """DIAMOnD 扩张(Ghiassian 2015,严格按原文逐步度校正超几何):
    每步对全部候选重算超几何尾部 p(陷阱#1:先正确后优化,不做增量优化),
    取最小者入模。tie-break 为原文未定义的复现性微决策,钉死:
    p 升序 → k_out 降序(连接更多者优先,贴近原文意图)→ 基因字典序
    (可复现红线 §1.2-8 同款)。种子先交集(陷阱#3)。"""
    nodes, idx, _Mt, degree_map = ctx
    import networkx as nx

    g = _graph_from_ctx(ctx)
    adjacency = {v: set(g[v]) for v in nodes}
    seed_in, seed_out = _seeds_intersect(ctx, seeds)
    n_target = int(_module_size(n))
    module = list(seed_in)          # 种子在前(基因序)
    mod_set = set(module)
    join_order: list[str] = []
    candidates = [v for v in nodes if v not in mod_set]
    while len(module) < n_target and candidates:
        M = len(nodes)
        n_mod = len(mod_set)
        best_key, best_v = None, None
        for v in candidates:
            k_out = len(adjacency[v] & mod_set)
            key = (_diamond_pvalue(k_out, n_mod, int(degree_map[v]), M),
                   -k_out, v)
            if best_key is None or key < best_key:
                best_key, best_v = key, v
        module.append(best_v)
        mod_set.add(best_v)
        join_order.append(best_v)
        candidates.remove(best_v)
    return {"module": module, "join_order": join_order,
            "n_added": len(join_order),
            "coverage": {"seeds_in_lcc": len(seed_in),
                         "seeds_total": len(seed_in) + len(seed_out),
                         "seeds_out_of_lcc": seed_out}}


def neighbor_expansion(ctx, seeds, n=None) -> dict:
    """近邻度扩张(方案 §4.1 最朴素对照):按对种子集的直接连接数降序扩张。
    兜底规则(钉注释):连接数 ≥1 的候选不足时,余量按全图度数降序补齐
    ——为保三方法同一个 N 的可比性(陷阱#4),不是算法语义的一部分。"""
    nodes, idx, _Mt, degree_map = ctx
    import networkx as nx

    g = _graph_from_ctx(ctx)
    adjacency = {v: set(g[v]) for v in nodes}
    seed_in, seed_out = _seeds_intersect(ctx, seeds)
    n_target = int(_module_size(n))
    module = list(seed_in)
    mod_set = set(module)
    seed_set = set(seed_in)
    candidates = [v for v in nodes if v not in mod_set]
    scored = sorted(candidates,
                    key=lambda v: (-len(adjacency[v] & seed_set), -degree_map[v], v))
    need = max(0, n_target - len(module))
    join_order = [v for v in scored if len(adjacency[v] & seed_set) >= 1][:need]
    if len(join_order) < need:  # 兜底补齐(保同 N 可比,见 docstring)
        picked = set(join_order)
        join_order += [v for v in scored if v not in picked][:need - len(join_order)]
    module += join_order
    return {"module": module, "join_order": join_order,
            "n_added": len(join_order),
            "coverage": {"seeds_in_lcc": len(seed_in),
                         "seeds_total": len(seed_in) + len(seed_out),
                         "seeds_out_of_lcc": seed_out}}


def rwr_module(ctx, seeds, n=None) -> dict:
    """RWR 取 top-N 为模块(方案 §4.1:一期 network_algorithms 零改动复用,
    陷阱:禁止复制实现——本函数只做编排)。alpha/tol/max_iter 读一期
    config.RWR_*;秩 = RWR 得分降序(种子行在前)。"""
    from WY import config
    from .network_algorithms import _rwr_from_context

    res = _rwr_from_context(ctx, seeds, config.RWR_ALPHA, config.RWR_TOL,
                            config.RWR_MAX_ITER, top_nonseed=None)
    seed_in, seed_out = _seeds_intersect(ctx, seeds)
    if not res.get("ok"):
        return {"ok": False, "error": res.get("error"), "module": [],
                "join_order": [], "n_added": 0,
                "coverage": {"seeds_in_lcc": len(seed_in),
                             "seeds_total": len(seeds),
                             "seeds_out_of_lcc": seed_out}}
    n_target = int(_module_size(n))
    seed_in = [row["gene"] for row in res["seeds"]]
    take = max(0, n_target - len(seed_in))
    module = seed_in + [row["gene"] for row in res["top_nonseed"][:take]]
    return {"ok": True, "module": module,
            "join_order": [row["gene"] for row in res["top_nonseed"][:take]],
            "n_added": min(take, len(res["top_nonseed"])),
            "coverage": {"seeds_in_lcc": int(res["coverage"]["seed_in_lcc"]),
                         "seeds_total": int(res["coverage"]["seed_total"]),
                         "seeds_out_of_lcc": seed_out}}


def _method_ranks(module: list[str], join_order: list[str]) -> dict:
    """模块成员 → 秩(陷阱#5:秩从 1 开始):种子在前(基因序),扩张按
    加入/得分序——三方法共用同一秩语义,共识才可比。"""
    head = sorted(set(module) - set(join_order))
    order = head + [g for g in module if g in set(join_order)]
    return {g: i + 1 for i, g in enumerate(dict.fromkeys(order))}


def consensus_targets(ctx, seeds, *, n=None, k=None) -> dict:
    """三方法 + Borda 共识(§4.2)。返回:
      schema / modules(三方法同 N 同输入,陷阱#4)/ consensus 表
      (gene, rwr_rank, diamond_rank, neighbor_rank, n_methods,
       consensus_score, consistent, is_seed, 各方法成员布尔列)/ coverage。
    秩从 1 开始(陷阱#5);非成员 = NaN;consensus_score = 所属方法秩的
    平均(Borda);consistent = n_methods ≥ k(默认 2)。
    consensus 表序:consensus_score 升序(更优在前)→ n_methods 降序 → 基因序。"""
    import pandas as pd

    dm = diamond_expand(ctx, seeds, n)
    nb = neighbor_expansion(ctx, seeds, n)
    rw = rwr_module(ctx, seeds, n)
    seed_in, seed_out = _seeds_intersect(ctx, seeds)
    if not rw.get("ok"):
        # 错误值走稳定枚举(一期 empty_after_filter 同款),文案归 API 层按
        # 枚举映射 i18n——core 词典纪律:test_architecture 扫描 "error": 中文字面量
        logger.error("疾病模块 RWR 失败: %s", rw.get("error"))
        return {"schema": CONSENSUS_SCHEMA, "ok": False,
                "error": rw.get("error") or "rwr_failed", "consensus": None,
                "modules": {"diamond": dm["module"], "neighbor": nb["module"],
                            "rwr": []},
                "coverage": dm["coverage"]}
    ranks = {"rwr": _method_ranks(rw["module"], rw["join_order"]),
             "diamond": _method_ranks(dm["module"], dm["join_order"]),
             "neighbor": _method_ranks(nb["module"], nb["join_order"])}
    in_mod = {m: set(ranks[m]) for m in ranks}
    kk = int(_consensus_k(k))
    union = sorted(set().union(*in_mod.values()))
    seed_set = set(seed_in)
    rows = []
    for gene in union:
        present = [ranks[m][gene] for m in ("rwr", "diamond", "neighbor")
                   if gene in in_mod[m]]
        rows.append({
            "gene": gene,
            "rwr_rank": ranks["rwr"].get(gene),
            "diamond_rank": ranks["diamond"].get(gene),
            "neighbor_rank": ranks["neighbor"].get(gene),
            "n_methods": len(present),
            "consensus_score": (float(sum(present) / len(present))
                                if present else None),
            "consistent": len(present) >= kk,
            "is_seed": gene in seed_set,
        })
    df = pd.DataFrame(rows).sort_values(
        ["consensus_score", "n_methods", "gene"],
        ascending=[True, False, True], na_position="last").reset_index(drop=True)
    coverage = dict(dm["coverage"])
    coverage["seeds_out_of_lcc"] = sorted(
        set(dm["coverage"]["seeds_out_of_lcc"])
        | set(nb["coverage"]["seeds_out_of_lcc"]))
    return {"schema": CONSENSUS_SCHEMA, "ok": True, "error": None,
            "modules": {"rwr": rw["module"], "diamond": dm["module"],
                        "neighbor": nb["module"]},
            "module_sizes": {"rwr": len(rw["module"]),
                             "diamond": len(dm["module"]),
                             "neighbor": len(nb["module"])},
            "consensus": df, "coverage": coverage, "k": kk}


def module_significance(ctx, seeds, *, n=None, n_perms=None,
                        null_seed=None) -> dict:
    """三方法模块显著性(§4.2,度匹配零模型,一期框架复用):
    观测统计 = 模块内部边数(cohesion);零分布 = degree_matched_sample
    重抽种子集 → 同一方法扩张 → 内部边数;经验单侧 p = (1+#{null>=obs})
    /(perms+1)(富集方向,一期邻近度经验 p 同式)。
    n_perms 默认 RWR_NULL_PERMS,null_seed 默认 RWR_NULL_SEED(一期同源)。
    成本注记(§12.1-E4 同款):DIAMOnD 零分布 = perms × 逐步扩张,真实图
    分钟级起,批量前先小 perms 试跑(陷阱#1 先正确后优化)。"""
    import numpy as np
    from .network_algorithms import _build_degree_bins, degree_matched_sample

    nodes, idx, Mt, degree_map = ctx
    seed_in, _out = _seeds_intersect(ctx, seeds)
    n_target = int(_module_size(n))
    from WY import config

    perms = int(n_perms) if n_perms is not None else int(config.RWR_NULL_PERMS)
    rng = np.random.default_rng(int(null_seed) if null_seed is not None
                                else int(config.RWR_NULL_SEED))
    g = _graph_from_ctx(ctx)
    bins = _build_degree_bins(degree_map, 100)

    def cohesion(members: set) -> int:
        return int(sum(1 for u, v in g.edges() if u in members and v in members))

    expanders = {
        "rwr": lambda ss: rwr_module(ctx, ss, n_target)["module"],
        "diamond": lambda ss: diamond_expand(ctx, ss, n_target)["module"],
        "neighbor": lambda ss: neighbor_expansion(ctx, ss, n_target)["module"],
    }
    out = {}
    for name, expand in expanders.items():
        module = expand(seed_in)
        obs = cohesion(set(module))
        nulls = np.empty(perms, dtype=float)
        for i in range(perms):
            sample = degree_matched_sample(degree_map, seed_in, 100, rng, bins=bins)
            nulls[i] = cohesion(set(expand(sorted(sample))))
        p = float((1 + int((nulls >= obs).sum())) / (perms + 1))
        out[name] = {"observed_cohesion": obs, "null_mean": float(nulls.mean()),
                     "null_sd": float(nulls.std()), "n_perms": perms,
                     "p": p}
    logger.info("疾病模块显著性: %s", {m: v["p"] for m, v in out.items()})
    return {"schema": CONSENSUS_SCHEMA, "n": n_target, "n_perms": perms,
            "methods": out, "seed_in_lcc": len(seed_in)}


def export_consensus_csv(df, path: str) -> str:
    """consensus_targets.csv 落盘(§4.2 输出物;原子写,照 caching 模式)。"""
    from .caching import save_csv_cache

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    save_csv_cache(df, path, skip_empty=False)
    return path
