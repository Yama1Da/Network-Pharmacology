"""分离度协同预测数学核心(二期方案 §2.1;契约 synergy.v1,2026-11-17 定稿)。

分工红线(文件清单 §1.3):本文件只做数学核心——s_AB 分离度、s_AB 度匹配
零分布、E3 预检、z_A/z_B/z_AB 对一期邻近度检验的编排复用;E1 分层排序、
E6 联合排序、输出格式化归队友(静默风险点强制本人审查)。

2026-11-17 定案(与清单文档同步):
  - z 路线照抄一期度匹配置换(Guney 框架),"解析近似"否决——无可引公式,
    一期 degree_matched_sample/proximity_zscore 现成且文献同源;
  - API 走 jobs 提交-轮询(三层置换 z 单对分钟级,同步请求否决);
  - 不做结果缓存(单对分钟级可接受;距离矩阵函数保持独立,将来加进程内
    memo 是 5 行改动,契约/API 不动);
  - geometry 标签只述 Menche 几何事实(overlap/separated/tie),不做
    协同/拮抗断言——"协同"声明只出自 DrugComb 基准(队友,§1.8)。

复用纪律(陷阱 #4/#7):z_A/z_B/z_AB 一律调一期 proximity_zscore(Guney
度匹配置换,自带经验 p),s_AB 零分布复用 degree_matched_sample,禁止复制
实现;LCC 上下文/_multi_source_bfs 借用一期同包私有函数(公共门面缺位,
若一期后续提供再迁移)。
不可达约定(陷阱 #2):沿用一期——距离一律在背景互作组 LCC 内算,LCC 内
必然可达;LCC 外靶点不进距离但计入 coverage 报告,交集外不静默丢(陷阱 #8)。
E3 预检(陷阱 #5):入 LCC 靶点 <5 抛 SynergyPrecheckError 拒绝(不是
return None);任一侧 <15 → 结果带 power_limited=true。门檻按入 LCC 数
(网络有效规模),LCC 外靶点算不得距离,按总数放行是答辩雷。
"""

import logging
import os
from datetime import datetime, timezone

from WY.core.network_algorithms import (
    _build_degree_bins,
    _lcc_context,
    _multi_source_bfs,
    degree_matched_sample,
    normalize_genes,
    proximity_zscore,
)

logger = logging.getLogger(__name__)

# 契约版本锚点:字段变更必须升版本并在清单文档 §1.3 留纪要(冻结纪律 §4)。
# v1.1(2026-11-18):meta.run 记录页级运行参数(靶点集截断 N/排序规则 id/proximity
# 置换轮数)——切片在调用方完成(定案 #7),核心不截断只记账;progress 只发 fraction
# 不带文案(核心无文案纪律),阶段码由 API 层按区间标注。
SYNERGY_SCHEMA = "synergy.v1.1"


class SynergyPrecheckError(RuntimeError):
    """E3 预检拒绝(入 LCC 靶点 <5)。reason 为稳定枚举,API 层按它映射 i18n 文案:
    too_few_targets_a / too_few_targets_b。"""

    def __init__(self, reason: str, detail: str):
        super().__init__(f"{reason}: {detail}")
        self.reason = reason
        self.detail = detail


def _mean_pair_within(Mt, idx, xs: list[str]) -> tuple[float, int]:
    """集内平均对距离(C(n,2) 对,LCC 内必然可达)。返回 (mean, n_pairs)。"""
    import numpy as np

    xs_idx = [idx[g] for g in xs]
    rows = [_multi_source_bfs(Mt, np.array([i])) for i in xs_idx]
    total = pairs = 0
    for i in range(len(xs_idx)):
        for j in range(i + 1, len(xs_idx)):
            total += int(rows[i][xs_idx[j]])
            pairs += 1
    if not pairs:  # 预检保证 >=5,理论不可达;防御None返回类型
        return float("nan"), 0
    return total / pairs, pairs


def _mean_pair_across(Mt, idx, xs: list[str], ys: list[str]) -> float:
    """跨集平均对距离(n_x×n_y 对;同节点跨集出现距离 0 计入,与一期
    _multi_source_bfs 源点置 0 约定一致——重叠靶点把 s_AB 拉向 overlap 是定义内行为)。"""
    import numpy as np

    xs_idx = [idx[g] for g in xs]
    ys_idx = np.array([idx[g] for g in ys])
    total = 0
    for i in xs_idx:
        total += int(_multi_source_bfs(Mt, np.array([i]))[ys_idx].sum())
    return total / (len(xs_idx) * len(ys_idx))


def _separation(Mt, idx, a_lcc: list[str], b_lcc: list[str]) -> float:
    """Menche 2015 网络分离度:s_AB = d_AB − (d_AA + d_BB)/2。
    全项目最危险的一行代码(清单陷阱 #1):<0 = 两药靶点模块重叠(同区作用),
    >0 = 分离(互补作用)。写反符号则全体颠倒且不报错,改动必须有构造小图对照。"""
    d_aa, _ = _mean_pair_within(Mt, idx, a_lcc)
    d_bb, _ = _mean_pair_within(Mt, idx, b_lcc)
    d_ab = _mean_pair_across(Mt, idx, a_lcc, b_lcc)
    return d_ab - (d_aa + d_bb) / 2.0


def _precheck_e3(a_lcc: list[str], b_lcc: list[str]) -> bool:
    """E3 三档门槛(陷阱 #5):<5 抛业务异常拒绝;5–14 放行并标 power_limited;
    >=15 为 normal。返回 power_limited。"""
    for side, genes in (("a", a_lcc), ("b", b_lcc)):
        if len(genes) < 5:
            raise SynergyPrecheckError(
                f"too_few_targets_{side}",
                f"入网靶点 {len(genes)} < 5,统计不可信,拒绝计算",
            )
    return not (len(a_lcc) >= 15 and len(b_lcc) >= 15)


def _analyze_pair_from_context(
    ctx: tuple, a_genes: list, b_genes: list, d_genes: list | None, *,
    drug_a: str = "", drug_b: str = "", disease_name: str = "",
    disease_source: str = "", species: str = "Human",
    interactome_meta: dict | None = None,
    n_perm: int = 200, seed: int = 0, bin_size: int = 100,
    top_n: int | None = None, sort_rule: str | None = None,
    proximity_n_perm: int | None = None, progress=None,
) -> dict:
    """契约 synergy.v1 组装(纯计算,LCC 上下文注入,便于小图自检——一期同模式)。
    d_genes 为空/None 时疾病为 None,z_* 三字段全 None(ML/CLI 批量可无疾病)。"""
    import numpy as np

    nodes, idx, Mt, degree_map = ctx
    a_norm = normalize_genes(a_genes)
    b_norm = normalize_genes(b_genes)
    d_norm = normalize_genes(d_genes) if d_genes else []
    # 入网映射:交集外计数进 coverage,不静默丢(陷阱 #8;一期 coverage 同约定)
    a_lcc = sorted(g for g in a_norm if g in idx)
    b_lcc = sorted(g for g in b_norm if g in idx)
    d_lcc = sorted(g for g in d_norm if g in idx)
    power_limited = _precheck_e3(a_lcc, b_lcc)
    if len(nodes) and (set(a_lcc) | set(b_lcc)) - set(nodes):
        raise AssertionError("入 LCC 靶点越界(上下文与映射不一致)")  # 防御,理论不可达

    # ---- s_AB 观测值 + 度匹配零分布(方案 §2.1:对 T_A/T_B 重采样) ----
    d_aa, _ = _mean_pair_within(Mt, idx, a_lcc)
    d_bb, _ = _mean_pair_within(Mt, idx, b_lcc)
    d_ab = _mean_pair_across(Mt, idx, a_lcc, b_lcc)
    s_obs = d_ab - (d_aa + d_bb) / 2.0
    rng = np.random.default_rng(seed)
    bins = _build_degree_bins(degree_map, bin_size)  # 分箱建一次循环复用(一期同优化)
    s_null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        ra = degree_matched_sample(degree_map, a_lcc, bin_size, rng, bins=bins)
        rb = degree_matched_sample(degree_map, b_lcc, bin_size, rng, bins=bins)
        s_null[i] = _separation(Mt, idx, ra, rb)
        if progress:
            progress(fraction=0.1 + 0.45 * (i + 1) / n_perm)  # 仅分数,文案归 API 层
    mu, sd = float(s_null.mean()), float(s_null.std())
    if sd > 0:
        z_sep = float((s_obs - mu) / sd)
    else:
        # 零分布退化:z 输出 None 不伪装 z=0(一期邻近度同约定)
        z_sep = None
        logger.info("synergy s_AB 零分布退化(sd=0),z_sep=None")
    p_lower = (1 + int((s_null <= s_obs).sum())) / (n_perm + 1)  # 重叠方向单尾,一期经验 p 同式

    # ---- z_A/z_B/z_AB:一期邻近度检验编排复用(陷阱 #7,禁止复制实现) ----
    # proximity_zscore 自带覆盖率保护/内部置换;传 None 让它读 PROXIMITY_* 配置。
    if d_norm:
        if progress:
            progress(fraction=0.6)
        prox: dict = {}
        for key, genes in (("z_a", a_norm), ("z_b", b_norm),
                           ("z_ab", list(dict.fromkeys(a_norm + b_norm)))):
            r = proximity_zscore(genes, d_norm, species=species,
                                 n_perm=proximity_n_perm)
            prox[key] = {"ok": r.get("ok"), "error": r.get("error"),
                         "d_c": r.get("d_c"), "z": r.get("z"),
                         "p": r.get("p_empirical")}
        if progress:
            progress(fraction=0.95)
    else:
        prox = {"z_a": None, "z_b": None, "z_ab": None}

    # ---- geometry 标签:只述 Menche 几何事实,不做协同/拮抗断言(2026-11-17 定案) ----
    if s_obs < 0:
        label, basis = "overlap", "s_ab < 0"
    elif s_obs > 0:
        label, basis = "separated", "s_ab > 0"
    else:
        label, basis = "tie", "s_ab == 0"

    return {
        "schema": SYNERGY_SCHEMA,
        "pair": {"drug_a": drug_a, "drug_b": drug_b},
        "disease": ({"name": disease_name, "n_targets": len(d_norm),
                     "source": disease_source} if d_norm else None),
        "inputs": {
            "coverage": {"a_total": len(a_norm), "a_in_lcc": len(a_lcc),
                         "b_total": len(b_norm), "b_in_lcc": len(b_lcc),
                         "d_total": len(d_norm), "d_in_lcc": len(d_lcc)},
            "n_overlap_ab": len(set(a_norm) & set(b_norm)),
            "power_limited": power_limited,
            "interactome": interactome_meta or {},
        },
        "separation": {
            "s_ab": s_obs, "d_aa": d_aa, "d_bb": d_bb, "d_ab": d_ab,
            "null": {"mean": mu, "sd": sd, "z": z_sep,
                     "p_lower": p_lower, "n_perm": int(n_perm)},
        },
        "proximity": prox,
        "geometry_label": label,
        "label_basis": basis,
        "ranks": None,       # E1 分层排序(队友):重叠型按 z_ab、互补型按分离强度各排 Top-K
        "reachability": None,  # 1.2 浓度可达性模块挂载位
        "meta": {"computed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                 "engine": "synergy.py",
                 "run": {"top_n": top_n, "sort": sort_rule,
                         "proximity_n_perm": proximity_n_perm}},
    }


def results_root() -> str:
    """协同运行归档根(结果管理页/归档写入共用):DOWNLOADS_DIR/synergy-results/。"""
    from WY import config

    return os.path.join(config.DOWNLOADS_DIR, "synergy-results")


def list_results() -> list[dict]:
    """运行清单(结果管理页,新→旧):每 run 读 result.json 出行级摘要;损坏或缺
    result.json 的条目不丢弃——带 ok=false 与 reason 出行,让用户能看到坏档并手动清。"""
    root = results_root()
    if not os.path.isdir(root):
        return []
    rows: list[dict] = []
    for name in sorted(os.listdir(root), reverse=True):
        path = os.path.join(root, name)
        if not os.path.isdir(path) or name.startswith("."):
            continue
        mtime = datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M:%S")
        row = {"run": name, "mtime": mtime, "ok": False, "reason": None,
               "drug_a": "", "drug_b": "", "disease": "", "geometry_label": "",
               "s_ab": None, "power_limited": None, "computed_at": ""}
        rf = os.path.join(path, "result.json")
        if not os.path.isfile(rf):
            row["reason"] = "缺 result.json"
            rows.append(row)
            continue
        try:
            import json

            with open(rf, encoding="utf-8") as f:
                data = json.load(f)
            row.update(ok=True,
                       drug_a=str(data.get("pair", {}).get("drug_a", "")),
                       drug_b=str(data.get("pair", {}).get("drug_b", "")),
                       disease=str((data.get("disease") or {}).get("name", "")),
                       geometry_label=str(data.get("geometry_label", "")),
                       s_ab=data.get("separation", {}).get("s_ab"),
                       power_limited=bool(data.get("inputs", {}).get("power_limited")),
                       computed_at=str(data.get("meta", {}).get("computed_at", "")))
        except (OSError, ValueError):
            row["reason"] = "result.json 解析失败"
            logging.exception("synergy 归档损坏: %s", name)
        rows.append(row)
    return rows


def resolve_run(name: str) -> str:
    """校验运行名(纯文件名、synergy_ 前缀、在清单内)并返回目录绝对路径;
    越权/不存在抛 ValueError,API 层转 400。"""
    if not name or name != os.path.basename(name) or not name.startswith("synergy_"):
        raise ValueError(f"非法运行名: {name!r}")
    path = os.path.join(results_root(), name)
    if not os.path.isdir(path):
        raise ValueError(f"运行不存在: {name}")
    return path


def read_run(name: str) -> dict:
    """读单次运行的契约 JSON(结果详情页/复现回填用)。"""
    import json

    path = os.path.join(resolve_run(name), "result.json")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def delete_runs(names: list[str]) -> dict:
    """删除勾选的运行目录:逐个校验,被拒/IO 失败的进 skipped(不中断整批)。"""
    import shutil

    deleted: list[str] = []
    skipped: list[dict] = []
    for n in names:
        try:
            path = resolve_run(n)
        except ValueError as exc:
            skipped.append({"name": n, "reason": str(exc)})
            continue
        try:
            shutil.rmtree(path)
            deleted.append(n)
        except OSError as exc:
            skipped.append({"name": n, "reason": str(exc)})
    return {"deleted": deleted, "skipped": skipped}


def zip_runs(names: list[str]) -> bytes:
    """把勾选的运行目录打包成 ZIP 字节流(批量下载);越权名抛 ValueError 上抛转 400。"""
    import io
    import zipfile

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for n in names:
            root = resolve_run(n)
            for base, _dirs, files in os.walk(root):
                for f in files:
                    full = os.path.join(base, f)
                    zf.write(full, os.path.join(n, os.path.relpath(full, root)))
    return buf.getvalue()


def analyze_pair(
    drug_a_targets: list, drug_b_targets: list, disease_targets: list | None = None, *,
    drug_a: str = "", drug_b: str = "", disease_name: str = "",
    disease_source: str = "", species: str = "Human",
    n_perm: int | None = None, seed: int | None = None,
    top_n: int | None = None, sort_rule: str | None = None,
    proximity_n_perm: int | None = None, progress=None,
) -> dict:
    """分离度协同预测主入口(API/CLI 共用;CLI 批量即对本函数按 pairs 循环)。
    top_n/sort_rule 仅记账进 meta.run——切片在调用方完成(定案 #7),核心不截断。
    proximity_n_perm=None 时一期邻近度读 PROXIMITY_N_PERM 配置默认。
    背景互作组未构建时抛 InteractomeNotBuiltError,由调用方转用户提示。"""
    from WY import config
    from WY.core import interactome as itx

    n_perm = int(n_perm if n_perm is not None else config.SYNERGY_N_PERM)
    seed = int(seed if seed is not None else config.SYNERGY_SEED)
    meta = itx.require_meta(species)  # 未构建 → InteractomeNotBuiltError 透传
    interactome_meta = {
        "species": species,
        "string_version": str(meta.get("string_version", "")),
        "confidence": int(meta.get("confidence", 0)),
        "built_at": str(meta.get("built_at", "")),
    }
    result = _analyze_pair_from_context(
        _lcc_context(species), drug_a_targets, drug_b_targets, disease_targets,
        drug_a=drug_a, drug_b=drug_b, disease_name=disease_name,
        disease_source=disease_source, species=species,
        interactome_meta=interactome_meta, n_perm=n_perm, seed=seed,
        bin_size=int(config.PROXIMITY_DEGREE_BIN_SIZE),
        top_n=top_n, sort_rule=sort_rule, proximity_n_perm=proximity_n_perm,
        progress=progress,
    )
    logger.info("synergy %s×%s: s_ab=%.4f(%s, p_lower=%.3f, power_limited=%s)",
                drug_a or "A", drug_b or "B", result["separation"]["s_ab"],
                result["geometry_label"], result["separation"]["null"]["p_lower"],
                result["inputs"]["power_limited"])
    return result
