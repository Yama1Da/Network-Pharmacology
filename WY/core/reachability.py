"""靶点浓度可达性:归一化 + R 值 + 加权排序(二期方案 §3.1 第 2/3 层)。

契约 reachability.v1.1(2026-11-22 冻结;2026-11-25 升版——exposure 增
synonym_fallback_sid 字段 + 名字解析改身份优先 + 冲突台账,纪要见清单 §1.2;
字段/枚举变更必须升版本并在清单文档留纪要,同 synergy.py SYNERGY_SCHEMA 纪律):
  - 分档枚举(§12.2-B2 定稿;unknown 为独立档,永不与 unreachable 合并——
    B3 给 unknown 0.5 > unreachable 0.2 的权重差,前端"未知"独立状态由此保证):
    reachable(R >= 1) / partial(0.1 <= R < 1) / unreachable(R < 0.1) / unknown
  - unknown_reason ∈ {no_exposure, no_activity, no_molweight}
    (无 Cmax 记录 / 无活性记录 / 缺分子量无法换算——清单 §1.2 陷阱 1/5)
  - E6 消费(§12.1-E6 定稿):协同结果表加 reachability_a/reachability_b 两列
    (取值为分档枚举);exposure_confident = 两药均 ∈ {reachable, partial}
  - 单药聚合口径(2026-11-22 随合同冻结):任一可计算 R 的靶点达档即取其
    最高权档(任一靶点可达 → 药可达);全部不可计算 → unknown + 原因
  - 排序惩罚(2026-11-22 随合同冻结,E1/E6 外圈实现仍属强制审查点):
    min(w_a, w_b) 乘入排序分(瓶颈逻辑,与高置信规则单调一致)
  - w(R) 权重(§12.2-B3 定稿初值):1.0/0.6/0.5/0.2——数值经 config
    REACH_TIER_WEIGHTS 供给(常量草案待负责人合入),1.7 基准网格微调后
    须负责人复核签字才改,不得手调
  - fub 缺失 → 保守上界 1.0 +"可能高估"标记,绝不当 0(除零 → inf)
  - 活性(摩尔) = 10 ** (-pChEMBL);方向写反则 R 全体倒置且不报错
    (清单 §1.2 陷阱 2,靠手算对照兜住——验收三案例见清单 §5.6)
  - v1.1 增补(2026-11-25,负责人定稿):exposure 增 synonym_fallback_sid
    (null=主 sid 直接命中;否则=实际取数所用的同义兄弟 sid,显式标记不静默)。
    名字解析改**身份优先**(名字==自身 sid/正式名 无条件胜过他物质同义词声明,
    与节点序无关;2026-11-25 实测 12/770 行经争议名解析、其中 repaglinide 2 行
    误挂致幻剂 5-meo-dalt);查名冲突台账落 _conflict_ledger.json(构建期写入),
    查询侧同义兄弟兜底(compute_reachability 的 synonym_fallback 参数,默认开:
    主 sid 无暴露记录才触发,按台账竞争 sid 依序取数)——归属正确性与数据
    覆盖率两全,兜底显式可见

分工红线(文件清单 §1.2):本文件=队友实现;w(R) 标定产出与验收自检
(二甲双胍/雷米普利手算对照、缺失标记两案例)须工程负责人复核签字后
方可合入(§5.4-9)。

分层纪律(架构红线,清单 §0.3):本文件属归一化层——读 pkdb_raw/ 镜像
(2026-11-17 探针定案"数字化源表直下":<sid>.json 的 tsv[].rows 为逐字
字符串二维表,Tab*.tsv 含 measurement_type=cmax/auc_end/clearance 行 +
intervention 列回链给药方案),换算但不联网;单位换算只发生在本文件;
下游只读本文件产出的版本化 Parquet,永不碰原始 JSON/TSV。

输出(提案 C,2026-11-22 送审;照 interactome.py 模式:文件名带版本 +
to_parquet + tmp + os.replace):
  substance_map_v{PK_NORM_VERSION}.parquet / cmax_norm_v{...} / fub_norm_v{...}

开口待负责人定(2026-11-22 提案三开口,相应逻辑在定案前不动工):
  1) ChEMBL 缓存无 standard_type(chembl.py 只存 assay_type+pchembl_value),
     B1"Ki/Kd 优先于 IC50/EC50"暂无法实现——扩缓存 vs 记录偏离,二选一
  2) 分子量来源:2026-11-22 实测 PK-DB substance.mass 与 xrefs(chembl) 自带
     质量与 ChEMBL ID——候选解=镜像数据直取(已按此实现,待负责人复核
     覆盖面);若覆盖不足再议镜像层补抓(归一化层不联网红线不动)
  3) assay_type F/B → 细胞/生化 映射为代理,待生物侧确认(§5.5-5 TODO)
"""

import logging
import os

logger = logging.getLogger(__name__)

# 契约版本锚点:枚举/字段变更必须升版本并在清单文档留纪要(冻结纪律)
# v1.1(2026-11-25):exposure 增 synonym_fallback_sid(见模块 docstring 契约块)
REACHABILITY_SCHEMA = "reachability.v1.1"

# 分档枚举(契约冻结值;B2 定稿三档 + unknown 显式档)
REACH_TIERS = ("reachable", "partial", "unreachable", "unknown")
# E6 高置信档集合:两药均落入此集合才打高置信标记(§12.1-E6 定稿)
_CONFIDENT_TIERS = frozenset(("reachable", "partial"))
# unknown 原因枚举(契约冻结值)
REACH_UNKNOWN_REASONS = ("no_exposure", "no_activity", "no_molweight")


def _norm_dir() -> str:
    """归一化 Parquet 落盘目录:CACHE_DIR/reachability/(沙箱随 CACHE_DIR 走,
    与 pk_mirror._raw_dir 同模式;目录名属提案 C 送审项,冻结前不改)。"""
    from WY import config

    return os.path.join(config.CACHE_DIR, "reachability")


def _versioned_path(stem: str) -> str:
    """版本化文件路径:<stem>_v{PK_NORM_VERSION}.parquet(TARGET_FILTER_VERSION
    同款逻辑:改归一化逻辑必递增版本,旧缓存自动失效)。PK_NORM_VERSION 由
    config 供给(常量草案待负责人合入);缺失时报错拒绝而非静默默认。"""
    from WY import config

    version = getattr(config, "PK_NORM_VERSION", None)
    if version is None:
        raise RuntimeError(
            "config.PK_NORM_VERSION 未合入(常量草案待负责人),归一化缓存不可构建")
    return os.path.join(_norm_dir(), f"{stem}_v{int(version)}.parquet")


# 方案定稿值与 config 的统一口径:config 未合入(提案 B)按方案定稿值兜底并
# 告警一次;config 已合入但与方案定稿值不一致 → 拒绝(防口径漂移;
# §3.2 敏感性对照走显式参数,不改 config,故分箱阈值不做此强校验)
_WARNED_ONCE = set()


def _warn_once(key: str, msg: str) -> None:
    if key not in _WARNED_ONCE:
        _WARNED_ONCE.add(key)
        logger.warning(msg)


def _plan_or_config(attr: str, plan_value, key: str) -> float:
    from WY import config

    if hasattr(config, attr):
        v = float(getattr(config, attr))
        if v != float(plan_value):
            raise RuntimeError(
                f"config.{attr}={v} 与方案定稿值 {plan_value} 不一致,拒绝({key})")
        return v
    _warn_once(f"cfg:{attr}", f"config.{attr} 未合入(提案 B),按方案定稿值 "
                             f"{plan_value} 运行({key})")
    return float(plan_value)


def exposure_confident(tier_a: str, tier_b: str) -> bool:
    """E6 高置信标记(§12.1-E6 定稿):两药均 ≥ 部分可达。契约冻结件,纯函数;
    非法枚举值显式报错(静默 False 会把打错的档伪装成"不高置信")。"""
    for t in (tier_a, tier_b):
        if t not in REACH_TIERS:
            raise ValueError(f"非法分档枚举 {t!r}(契约 {REACHABILITY_SCHEMA})")
    return tier_a in _CONFIDENT_TIERS and tier_b in _CONFIDENT_TIERS


def _read_json(path: str) -> dict:
    """显式 utf-8 读 JSON(清单 §0.1)。"""
    import json

    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _load_substance_nodes() -> list[dict]:
    """读镜像层 substance_map.json(pk_mirror 落盘的 info_nodes/substance 原文,
    信封双层 data 已由镜像层剥好,nodes 为纯列表)。文件缺失或 nodes 为空时
    报错拒绝——空映射会让 cmax 全部静默变 unknown,宁可停也不空转。"""
    import json

    from WY import config

    path = os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR,
                        "substance_map.json")
    if not os.path.exists(path):
        raise RuntimeError(f"镜像层物质映射缺失: {path}(先跑 pk_mirror.build_mirror)")
    with open(path, encoding="utf-8") as f:  # 显式 utf-8:清单 §0.1,Windows 默认 GBK
        data = json.load(f)
    nodes = data.get("nodes") or []
    if not nodes:
        raise RuntimeError(f"substance_map.json 无节点(镜像异常): {path}")
    return nodes


def _substance_rows(nodes: list[dict]) -> "pd.DataFrame":
    """info_nodes(substance) 原文 → substance_map 表(提案 C schema + 实数据
    修订列 dtype/deprecated/parent_sids,2026-11-22)。纯函数,不读盘不联网。
    字段映射(对 2026-11-22 实抓节点核对,非猜测):
      sid/name/dtype/deprecated 直取;synonyms 为字符串列表;
      annotations[].term 形如 "CHEBI:6073"(collection=chebi)→ chebi_ids;
      xrefs[].name=="chembl" 的 accession(如 CHEMBL113)→ chembl_id;
      parents[].sid → parent_sids(复合物质=代谢物总和,母药判定第 5 步用);
      substance.mass → mol_weight_g_mol(PK-DB 自带质量;2026-11-22 实测
      metformin/ramipril 等均带值——开口 2 由数据自解,无需改镜像层,
      归一化不联网红线不破;缺失按 None 走"无法换算",不近似硬算 §1.2-1)。
    NaN 纪律(§0.4):mol_weight_g_mol 进 DataFrame 后缺失即 NaN——下游
    判定一律 pd.notna/pd.isna,禁止 is None;该列 NaN = "无法换算"显式分档,
    绝不参与任何数值运算兜底。
    """
    import pandas as pd

    rows: list[dict] = []
    seen: set[str] = set()
    dup = 0
    for n in nodes:
        sid = str(n.get("sid", "") or "")
        if not sid:
            continue
        if sid in seen:  # 重复 sid 保留首条并计数,不静默吞
            dup += 1
            continue
        seen.add(sid)
        anns = n.get("annotations") or []
        chebi = sorted({str(a.get("term", "")).strip() for a in anns
                        if str(a.get("term", "")).upper().startswith("CHEBI:")})
        xrefs = n.get("xrefs") or []
        chembl = sorted({str(x.get("accession", "")).strip() for x in xrefs
                         if str(x.get("name", "")).lower() == "chembl"
                         and x.get("accession")})
        mass = (n.get("substance") or {}).get("mass")
        try:
            mass_f = float(mass) if mass is not None else None
        except (TypeError, ValueError):
            mass_f = None  # 非数值质量按缺失走"无法换算",不近似硬算
        rows.append({
            "sid": sid,
            "name": str(n.get("name", "") or ""),
            "dtype": str(n.get("dtype", "") or ""),
            "deprecated": bool(n.get("deprecated", False)),
            "synonyms": ";".join(str(s) for s in (n.get("synonyms") or [])),
            "chebi_ids": ";".join(chebi),
            "chembl_id": ";".join(chembl),
            "parent_sids": ";".join(str(p.get("sid", "") or "")
                                    for p in (n.get("parents") or [])),
            "mol_weight_g_mol": mass_f,
        })
    if dup:
        logger.warning("substance_map 重复 sid %d 个(保留首条)", dup)
    return pd.DataFrame(rows, columns=[
        "sid", "name", "dtype", "deprecated", "synonyms",
        "chebi_ids", "chembl_id", "parent_sids", "mol_weight_g_mol"])


def _save_table_atomic(df: "pd.DataFrame", path: str) -> str:
    """版本化表落盘(照 interactome._save_edge_table:优先 parquet,to_parquet +
    tmp + os.replace 原子替换;失败回退 CSV 并记异常)。返回实际写出路径。"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    try:
        tmp = path + ".tmp"
        df.to_parquet(tmp, index=False)
        os.replace(tmp, path)
        return path
    except Exception:
        logger.exception("Parquet 写盘失败,回退 CSV: %s", path)
        try:
            os.remove(path + ".tmp")
        except OSError:
            pass
    csv_path = os.path.splitext(path)[0] + ".csv"
    tmp = csv_path + ".tmp"
    df.to_csv(tmp, index=False, encoding="utf-8")
    os.replace(tmp, csv_path)
    return csv_path


def build_substance_map(force: bool = False) -> str:
    """substance_map_v{PK_NORM_VERSION}.parquet 构建(第 4 步):读镜像层
    substance_map.json → 纯函数解析 → 版本化落盘。force=False 且目标文件
    已存在时跳过(幂等续跑)。返回实际路径。"""
    out = _versioned_path("substance_map")
    if not force and os.path.exists(out):
        logger.debug("substance_map 命中已有版本,跳过: %s", out)
        return out
    df = _substance_rows(_load_substance_nodes())
    path = _save_table_atomic(df, out)
    ledger_path = _write_conflict_ledger(_load_substance_nodes())
    logger.info("substance_map 归一化完成: %d 物质 → %s(冲突台账: %s)",
                len(df), path, ledger_path)
    return path


# ---- cmax 归一化(第 5 步):TSV 参数表 → 单位统一到摩尔 → 中位数聚合 ----

# 单位规范表(§1.2 陷阱 1:PK-DB 单位五花八门;2026-11-22 全量 722 研究实测分布
# ng/ml 2255 / nmol/l 202 / microg/ml 155 / µg/ml 62 / mg/l 50 / µm 47 …)。
# molar 族因子直接到 mol/L;mass 族因子到 g/L,再除摩尔质量(缺口=无法换算,
# 不许近似硬算)。新单位一律先进表再算,禁止行内硬编码换算。
_MOLAR_FACTOR = {
    "mol/l": 1.0, "mmol/l": 1e-3, "mm": 1e-3,
    "µmol/l": 1e-6, "µm": 1e-6, "µmol/ml": 1e-3,
    "nmol/l": 1e-9, "nm": 1e-9, "nmol/ml": 1e-6,
    "pmol/l": 1e-12, "pmol/ml": 1e-9,
}
_MASS_FACTOR = {  # 因子 → g/L
    "g/l": 1.0, "mg/l": 1e-3, "mg/ml": 1.0, "kg/l": 1000.0,
    "µg/l": 1e-9, "ng/l": 1e-12, "pg/l": 1e-15,
    "µg/ml": 1e-6, "ng/ml": 1e-6, "pg/ml": 1e-9,
}
_UNIT_ALIASES = {  # 拼写归一(先做别名映射再查因子表)
    "umol/l": "µmol/l", "µmole/l": "µmol/l", "umole/l": "µmol/l",
    "micromol/l": "µmol/l", "um": "µm", "umol/ml": "µmol/ml",
    "ug/l": "µg/l", "microg/l": "µg/l", "microgr/l": "µg/l",
    "ug/ml": "µg/ml", "microg/ml": "µg/ml", "microgr/ml": "µg/ml",
    "umol": "µmol/l", "nmole/l": "nmol/l",
}
_TRACER_MARKERS = ("c14", "tracer", "microdose", "scintillation")
_RADIO_DOSE_UNITS = ("bq", "kbq", "mbq", "gbq", "uci", "mci", "ci")


def _canon_unit(raw: str) -> str:
    """单位拼写归一:小写、去空格、希腊 mu 统一为 µ、别名收敛(只做拼写,
    不做数值换算——换算只发生在因子表,防千倍错 §1.2-1)。"""
    u = str(raw or "").strip().lower().replace("\u03bc", "\u00b5")
    u = "".join(u.split())
    return _UNIT_ALIASES.get(u, u)


def _parse_num(raw) -> float | None:
    """TSV 数值解析:'NA'/'na'/空 → None;其余 float()。"""
    s = str(raw or "").strip()
    if not s or s.lower() in ("na", "nan", "none"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _build_name_index(nodes: list[dict]) -> tuple[dict, dict[str, list[str]], int]:
    """物质名索引(身份优先两遍构建,2026-11-25 定稿):
    第一遍——身份声明:名字精确等于自身 sid 或正式名的,无条件占位;
    第二遍——名字+同义词补位:已被他 sid 占用的名字不覆盖,记入台账
    (key → 竞争 sid 列表,依遭遇序)并计数。
    返回 (lookup: 名字 → (sid, mol_weight_g_mol), ledger: 争议名台账, 冲突计数)。
    为什么不是单纯"保留先到":节点序=镜像构建时上游返回序,既无语义优先级
    又随重拉漂移——2026-11-25 实测 12/770 行经争议名解析,其中 repaglinide
    2 行误挂 5-meo-dalt(上游脏同义词)。身份声明与顺序无关,修此错类;
    剩余争议(真不同物质共享同义词,如 oxazepam-glucuronide 声称 diazepam)
    本就不该消化,留台账供查询侧兜底与审计。"""
    lookup: dict[str, tuple[str, float | None]] = {}
    mw_by_sid: dict[str, float | None] = {}
    for n in nodes:
        sid = str(n.get("sid", "") or "")
        if not sid:
            continue
        mass = (n.get("substance") or {}).get("mass")
        try:
            mw_by_sid[sid] = float(mass) if mass is not None else None
        except (TypeError, ValueError):
            mw_by_sid[sid] = None
        own = str(n.get("name", "") or "").strip().lower()
        for key in {sid.lower(), own} - {""}:
            if key not in lookup:
                lookup[key] = (sid, mw_by_sid[sid])
    ledger: dict[str, list[str]] = {}
    conflicts = 0
    for n in nodes:
        sid = str(n.get("sid", "") or "")
        if not sid:
            continue
        names = [str(n.get("name", "") or "")] + [str(s) for s in (n.get("synonyms") or [])]
        for name in names:
            key = name.strip().lower()
            if not key:
                continue
            holder = lookup.get(key)
            if holder is not None:
                if holder[0] != sid:
                    conflicts += 1
                    ledger.setdefault(key, [])
                    if sid not in ledger[key]:
                        ledger[key].append(sid)
                continue
            lookup[key] = (sid, mw_by_sid[sid])
    if conflicts:
        # 告警放原语内:构建期与运行期(compute 直调本函数)两条路都要响,
        # 每进程一次(2026-11-25 教训:告警留在包装里曾致运行路径静默)
        _warn_once("substance_name_conflict",
                   "substance 名称映射冲突 %d 处(身份优先已定归属;竞争者明细见"
                   "归一化缓存目录 _conflict_ledger.json)" % conflicts)
    return lookup, ledger, conflicts


def _substance_lookup(nodes: list[dict]) -> dict:
    """物质名/同义词 → (sid, mol_weight_g_mol) 查找表(_build_name_index 的
    薄包装;构建侧只需 lookup,查询侧兜底用 _build_name_index 拿台账)。"""
    return _build_name_index(nodes)[0]


def _write_conflict_ledger(nodes: list[dict]) -> str:
    """查名冲突台账落盘(构建期,substance_map 产物随行;审计与查询侧兜底共用):
    _norm_dir()/_conflict_ledger.json,含 winner/claimants 与版本戳。"""
    from datetime import datetime, timezone

    from WY import config
    from .caching import save_json_file

    _lookup, ledger, conflicts = _build_name_index(nodes)
    payload = {
        "pk_norm_version": int(getattr(config, "PK_NORM_VERSION", 0) or 0),
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "conflicts_total": conflicts,
        "keys": {k: {"winner_sid": _lookup.get(k, (None,))[0], "claimants": v}
                 for k, v in sorted(ledger.items())},
    }
    path = os.path.join(_norm_dir(), "_conflict_ledger.json")
    save_json_file(path, payload, indent=1)
    return path


def _collect_dosing(study: dict) -> dict:
    """单研究给药方案表:TabInterventions 类表(name 列) → {name: 剂量元组}。
    供治疗剂量过滤做 intervention 关联(2026-11-22 实测:125 研究可关联,
    17 研究组名不命中——不命中者保留并打标,不剔除)。"""
    dosing: dict[str, dict] = {}
    for t in study.get("tsv", []):
        rows = t.get("rows") or []
        if len(rows) < 2:
            continue
        header = rows[0]
        if "measurement_type" not in header or "route" not in header or "name" not in header:
            continue
        mi, ni = header.index("measurement_type"), header.index("name")
        for r in rows[1:]:
            if len(r) <= mi or r[mi].strip().lower() != "dosing":
                continue

            def col(c):
                return r[header.index(c)] if c in header and len(r) > header.index(c) else ""
            name = r[ni].strip()
            if name and name not in dosing:
                dosing[name] = {
                    "dose_value": _parse_num(col("value")),
                    "dose_unit": col("unit").strip().lower(),
                    "route": col("route").strip().lower(),
                    "form": col("form").strip().lower(),
                }
    return dosing


def _is_tracer(substance_raw: str, intervention: str, method: str,
               dose_unit: str) -> bool:
    """示踪剂/微剂量判定(结构性规则,非业务阈值拍板):放射性剂量单位,
    或物质/干预/方法字段带 C14/tracer/microdose/示踪计数标记。
    更严格的"治疗剂量阈值"属方案级决定,本层不自行拍板(§5.4-7)。"""
    if _canon_unit(dose_unit) in _RADIO_DOSE_UNITS:
        return True
    blob = " ".join((substance_raw, intervention, method)).lower()
    return any(m in blob for m in _TRACER_MARKERS)


def _cmax_observations() -> tuple["pd.DataFrame", dict]:
    """遍历 pkdb_raw 全部研究,抽取 cmax 参数行并归一(纯读取,不落盘不联网)。
    返回 (聚合表, 统计 dict)。数值列择优:value > mean > median(实测
    Anderson1961 个体行 value 有值而 mean=NA;组行反之)。聚合键 =
    (substance_key, intervention),组内取中位数(§3.1 多观测取中位数),
    dropna=False 保证未匹配物质行不静默丢(§1.2-5/7)。"""
    import glob as _glob

    import numpy as np
    import pandas as pd

    lookup = _substance_lookup(_load_substance_nodes())
    from WY import config

    raw_dir = os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR)
    rows: list[dict] = []
    stats = {"studies": 0, "cmax_rows": 0, "tracer_excluded": 0, "no_value": 0,
             "unknown_unit": 0, "no_molweight": 0, "converted": 0,
             "dose_unmatched": 0}
    for path in sorted(_glob.glob(os.path.join(raw_dir, "*.json"))):
        fname = os.path.basename(path)
        # sid 命名非统一前缀(实测 645 个 PKDB* + 157 个数字 sid):只排除
        # 镜像层元文件,其余全部纳入,不按前缀猜(2026-11-22 修正,防静默漏数据)
        if fname in ("_meta.json", "substance_map.json"):
            continue
        try:
            study = _read_json(path)
        except Exception:
            logger.warning("研究文件读取失败,跳过: %s", fname)
            continue
        stats["studies"] += 1
        sid = fname[:-5]
        dosing = _collect_dosing(study)
        for t in study.get("tsv", []):
            trows = t.get("rows") or []
            if len(trows) < 2 or "measurement_type" not in trows[0]:
                continue
            header = trows[0]

            def col(r, c):
                return r[header.index(c)] if c in header and len(r) > header.index(c) else ""
            for r in trows[1:]:
                if len(r) <= header.index("measurement_type") \
                        or r[header.index("measurement_type")].strip().lower() != "cmax":
                    continue
                stats["cmax_rows"] += 1
                substance_raw = str(col(r, "substance") or "").strip()
                intervention = str(col(r, "intervention") or "").strip()
                method = str(col(r, "method") or "").strip()
                dose = dosing.get(intervention)
                if not dose:
                    stats["dose_unmatched"] += 1
                    dose = {"dose_value": None, "dose_unit": "", "route": "", "form": ""}
                if _is_tracer(substance_raw, intervention, method, dose.get("dose_unit", "")):
                    stats["tracer_excluded"] += 1
                    continue
                value = _parse_num(col(r, "value"))
                if value is None:
                    value = _parse_num(col(r, "mean"))
                if value is None:
                    value = _parse_num(col(r, "median"))
                if value is None:
                    stats["no_value"] += 1
                    continue
                unit_raw = str(col(r, "unit") or "").strip()
                unit = _canon_unit(unit_raw)
                sub_sid, matched, mw = None, "unmatched", None
                key = substance_raw.lower()
                if key in lookup:
                    sub_sid, mw = lookup[key]
                    matched = "name"
                u_conv = "unknown_unit"
                molar = None
                if unit in ("na", ""):
                    u_conv = "na_unit"
                    stats["unknown_unit"] += 1
                elif unit in _MOLAR_FACTOR:
                    molar = value * _MOLAR_FACTOR[unit]
                    u_conv = "molar_direct"
                    stats["converted"] += 1
                elif unit in _MASS_FACTOR:
                    if mw:
                        molar = value * _MASS_FACTOR[unit] / mw
                        u_conv = "mass_via_mw"
                        stats["converted"] += 1
                    else:
                        u_conv = "no_molweight"
                        stats["no_molweight"] += 1
                else:
                    stats["unknown_unit"] += 1
                disp = "none"
                for c in ("sd", "se", "sem", "min", "max"):
                    if _parse_num(col(r, c)) is not None:
                        disp = "reported"
                        break
                rows.append({
                    "substance_sid": sub_sid,
                    "substance_raw": substance_raw,
                    "match_method": matched,
                    "intervention": intervention,
                    "group": str(col(r, "group") or "").strip(),
                    "tissue": str(col(r, "tissue") or "").strip(),
                    "value_raw": value,
                    "unit_raw": unit_raw,
                    "cmax_molar_M": molar,
                    "unit_conv_status": u_conv,
                    "dispersion_flag": disp,
                    "dose_value": dose["dose_value"],
                    "dose_unit": dose["dose_unit"],
                    "route": dose["route"],
                    "exclusion_flag": "therapeutic_kept",
                    "licence": str(study.get("licence", "") or ""),
                    "pkdb_sid": sid,
                })
    df = pd.DataFrame(rows, columns=[
        "substance_sid", "substance_raw", "match_method", "intervention", "group",
        "tissue", "value_raw", "unit_raw", "cmax_molar_M", "unit_conv_status",
        "dispersion_flag", "dose_value", "dose_unit", "route", "exclusion_flag",
        "licence", "pkdb_sid"])
    # 中位数聚合(§3.1):键 = 物质 × 干预(不同剂量不能混聚);
    # dropna=False 防未匹配物质键被 groupby 静默丢;回填 reset_index(§1.2-7)
    df["substance_key"] = df["substance_sid"].fillna(
        "raw:" + df["substance_raw"].str.lower())
    agg = (df[df["cmax_molar_M"].notna()]
           .groupby(["substance_key", "intervention"], dropna=False)
           .agg(cmax_molar_M=("cmax_molar_M", "median"),
                n_obs=("cmax_molar_M", "size"),
                cmax_min=("cmax_molar_M", "min"),
                cmax_max=("cmax_molar_M", "max"),
                substance_raw=("substance_raw", "first"),
                substance_sid=("substance_sid", "first"),
                match_method=("match_method", "first"),
                group=("group", "first"),
                tissue=("tissue", "first"),
                dispersion_flag=("dispersion_flag", "max"),
                dose_value=("dose_value", "first"),
                dose_unit=("dose_unit", "first"),
                route=("route", "first"),
                licence=("licence", "first"),
                source_sids=("pkdb_sid", lambda s: ";".join(sorted(set(s)))))
           .reset_index())
    logger.info("cmax 归一化统计: %s", stats)
    return agg, stats


def build_cmax_norm(force: bool = False) -> str:
    """cmax_norm_v{PK_NORM_VERSION}.parquet 构建(第 5 步):全量研究 cmax 行
    → 单位统一摩尔 → (物质×干预) 中位数聚合 → 版本化落盘。"""
    out = _versioned_path("cmax_norm")
    if not force and os.path.exists(out):
        logger.debug("cmax_norm 命中已有版本,跳过: %s", out)
        return out
    agg, _stats = _cmax_observations()
    if agg.empty:
        raise RuntimeError("cmax_norm 聚合结果为空(镜像数据异常,拒绝写空表)")
    path = _save_table_atomic(agg, out)
    logger.info("cmax_norm 归一化完成: %d 行 → %s", len(agg), path)
    return path


# ---- fub 归一化(第 6 步):PK-DB fu 采样 → 缺失保守上界 1.0 + 可能高估标记 ----

def _fub_observations() -> tuple["pd.DataFrame", dict]:
    """遍历 pkdb_raw,抽取 fu 相关参数行并换算成游离分数(纯读取,不落盘)。
    2026-11-22 全量实测:PK-DB 的 fu 极稀缺——fraction_unbound 13 行(约 5 个
    物质,直接小数)+ plasma_binding 19 行(百分数,需 1 − 结合率/100);
    方案 §3.1"PK-DB 为主"在数据层面不成立(已报负责人),本层照常实现,
    缺口由 FDA/DailyMed 线与保守上界承接。
    返回 (按物质聚合表, 统计 dict);来源优先级:直接 fraction_unbound >
    plasma_binding 派生(混合时只取前者,注释钉死)。"""
    import glob as _glob

    import pandas as pd

    lookup = _substance_lookup(_load_substance_nodes())
    from WY import config

    raw_dir = os.path.join(config.CACHE_DIR, config.PKDB_CACHE_SUBDIR)
    rows: list[dict] = []
    stats = {"fraction_unbound": 0, "plasma_binding": 0, "no_value": 0,
             "unknown_unit": 0, "out_of_range": 0}
    for path in sorted(_glob.glob(os.path.join(raw_dir, "*.json"))):
        fname = os.path.basename(path)
        if fname in ("_meta.json", "substance_map.json"):
            continue
        try:
            study = _read_json(path)
        except Exception:
            logger.warning("研究文件读取失败,跳过: %s", fname)
            continue
        for t in study.get("tsv", []):
            trows = t.get("rows") or []
            if len(trows) < 2 or "measurement_type" not in trows[0]:
                continue
            header = trows[0]

            def col(r, c):
                return r[header.index(c)] if c in header and len(r) > header.index(c) else ""
            for r in trows[1:]:
                if len(r) <= header.index("measurement_type"):
                    continue
                mt = r[header.index("measurement_type")].strip().lower()
                if mt not in ("fraction_unbound", "plasma_binding"):
                    continue
                stats[mt] += 1
                value = _parse_num(col(r, "value"))
                if value is None:
                    value = _parse_num(col(r, "mean"))
                if value is None:
                    value = _parse_num(col(r, "median"))
                if value is None:
                    stats["no_value"] += 1
                    continue
                substance_raw = str(col(r, "substance") or "").strip()
                key = substance_raw.lower()
                sub_sid, matched = None, "unmatched"
                if key in lookup:
                    sub_sid, matched = lookup[key][0], "name"
                unit = _canon_unit(col(r, "unit"))
                if mt == "fraction_unbound":
                    fub, u_conv = value, "fraction_direct"
                else:  # plasma_binding:结合率 → 游离分数
                    if unit in ("percent", "%"):
                        fub, u_conv = 1.0 - value / 100.0, "percent_derived"
                    elif unit in ("", "fraction"):
                        fub, u_conv = 1.0 - value, "fraction_derived"
                    else:
                        u_conv = "unknown_unit"
                        stats["unknown_unit"] += 1
                        fub = None
                if fub is not None and not (0.0 < fub <= 1.0):
                    u_conv = "out_of_range"
                    stats["out_of_range"] += 1
                    fub = None
                rows.append({
                    "substance_sid": sub_sid, "substance_raw": substance_raw,
                    "match_method": matched, "method": mt,
                    "value_raw": value, "unit_raw": str(col(r, "unit") or "").strip(),
                    "fub_value": fub, "conv_status": u_conv,
                    "group": str(col(r, "group") or "").strip(),
                    "tissue": str(col(r, "tissue") or "").strip(),
                    "pkdb_sid": fname[:-5],
                })
    df = pd.DataFrame(rows, columns=[
        "substance_sid", "substance_raw", "match_method", "method", "value_raw",
        "unit_raw", "fub_value", "conv_status", "group", "tissue", "pkdb_sid"])
    ok = df[df["fub_value"].notna()]
    direct = ok[ok["method"] == "fraction_unbound"]
    derived = ok[ok["method"] == "plasma_binding"]
    direct_sids = set(direct["substance_sid"].dropna())
    # 优先级钉死:有直接 fraction_unbound 的物质不混入派生值(两类口径不同)
    derived = derived[~derived["substance_sid"].isin(direct_sids)]
    use = pd.concat([direct, derived])

    def _agg(g: "pd.DataFrame") -> "pd.Series":
        # 分组键 substance_sid 由 reset_index() 从索引回填(§1.2-7 同款);
        # include_groups=False 下分组列不在 g 内,勿从 g 取键
        return pd.Series({
            "fub_value": float(g["fub_value"].median()),
            "n_obs": int(len(g)),
            "method": ("fraction_unbound"
                       if (g["method"] == "fraction_unbound").any() else "plasma_binding"),
            "conv_status": ";".join(sorted(set(g["conv_status"]))),
            "substance_raw": str(g["substance_raw"].iloc[0]),
            "source_sids": ";".join(sorted(set(g["pkdb_sid"]))),
        })

    agg = (use.groupby("substance_sid", dropna=False)
           .apply(_agg, include_groups=False)
           .reset_index(drop=False))
    logger.info("fub 采样统计: %s", stats)
    return agg, stats


def build_fub_norm(force: bool = False) -> str:
    """fub_norm_v{PK_NORM_VERSION}.parquet 构建(第 6 步):物质全集 = substance_map;
    PK-DB 采到的用真值(中位数),其余缺失 → 保守上界 1.0 + "可能高估"标记
    (§3.1 定稿;fub 绝不能当 0,除零 → inf,清单 §1.2-4)。"""
    out = _versioned_path("fub_norm")
    if not force and os.path.exists(out):
        logger.debug("fub_norm 命中已有版本,跳过: %s", out)
        return out
    import pandas as pd

    from WY import config

    agg, _stats = _fub_observations()
    substances = _substance_rows(_load_substance_nodes())[["sid"]]
    fub_by_sid = {r["substance_sid"]: r for r in agg.to_dict("records")
                  if r["substance_sid"] is not None}
    default_fub = _plan_or_config("FUB_DEFAULT", 1.0, "fub 保守上界")
    rows = []
    for sid in substances["sid"]:
        obs = fub_by_sid.get(sid)
        if obs is not None:
            rows.append({"substance_sid": sid, "fub_value": obs["fub_value"],
                         "source": "pkdb", "imputed_flag": False,
                         "method": obs["method"], "n_obs": obs["n_obs"],
                         "source_sids": obs["source_sids"]})
        else:
            rows.append({"substance_sid": sid, "fub_value": default_fub,
                         "source": "imputed", "imputed_flag": True,
                         "method": None, "n_obs": 0, "source_sids": ""})
    df = pd.DataFrame(rows, columns=[
        "substance_sid", "fub_value", "source", "imputed_flag", "method",
        "n_obs", "source_sids"])
    # 硬断言(除零防线):fub 必须落在 (0, 1],无 NaN——fub=0 会让 Cmax_free=0/R=0
    # 与 fub 缺失不可区分(§0.4/§1.2-4)
    if df["fub_value"].isna().any() or not ((df["fub_value"] > 0) & (df["fub_value"] <= 1)).all():
        raise RuntimeError("fub_norm 出现 0/越界/NaN 值,违反保守上界约束,拒绝写盘")
    path = _save_table_atomic(df, out)
    logger.info("fub_norm 归一化完成: %d 物质(实测 %d, 保守上界 %d) → %s",
                len(df), int((df["source"] == "pkdb").sum()),
                int((df["source"] == "imputed").sum()), path)
    return path


def build_normalized_cache(force: bool = False, progress=None,
                           max_studies: int | None = None) -> dict:
    """归一化主入口(第 2 层,方案 §3.1):读 pkdb_raw/ → 产出三个版本化
    Parquet(substance_map / cmax_norm / fub_norm),返回各产物实际路径。
    force=False 且同名版本文件已存在时逐件跳过(幂等续跑);max_studies 为
    镜像侧参数,归一化恒全量读(截断会在归一化层丢数据,不实现)。"""
    paths = {
        "substance_map": build_substance_map(force=force),
        "cmax_norm": build_cmax_norm(force=force),
        "fub_norm": build_fub_norm(force=force),
    }
    logger.info("归一化三件套完成: %s", paths)
    return paths


# ---- R 值计算(第 7 步,第 3 层:只读归一化 Parquet 与既有靶点缓存,零联网) ----

# B1 活性类别映射(§12.2-B1 定稿"细胞实验优先"的落地):ChEMBL assay_type
# F(functional)→细胞、B(binding)→生化 为代理映射——开口3 待生物侧确认
# (§5.5-5 TODO,不静默断言);A(ADME)/T(毒性)/P(理化)/U(未分类)非靶点
# 活性,不进主 R 值。Ki/Kd 优先于 IC50/EC50(§12.2-B1)当前无法实现:
# ChEMBL 缓存无 standard_type 字段(开口1,已报负责人二选一)
_ACTIVITY_CLASS = {"F": "cell", "B": "biochemical"}


def r_value(cmax_total_M, fub, activity_M) -> float | None:
    """R = (Cmax × fub) ÷ 活性(均为摩尔,§3.2 定稿)。任一缺失或非正 → None
    (显式缺失不硬算——活性≤0 会造 inf/负 R,fub>0 已由 fub_norm 硬断言
    保证,此处双保险)。"""
    if cmax_total_M is None or fub is None or activity_M is None:
        return None
    if activity_M <= 0 or fub <= 0 or cmax_total_M <= 0:
        return None
    return (cmax_total_M * fub) / activity_M


def _tier(r: float, *, reachable_min=None, partial_min=None) -> str:
    """三档分箱(§12.2-B2 定稿:主定义不加安全系数 k;R=1⇔50% 占据,
    R=3⇔75%,R=9⇔90%)。显式传参供 §3.2 边界敏感性对照(R≥3/≥10 两档
    另跑报告排序位移),不传则读 config,未合入按定稿值兜底+告警。"""
    rmin = (float(reachable_min) if reachable_min is not None
            else _plan_or_config("R_REACHABLE_MIN", 1.0, "R 分档阈值"))
    pmin = (float(partial_min) if partial_min is not None
            else _plan_or_config("R_PARTIAL_MIN", 0.1, "R 分档阈值"))
    if r >= rmin:
        return "reachable"
    if r >= pmin:
        return "partial"
    return "unreachable"


def _expand_merged_activity(merged) -> "pd.DataFrame":
    """merged 缓存格式 → B1 记录长表(TARGET_FILTER_VERSION ≥ 6 起药物 merged
    带 chmbl_assay_types/chmbl_pchembl_values 明细列,分号串联、位置对齐)。
    两列长度不一致属缓存损坏——报错拒绝,不静默截断;原始格式(直接含
    pchembl_value 列)不经此函数,行为不变。"""
    import pandas as pd

    rows = []
    for rec in merged.to_dict("records"):
        t_raw, v_raw = rec.get("chmbl_assay_types"), rec.get("chmbl_pchembl_values")
        if pd.isna(t_raw) or pd.isna(v_raw):
            continue  # 无 ChEMBL 定量明细的基因(纯 DGIdb 佐证等)→ 无活性
        types = [s.strip() for s in str(t_raw).split(";")]
        vals = [s.strip() for s in str(v_raw).split(";")]
        if len(types) != len(vals):
            raise RuntimeError(
                f"merged 缓存活性明细列错位({rec.get('approvedSymbol')}: "
                f"{len(types)} 型 / {len(vals)} 值)——缓存损坏,请按现行 "
                f"TARGET_FILTER_VERSION 重建 merged 缓存")
        for t, v in zip(types, vals):
            rows.append({"approvedSymbol": rec.get("approvedSymbol"),
                         "assay_type": t, "pchembl_value": v})
    if not rows:
        return pd.DataFrame(columns=["approvedSymbol", "assay_type", "pchembl_value"])
    out = pd.DataFrame(rows)
    out["pchembl_value"] = pd.to_numeric(out["pchembl_value"], errors="coerce")
    return out


def _b1_activity(targets) -> "pd.DataFrame":
    """B1 活性聚合(§12.2-B1 定稿):同基因多 assay——细胞(F)优先取中位数,
    无细胞值才回退生化(B)并打回退标记;另一类中位数保留为敏感性对照列。
    活性(摩尔) = 10 ** (-pChEMBL)(§3.2;方向写反则 R 全体倒置且不报错,
    清单 §1.2-2,靠 §5.6 手算对照兜住)。Ki/Kd 优先规则暂缺(开口1)。
    纯函数,不读盘不联网;gene 列取一期缓存列名 approvedSymbol。
    输入双格式(2026-11-23):原始 ChEMBL 记录表,或 TARGET_FILTER_VERSION≥6
    的 merged 缓存(明细列经 _expand_merged_activity 展开后再聚合)。"""
    import pandas as pd

    cols = ["gene", "activity_M", "activity_class", "n_assays",
            "activity_alt_M", "fallback_flag"]
    if targets is None or getattr(targets, "empty", True) \
            or "approvedSymbol" not in targets.columns:
        return pd.DataFrame(columns=cols)
    if "pchembl_value" not in targets.columns:
        if "chmbl_pchembl_values" not in targets.columns:
            return pd.DataFrame(columns=cols)
        targets = _expand_merged_activity(targets)  # merged 明细列 → 记录长表
        if targets.empty:
            return pd.DataFrame(columns=cols)
    df = targets.copy()
    df["_p"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    at = df["assay_type"] if "assay_type" in df.columns else ""
    df["_cls"] = (at.map(_ACTIVITY_CLASS) if hasattr(at, "map") else None)
    df = df[df["_p"].notna() & df["_cls"].notna() & df["approvedSymbol"].notna()]
    if df.empty:
        return pd.DataFrame(columns=cols)
    rows = []
    for gene, g in df.groupby("approvedSymbol"):
        cell_p = g.loc[g["_cls"] == "cell", "_p"]
        bio_p = g.loc[g["_cls"] == "biochemical", "_p"]
        # 中位数取在 pChEMBL(log 空间)再换算摩尔——记录的原生尺度,效价聚合
        # 标准做法;若取在摩尔线性值,偶数条记录的中位数会落成"没有任何一条
        # 记录的值"(如 [1e-6,1e-7] → 5.5e-7)。此口径为静默风险点,负责人
        # 强制审查项(2026-11-22 定,偏差须改此处并递增 PK_NORM_VERSION)
        if len(cell_p):  # 细胞优先用于主 R 值(B1 定稿),生化留作敏感性对照
            main, cls, fb, n = (10.0 ** -float(cell_p.median()), "cell",
                                False, int(len(cell_p)))
            alt = 10.0 ** -float(bio_p.median()) if len(bio_p) else None
        else:  # 细胞值缺失才回退生化(回退留标记列,§1.2-3)
            main, cls, fb, n = (10.0 ** -float(bio_p.median()), "biochemical",
                                True, int(len(bio_p)))
            alt = None
        rows.append({"gene": str(gene), "activity_M": main, "activity_class": cls,
                     "n_assays": n, "activity_alt_M": alt, "fallback_flag": fb})
    return pd.DataFrame(rows, columns=cols)


def _read_norm(stem: str) -> "pd.DataFrame":
    """读归一化版本化表(parquet 优先,CSV 兜底,与 _save_table_atomic 对偶)。
    缺失时报错提示先跑 build_normalized_cache——计算层零联网,缺数据是响亮
    错误,不是静默空结果。"""
    import pandas as pd

    path = _versioned_path(stem)
    if not os.path.exists(path):
        raise RuntimeError(f"归一化缓存缺失: {path}(先跑 build_normalized_cache)")
    if path.endswith(".parquet"):
        return pd.read_parquet(path)
    return pd.read_csv(path, encoding="utf-8")


def compute_reachability(drug_targets, drug_name: str = "", *,
                         synonym_fallback: bool = True,
                         reachable_min=None, partial_min=None) -> dict:
    """单药可达性(第 3 层计算,§3.2;契约 reachability.v1.1 单药块):
    B1 活性聚合 → R = (Cmax 中位数 × fub) ÷ 活性 → 三档分箱 → 逐靶点明细
    + 药级档(任一靶点达档即取最高权档,合同 2026-11-22 冻结)。
    drug_targets:一期 ChEMBL 靶点缓存 DataFrame(调用方保证已取——本函数
    不联网,缓存未就绪是调用方职责);reachable_min/partial_min 显式传参
    供 §3.2 边界敏感性对照。Cmax 取该物质全部干预中位数(§3.1 多观测取
    中位数;跨干预聚合口径属静默风险点,负责人强制审查项)。
    synonym_fallback(2026-11-25 负责人定稿,默认开):主 sid 无暴露记录时,
    按查名冲突台账(_build_name_index 的争议 key→竞争 sid)依序取数,命中即
    用并在 exposure.synonym_fallback_sid 显式标记;主 sid 有数据永不兜底
    (身份归属不受污染)。关闭则严格只认身份归属的 sid。"""
    import pandas as pd

    unknown_reason = None
    sid = None
    cmax_total = fub = None
    fub_imputed = None
    fallback_sid = None
    cmax_meta = {"cmax_molar_M": None, "cmax_interventions": "", "n_cmax_obs": 0,
                 "fub_value": None, "fub_imputed": None,
                 "synonym_fallback_sid": None}
    key = str(drug_name or "").strip().lower()
    lookup, ledger, _conflicts = _build_name_index(_load_substance_nodes())
    if key:
        sid = lookup.get(key, (None, None))[0]
    if sid is not None:
        def _fetch(sid_try: str):
            crows = _read_norm("cmax_norm")
            crows = crows[crows["substance_sid"] == sid_try]
            frow = _read_norm("fub_norm")
            frow = frow[frow["substance_sid"] == sid_try]
            fub_v = fub_imp = None
            if not frow.empty:
                fub_v = float(frow["fub_value"].iloc[0])
                fub_imp = bool(frow["imputed_flag"].iloc[0])
            return crows, fub_v, fub_imp

        crows, fub, fub_imputed = _fetch(sid)
        if (crows.empty or fub is None) and synonym_fallback and key:
            for alt in ledger.get(key, ()):
                c2, f2, fi2 = _fetch(alt)
                if not c2.empty and f2 is not None:
                    crows, fub, fub_imputed = c2, f2, fi2
                    fallback_sid = alt
                    logger.info("reachability %s: 主 sid 无暴露记录,"
                                "同义兄弟兜底取数 sid=%s", key, alt)
                    break
        if not crows.empty and fub is not None:
            # 多干预中位数(§3.1);Cmax_free 在 r_value 内乘 fub
            cmax_total = float(crows["cmax_molar_M"].median())
            cmax_meta = {
                "cmax_molar_M": cmax_total,
                "cmax_interventions": ";".join(sorted(set(crows["intervention"]))),
                "n_cmax_obs": int(crows["n_obs"].sum()),
                "fub_value": fub, "fub_imputed": fub_imputed,
                "synonym_fallback_sid": fallback_sid,
            }
        else:
            unknown_reason = "no_exposure"
    else:
        unknown_reason = "no_exposure"

    b1 = _b1_activity(drug_targets)
    act_by_gene = {r["gene"]: r for r in b1.to_dict("records")}
    genes = sorted({str(g) for g in (drug_targets["approvedSymbol"].dropna()
                                     if drug_targets is not None
                                     and "approvedSymbol" in getattr(drug_targets, "columns", [])
                                     else [])})
    weights = {"reachable": 1.0, "partial": 0.6, "unknown": 0.5, "unreachable": 0.2}
    from WY import config

    if hasattr(config, "REACH_TIER_WEIGHTS"):
        weights = dict(getattr(config, "REACH_TIER_WEIGHTS"))
    else:
        _warn_once("cfg:REACH_TIER_WEIGHTS",
                   "config.REACH_TIER_WEIGHTS 未合入(提案 B),按方案 §12.2-B3 "
                   "定稿值 1.0/0.6/0.5/0.2 运行")
    rows = []
    for gene in genes:
        a = act_by_gene.get(gene)
        r = r_value(cmax_total, fub, a["activity_M"] if a else None)
        if r is not None:
            tier = _tier(r, reachable_min=reachable_min, partial_min=partial_min)
        else:
            tier = "unknown"  # 活性未知独立档,永不并入不可达(§1.2-5)
        rows.append({
            "gene": gene,
            "r_value": r,
            "tier": tier,
            "activity_M": a["activity_M"] if a else None,
            "activity_class": a["activity_class"] if a else None,
            "n_assays": a["n_assays"] if a else 0,
            "activity_alt_M": a["activity_alt_M"] if a else None,
            "fallback_flag": a["fallback_flag"] if a else None,
            "unknown_reason": None if r is not None else "no_activity",
        })
    detail = pd.DataFrame(rows, columns=[
        "gene", "r_value", "tier", "activity_M", "activity_class", "n_assays",
        "activity_alt_M", "fallback_flag", "unknown_reason"])
    # 回退标记用可空布尔(None/True/False 混存 object 会触发 pandas 降级警告,
    # 且下游 fillna 语义含糊)——nullable BooleanDtype 一次定清
    detail["fallback_flag"] = detail["fallback_flag"].astype("boolean")
    # 查看序:R 降序 + 基因名字典序 tie-break(可复现红线 §1.2-8;
    # 注意这是明细查看序,E1/E6 联合排序契约在第 8 步另实现)
    detail = detail.sort_values(["r_value", "gene"], ascending=[False, True],
                                na_position="last").reset_index(drop=True)
    computed = detail[detail["r_value"].notna()]
    if not computed.empty:
        drug_tier = max(computed["tier"], key=lambda t: weights.get(t, 0.0))
        top_r = computed.loc[computed["tier"] == drug_tier, "r_value"].max()
        return {
            "schema": REACHABILITY_SCHEMA,
            "drug": drug_name, "substance_sid": sid,
            "tier": drug_tier, "r_primary": float(top_r),
            "unknown_reason": None,
            "exposure": cmax_meta,
            "n_targets": int(len(detail)), "n_computed": int(len(computed)),
            "targets": detail,
        }
    return {
        "schema": REACHABILITY_SCHEMA,
        "drug": drug_name, "substance_sid": sid,
        "tier": "unknown", "r_primary": None,
        "unknown_reason": unknown_reason or "no_activity",
        "exposure": cmax_meta,
        "n_targets": int(len(detail)), "n_computed": 0,
        "targets": detail,
    }


# ---- w(R) 加权排序与 E1/E6 联合排序(第 8 步,工程外圈,排序键方向=强制审查点) ----

_PLAN_WEIGHTS = {"reachable": 1.0, "partial": 0.6, "unknown": 0.5, "unreachable": 0.2}


def reach_weights() -> dict:
    """w(R) 分档常数权重(§12.2-B3 定稿初值,可调参数总量≤6 之四;1.7 基准
    网格微调后须负责人复核签字才改——本函数不提供手调入口)。"""
    from WY import config

    if hasattr(config, "REACH_TIER_WEIGHTS"):
        w = dict(getattr(config, "REACH_TIER_WEIGHTS"))
    else:
        _warn_once("cfg:REACH_TIER_WEIGHTS",
                   "config.REACH_TIER_WEIGHTS 未合入(提案 B),按方案 §12.2-B3 "
                   "定稿值 1.0/0.6/0.5/0.2 运行")
        w = dict(_PLAN_WEIGHTS)
    missing = [k for k in REACH_TIERS if k not in w]
    if missing:
        raise RuntimeError(f"REACH_TIER_WEIGHTS 缺分档键: {missing}"
                           f"(契约 {REACHABILITY_SCHEMA})")
    return w


def reach_weight(tier: str, weights: dict | None = None) -> float:
    """单档权重;非法枚举与残缺权重表都显式报错(静默回 0 会把打错的档
    伪装成最低优先)。"""
    if tier not in REACH_TIERS:
        raise ValueError(f"非法分档枚举 {tier!r}(契约 {REACHABILITY_SCHEMA})")
    w = weights if weights is not None else reach_weights()
    missing = [k for k in REACH_TIERS if k not in w]
    if missing:
        raise RuntimeError(f"权重表缺分档键: {missing}(契约 {REACHABILITY_SCHEMA})")
    return float(w[tier])


def exposure_penalty(w_a: float, w_b: float) -> float:
    """E6 可达性惩罚:两药权重取 min(瓶颈逻辑,合同 2026-11-22 随提案 A 冻结;
    与 exposure_confident 单调一致——两药均≥部分可达才不被惩罚封顶)。"""
    return float(min(w_a, w_b))


def apply_e1_sort(df, *, s_ab: str = "s_ab", z_ab: str = "z_ab",
                  pair_cols: tuple = ("drug_a", "drug_b"),
                  top_k: int | None = None):
    """E1 分层排序(§12.1-E1 定稿):重叠型(s_AB<0)按联合邻近度 z_AB 降序、
    互补型(s_AB>0)按分离强度 s_AB 降序,各自取 Top-K,类型为标注列,不做加权
    合成(分层依据引 Cheng 2019)。

    键方向为全项目强制审查点(清单 §1.3:取反即静默错)——两处降序逐字照
    定稿实现,后续修改必须同步 §1.3 纪要。s_AB==0(tie)归互补侧处理,为
    定稿未覆盖的微决策,已在此钉注释待负责人复核。
    tie-break:药对名(药 A,药 B)字典序,可复现红线(§1.2-8)。"""
    import pandas as pd

    out = df.copy()
    out["geometry_label"] = "separated"
    out.loc[out[s_ab] < 0, "geometry_label"] = "overlap"
    ov = out[out[s_ab] < 0].sort_values(
        [z_ab, *pair_cols], ascending=[False] + [True] * len(pair_cols))
    sep = out[out[s_ab] >= 0].sort_values(
        [s_ab, *pair_cols], ascending=[False] + [True] * len(pair_cols))
    if top_k is not None:
        ov, sep = ov.head(int(top_k)), sep.head(int(top_k))
    return pd.concat([ov, sep]).reset_index(drop=True)


def apply_e6_sort(df, *, score_col: str, tier_a_col: str, tier_b_col: str,
                  pair_cols: tuple = ("drug_a", "drug_b"),
                  weights: dict | None = None):
    """E6 A×B 联合排序(§12.1-E6 定稿):打 exposure_confident 标记
    (两药均≥部分可达)与 reach_penalty(=min(w_a,w_b),合同冻结口径),
    整表按 协同置信×可达性惩罚 降序,tie-break 药名字典序。
    score_col=协同置信数值列(口径归集成层定义,本函数只消费);分层报告
    (重叠/互补)由调用方先按 apply_e1_sort 分块。非法分档枚举直接报错。"""
    w = weights or reach_weights()
    out = df.copy()
    out["reach_weight_a"] = [reach_weight(t, w) for t in out[tier_a_col]]
    out["reach_weight_b"] = [reach_weight(t, w) for t in out[tier_b_col]]
    out["exposure_confident"] = [exposure_confident(a, b) for a, b in
                                 zip(out[tier_a_col], out[tier_b_col])]
    out["reach_penalty"] = [exposure_penalty(a, b) for a, b in
                            zip(out["reach_weight_a"], out["reach_weight_b"])]
    out["_e6_score"] = out[score_col] * out["reach_penalty"]
    out = out.sort_values(["_e6_score", *pair_cols],
                          ascending=[False] + [True] * len(pair_cols))
    return out.drop(columns="_e6_score").reset_index(drop=True)


def exposure_compare(baseline, reach_result) -> "pd.DataFrame":
    """"暴露约束前 vs 后"靶点优先级对比表(§3.2 招牌图数据侧):
    baseline = 一期排序(邻近度/RWR 产出的有序基因列表);加权口径 = 按 w(R)
    降序分桶,桶内保持一期相对次序(baseline_rank 升序),最终 tie-break
    基因名字典序(可复现红线)。分数版(网络分×w 重排)待集成层对齐邻近度
    分数列后实现——当前只承诺秩口径,不臆造分数。delta>0 = 加权后前移。"""
    import pandas as pd

    w = reach_weights()
    detail = (reach_result or {}).get("targets")
    tmap = {r["gene"]: (r["tier"], r["r_value"])
            for r in (detail.to_dict("records") if detail is not None else [])}
    recs = []
    for i, gene in enumerate(baseline, start=1):
        tier, rv = tmap.get(str(gene), ("unknown", None))
        recs.append({"gene": str(gene), "baseline_rank": int(i), "tier": tier,
                     "r_value": rv, "w": float(w[tier])})
    df = pd.DataFrame(recs, columns=["gene", "baseline_rank", "tier",
                                     "r_value", "w"])
    ordered = df.sort_values(["w", "baseline_rank", "gene"],
                             ascending=[False, True, True]).reset_index(drop=True)
    ordered["new_rank"] = ordered.index + 1
    out = ordered.sort_values("baseline_rank").reset_index(drop=True)
    out["delta"] = out["baseline_rank"] - out["new_rank"]
    return out[["gene", "baseline_rank", "new_rank", "delta", "tier",
                "w", "r_value"]]
