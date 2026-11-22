"""靶点浓度可达性:归一化 + R 值 + 加权排序(二期方案 §3.1 第 2/3 层)。

契约 reachability.v1(2026-11-22 冻结,接口合同随 E6 消费端先行定稿;
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
REACHABILITY_SCHEMA = "reachability.v1"

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
    logger.info("substance_map 归一化完成: %d 物质 → %s", len(df), path)
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


def _substance_lookup(nodes: list[dict]) -> dict:
    """物质名/同义词 → (sid, mol_weight_g_mol) 查找表(cmax 行的 substance 是
    名义字符串,如 'atorvastatin acid',需经 substance_map 桥到 sid 与质量)。
    全部小写化;同条映射冲突时保留先到并计数(不静默覆盖)。"""
    lookup: dict[str, tuple[str, float | None]] = {}
    conflicts = 0
    for n in nodes:
        sid = str(n.get("sid", "") or "")
        if not sid:
            continue
        mass = (n.get("substance") or {}).get("mass")
        try:
            mw = float(mass) if mass is not None else None
        except (TypeError, ValueError):
            mw = None
        names = [str(n.get("name", "") or "")] + [str(s) for s in (n.get("synonyms") or [])]
        for name in names:
            key = name.strip().lower()
            if not key:
                continue
            if key in lookup and lookup[key][0] != sid:
                conflicts += 1
                continue
            lookup[key] = (sid, mw)
    if conflicts:
        logger.warning("substance 名称映射冲突 %d 个(保留先到)", conflicts)
    return lookup


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
    default_fub = float(getattr(config, "FUB_DEFAULT", 1.0))
    if default_fub != 1.0:
        raise RuntimeError("FUB_DEFAULT 与方案 §3.1 定稿的保守上界 1.0 不一致,拒绝构建")
    if not hasattr(config, "FUB_DEFAULT"):
        logger.warning("config.FUB_DEFAULT 未合入(提案 B),按方案 §3.1 定稿值 1.0 运行")
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


def compute_reachability(drug_name: str, drug_targets) -> dict:
    """单药可达性(第 3 层计算,只读归一化 Parquet):B1 活性聚合(中位数 +
    细胞优先,回退留标记)→ R = Cmax_free ÷ 活性(摩尔)→ 三档分箱 +
    unknown 显式标注 → 逐靶点明细与药级档。
    返回契约 reachability.v1 单药块(供 synergy E6 挂载)与逐靶点明细表。
    第 7-8 步填充(三开口定案前 Ki/Kd 优先规则不动工)。"""
    raise NotImplementedError(
        "B1 聚合 + R 值 + 分箱在第 7 步实现,加权排序与对比表在第 8 步实现")
