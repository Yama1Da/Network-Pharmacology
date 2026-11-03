"""富集分析的缓存与状态公共件。

缓存键生成、缓存文件命名(含旧版哈希命名的解析与迁移)、同名覆盖式缓存的命中判定、
缓存目录限频清理,以及 status 写入小工具。ORA(g:Profiler/Enrichr)与 GSEA 两条取数路径共用。
"""

import hashlib
import logging
import os
import re

import pandas as pd

from .caching import (clean_expired_cache_files, is_cache_fresh, read_cache_key)


def _enrichment_cache_key(genes: list[str], organism: str, sources, pvalue: float) -> str:
    """根据基因列表与查询参数生成富集缓存文件名（md5）：相同输入命中同一缓存文件。"""
    payload = (
        "\r".join(sorted(set(genes)))
        + f"|{organism}|{','.join(sources)}|{pvalue}"
    )
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


# 缓存文件名中哈希的保留长度:完整 md5 的前 12 位(48 位区分度,同一标识下几乎不会碰撞),
# 避免整串 32 位哈希把文件名拖得又长又难认。
_HASH_LEN = 12
# 文件名解析时兼容的哈希长度:新文件 12 位,旧文件完整 32 位
_HASH_RE = r"[0-9a-f]{8,32}"


def enrichment_cache_filename(source_name: str, label: str | None, md5: str) -> str:
    """构造富集缓存文件名:数据源_标识.csv。同名覆盖式缓存:同一标识只保留一份,是否命中由文件头部
    缓存键(save_csv_cache 的 header_lines 写入)校验,标识同但基因集/参数变了会重算后覆盖,不返回过期结果;
    label 为空(裸查询)退回短哈希命名避免互相覆盖。network_algorithms 推算模块缓存也用它,命名规则只有一份。"""
    if label:
        return f"{source_name}_{label}.csv"
    return f"{source_name}_{md5[:_HASH_LEN]}.csv"


# 富集缓存目录过期清理的上次执行时刻(monotonic):一次分析会多次进入取数路径,
# 限频避免反复全目录扫描;间隔取缓存有效期的 1/4(不少于 5 分钟)。
_cache_cleaned_at = -float("inf")


def _clean_cache_dir_limited() -> None:
    """清理富集缓存目录里超过 TTL 的文件。TTL 在命中判定里只挡复用、不删文件,
    不清理则目录长期只增不减;结果已永久保存在分析结果目录,删除无损失。"""
    global _cache_cleaned_at
    import time

    from WY import config

    now = time.monotonic()
    if now - _cache_cleaned_at < max(config.CACHE_TTL_SECONDS / 4, 300):
        return
    _cache_cleaned_at = now
    clean_expired_cache_files(config.ENRICHMENT_CACHE_DIR)


# ---------- 缓存文件名解析(页面友好展示与分子对接选择共用) ----------
# 新版命名不带哈希(gProfiler_标识.csv);旧版带 8~32 位哈希的文件名同样兼容解析。

_ORA_NAME_RE = re.compile(
    r"^(?P<source>gProfiler|Enrichr)_(?P<label>.+?)(?:_(?P<hash>[0-9a-f]{8,32}))?\.csv$")
_ORA_NAME_NO_LABEL_RE = re.compile(
    r"^(?P<source>gProfiler|Enrichr)_(?P<hash>[0-9a-f]{8,32})\.csv$")
_GSEA_NAME_RE = re.compile(r"^GSEA_(?P<rest>.+?)(?:_(?P<hash>[0-9a-f]{8,32}))?\.csv$")
_GSEA_NAME_NO_LABEL_RE = re.compile(r"^GSEA_(?P<hash>[0-9a-f]{8,32})\.csv$")


def parse_enrichment_filename(file_name: str) -> dict | None:
    """解析富集缓存文件名;不是富集缓存命名规则时返回 None,命中时含 kind/source/label/side 字段。
    新命名为 GSEA_侧_标识(无哈希);旧命名 GSEA_侧_标识_哈希 与 GSEA_标识_侧_哈希 也兼容解析。
    """
    if not isinstance(file_name, str):
        return None
    # 无标识的裸哈希命名先判(更具体),避免被通用规则当作带该哈希的 label
    m = _ORA_NAME_NO_LABEL_RE.match(file_name)
    if m:
        return {"kind": "ora", "source": m.group("source"), "label": "", "side": None}
    m = _ORA_NAME_RE.match(file_name)
    if m:
        return {"kind": "ora", "source": m.group("source"), "label": m.group("label"),
                "side": None}
    m = _GSEA_NAME_NO_LABEL_RE.match(file_name)
    if m:
        return {"kind": "gsea", "source": "GSEA", "label": "", "side": None}
    m = _GSEA_NAME_RE.match(file_name)
    if m:
        rest = m.group("rest")
        sm = re.match(r"^(?P<side>drug|disease)_(?P<prefix>.+)$", rest)      # 新命名:侧别在最前
        if not sm:
            sm = re.match(r"^(?P<prefix>.+)_(?P<side>drug|disease)$", rest)  # 旧命名:侧别在最后
        if sm:
            return {"kind": "gsea", "source": "GSEA",
                    "label": sm.group("prefix"), "side": sm.group("side")}
        return {"kind": "gsea", "source": "GSEA", "label": rest, "side": None}
    return None


def enrichment_filename_label(file_name: str) -> str:
    """从富集缓存文件名解析"药物_疾病"标识,用于命名后续结果;无法识别命名时返回空串。"""
    info = parse_enrichment_filename(file_name)
    return info["label"] if info else ""


def migrate_legacy_cache_names() -> int:
    """把旧版带哈希命名的富集缓存改为"数据源_标识"同名覆盖式命名,返回处理(改名/去重)的文件数。
    同一标识多变体时保留最新一份、删除其余,同目标的新格式文件同样参与"最新者胜";幂等,重复调用无副作用。
    """
    from WY import config

    try:
        names = sorted(os.listdir(config.ENRICHMENT_CACHE_DIR))
    except OSError:
        return 0
    # 目标名 -> [(修改时间, 文件名)];新格式文件本身也是"最新者胜"的候选之一
    groups: dict[str, list[tuple[float, str]]] = {}
    for name in names:
        info = parse_enrichment_filename(name)
        if info is None:
            continue
        m = re.search(rf"_({_HASH_RE})\.csv$", name)
        if m:
            label = info["label"]
            if info["kind"] == "gsea" and info["side"]:
                label = f"{info['side']}_{label}"  # GSEA 统一为侧别在最前
            target = enrichment_cache_filename(info["source"], label, m.group(1))
        else:
            target = name
        path = os.path.join(config.ENRICHMENT_CACHE_DIR, name)
        try:
            groups.setdefault(target, []).append((os.path.getmtime(path), name))
        except OSError:
            continue
    moved = 0
    for target, cands in groups.items():
        if len(cands) == 1 and cands[0][1] == target:
            continue
        cands.sort(reverse=True)
        newest = cands[0]
        for _, name in cands[1:]:
            try:
                os.remove(os.path.join(config.ENRICHMENT_CACHE_DIR, name))
                moved += 1
            except OSError:
                logging.exception("Failed to remove duplicate enrichment cache %s", name)
        if newest[1] != target:
            try:
                os.replace(os.path.join(config.ENRICHMENT_CACHE_DIR, newest[1]),
                           os.path.join(config.ENRICHMENT_CACHE_DIR, target))
                moved += 1
            except OSError:
                logging.exception("Failed to rename enrichment cache %s", newest[1])
    if moved:
        logging.info("Enrichment cache migrated to label-keyed naming: %d ops", moved)
    return moved


def _set_status(status: dict | None, ok: bool, reason: str) -> None:
    """把本次调用状态写入 status(若调用方提供)，供页面区分"真无显著项"与"服务失败"。"""
    if status is not None:
        status.update(ok=ok, reason=reason)


def _hit_csv_cache(cache_file: str, md5: str, label: str) -> pd.DataFrame | None:
    """同名覆盖式缓存的命中判定:文件存在、有效期内、且头部键与当前查询一致。
    键不一致(同一标识但基因集/参数变了)或旧格式无键文件按未命中处理,重算后覆盖,避免旧输入结果被误当新分析结果。
    """
    if not is_cache_fresh(cache_file):
        return None
    if read_cache_key(cache_file) != md5:
        logging.info("%s: cache key changed, recompute %s", label, os.path.basename(cache_file))
        return None
    logging.info("Loading cached %s", label)
    return pd.read_csv(cache_file, encoding="utf-8-sig", comment="#")
