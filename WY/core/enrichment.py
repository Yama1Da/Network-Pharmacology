"""GO/KEGG 功能富集分析：主用 g:Profiler，连不上时自动改用 Enrichr，并绘制结果图。"""

import hashlib
import json
import logging
import os
import re
import socket

import pandas as pd

from .caching import (clean_expired_cache_files, is_cache_fresh, read_cache_key,
                      save_csv_cache, save_json_file)

# 两个数据源的基因集库名对应关系（GO 生物过程/分子功能/细胞组分、KEGG 通路）。
ENRICHR_LIB_MAP = {
    "GO:BP": "GO_Biological_Process_2021",
    "GO:MF": "GO_Molecular_Function_2021",
    "GO:CC": "GO_Cellular_Component_2021",
    "KEGG": "KEGG_2021_Human",
}

# g:Profiler 物种名 → Enrichr 物种名。Enrichr 只支持人/小鼠等，大鼠没有对应库，
# 因此大鼠场景下 Enrichr 降级会被跳过，只能用 g:Profiler。
_ENRICHR_ORGANISM_MAP = {
    "hsapiens": "human",
    "mmusculus": "mouse",
}

# fetch_enrichment_* 返回状态中代表"调用成功"的取值;"circuit_open"/"probe_failed"/"request_failed"
# 等服务侧失败均视为不可信。module 富集与主富集的"真无显著项 vs 服务失败"判定共用这一份定义。
TRUSTED_REASONS = ("success", "no_significant", "cache")


def _enrichr_libs(sources) -> list[str]:
    """把 g:Profiler 数据源名映射为对应的 Enrichr 库名，过滤掉 Enrichr 没有的数据源。"""
    return [ENRICHR_LIB_MAP[s] for s in sources if s in ENRICHR_LIB_MAP]


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

    import config

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


def migrate_legacy_cache_names() -> int:
    """把旧版带哈希命名的富集缓存改为"数据源_标识"同名覆盖式命名,返回处理(改名/去重)的文件数。
    同一标识多变体时保留最新一份、删除其余,同目标的新格式文件同样参与"最新者胜";幂等,重复调用无副作用。
    """
    import config

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


# ---------- 富集数据源熔断（g:Profiler 主源 / Enrichr 兜底源共用状态机） ----------
# 数据源暂时连不上时自动改用另一个，避免反复苦等超时：
# - g:Profiler（主源）：单次调用重试耗尽（连续 ENRICHMENT_MAX_RETRIES 次异常）即"跳闸"；
# - Enrichr（兜底源）：每个请求只尝试 1 次，失败即"跳闸"——兜底失败时本次分析已注定告警，
#   重试上移到熔断层（冷却后的探测即重试），不在请求里白耗；
# - 冷却时长按阶梯逐次加长（各自的 *_CIRCUIT_COOLDOWN_LADDER）：每次"恢复失败"（冷却到期后的
#   探测失败）推进一档，最后一档封顶；任何有效响应（含"无显著条目"的空结果）即完全复位；
# - 冷却到期后的下一次请求作为"探测"：单次尝试 + 短超时（*_PROBE_TIMEOUT_SECONDS），
#   快速判定服务是否恢复，避免触发探测的分析苦等整个默认超时。
# 状态存于各自标记文件（失败时刻 + 阶梯档位），跨页面刷新与进程重启仍生效。


class _FileCircuitBreaker:
    """基于标记文件的三态熔断器：open（冷却期内,请求跳过）→ probe（冷却到期,下一请求试探）
    → closed（探测到有效响应,复位）。两个数据源只是标记文件与阶梯/探测参数不同。"""

    def __init__(self, marker_name: str, ladder_key: str, probe_timeout_key: str):
        self._marker_name = marker_name
        self._ladder_key = ladder_key
        self._probe_timeout_key = probe_timeout_key

    def _file(self) -> str:
        """熔断标记文件路径（记录失败时刻与冷却阶梯档位）。"""
        import config

        return os.path.join(config.ENRICHMENT_CACHE_DIR, self._marker_name)

    def read(self) -> dict | None:
        """读熔断标记，返回 {"failed_at", "stage"}；无标记或损坏返回 None（损坏按未熔断处理，
        不会因坏文件永久禁用数据源）。兼容旧版"纯时间戳"标记（仅 g:Profiler 有），视作第 1 档。"""
        path = self._file()
        if not os.path.exists(path):
            return None
        raw = ""
        try:
            with open(path, "r", encoding="utf-8") as f:
                raw = f.read().strip()
        except OSError:
            return None
        try:
            data = json.loads(raw)
            return {"failed_at": float(data["failed_at"]),
                    "stage": max(1, int(data.get("stage", 1)))}
        except Exception:
            try:  # 旧版格式：文件里只有失败时刻的 Unix 时间戳
                return {"failed_at": float(raw), "stage": 1}
            except Exception:
                return None

    def cooldown(self, stage: int) -> float:
        """第 stage 档对应的冷却时长（秒），超出阶梯按最后一档封顶。"""
        import config

        ladder = getattr(config, self._ladder_key)
        return ladder[min(stage, len(ladder)) - 1]

    def state(self) -> str:
        """当前熔断状态："open"（冷却期内,请求跳过）/ "probe"（冷却到期待探测）/ "closed"。"""
        import time

        mark = self.read()
        if mark is None:
            return "closed"
        if time.time() - mark["failed_at"] < self.cooldown(mark["stage"]):
            return "open"
        return "probe"

    def remaining_seconds(self) -> float:
        """冷却剩余秒数（probe 状态已到期待探测,返回接近 0）,
        供页面提示"约 X 分钟后自动重试"。"""
        import time

        mark = self.read()
        if mark is None:
            return 0.0
        return max(0.0, self.cooldown(mark["stage"]) - (time.time() - mark["failed_at"]))

    def trip(self, advance: bool) -> None:
        """记录失败时刻并开启熔断。advance=True 用于冷却到期后的探测失败：仅当上一条失败记录
        确实已过冷却期才推进一档（封顶最后一档），并发探测同时失败只有第一个推进；
        advance=False 用于常规请求失败：保持现有档位（无标记则为第 1 档）只刷新时间戳，
        避免同一次故障里并发线程的逐次失败把冷却阶梯一次推满。"""
        import time

        import config

        prev = self.read()
        stage = int((prev or {}).get("stage", 0))
        if advance:
            if prev and (time.time() - prev["failed_at"]) < self.cooldown(stage):
                stage = max(1, stage)  # 并发探测撞车：只刷新时间戳，不推进
            else:
                stage = min(stage + 1, len(getattr(config, self._ladder_key)))
        else:
            stage = max(1, stage)
        try:
            save_json_file(self._file(), {"failed_at": time.time(), "stage": stage})
        except Exception:
            logging.exception("Failed to write circuit breaker file: %s", self._marker_name)

    def trip_normal(self) -> None:
        """常规请求失败：保持现有档位（无标记则第 1 档），只刷新失败时刻。"""
        self.trip(advance=False)

    def trip_probe(self) -> None:
        """冷却到期后的探测失败：属于"恢复失败"，推进冷却阶梯。"""
        self.trip(advance=True)

    def close(self) -> None:
        """请求成功（任何有效响应）时清除熔断标记，完全复位。"""
        path = self._file()
        try:
            if os.path.exists(path):
                os.remove(path)
        except Exception:
            logging.exception("Failed to remove circuit breaker file: %s", self._marker_name)

    def skip_reason(self) -> str | None:
        """冷却期内返回跳过原因；冷却到期返回 None 放行（该次请求将按探测预算执行）。"""
        return "circuit_open" if self.state() == "open" else None

    def request_plan(self, normal_attempts: int | None = None) -> dict:
        """本次请求的执行预算：冷却到期（probe）按单次尝试 + 短超时探测，失败按"恢复失败"
        推进阶梯；closed 按常规尝试次数（normal_attempts 为 None 时取 ENRICHMENT_MAX_RETRIES）。
        open 不会走到（skip_reason 已拦截），防御性按常规预算返回。"""
        import config

        plan: dict = {"attempts": normal_attempts or config.ENRICHMENT_MAX_RETRIES}
        if self.state() == "probe":
            plan.update(attempts=1,
                        hard_deadline=getattr(config, self._probe_timeout_key),
                        on_fail=self.trip_probe)
        return plan


_gp_circuit = _FileCircuitBreaker(".gprofiler_unreachable",
                                  "GPROFILER_CIRCUIT_COOLDOWN_LADDER",
                                  "GPROFILER_PROBE_TIMEOUT_SECONDS")
_enrichr_circuit = _FileCircuitBreaker(".enrichr_unreachable",
                                       "ENRICHR_CIRCUIT_COOLDOWN_LADDER",
                                       "ENRICHR_PROBE_TIMEOUT_SECONDS")


def enrichment_retry_hint_minutes() -> int | None:
    """富集数据源熔断后距下次自动重试的大致分钟数（取两个数据源中较近的重试时机,
    向上取整,不足 1 分钟按 1 分钟）,供页面横幅提示"约 X 分钟后自动重试";
    两个数据源都不在熔断冷却期时返回 None。"""
    import math

    remaining = [c.remaining_seconds() for c in (_gp_circuit, _enrichr_circuit)
                 if c.state() != "closed"]
    if not remaining:
        return None
    return max(1, math.ceil(min(remaining) / 60))


def _run_with_deadline(fn, timeout: float):
    """在独立线程里限时执行一次调用，超时返回 None（由调用方按失败处理）。
    不能用 with 管理执行器——shutdown(wait=True) 会等孤儿线程跑完，限时就失去意义；
    孤儿线程最迟被进程级默认 socket 超时结束，不会累积。"""
    import concurrent.futures

    executor = concurrent.futures.ThreadPoolExecutor(max_workers=1,
                                                     thread_name_prefix="enrichcircuit")
    try:
        return executor.submit(fn).result(timeout=timeout)
    except concurrent.futures.TimeoutError:
        return None
    finally:
        executor.shutdown(wait=False)


_socket_timeout_set = False


def _ensure_default_socket_timeout() -> None:
    """给不接收超时参数的第三方富集库(g:Profiler/gseapy)兜底进程级默认 socket 超时。
    只设置一次、之后不再复位:socket.setdefaulttimeout 是进程级全局状态,多线程里
    逐请求 set/reset 会互相覆盖(线程 A 刚设好超时进入请求,线程 B 的 finally 复位
    None,A 之后建立的连接将无超时);一次性设为 API_TIMEOUT 对显式传超时的
    requests 链路无影响,只是让不带超时的库调用不再可能无限挂起。"""
    global _socket_timeout_set
    if not _socket_timeout_set:
        import config

        socket.setdefaulttimeout(config.API_TIMEOUT)
        _socket_timeout_set = True


def _call_with_request_guard(fn, label: str, on_fail=None,
                             attempts: int | None = None,
                             hard_deadline: float | None = None):
    """给第三方富集查询库统一补上请求保护：请求间隔、默认超时、失败重试;全部尝试失败返回 None。
    on_fail 在每次失败后被调用（g:Profiler 用它累计熔断）。attempts/hard_deadline 供熔断探测用：
    限定尝试次数并在独立线程限时执行（超时按失败处理），让冷却到期后的探测快速给出结论。"""
    import time

    import config

    _ensure_default_socket_timeout()
    if attempts is None:
        attempts = config.ENRICHMENT_MAX_RETRIES
    res = None
    for attempt in range(attempts):
        time.sleep(config.API_REQUEST_DELAY)  # 每次请求前稍等，避免请求过快被服务器拒绝
        try:
            if hard_deadline is None:
                res = fn()
            else:
                res = _run_with_deadline(fn, hard_deadline)
                if res is None:
                    raise TimeoutError(f"no result within {hard_deadline:g}s")
            break
        except Exception:
            logging.exception(
                "%s enrichment failed (attempt %d/%d)",
                label, attempt + 1, attempts,
            )
            res = None
            if on_fail is not None:
                on_fail()
            if attempt < attempts - 1:
                time.sleep(config.API_RETRY_DELAY)
    return res


def _set_status(status: dict | None, ok: bool, reason: str) -> None:
    """把本次调用状态写入 status(若调用方提供)，供页面区分"真无显著项"与"服务失败"。"""
    if status is not None:
        status.update(ok=ok, reason=reason)


def _profile_gprofiler(genes: list[str], organism: str, sources, pvalue: float):
    """单次 g:Profiler 查询（库较重,延迟导入）。"""
    from gprofiler import GProfiler

    return GProfiler(return_dataframe=True).profile(
        query=genes,
        organism=organism,
        sources=list(sources),
        user_threshold=pvalue,
        no_evidences=False,  # 需要返回每个通路命中的基因列表
    )


def _query_enrichr(genes: list[str], libs, organism: str, pvalue: float):
    """单次 Enrichr 查询（经 gseapy;不让库自行写盘,结果由本程序统一缓存）。"""
    import gseapy

    return gseapy.enrichr(
        gene_list=genes,
        gene_sets=list(libs),
        organism=organism,
        outdir=None,
        cutoff=pvalue,
    ).results


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


def _fetch_enrichment_cached(
    genes: list[str],
    organism: str,
    pvalue: float,
    label: str | None,
    status: dict | None,
    *,
    source_name: str,
    cache_key_extra,
    input_ok: bool,
    request,
    normalize,
    skip_if=None,
    on_fail=None,
    on_result_ok=None,
    request_plan=None,
) -> pd.DataFrame:
    """g:Profiler / Enrichr 两个富集数据源共用的"清洗输入→查缓存→带保护请求→整理→缓存"流程;
    网络失败时由请求保护返回 None。skip_if 返回非 None 时按该原因跳过请求(g:Profiler 的熔断),
    on_fail 在每次失败后被调用(用于累计熔断),on_result_ok 在拿到有效响应后调用(用于关闭熔断;
    空结果也算有效——服务可用性与结果多少是两回事,先复位再判空),request_plan 返回探测预算时
    按单次限时执行(g:Profiler 冷却到期后的探测)。"""
    import config

    if not input_ok:
        _set_status(status, False, "empty_input")
        return pd.DataFrame()
    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    _clean_cache_dir_limited()
    md5 = _enrichment_cache_key(genes, organism, cache_key_extra, pvalue)
    cache_file = os.path.join(config.ENRICHMENT_CACHE_DIR,
                              enrichment_cache_filename(source_name, label, md5))
    cached = _hit_csv_cache(cache_file, md5,
                            f"{source_name} enrichment for {len(genes)} genes")
    if cached is not None:
        _set_status(status, True, "cache")
        return cached

    skip = skip_if() if skip_if is not None else None
    if skip:
        logging.info("%s request skipped: %s", source_name, skip)
        _set_status(status, False, skip)
        return pd.DataFrame()

    attempts = hard_deadline = None
    if request_plan is not None:
        plan = request_plan() or {}
        attempts = plan.get("attempts")
        hard_deadline = plan.get("hard_deadline")
        if plan.get("on_fail") is not None:
            # 探测模式下的失败按"恢复失败"计（推进冷却阶梯），替换常规失败回调
            on_fail = plan["on_fail"]

    res = _call_with_request_guard(request, source_name, on_fail=on_fail,
                                   attempts=attempts, hard_deadline=hard_deadline)
    if res is None:
        # 重试(或探测)后仍失败（网络/服务不可用）
        reason = "probe_failed" if hard_deadline is not None else "request_failed"
        logging.warning("%s unreachable after %d attempt(s) for %d genes",
                        source_name, attempts or config.ENRICHMENT_MAX_RETRIES, len(genes))
        _set_status(status, False, reason)
        return pd.DataFrame()
    if on_result_ok is not None:
        # 任何有效响应都证明服务可用，先复位熔断再判定结果是否为空
        on_result_ok()
    if len(res) == 0:
        # 请求成功但无显著条目（真"空结果"，非失败）
        logging.info("%s: no significant terms (p < %s) for %d genes", source_name, pvalue, len(genes))
        _set_status(status, True, "no_significant")
        return pd.DataFrame()
    df = normalize(res)
    if df.empty:
        # 结果整理后为空时不写缓存，避免下次把异常空表当成有效结果
        logging.warning("%s result normalized to empty — skip caching (%d genes)", source_name, len(genes))
        _set_status(status, False, "normalize_empty")
        return df
    _set_status(status, True, "success")
    save_csv_cache(df, cache_file, header_lines=[md5])
    logging.info("%s enrichment for %d genes: %d terms cached", source_name, len(genes), len(df))
    return df


def fetch_enrichment_gprofiler(
    genes: list[str],
    organism: str = "hsapiens",
    sources: tuple = ("GO:BP", "GO:MF", "GO:CC", "KEGG"),
    pvalue: float = 0.05,
    label: str | None = None,
    status: dict | None = None,
) -> pd.DataFrame:
    """用 g:Profiler 查询基因列表的功能富集（GO/KEGG）,结果整理为统一列后缓存,相同输入直接读缓存。
    失败或无显著条目返回空表,status（若提供）区分"真无显著项"与"服务失败";熔断冷却期内跳过请求
    (reason="circuit_open"),冷却到期的首次请求按单次短超时探测,探测失败 reason="probe_failed";
    两种状态均由调用方降级走 Enrichr 兜底。
    """
    genes = sorted({g.strip() for g in genes if g and g.strip()})
    return _fetch_enrichment_cached(
        genes, organism, pvalue, label, status,
        source_name="gProfiler",
        cache_key_extra=sources,
        input_ok=bool(genes),
        skip_if=_gp_circuit.skip_reason,
        request=lambda: _profile_gprofiler(genes, organism, sources, pvalue),
        normalize=_normalize_enrichment_gprofiler,
        on_fail=_gp_circuit.trip_normal,
        on_result_ok=_gp_circuit.close,
        request_plan=_gp_circuit.request_plan,
    )


def _normalize_enrichment_gprofiler(df: pd.DataFrame) -> pd.DataFrame:
    """把 g:Profiler 返回的结果整理为统一列;其返回的 p 值已做多重检验校正,因此校正前后 p 值相同。
    同时兼容不同版本返回的列名差异。
    """
    df = df.reset_index()  # 通路 ID 有时在行索引里，先还原成列

    def _pick(*names) -> pd.Series:
        """按优先级取第一个存在的列，都没有时返回空值。"""
        for n in names:
            if n in df.columns:
                return df[n]
        return pd.Series(pd.NA, index=df.index, dtype=object)

    out = pd.DataFrame({
        "source": df["source"].astype(str),
        "term_id": _pick("term_id", "native"),
        "term_name": _pick("term_name", "name"),
        "p_value": pd.to_numeric(df["p_value"], errors="coerce"),
        "adjusted_p_value": pd.to_numeric(df["p_value"], errors="coerce"),
        "intersection_size": pd.to_numeric(df["intersection_size"], errors="coerce"),
        "query_size": pd.to_numeric(df["query_size"], errors="coerce"),
        "genes": _pick("intersection", "intersections").fillna(""),
    })
    out = out.dropna(subset=["term_id", "p_value"]).sort_values("p_value")
    out["term_id"] = out["term_id"].astype(str)
    out["term_name"] = out["term_name"].astype(str)
    return out.reset_index(drop=True)


def fetch_enrichment_enrichr(
    genes: list[str],
    libs: list[str],
    organism: str = "human",
    pvalue: float = 0.05,
    label: str | None = None,
    status: dict | None = None,
) -> pd.DataFrame:
    """用 Enrichr 查询功能富集,作为 g:Profiler 连不上时的备选方案;结果同样缓存,相同输入直接读缓存。
    兜底源自带熔断:每个请求只尝试 1 次,失败即跳闸,冷却期内跳过请求(reason="circuit_open"),
    冷却到期的首次请求按单次短超时探测,探测失败 reason="probe_failed"。
    请求失败或没有显著条目时返回空表,status 含义与 g:Profiler 版本一致。
    """
    genes = sorted({g.strip() for g in genes if g and g.strip()})
    n_query = len(genes)

    def _normalize_with_query_size(df: pd.DataFrame) -> pd.DataFrame:
        out = _normalize_enrichment_enrichr(df)
        # Enrichr 的 Overlap "a/b" 中 b 是通路基因集大小,不是查询集大小;
        # 统一列的 query_size 语义(g:Profiler=输入查询基因数)用真实查询数覆盖,
        # 否则兜底路径下页面表格的 query_size 会显示成几百的库大小,误导命中比例
        out["query_size"] = n_query
        return out

    return _fetch_enrichment_cached(
        genes, organism, pvalue, label, status,
        source_name="Enrichr",
        cache_key_extra=libs,
        input_ok=bool(genes and libs),
        skip_if=_enrichr_circuit.skip_reason,
        request=lambda: _query_enrichr(genes, libs, organism, pvalue),
        normalize=_normalize_with_query_size,
        on_fail=_enrichr_circuit.trip_normal,
        on_result_ok=_enrichr_circuit.close,
        request_plan=lambda: _enrichr_circuit.request_plan(normal_attempts=1),
    )


def _normalize_enrichment_enrichr(df: pd.DataFrame) -> pd.DataFrame:
    """把 Enrichr 返回的结果整理为统一列，库名换回与 g:Profiler 一致的写法。"""
    import re

    def _parse_overlap(s) -> tuple[int, int]:
        try:
            a, b = str(s).split("/")
            return int(a), int(b)
        except (ValueError, TypeError):
            return 0, 0

    def _extract_term_id(term) -> str:
        # 通路名形如 "XXX (GO:0018212)"，提取末尾括号内的通路 ID
        m = re.search(r"\(([^()]+)\)\s*$", str(term))
        return m.group(1) if m else str(term)

    def _strip_suffix(term) -> str:
        return re.sub(r"\s*\([^()]*\)\s*$", "", str(term)).strip() or str(term)

    parsed = df["Overlap"].map(_parse_overlap)
    source_rev = {v: k for k, v in ENRICHR_LIB_MAP.items()}  # 库名 → g:Profiler source 命名
    out = pd.DataFrame({
        "source": df["Gene_set"].map(source_rev).fillna(df["Gene_set"]),
        "term_id": df["Term"].map(_extract_term_id),
        "term_name": df["Term"].map(_strip_suffix),
        "p_value": pd.to_numeric(df["P-value"], errors="coerce"),
        "adjusted_p_value": pd.to_numeric(df["Adjusted P-value"], errors="coerce"),
        "intersection_size": parsed.map(lambda t: t[0]),
        "query_size": parsed.map(lambda t: t[1]),
        "genes": df["Genes"].fillna(""),
    })
    out = out.dropna(subset=["term_id", "p_value"]).sort_values("p_value")
    return out.reset_index(drop=True)


def run_enrichment(
    genes: list[str],
    organism: str = "hsapiens",
    pvalue: float = 0.05,
    sources: tuple | None = None,
    label: str | None = None,
    status: dict | None = None,
) -> tuple[pd.DataFrame, str | None]:
    """富集分析主入口：输入基因列表（通常是药物与疾病的交集靶点），先请求 g:Profiler，失败自动改用 Enrichr。
    两条链路都拿不到结果时返回空表;status（若提供）区分"真无显著项"（reason="no_significant"）
    与"服务不可用"（reason="both_sources_down"/"gprofiler_down_no_fallback"），供页面显著告警。"""
    import config

    sources = sources or config.ENRICHMENT_SOURCES

    st_gp: dict = {}
    df = fetch_enrichment_gprofiler(genes, organism=organism, sources=sources,
                                    pvalue=pvalue, label=label, status=st_gp)
    if not df.empty:
        return df, "gProfiler"

    # 备：Enrichr（不支持大鼠，跳过）
    enrichr_organism = _ENRICHR_ORGANISM_MAP.get(organism)
    st_enr: dict = {}
    if enrichr_organism:
        libs = _enrichr_libs(sources)
        df = fetch_enrichment_enrichr(genes, libs=libs, organism=enrichr_organism,
                                      pvalue=pvalue, label=label, status=st_enr)
        if not df.empty:
            return df, "Enrichr"

    # 区分"真无显著项"与"服务不可用"：任一已启用链路给出可信空结果即算 no_significant;
    # 否则按服务失败标注,页面据此显著告警而非当作普通空结果
    if st_gp.get("reason") in TRUSTED_REASONS or st_enr.get("reason") in TRUSTED_REASONS:
        _set_status(status, True, "no_significant")
    elif enrichr_organism:
        _set_status(status, False, "both_sources_down")
    else:
        _set_status(status, False, "gprofiler_down_no_fallback")
    return pd.DataFrame(), None


def render_enrichment_plots(df: pd.DataFrame) -> tuple:
    """绘制富集结果图：柱状图（每类数据源各取前 N 条,横轴为 -log10 校正 p 值）和气泡图（全局前 N 条,
    点越大命中基因越多）;没有有效结果时返回 (None, None)。
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import config

    df = df[df["adjusted_p_value"].notna() & (df["adjusted_p_value"] > 0)].copy()
    if df.empty:
        return None, None
    df["neg_log10_p"] = -np.log10(df["adjusted_p_value"])

    source_colors = {
        "GO:BP": "#4C72B0", "GO:MF": "#DD8452",
        "GO:CC": "#55A868", "KEGG": "#C44E52",
    }

    top = (
        df.sort_values("adjusted_p_value")
        .groupby("source", as_index=False)
        .head(config.ENRICHMENT_TOP_PER_SOURCE)
    )
    fig_bar, ax_bar = plt.subplots(figsize=(9, max(5, 0.28 * len(top))))
    colors = [source_colors.get(s, "#7f7f7f") for s in top["source"]]
    ax_bar.barh(top["term_name"], top["neg_log10_p"], color=colors, height=0.6)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("-log10(adjusted p-value)")
    ax_bar.set_title("Functional Enrichment (top terms per source)")
    ax_bar.tick_params(axis="y", labelsize=7)
    fig_bar.tight_layout()

    bubble = df.sort_values("adjusted_p_value").head(config.PATHWAY_DISPLAY_NUM)
    sizes = 40 + 120 * (bubble["intersection_size"] / max(bubble["intersection_size"].max(), 1))
    fig_dot, ax_dot = plt.subplots(figsize=(10, max(6, 0.3 * len(bubble))))
    sc = ax_dot.scatter(
        bubble["neg_log10_p"], bubble["term_name"],
        s=sizes, c=bubble["adjusted_p_value"], cmap="viridis_r",
        edgecolors="grey", alpha=0.85,
    )
    cb = fig_dot.colorbar(sc, ax=ax_dot)
    cb.set_label("adjusted p-value")
    ax_dot.set_xlabel("-log10(adjusted p-value)")
    ax_dot.set_title("Enrichment bubble plot")
    ax_dot.tick_params(axis="y", labelsize=7)
    fig_dot.tight_layout()

    return fig_bar, fig_dot


def _parse_hit_genes(raw) -> set[str]:
    """把富集结果里的命中基因列解析成大写基因集合;兼容 g:Profiler 列表样式字符串、Enrichr 分号/逗号分隔三种格式。"""
    import ast

    s = str(raw).strip()
    if not s or s.lower() == "nan":
        return set()
    if s.startswith("["):
        try:
            val = ast.literal_eval(s)
            if isinstance(val, (list, tuple, set)):
                return {str(g).strip().upper() for g in val if str(g).strip()}
        except (ValueError, SyntaxError):
            pass  # 不是合法的列表字面量,退回分隔符解析
    return {g.strip().strip("'\"").upper()
            for g in re.split(r"[;,]", s) if g.strip().strip("'\"")}


def render_drug_target_pathway_plot(
    drug_name: str,
    targets: list[str],
    enrichment_df: pd.DataFrame | None,
    top_n_pathways: int | None = None,
) -> "matplotlib.figure.Figure | None":
    """绘制"药物-靶点-通路"网络图（论文用）：红色三角为药物、黄色方块为靶点、蓝色 V 为 KEGG 通路。
    通路取 KEGG 来源、校正 p 值最小的前 N 条;只画命中至少一条所示通路的靶点,手动基因模式无药物名时省略药物节点。
    无 KEGG 通路或无命中靶点时返回 None。"""
    import textwrap

    import numpy as np
    import matplotlib.pyplot as plt
    import config

    if enrichment_df is None or enrichment_df.empty:
        return None
    top_n = top_n_pathways or config.PATHWAY_DISPLAY_NUM
    kegg = (enrichment_df[enrichment_df["source"] == "KEGG"]
            .sort_values("adjusted_p_value").head(top_n))
    if kegg.empty:
        return None

    # 每条通路的命中基因（见 _parse_hit_genes,兼容列表/分号/逗号三种格式），只保留输入靶点
    target_set = {str(tg).strip().upper() for tg in targets}
    pathway_hits: list[tuple[str, list[str]]] = []
    hit_targets: set[str] = set()
    for _, row in kegg.iterrows():
        hits = _parse_hit_genes(row.get("genes", "")) & target_set
        if hits:
            pathway_hits.append((str(row["term_name"]), sorted(hits)))
            hit_targets.update(hits)
    if not pathway_hits or not hit_targets:
        return None

    shown = sorted(hit_targets)
    n_t = len(shown)
    n_cols = max(6, min(12, int(np.ceil(np.sqrt(n_t * 1.5)))))
    n_rows = int(np.ceil(n_t / n_cols))
    grid_pos: dict[str, tuple[float, float]] = {}
    for i, gene in enumerate(shown):
        r, c = divmod(i, n_cols)
        row_len = min(n_cols, n_t - r * n_cols)
        x = c - (n_cols - 1) / 2 + (n_cols - row_len) / 2  # 末行不满时居中
        grid_pos[gene] = (x, (n_rows - 1) / 2 - r)

    k = len(pathway_hits)
    rx = n_cols / 2 + 2.6 + (0.8 if k > 24 else 0.0)  # 通路多时加大椭圆,给标签留出间距
    ry = max(n_rows / 2 + 1.4, 2.4) + (0.5 if k > 24 else 0.0)
    # 各节点整体偏移半个步长，避免正好落在正上/正下方与标题、药物挤在一起
    pathway_pos = {
        name: (rx * np.cos(ang), ry * np.sin(ang))
        for name, ang in (
            (name, np.pi / 2 + (i + 0.5) * 2 * np.pi / k)
            for i, (name, _) in enumerate(pathway_hits)
        )
    }
    drug_label = drug_name.strip().title() if drug_name and drug_name.strip() else ""
    drug_xy = (0.0, -(ry + 1.5))

    # ---- 画布与边（先画边再画节点，节点覆盖在线上层） ----
    fig_w = max(12.0, 1.05 * n_cols + 7.0)
    fig_h = max(8.0, 0.75 * n_rows + 4.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    if drug_label:
        for x, y in grid_pos.values():
            ax.plot([drug_xy[0], x], [drug_xy[1], y], color="#9a9a9a",
                    lw=0.7, alpha=0.45, zorder=1)
    for name, hits in pathway_hits:
        px, py = pathway_pos[name]
        for gene in hits:
            x, y = grid_pos[gene]
            ax.plot([px, x], [py, y], color="#9a9a9a", lw=0.6, alpha=0.25, zorder=1)

    ax.scatter(*zip(*grid_pos.values()), marker="s", s=560, c="#F3C63F",
               edgecolors="#8a7508", linewidths=0.6, zorder=3)
    gene_font = 6.0 if n_t > 100 else (6.5 if n_t > 60 else 7.0)
    for gene, (x, y) in grid_pos.items():
        ax.annotate(gene, (x, y), ha="center", va="center",
                    fontsize=gene_font, color="#222222", zorder=4)

    ax.scatter(*zip(*pathway_pos.values()), marker="v", s=300, c="#8FD4E8",
               edgecolors="#33708c", linewidths=0.6, zorder=3)
    for i, (name, (x, y)) in enumerate(pathway_pos.items()):
        # 标签沿椭圆法线方向外推;相邻节点按奇偶交错深浅两级偏移,
        # 斜向区与顶/底帽区的标签因此错层,不再互相压叠
        r = float(np.hypot(x, y))
        ux, uy = (x / r, y / r) if r else (0.0, 1.0)
        off = 0.85 if i % 2 else 0.3
        anchor = (x + ux * off, y + uy * off)
        if uy > 0.5:
            va = "bottom"
        elif uy < -0.5:
            va = "top"
        else:
            va = "center"
        ha = ("center" if abs(ux) < 0.35 else ("left" if ux > 0 else "right"))
        ax.annotate("\n".join(textwrap.wrap(name, width=24, max_lines=2, placeholder="…")),
                    anchor, ha=ha, va=va, fontsize=7, color="#1f4e5f", zorder=4)

    if drug_label:
        ax.scatter(*drug_xy, marker="^", s=1100, c="#E23B3B",
                   edgecolors="#8f1d1d", linewidths=0.8, zorder=3)
        ax.annotate(drug_label, drug_xy, xytext=(0, -16), textcoords="offset points",
                    ha="center", va="top", fontsize=11, fontweight="bold",
                    color="#8f1d1d", zorder=4)

    title = f"{drug_label} – target – pathway network" if drug_label \
        else "Target – pathway network"
    ax.set_title(f"{title} · top {k} KEGG pathways", fontsize=13)
    ax.set_xlim(-rx - 7.5, rx + 7.5)
    ax.set_ylim(-(ry + 3.0), ry + 2.8)
    ax.axis("off")
    fig.tight_layout()
    return fig


# ---------- 排序富集（GSEA prerank，仅智人） ----------
# 与普通富集的区别：输入是"基因|分数"排序列表（如药物亲和力、疾病关联分），
# 而非单纯的基因列表，用于检验高分靶点整体偏向哪些通路。
# 基因集来自 Enrichr 库的本地文件（首次使用时联网下载一次，之后不再重复下载）；
# 失败时返回空表，且失败结果不会写进缓存。


def _gsea_gmt_path(lib: str) -> str | None:
    """返回 Enrichr 库对应的本地基因集文件（GMT）路径；本地没有时联网下载一次，
    之后一直复用，下载失败返回 None。删除 GSEA_GMT_DIR 目录可强制重新下载。
    """
    import glob
    import config

    os.makedirs(config.GSEA_GMT_DIR, exist_ok=True)
    snaps = sorted(glob.glob(os.path.join(config.GSEA_GMT_DIR, f"{lib}_Human_*.gmt")))
    if snaps:
        return snaps[-1]

    import gseapy

    try:
        lib_dict = gseapy.get_library(lib, organism="Human")
    except Exception:
        logging.exception("Failed to download Enrichr library %s for GSEA", lib)
        return None
    if not lib_dict:
        logging.warning("Enrichr library %s downloaded empty", lib)
        return None
    from datetime import date

    path = os.path.join(config.GSEA_GMT_DIR, f"{lib}_Human_{date.today().isoformat()}.gmt")
    try:
        with open(path, "w", encoding="utf-8") as f:
            for term, gene_set in lib_dict.items():
                f.write("\t".join([str(term), "na", *[str(g) for g in gene_set]]) + "\n")
    except Exception:
        logging.exception("Failed to write GMT snapshot for %s", lib)
        return None
    logging.info("GSEA GMT snapshot saved: %s (%d gene sets)",
                 os.path.basename(path), len(lib_dict))
    return path


def gsea_gene_sets() -> dict[str, str]:
    """返回可用的基因集库（库名 → 文件路径）；单个库缺失自动跳过，全部缺失返回空字典。"""
    import config

    out: dict[str, str] = {}
    for lib in config.GSEA_ENRICHR_LIBS:
        path = _gsea_gmt_path(lib)
        if path:
            out[lib] = path
        else:
            logging.warning("GSEA gene set library unavailable, skipped: %s", lib)
    return out


def build_ranked_list(targets_df: pd.DataFrame, gene_col: str, score_col: str) -> tuple[pd.DataFrame, dict]:
    """把带打分的靶点表转为排序富集的输入（基因 | 分数，按分数降序）。
    基因名统一大写、同一基因取最高分;没有分数的行直接剔除不补零（补零会制造假信号）,剔除比例记入 stats。
    """
    if targets_df is None or targets_df.empty:
        return pd.DataFrame(columns=["gene", "score"]), {
            "n_input": 0, "n_kept": 0, "n_dropped": 0, "drop_ratio": 0.0,
            "missing_cols": [c for c in (gene_col, score_col)],
        }
    missing_cols = [c for c in (gene_col, score_col) if c not in targets_df.columns]
    if missing_cols:
        return pd.DataFrame(columns=["gene", "score"]), {
            "n_input": len(targets_df), "n_kept": 0, "n_dropped": len(targets_df),
            "drop_ratio": 1.0, "missing_cols": missing_cols,
        }

    df = targets_df[[gene_col, score_col]].rename(columns={gene_col: "gene", score_col: "score"})
    n_input = len(df)
    df["gene"] = df["gene"].astype(str).str.strip().str.upper()
    df = df[df["gene"] != ""]
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"])
    df = df.groupby("gene", as_index=False)["score"].max()
    df = df.sort_values("score", ascending=False).reset_index(drop=True)
    n_kept = len(df)
    stats = {
        "n_input": n_input,
        "n_kept": n_kept,
        "n_dropped": n_input - n_kept,
        "drop_ratio": (n_input - n_kept) / n_input if n_input else 0.0,
        "missing_cols": [],
    }
    logging.info("Ranked list built: %d inputs -> %d scored genes (drop_ratio=%.2f)",
                 n_input, n_kept, stats["drop_ratio"])
    return df, stats


def _normalize_enrichment_gsea(res2d: pd.DataFrame) -> pd.DataFrame:
    """把排序富集的原始结果整理为统一列；兼容不同版本软件的列名差异。"""
    if res2d is None or res2d.empty:
        return pd.DataFrame()
    cols = list(res2d.columns)

    def pick(*cands: str) -> str | None:
        for c in cands:  # 先精确（大小写不敏感）
            for col in cols:
                if col.lower() == c.lower():
                    return col
        for c in cands:  # 再子串
            for col in cols:
                if c.lower() in col.lower():
                    return col
        return None

    name_of = {
        "term_name": pick("Term"),
        "es": pick("ES"),
        "nes": pick("NES"),
        "nominal_p_value": pick("NOM p-val", "NOM p"),
        "fdr_q_value": pick("FDR q-val", "FDR q", "FDR"),
        "tag_pct": pick("Tag %"),
        "gene_pct": pick("Gene %"),
        "lead_genes": pick("Lead_genes", "Lead genes"),
    }
    if name_of["term_name"] is None or name_of["nes"] is None:
        logging.warning("GSEA res2d missing expected columns: %s", cols)
        return pd.DataFrame()

    out = pd.DataFrame({
        "term_name": res2d[name_of["term_name"]].astype(str),
        "es": pd.to_numeric(res2d[name_of["es"]], errors="coerce") if name_of["es"] else pd.NA,
        "nes": pd.to_numeric(res2d[name_of["nes"]], errors="coerce"),
        "nominal_p_value": (pd.to_numeric(res2d[name_of["nominal_p_value"]], errors="coerce")
                            if name_of["nominal_p_value"] else pd.NA),
        "fdr_q_value": (pd.to_numeric(res2d[name_of["fdr_q_value"]], errors="coerce")
                        if name_of["fdr_q_value"] else pd.NA),
        "tag_pct": res2d[name_of["tag_pct"]] if name_of["tag_pct"] else "",
        "gene_pct": res2d[name_of["gene_pct"]] if name_of["gene_pct"] else "",
        "lead_genes": res2d[name_of["lead_genes"]].fillna("") if name_of["lead_genes"] else "",
        "source": res2d["source"],
    })
    out = out.dropna(subset=["term_name", "nes"])
    # 统一按 |NES| 降序（页面再按 FDR 筛选展示）
    out = out.reindex(out["nes"].abs().sort_values(ascending=False).index)
    return out.reset_index(drop=True)


def _gsea_cache_md5(rnk: pd.DataFrame, libs: dict[str, str]) -> str:
    """GSEA 缓存键:由基因排序与分数、基因集库、计算参数共同决定,相同输入命中同一缓存。"""
    import config

    gene_part = "\r".join(rnk["gene"].astype(str))
    score_part = "\r".join(f"{v:.6g}" for v in rnk["score"].astype(float))
    libs_part = "|".join(f"{k}:{os.path.basename(v)}" for k, v in sorted(libs.items()))
    payload = (gene_part + "||" + score_part + "||" + libs_part
               + f"|{config.GSEA_SET_MIN_SIZE}|{config.GSEA_SET_MAX_SIZE}"
               + f"|{config.GSEA_PERMUTATION_NUM}|{config.GSEA_SEED}")
    return hashlib.md5(payload.encode("utf-8")).hexdigest()


def _gsea_cache_file(rnk: pd.DataFrame, libs: dict[str, str], label: str | None,
                     md5: str | None = None) -> str:
    """GSEA 结果缓存路径:同名覆盖式命名 GSEA_侧_标识.csv(键记录在文件头部,见 _hit_csv_cache)。
    label 约定"侧别_分析标识",侧别放在最前,避免与疾病名里可能含的"disease"字样混淆。
    """
    import config

    if md5 is None:
        md5 = _gsea_cache_md5(rnk, libs)
    return os.path.join(config.ENRICHMENT_CACHE_DIR,
                        enrichment_cache_filename("GSEA", label, md5))


def fetch_enrichment_gsea_prerank(
    rnk: pd.DataFrame,
    libs: dict[str, str] | None = None,
    label: str | None = None,
    status: dict | None = None,
) -> pd.DataFrame:
    """排序富集（GSEA prerank）主入口：在本地计算,逐个基因集库运行后合并,结果按 |NES| 降序。
    输入为 build_ranked_list 输出的"基因|分数"表;失败或无显著结果返回空表（失败不写缓存）,status 含义与 g:Profiler 版本一致。
    """
    import config

    if libs is None:
        libs = gsea_gene_sets()
    if rnk is None or len(rnk) == 0:
        _set_status(status, False, "empty_input")
        return pd.DataFrame()
    if len(rnk) < config.GSEA_MIN_RANKED_GENES:
        logging.info("GSEA skipped: only %d scored genes (< %d)",
                     len(rnk), config.GSEA_MIN_RANKED_GENES)
        _set_status(status, False, "insufficient_input")
        return pd.DataFrame()
    if not libs:
        logging.warning("GSEA skipped: no gene set library available")
        _set_status(status, False, "request_failed")
        return pd.DataFrame()

    os.makedirs(config.ENRICHMENT_CACHE_DIR, exist_ok=True)
    _clean_cache_dir_limited()
    gsea_md5 = _gsea_cache_md5(rnk, libs)
    cache_file = _gsea_cache_file(rnk, libs, label, md5=gsea_md5)
    cached = _hit_csv_cache(cache_file, gsea_md5, f"GSEA prerank result for {len(rnk)} genes")
    if cached is not None:
        _set_status(status, True, "cache")
        return cached

    import gseapy

    frames: list[pd.DataFrame] = []
    n_error = 0
    for lib, gmt_path in libs.items():
        try:
            res = gseapy.prerank(
                rnk=rnk, gene_sets=gmt_path, outdir=None, verbose=False,
                min_size=config.GSEA_SET_MIN_SIZE, max_size=config.GSEA_SET_MAX_SIZE,
                permutation_num=config.GSEA_PERMUTATION_NUM, seed=config.GSEA_SEED,
            )
        except LookupError:
            # 该库与本次基因列表交集过小，没有基因集通过 min/max 过滤：
            # 属数据特征而非服务失败，跳过该库，让其余库继续计算
            logging.info("GSEA prerank: no gene sets in %s matched %d genes, skipped",
                         lib, len(rnk))
            continue
        except Exception:
            logging.exception("GSEA prerank failed for %s (%d genes)", lib, len(rnk))
            n_error += 1
            continue
        r2d = res.res2d
        if r2d is not None and len(r2d) > 0:
            frames.append(r2d.assign(source=lib))

    if not frames:
        if n_error:
            _set_status(status, False, "request_failed")
            return pd.DataFrame()
        # 各基因集库都没有可用的基因集或没有显著富集的通路——真"无显著结果"，非失败
        logging.info("GSEA prerank: no gene sets passed for %d genes", len(rnk))
        _set_status(status, True, "no_significant")
        return pd.DataFrame()

    raw = pd.concat(frames, ignore_index=True)
    df = _normalize_enrichment_gsea(raw)
    if df.empty:
        # 结果整理后为空时不写缓存
        logging.warning("GSEA result normalized to empty — skip caching (%d genes)", len(rnk))
        _set_status(status, False, "normalize_empty")
        return df
    _set_status(status, True, "success")
    save_csv_cache(df, cache_file, header_lines=[gsea_md5])
    logging.info("GSEA prerank for %d genes: %d terms cached", len(rnk), len(df))
    return df


def render_gsea_plot(df: pd.DataFrame) -> "matplotlib.figure.Figure | None":
    """绘制 GSEA 结果条形图：按 |NES| 取前 N 条通路,NES > 0 红色（在高分靶点中富集/激活）、NES < 0 蓝色（抑制）。
    没有有效结果时返回 None。
    """
    import matplotlib.pyplot as plt
    import config

    df = df[df["nes"].notna()].copy()
    if df.empty:
        return None
    top = df.reindex(df["nes"].abs().sort_values(ascending=False).index)
    top = top.head(config.GSEA_TOP_DISPLAY)
    colors = ["#C44E52" if v > 0 else "#4C72B0" for v in top["nes"]]
    fig, ax = plt.subplots(figsize=(9, max(4, 0.3 * len(top))))
    ax.barh(top["term_name"], top["nes"], color=colors, height=0.6)
    ax.axvline(0, color="#888888", linewidth=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("NES (normalized enrichment score)")
    ax.set_title("GSEA prerank (red = enriched in high-score targets)")
    ax.tick_params(axis="y", labelsize=7)
    fig.tight_layout()
    return fig
