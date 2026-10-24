"""缓存/结果文件的管理域逻辑:来源文件与整合文件的刷新编排、设置持久化、
缓存文件名的用户可读化。供药物/疾病管理页(经 ui.source_management)与
设置面板等界面调用;本模块不依赖任何 UI 框架,可离线使用。

页面通过 spec 字典描述各自的数据源差异:
    cache_dir        来源文件所在目录
    merged_dir       merged 整合文件所在目录
    merged_prefix    整合文件前缀("merged_drug_v1_" 等,带过滤规则版本号,取自 config)
    source_prefixes  来源文件前缀元组,如 ("ChEMBL_", "DGIdb_")
    fetch_source     {前缀: fn(name, cache_hours) -> (表, 是否成功)}:重新抓取一个来源文件
    fetch_merged     fn(name) -> 表:重建整合文件
    invalidate_merged fn(name):删除对应整合缓存
"""

import logging
import os
import threading
from collections.abc import Callable

import pandas as pd

from WY.i18n import t
from WY.core.enrichment import parse_enrichment_filename
from WY.core.http_client import clear_http_cache
from WY.core.text_utils import resolve_original_name, safe_name

# ui_settings.json 的进程内读写锁:齿轮面板与对接页都会高频保存设置
_ui_settings_lock = threading.Lock()


def save_ui_settings(overrides: dict) -> None:
    """把齿轮面板/对接页/日志页修改的分析参数持久化到 ui_settings.json。

    config.py 启动时会读取该文件套用覆盖值,因此设置重启后仍生效;
    只写白名单内的键(数值/布尔参数、ChEMBL 黑名单与日志级别),其余忽略。
    与文件中已有内容按键合并(不同入口各自改各自的参数,互不覆盖)。
    """
    import json

    from WY import config

    clean: dict = {}
    for k, v in overrides.items():
        if k == "CHEMBL_TARGET_BLACKLIST":
            if isinstance(v, (list, tuple)):
                clean[k] = sorted({str(x).strip().upper() for x in v if str(x).strip()})
        elif k == "STRING_VERSION":
            # 版本号统一校验归一(保小写,供 URL 派生),非法格式/非字符串一律不写盘;
            # 必须先于通用数值/布尔分支,避免数字被原样落盘
            sv = config.normalized_string_version(v)
            if sv:
                clean[k] = sv
        elif k in config.UI_SETTINGS_KEYS and isinstance(v, bool):
            clean[k] = v
        elif k in config.UI_SETTINGS_KEYS and isinstance(v, int):
            clean[k] = v  # 步进器参数按整数保存,避免 +- 调节后出现小数
        elif k in config.UI_SETTINGS_KEYS and isinstance(v, (int, float)) \
                and not isinstance(v, bool):
            clean[k] = float(v)
        elif k in config.UI_SETTINGS_KEYS and isinstance(v, str):
            clean[k] = v.strip().upper()  # 日志级别等枚举字符串,统一大写
    if not clean:
        return
    with _ui_settings_lock:  # 读-改-写整段互斥:不同入口各自改各自的参数,并发时不能互相覆盖
        merged: dict = {}
        if os.path.exists(config.UI_SETTINGS_FILE):
            try:
                with open(config.UI_SETTINGS_FILE, "r", encoding="utf-8") as f:
                    merged = json.load(f) or {}
            except Exception:
                merged = {}  # 旧文件损坏时按空文件处理,以本次写入为准
        merged.update(clean)
        os.makedirs(os.path.dirname(config.UI_SETTINGS_FILE), exist_ok=True)
        # 原子写(临时文件+替换),中断不会留下半截 JSON
        from WY.core.caching import save_json_file
        save_json_file(config.UI_SETTINGS_FILE, merged)


def purge_merged_cache() -> int:
    """清空 merged 整合缓存目录,返回删除的文件数。

    降噪阈值等参数变更后,旧参数算出的整合结果不再可信,统一清除,
    下次分析按新参数从原始缓存重建(原始缓存不受影响,无需重新联网)。
    """
    import glob

    from WY import config

    n = 0
    for path in glob.glob(os.path.join(config.MERGED_CACHE_DIR, "*.csv")):
        try:
            os.remove(path)
            n += 1
        except OSError:
            pass
    if n:
        logging.info("Purged %d merged cache files (settings changed)", n)
    return n


# ---------- 缓存文件名的友好命名展示 ----------
# 缓存文件名(数据源_标识_短哈希.csv)面向程序;下面两个函数把它翻译成用户可读的名称,
# 富集信息页与分子对接页的文件选择共用,保证两处看到的名字一致。


def humanize_analysis_label(label: str) -> str:
    """"药物_疾病"分析标识转成人读文本:首个下划线换成"→",其余下划线换回空格。"""
    label = (label or "").strip()
    if not label:
        return ""
    if "_" in label:
        drug, disease = label.split("_", 1)
        return f"{drug} → {disease.replace('_', ' ')}"
    return label.replace("_", " ")


def enrichment_display_name(file_name: str) -> str:
    """富集缓存文件名 → 用户友好名称,如 gProfiler_aspirin_cancer_<哈希>.csv
    → "aspirin → cancer · gProfiler";模块缓存的 `_moduleN` 后缀单独译成"模块 N",
    避免被埋进疾病名文本里难以辨认。无法识别的文件名原样返回。"""
    info = parse_enrichment_filename(file_name)
    if info is None:
        return file_name
    label = info["label"]
    module_desc = ""
    if "_module" in label:
        head, _, num = label.rpartition("_module")
        if num.isdigit():
            module_desc = t("enrich_file.module_n", n=int(num))
            label = head
    desc = humanize_analysis_label(label)
    if info["kind"] == "gsea":
        side = ""
        if info["side"] == "drug":
            side = t("enrich_file.side_drug")
        elif info["side"] == "disease":
            side = t("enrich_file.side_disease")
        return " · ".join(x for x in (desc, "GSEA", side) if x)
    return " · ".join(x for x in (desc, info["source"], module_desc) if x)


# ---------- 数据管理页公共逻辑(药物/疾病两个管理页共用) ----------


def _mgmt_file_path(filename: str, spec: dict) -> str:
    """整合前缀的文件位于 merged 目录,其余位于来源缓存目录。"""
    if filename.startswith(spec["merged_prefix"]):
        return os.path.join(spec["merged_dir"], filename)
    return os.path.join(spec["cache_dir"], filename)


def refresh_cache_file(filename: str, spec: dict, cache_hours: int | None = None,
                       clear_http: bool = True) -> dict:
    """按文件名重新抓取或整合数据,更新指定文件(药物/疾病管理页共用)。

    来源文件先备份旧文件,抓取成功才替换,失败则恢复旧文件,避免数据丢失;
    merged 整合文件先更新两个来源文件再重建,任一来源失败则保留旧的整合结果。
    fetch_source 返回 (表, 状态),状态为 "ok"(有数据)/"empty"(抓取成功但 0 靶点)/
    "failed"(请求失败);empty 不算失败,返回 Status=EMPTY,由页面以警告提示检查查询物。
    cache_hours 透传给 fetch_source(ChEMBL 更新需区分单文件/全部更新的 TTL)。
    clear_http=False 供并发批处理使用:由 update_all_cache_files 在入口统一清一次,
    避免批量中途清缓存干扰其他线程的在途请求。

    文件名里的名字是 safe_name 有损转换的结果(如 Alzheimer_s_Disease),
    先按映射还原成原始查询名(Alzheimer's Disease)再抓取,否则查不到数据。
    """
    # merged 整合文件:先更新两个来源文件,再重建
    if filename.startswith(spec["merged_prefix"]):
        stem = filename[len(spec["merged_prefix"]):-len(".csv")]
        name = resolve_original_name(stem) or stem
        for prefix in spec["source_prefixes"]:
            # 来源文件名用 safe 名(name 已还原为原始名,含空格/撇号等)
            source = f"{prefix}{safe_name(name)}.csv"
            result = refresh_cache_file(source, spec, cache_hours)
            if result["Status"] == "FAILED":  # 来源 0 靶点(EMPTY)不算失败,继续重建
                return {"Status": "FAILED",
                        "Message": f"source '{source}' update failed: {result['Message']}"}
        merged_path = _mgmt_file_path(filename, spec)
        if os.path.exists(merged_path):
            os.remove(merged_path)  # 可能已在上一步删除,这里兜底再删一次
        df = spec["fetch_merged"](name)
        # 两来源都为 0 靶点时整合结果也为空:不算失败,由页面以警告提示检查查询物
        return {"Status": "OK" if not df.empty else "EMPTY",
                "Message": f"{len(df)} targets"}

    # 来源文件:先备份,抓取成功才替换,失败则恢复旧文件
    path = _mgmt_file_path(filename, spec)
    backup = path + ".bak"
    if os.path.exists(path):
        os.replace(path, backup)

    # 清空 HTTP 响应缓存,避免强制更新时命中旧响应导致更新无效;
    # 并发批处理改为入口统一清一次,见 clear_http 参数说明
    if clear_http:
        clear_http_cache()

    status = "failed"
    df = pd.DataFrame()
    try:
        for prefix in spec["source_prefixes"]:
            if filename.startswith(prefix):
                stem = filename[len(prefix):-len(".csv")]
                name = resolve_original_name(stem) or stem
                df, status = spec["fetch_source"][prefix](name, cache_hours)
                break
        else:
            raise ValueError(f"unknown file prefix: {filename}")
    except Exception:
        logging.exception("refresh_file failed for '%s'", filename)
        status = "failed"

    if status in ("ok", "empty"):
        # 抓取方已写出新文件(ok 必有;empty 时仅部分来源写空缓存)则丢弃备份;
        # 没写(如 ChEMBL 查无此药不落盘)则恢复旧文件占位,避免条目从列表里消失
        if os.path.exists(path):
            if os.path.exists(backup):
                os.remove(backup)
        elif os.path.exists(backup):
            os.replace(backup, path)
        spec["invalidate_merged"](name)
        return {"Status": "OK" if status == "ok" else "EMPTY",
                "Message": f"{len(df)} targets"}
    # 失败:恢复旧文件;状态为 "failed:原因" 时把原因并入 Message 透传给页面,
    # 供界面区分配置问题(如未配置 DisGeNET 密钥)与一般请求失败
    if os.path.exists(backup):
        os.replace(backup, path)
    reason = status[len("failed:"):] if status.startswith("failed:") else ""
    return {"Status": "FAILED",
            "Message": f"fetch failed: {reason}" if reason else "fetch failed"}


def update_all_cache_files(spec: dict, cache_hours: int | None = None, *,
                           on_progress: Callable[[float, str], None] | None = None,
                           max_workers: int = 4,
                           filenames: list[str] | None = None):
    """更新某类别全部来源文件并重建全部整合文件,返回结果表(File/Status/Message)。

    filenames 提供时只更新这些来源文件,并只重建它们对应的整合文件(管理页"更新选中"
    的批量入口);选中清单里的 merged 整合文件按当前来源直接重建,不重抓来源。
    来源文件互相独立,线程池并发抓取提速;抓取链路自带按源限流(throttle 与
    DisGeNET 配额均有锁保护),并发安全。HTTP 响应缓存在入口只清一次:
    每文件中途清理会干扰其他线程的在途请求,且批量开始后本批抓取是唯一的
    缓存写入方,清一次即可保证全部抓到新数据。
    merged 整合依赖全部来源的最新结果,保持串行收尾。
    on_progress 为可选进度回调 on_progress(frac, text):frac∈[0,1],text 为已翻译的
    进度文案,由界面层决定呈现方式(如 st.progress,仅限页面线程传入);本模块不碰 UI。
    max_workers 为并发上限。
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    def _tick(frac: float, key: str, **fmt) -> None:
        if on_progress is not None:
            on_progress(frac, t(key, **fmt))

    def _short(name: str) -> str:
        """进度文案用的短文件名:去 .csv 后缀,超长截断加省略号。"""
        short = name.removesuffix(".csv")
        return short[:31] + "…" if len(short) > 32 else short

    # 第一步:先并发更新所有来源文件(filenames 提供时只更新选中的)
    jobs: list[str] = []
    names = set()
    # 全新安装、从未抓取过数据时缓存目录尚不存在:按无来源文件处理,直接跳到 merged 重建
    try:
        cache_entries = sorted(os.listdir(spec["cache_dir"]))
    except OSError:
        cache_entries = []
    for filename in cache_entries:
        if not filename.endswith(".csv"):
            continue
        if filenames is not None and filename not in filenames:
            continue
        for prefix in spec["source_prefixes"]:
            if filename.startswith(prefix):
                stem = filename[len(prefix):-len(".csv")]
                names.add(resolve_original_name(stem) or stem)
                jobs.append(filename)
                break

    rows: list[dict] = []
    if jobs:
        clear_http_cache()

        def _run_one(filename: str) -> dict:
            return refresh_cache_file(filename, spec, cache_hours, clear_http=False)

        results: dict[str, dict] = {}
        workers = max(1, min(max_workers, len(jobs)))
        if workers == 1:
            for done, filename in enumerate(jobs, 1):
                try:
                    results[filename] = _run_one(filename)
                except Exception:
                    results[filename] = {"Status": "FAILED", "Message": "update error"}
                _tick(0.8 * done / len(jobs), "mgmt.prog.sources",
                      done=done, total=len(jobs), latest=_short(filename))
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(_run_one, f): f for f in jobs}
                for done, fut in enumerate(as_completed(futures), 1):
                    filename = futures[fut]
                    try:
                        results[filename] = fut.result()
                    except Exception:
                        results[filename] = {"Status": "FAILED", "Message": "update error"}
                    _tick(0.8 * done / len(jobs), "mgmt.prog.sources",
                          done=done, total=len(jobs), latest=_short(filename))
        # 按文件名排序输出,与并发完成顺序无关
        rows = [{"File": f, **results[f]} for f in sorted(results)]

    # 第二步:再重建所有 merged 整合文件(name 已还原为原始名,文件名需再走 safe_name);
    # 选中清单里的整合文件也加入重建(按当前来源重建,不重抓来源)
    if filenames is not None:
        for filename in filenames:
            if filename.startswith(spec["merged_prefix"]):
                stem = filename[len(spec["merged_prefix"]):-len(".csv")]
                names.add(resolve_original_name(stem) or stem)
    merged_names = sorted(names)
    for idx, name in enumerate(merged_names, 1):
        merged_file = f"{spec['merged_prefix']}{safe_name(name)}.csv"
        merged_path = os.path.join(spec["merged_dir"], merged_file)
        if os.path.exists(merged_path):
            os.remove(merged_path)
        try:
            df = spec["fetch_merged"](name)
            if not df.empty:
                rows.append({"File": merged_file, "Status": "OK", "Message": f"{len(df)} targets"})
            else:
                # 两来源都 0 靶点:不算失败,由汇总警告提示检查查询物
                rows.append({"File": merged_file, "Status": "EMPTY", "Message": "0 targets"})
        except Exception:
            rows.append({"File": merged_file, "Status": "FAILED", "Message": "update error"})
        _tick(0.8 + 0.2 * idx / len(merged_names), "mgmt.prog.merged",
              done=idx, total=len(merged_names))
    _tick(1.0, "mgmt.prog.done")
    return pd.DataFrame(rows, columns=["File", "Status", "Message"])


def save_api_env(changes: dict) -> list[str]:
    """把齿轮面板提交的 API.env 参数写入 API.env 并回填运行时(config.apply_api_env)。

    只接受 config.API_ENV_KEYS 白名单内的键;值先消毒(去首尾空白、剥掉一层成对的
    首尾引号——从文件属性等处复制路径时常把引号一起带上来,留着会被 dotenv 读回
    成值的一部分;含换行/井号的值直接丢弃,堵 env 文件的注入面);空值一律跳过——
    密钥框不回显当前值,空输入天然等于"不修改"。返回实际生效的键名;写盘失败时
    异常原样抛出由调用方提示。
    落盘格式由 set_key 决定(单引号包裹、反斜杠转义为 \\\\):这是设计内形态,
    dotenv 读回时无损还原,与手写单反斜杠的旧行为等价。
    CLI 与 Web 进程共享 API.env:落盘对 CLI 在其下次启动时生效,本进程即时生效。
    日志只记键名不记值,密钥与路径都不应进入日志文件。"""
    from WY import config
    from dotenv import set_key

    clean: dict = {}
    for k, v in changes.items():
        if k not in config.API_ENV_KEYS:
            continue
        v = str(v).strip()
        if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
            v = v[1:-1].strip()
        if not v or "\n" in v or "\r" in v or "#" in v:
            continue
        clean[k] = v
    if not clean:
        return []
    env_path = os.path.join(config.BASE_DIR, "API.env")
    for k, v in clean.items():
        set_key(env_path, k, v)  # 只改目标行,保留注释与其他行;文件不存在则创建
    applied = config.apply_api_env(clean)
    logging.info("API.env updated: %s", ", ".join(applied))
    return applied
