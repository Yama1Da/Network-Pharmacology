"""ML screen 结果目录的管理数据层(方案 §10 增补:结果保留语义与结果管理页)。

历史 run 目录只增不减(每次运行一个 <tag>_<时间戳> 目录),本模块是治理能力的
唯一数据出口 —— 结果管理页面只做渲染与调用,不自带文件逻辑;前端后续改写 TS
时以本模块为 API 边界:

- list_ml_runs:         双根扫描(CLI 的 results/mlscreen 与页面的 GEOresult)
                        列出全部可识别 run,含 summary 概要、陈旧原因、体积;
- delete_ml_runs:       带双重校验的删除(必须位于已登记根下、必须是含
                        summary.json 的 run 目录),逐条 INFO 留痕;
- zip_ml_runs:          把一个/多个 run 目录打包为 zip 字节(管理页下载用);
- cleanup_keep_per_gse: 写侧自动清理 —— 同一 tag(GSE/矩阵名)只保留最近
                        ML_RESULTS_KEEP_PER_GSE 份;新落盘的 run 必在其中。
"""

import json
import logging
import os
import re
import shutil
import time

# run 目录名以 _YYYYMMDD_HHMMSS 结尾;tag 自身可含下划线(本地矩阵文件名)
_RUN_TS = re.compile(r"_(\d{8})_(\d{6})$")


def _roots() -> list[str]:
    """登记的结果根:页面产物(GEOresult,与 MDresult 平级)+ CLI 产物(results/mlscreen)。"""
    from WY import config
    from .ml_pipeline import default_results_root
    return [config.GEO_RESULT_DIR, os.path.join(default_results_root(), "mlscreen")]


def run_tag(dir_name: str) -> str | None:
    """从 run 目录名提取 tag(去掉尾部时间戳);不含时间戳的目录名返回 None。"""
    m = _RUN_TS.search(dir_name)
    return dir_name[:m.start()] if m else None


def _load_summary(run_dir: str) -> dict | None:
    sp = os.path.join(run_dir, "summary.json")
    if not os.path.isfile(sp):
        return None
    try:
        with open(sp, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def _dir_size(path: str) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            try:
                total += os.path.getsize(os.path.join(root, name))
            except OSError:
                pass
    return total


def list_ml_runs(roots: list[str] | None = None) -> list[dict]:
    """列出全部可识别的 run 目录(新→旧),供结果管理页与页内"最近 N 次"使用。

    每行:{name, tag, root, path, summary, size_bytes, ts, stale}。
    summary 为 None 表示 summary.json 损坏/缺失(陈旧列会标注不可读);
    stale 为 stale_reasons 的机器可读原因列表(空列表 = 未陈旧)。
    """
    from .ml_pipeline import stale_reasons  # 运行时导入:管线载入后无循环问题

    rows: list[dict] = []
    for root in (roots or _roots()):
        if not os.path.isdir(root):
            continue
        for name in sorted(os.listdir(root)):
            run_dir = os.path.join(root, name)
            if not os.path.isdir(run_dir):
                continue
            summary = _load_summary(run_dir)
            if summary is None and not os.path.isfile(os.path.join(run_dir, "summary.json")):
                continue  # 无 summary 的目录不是可识别 run(如手工残骸),不列不删
            m = _RUN_TS.search(name)
            ts = (m.group(1) + m.group(2)) if m else "0"
            try:
                stale = stale_reasons(summary) if summary else [{"type": "unreadable"}]
            except Exception:
                stale = [{"type": "unreadable"}]
            meta = summary.get("ml_meta", {}) if summary else {}
            rows.append({
                "name": name, "tag": run_tag(name), "root": root, "path": run_dir,
                "summary": summary, "ts": ts,
                "size_bytes": _dir_size(run_dir),
                "mtime": os.path.getmtime(run_dir),
                "stale": stale,
                "gse": summary.get("gse") if summary else None,
                "source": summary.get("source") if summary else None,
                "n_samples": summary.get("n_samples") if summary else None,
                "auc": meta.get("cv_auc") if meta else None,
                "shortlist_count": summary.get("shortlist_count") if summary else None,
                "ml_version": summary.get("ml_version") if summary else None,
            })
    rows.sort(key=lambda r: r["ts"], reverse=True)
    return rows


def delete_ml_runs(paths: list[str], roots: list[str] | None = None) -> list[str]:
    """带双重校验的删除:路径必须位于已登记结果根下,且必须是含 summary.json
    的 run 目录。任一校验不过则跳过并 warning(防御性:绝不 rmtree 未知目录)。
    每个成功删除 INFO 留痕。返回实际删除的路径列表。"""
    registered = [os.path.abspath(r) for r in (roots or _roots())]
    deleted: list[str] = []
    for p in paths:
        ap = os.path.abspath(p)
        if not any(ap == r or ap.startswith(r + os.sep) for r in registered):
            logging.warning("拒绝删除登记结果根之外的路径: %s", p)
            continue
        if not os.path.isfile(os.path.join(ap, "summary.json")):
            logging.warning("拒绝删除非 run 目录(缺 summary.json): %s", p)
            continue
        shutil.rmtree(ap, ignore_errors=True)
        if not os.path.isdir(ap):
            deleted.append(ap)
            logging.info("ML 结果目录已删除: %s", ap)
    return deleted


def zip_ml_runs(paths: list[str]) -> bytes:
    """把一个/多个 run 目录打包为 zip 字节;单个用 make_dir_zip,多个用 make_dirs_zip
    (zip 内以目录名为一层文件夹)。"""
    from .files import make_dir_zip, make_dirs_zip
    dirs = [p for p in paths if os.path.isdir(p)]
    if not dirs:
        return b""
    if len(dirs) == 1:
        return make_dir_zip(dirs[0])
    return make_dirs_zip({os.path.basename(d): d for d in dirs})


def cleanup_keep_per_gse(keep: int | None = None,
                         roots: list[str] | None = None) -> list[str]:
    """写侧自动清理(方案 §10 增补):同一 tag 只保留最近 keep 份,其余删除。

    语义边界:按新近度计数、与时间无关(今天连跑 keep+1 次,第 1 份也会被删);
    跨 tag 的总量永不自动清 —— 不同数据集的结果是实验记录,治理走管理页手动勾选。
    每次删除 INFO 留痕。返回实际删除的路径列表。
    """
    from WY import config
    keep = int(keep or getattr(config, "ML_RESULTS_KEEP_PER_GSE", 6))
    by_tag: dict[str, list[dict]] = {}
    for r in list_ml_runs(roots):
        if r["tag"]:
            by_tag.setdefault(r["tag"], []).append(r)
    doomed: list[str] = []
    for tag, rs in by_tag.items():
        rs.sort(key=lambda r: r["ts"], reverse=True)
        doomed.extend(r["path"] for r in rs[keep:])
    if doomed:
        logging.info("结果自动清理(同 tag 保留最近 %d 份): 删除 %d 个目录: %s",
                     keep, len(doomed), [os.path.basename(d) for d in doomed])
    return delete_ml_runs(doomed, roots=roots)
