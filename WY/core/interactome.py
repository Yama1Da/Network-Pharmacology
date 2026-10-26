"""背景互作组模块：从 STRING 数据库构建"全蛋白互作背景网络"的本地缓存。

背景网络是网络邻近度、RWR 等靶点优先级计算的参照;输入基因之间的查询型
PPI 由 core.ppi 负责,两者互不替代。缓存缺失时会报错提示先构建,绝不会
拿查询基因的小网络冒充背景(那会让随机对照与显著性计算失真)。
"""

import functools
import json
import logging
import os

import numpy as np
import pandas as pd


class InteractomeNotBuiltError(RuntimeError):
    """背景网络缓存尚未构建时报错,提示先运行构建脚本。"""


def _taxon(species: str) -> str:
    """检查物种是否已启用,返回其 NCBI 物种编号(如 'Human' -> '9606')。"""
    from WY import config

    if species not in config.INTERACTOME_ENABLED_SPECIES:
        raise ValueError(
            f"species {species!r} not enabled for interactome; "
            f"enabled: {config.INTERACTOME_ENABLED_SPECIES}")
    return str(config.SPECIES_MAP[species]["ncbi_taxon"])


def _file_stem(species: str, confidence: int, string_version: str) -> str:
    return f"string_{_taxon(species)}_{confidence}_{string_version}"


def _meta_path(species: str, confidence: int, string_version: str) -> str:
    from WY import config

    return os.path.join(
        config.INTERACTOME_CACHE_DIR,
        _file_stem(species, confidence, string_version) + "_meta.json")


def _read_meta(species: str, confidence: int, string_version: str) -> dict | None:
    """读取构建记录;缺失或损坏时返回 None,视同缓存不存在。"""
    path = _meta_path(species, confidence, string_version)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            meta = json.load(f)
    except Exception:
        logging.exception("Interactome meta file unreadable: %s", path)
        return None
    # 构建记录里登记的边表/度数表文件必须都存在,否则视同缓存缺失
    for key in ("edges", "degrees"):
        name = meta.get("files", {}).get(key)
        if not name or not os.path.exists(os.path.join(_cache_dir(), name)):
            logging.warning("Interactome cache file missing for key %s: %r", key, name)
            return None
    return meta


def _cache_dir() -> str:
    from WY import config

    return config.INTERACTOME_CACHE_DIR


def is_built(species: str = "Human",
             confidence: int | None = None,
             string_version: str | None = None) -> bool:
    """背景网络缓存是否已构建且文件齐全可用。"""
    from WY import config

    confidence = int(confidence if confidence is not None else config.INTERACTOME_CONFIDENCE)
    string_version = string_version or config.STRING_VERSION
    return _read_meta(species, confidence, string_version) is not None


# ---------- 构建侧(命令行专用) ----------

# 分段并发下载的适用下限:文件太小不值得多连接;单段太小则握手/重试开销占主导
_RANGE_MIN_SIZE = 8 * 1024 * 1024
_RANGE_SEG_MIN_SIZE = 8 * 1024 * 1024
# 分段整体重试轮数:每轮只重跑未完成的段(段内从盘上续传),避免个别掉队段
# 连累全部已完成段一起作废重下
_RANGED_MAX_ROUNDS = 3


def _probe_range_support(url: str) -> int:
    """探测下载源是否支持 HTTP Range 分段:支持返回文件总字节数,否则返回 0。
    用 1 字节区间实测(部分服务器 HEAD 不回 Accept-Ranges,206 + Content-Range 才算数);
    stream=True 下立即关闭,不会把 200 的完整响应体拉下来。"""
    import requests

    from WY import config

    try:
        with requests.get(url, headers={"Range": "bytes=0-0"}, stream=True,
                          timeout=config.API_TIMEOUT) as resp:
            if resp.status_code != 206:
                return 0
            total = resp.headers.get("Content-Range", "").rsplit("/", 1)[-1]
            return int(total) if total.isdigit() else 0
    except Exception:
        logging.exception("Range probe failed: %s", url)
        return 0


def string_version_exists(version: str) -> bool:
    """轻量探测下载站是否已提供该版本的数据文件(protein.info 小文件,流式即开即关)。
    200/206 即认为版本存在(不要求 Range 支持);任一下载源可用即可。
    供设置面板在保存版本切换前拦截拼错/尚未发布的版本号。"""
    import requests

    from WY import config

    rel = f"/protein.info.{version}/{_taxon('Human')}.protein.info.{version}.txt.gz"
    for base in config.INTERACTOME_DOWNLOAD_SOURCES:
        try:
            with requests.get(base.rstrip("/") + rel, stream=True,
                              timeout=config.API_TIMEOUT) as resp:
                if resp.status_code in (200, 206):
                    return True
        except Exception:
            logging.exception("STRING version probe failed: %s", base + rel)
    return False


def _download_segment(url: str, start: int, end: int, path: str,
                      done: list, lock) -> None:
    """下载 [start, end] 字节区间到 path(工作线程内执行,不碰进度回调);
    整段重试耗尽才抛错,由调用方按轮次只重跑失败段。支持断点续传:重试只拉取
    未完成的尾部区间并追加写入;path 已有部分内容(上一轮的成果)时从盘上字节数
    续传。done[0] 只累计新落盘的字节,某次尝试失败时回滚该次计入的部分。"""
    import time

    import requests

    from WY import config

    expected = end - start + 1
    have = 0
    if os.path.exists(path):
        have = os.path.getsize(path)
        if have > expected:  # 残留文件不可信(如上次越界多写),删除重下
            try:
                os.remove(path)
            except OSError:
                pass
            have = 0
    if have == expected:
        return
    last_err: Exception | None = None
    for attempt in range(1, config.API_MAX_RETRIES + 1):
        base = have
        overrun = False
        try:
            with requests.get(url, headers={"Range": f"bytes={start + have}-{end}"},
                              stream=True, timeout=config.API_TIMEOUT) as resp:
                if resp.status_code != 206:
                    raise RuntimeError(f"http {resp.status_code}: range not honored")
                with open(path, "ab" if have else "wb") as f:
                    for chunk in resp.iter_content(chunk_size=1 << 20):
                        if not chunk:
                            continue
                        f.write(chunk)
                        have += len(chunk)
                        with lock:
                            done[0] += len(chunk)
                        if have > expected:
                            overrun = True
                            break
                if overrun:
                    raise RuntimeError(f"server overran range ({have} > {expected})")
                if have == expected:
                    return
                raise RuntimeError(f"incomplete segment: {have}/{expected} bytes")
        except Exception as exc:
            last_err = exc
            with lock:
                done[0] -= have - base
            if have > expected:
                # 越界多写后文件不可信,删除并从段首重来
                try:
                    os.remove(path)
                except OSError:
                    pass
                with lock:
                    done[0] -= have
                have = 0
            logging.warning("Segment %s failed (attempt %d/%d): %s",
                            os.path.basename(path), attempt,
                            config.API_MAX_RETRIES, exc)
            if attempt < config.API_MAX_RETRIES:
                time.sleep(config.API_RETRY_DELAY)
    raise RuntimeError(f"segment exhausted retries: {last_err}")


def _download_gz_ranged(url: str, dest: str, progress_cb,
                        lo: float, hi: float, progress_label: str) -> tuple[str, int]:
    """并发分段下载:源支持 HTTP Range 时按 config.INTERACTOME_DOWNLOAD_CONCURRENCY
    切段并行拉取;某段失败只重跑该段(段内续传,已完成段不动),最多
    _RANGED_MAX_ROUNDS 轮,全部成功后按序拼接(先写 .part 再原子替换,不留半截文件)。
    返回 (状态, 探测到的总字节数):"ok" 成功;"skip" 源不支持 Range 或文件太小
    (本就不该分段,不算失败);"failed" 重试轮次耗尽(半成品已清掉),由调用方回退单流。"""
    import shutil
    import threading
    import time
    from concurrent.futures import ThreadPoolExecutor

    from WY import config

    total = _probe_range_support(url)
    if total < _RANGE_MIN_SIZE:
        return "skip", total
    n_seg = min(int(getattr(config, "INTERACTOME_DOWNLOAD_CONCURRENCY", 1) or 1),
                total // _RANGE_SEG_MIN_SIZE)
    if n_seg <= 1:
        return "skip", total
    seg_size = total // n_seg
    # 末段到文件末尾,吸收整除余数
    bounds = [(i * seg_size, (i + 1) * seg_size - 1) for i in range(n_seg - 1)]
    bounds.append(((n_seg - 1) * seg_size, total - 1))
    parts = [f"{dest}.part{i}" for i in range(n_seg)]
    for p in parts:  # 清掉上次进程中断可能残留的段文件:本次调用从零开始
        try:
            os.remove(p)
        except OSError:
            pass
    done = [0]
    lock = threading.Lock()
    logging.info("Ranged download %s: %d segments x ~%.0f MB",
                 url, n_seg, seg_size / 1e6)
    try:
        pending = list(range(n_seg))
        for round_no in range(1, _RANGED_MAX_ROUNDS + 1):
            with ThreadPoolExecutor(max_workers=len(pending)) as pool:
                futures = {pool.submit(_download_segment, url, bounds[i][0],
                                       bounds[i][1], parts[i], done, lock): i
                           for i in pending}
                # 进度回调只在本线程发(回调可能操作 Streamlit 组件,不能跨线程调用):
                # 工作线程只更新计数字节,这里每半秒聚合一次
                while not all(f.done() for f in futures):
                    time.sleep(0.5)
                    if progress_cb:
                        with lock:
                            cur = done[0]
                        progress_cb(progress_label,
                                    lo + (hi - lo) * min(cur / total, 1.0),
                                    f"{cur / 1e6:.0f}/{total / 1e6:.0f} MB")
            pending = [i for f, i in futures.items() if f.exception()]
            if not pending:
                break
            logging.warning("Ranged download round %d: %d/%d segment(s) failed, "
                            "retrying only those", round_no, len(pending), n_seg)
            time.sleep(config.API_RETRY_DELAY)
        if pending:
            logging.warning("Ranged download gave up after %d round(s) (%d segment "
                            "unrecovered), falling back to single stream",
                            _RANGED_MAX_ROUNDS, len(pending))
            return "failed", total
        tmp = dest + ".part"
        with open(tmp, "wb") as out:
            for p in parts:
                with open(p, "rb") as f:
                    shutil.copyfileobj(f, out)
        if os.path.getsize(tmp) != total:
            logging.warning("Ranged download size mismatch (%d != %d), falling "
                            "back to single stream", os.path.getsize(tmp), total)
            return "failed", total
        os.replace(tmp, dest)
        logging.info("Ranged download finished: %s (%.1f MB)", dest, total / 1e6)
        return "ok", total
    finally:
        for p in parts:  # 成功时各段已拼入目标文件,失败时丢弃半成品
            try:
                os.remove(p)
            except OSError:
                pass
        try:
            os.remove(dest + ".part")
        except OSError:
            pass


def _download_gz(url_path: str, dest: str,
                 progress_cb=None, progress_range: tuple[float, float] = (0.0, 1.0),
                 progress_label: str = "download") -> str:
    """从 STRING 下载压缩文件到本地,返回文件路径;本地已有同名文件时直接跳过下载,
    多个下载源依次尝试,全部失败时报错。每个源先尝试并发分段(Range),未成功自动
    回退单流;大文件首次回退时经进度回调发一次用户提醒(单流约需 10 分钟以上),
    小文件单流本来就快,不提醒以免误报。"""
    import time

    import requests

    from WY import config

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    lo, hi = progress_range
    if os.path.exists(dest) and os.path.getsize(dest) > 1024:
        logging.info("Download skipped, file already present: %s", dest)
        if progress_cb:
            progress_cb(progress_label, hi, "cached")
        return dest

    bases = [b.rstrip("/") for b in config.INTERACTOME_DOWNLOAD_SOURCES]
    warned_fallback = False  # 提醒每次下载只发一次,两个源都回退时不重复打扰
    last_err: Exception | None = None
    for base in bases:
        url = base + url_path
        status, total = _download_gz_ranged(url, dest, progress_cb, lo, hi,
                                            progress_label)
        if status == "ok":
            return dest
        if not warned_fallback and total >= _RANGE_MIN_SIZE:
            warned_fallback = True
            host = base.split("//", 1)[-1]
            logging.warning("Concurrent download did not finish on %s (status %s), "
                            "falling back to single stream (may take 10+ minutes)",
                            host, status)
            if progress_cb:
                progress_cb("fallback", None, host)
        for attempt in range(1, config.API_MAX_RETRIES + 1):
            try:
                with requests.get(url, stream=True, timeout=config.API_TIMEOUT) as resp:
                    resp.raise_for_status()
                    total = int(resp.headers.get("Content-Length") or 0)
                    tmp = dest + ".part"
                    done = 0
                    with open(tmp, "wb") as f:
                        for chunk in resp.iter_content(chunk_size=1 << 20):
                            if not chunk:
                                continue
                            f.write(chunk)
                            done += len(chunk)
                            if progress_cb and total:
                                frac = lo + (hi - lo) * min(done / total, 1.0)
                                progress_cb(progress_label, frac,
                                            f"{done / 1e6:.0f}/{total / 1e6:.0f} MB")
                    os.replace(tmp, dest)
                logging.info("Downloaded %s -> %s", url, dest)
                return dest
            except requests.HTTPError as exc:
                last_err = exc
                status = getattr(exc.response, "status_code", None)
                if status is not None and 400 <= status < 500:
                    # 404/410 等 4xx 是确定性失败(实测镜像仅覆盖 v12.0,其余版本一律
                    # 404),重试只会白等,立即换下一个下载源,并在进度条上告知原因
                    logging.warning("Download source %s returned %d for %s, skipping "
                                    "retries and trying next source", base, status, url)
                    if not warned_fallback and progress_cb:
                        warned_fallback = True
                        progress_cb("mirror_skip", None, base.split("//", 1)[-1])
                    break
                logging.exception("Download failed (attempt %d/%d): %s",
                                  attempt, config.API_MAX_RETRIES, url)
                if attempt < config.API_MAX_RETRIES:
                    time.sleep(config.API_RETRY_DELAY)
            except Exception as exc:
                last_err = exc
                logging.exception("Download failed (attempt %d/%d): %s",
                                  attempt, config.API_MAX_RETRIES, url)
                if attempt < config.API_MAX_RETRIES:
                    time.sleep(config.API_RETRY_DELAY)
        logging.warning("Download source exhausted, trying next: %s failed", base)
    raise RuntimeError(f"download failed from all sources: {url_path} ({last_err})")


def _info_column_names(info_df: pd.DataFrame) -> tuple[str, str]:
    """从蛋白信息表头里找出蛋白编号列和基因符号列;找不到时报错。"""
    id_candidates = ("protein_external_id", "#protein_external_id",
                     "string_protein_id", "#string_protein_id", "protein_id")
    sym_candidates = ("preferredName", "preferred_name", "symbol")
    id_col = next((c for c in id_candidates if c in info_df.columns), None)
    sym_col = next((c for c in sym_candidates if c in info_df.columns), None)
    if not id_col or not sym_col:
        raise ValueError(f"unexpected protein.info columns: {list(info_df.columns)}")
    return id_col, sym_col


def _links_to_edges(links_df: pd.DataFrame, id2sym: dict, confidence: int) -> pd.DataFrame:
    """把 STRING 原始互作表整理成基因符号层面的边表:按置信度阈值过滤后映射为基因符号,去掉自身互作,
    同一对基因只保留最高分的一条边,避免基因连接数被高估;结果按分数降序。"""
    needed = {"protein1", "protein2", "combined_score"}
    if not needed.issubset(links_df.columns):
        raise ValueError(f"unexpected links columns: {list(links_df.columns)}")

    edges = links_df[links_df["combined_score"] >= confidence].copy()
    edges["protein_a"] = edges["protein1"].map(id2sym)
    edges["protein_b"] = edges["protein2"].map(id2sym)
    edges = edges.dropna(subset=["protein_a", "protein_b"])
    edges = edges[edges["protein_a"] != edges["protein_b"]]

    # 向量化配对规范化:lo/hi 使 (A,B) 与 (B,A) 归为同一条边;
    # 列此时已是 dropna 后的纯字符串,逐元素字符串比较与逐行 min/max 结果一致
    a = edges["protein_a"].to_numpy()
    b = edges["protein_b"].to_numpy()
    edges["_lo"] = pd.Series(np.where(a <= b, a, b), index=edges.index)
    edges["_hi"] = pd.Series(np.where(a <= b, b, a), index=edges.index)
    edges = (edges.sort_values("combined_score", ascending=False)
                  .drop_duplicates(subset=["_lo", "_hi"]))
    out = (edges[["_lo", "_hi", "combined_score"]]
           .rename(columns={"_lo": "protein_a", "_hi": "protein_b",
                            "combined_score": "score"})
           .reset_index(drop=True))
    logging.info("Interactome edges: %d (confidence>=%d)", len(out), confidence)
    return out


def _write_edges_auto(df: pd.DataFrame, stem: str) -> str:
    """把边表保存到磁盘,返回实际文件路径。优先 parquet;失败时回退 CSV——
    纯文本没有反序列化执行风险,不再像旧版那样回退 pickle(本地缓存一旦被篡改即任意代码执行)。
    先写 .tmp 再原子替换:构建中途被 Stop/中断不会留下半截边表顶替旧缓存。"""
    try:
        path = stem + ".parquet"
        tmp = path + ".tmp"
        df.to_parquet(tmp, index=False)
        os.replace(tmp, path)
        return path
    except Exception:
        logging.exception("Parquet write failed — falling back to CSV")
        path = stem + ".csv"
        tmp = path + ".tmp"
        df.to_csv(tmp, index=False, encoding="utf-8-sig")
        os.replace(tmp, path)
        return path


def _read_edges_auto(path: str) -> pd.DataFrame:
    if path.endswith(".parquet"):
        return pd.read_parquet(path)
    if path.endswith(".pkl"):
        # 仅兼容读取旧版产生的 .pkl 缓存;新构建不再写 pickle
        return pd.read_pickle(path)
    return pd.read_csv(path, encoding="utf-8-sig", comment="#")


def build_interactome(species: str = "Human",
                      confidence: int | None = None,
                      string_version: str | None = None,
                      keep_raw: bool = False,
                      progress_cb=None) -> dict:
    """下载 STRING 全互作数据并构建本地背景网络缓存;同参数重跑会直接覆盖原文件,默认构建成功后删除原始下载文件。
    progress_cb 为进度回调,供界面显示进度条;回调出错只记日志,不影响构建。"""
    import shutil
    from datetime import datetime

    from WY import config

    confidence = int(confidence if confidence is not None else config.INTERACTOME_CONFIDENCE)
    string_version = string_version or config.STRING_VERSION
    taxon = _taxon(species)
    os.makedirs(config.INTERACTOME_CACHE_DIR, exist_ok=True)
    raw_dir = os.path.join(config.INTERACTOME_CACHE_DIR, "raw")

    # 下载基址已含 /download 前缀,这里只拼相对路径
    links_rel = f"/protein.links.{string_version}/{taxon}.protein.links.{string_version}.txt.gz"
    info_rel = f"/protein.info.{string_version}/{taxon}.protein.info.{string_version}.txt.gz"

    def report(stage: str, fraction: float | None, detail: str = "") -> None:
        """进度回调的包装:回调出错只记日志,不会中断构建。"""
        if progress_cb is None:
            return
        try:
            progress_cb(stage, fraction, detail)
        except Exception:
            logging.exception("Interactome progress callback raised; ignoring")

    logging.info("Downloading STRING %s interactome (taxon %s, confidence>=%d)...",
                 string_version, taxon, confidence)
    # 进度条区间:大文件下载占 0.05~0.95,小文件占 0.95~1.0
    links_gz = _download_gz(links_rel, os.path.join(raw_dir, os.path.basename(links_rel)),
                            progress_cb=report, progress_range=(0.05, 0.95),
                            progress_label="download_links")
    info_gz = _download_gz(info_rel, os.path.join(raw_dir, os.path.basename(info_rel)),
                           progress_cb=report, progress_range=(0.95, 1.0),
                           progress_label="download_info")

    report("processing", None)
    info_df = pd.read_csv(info_gz, sep="\t", compression="gzip")
    id_col, sym_col = _info_column_names(info_df)
    id2sym = dict(zip(info_df[id_col].astype(str),
                      info_df[sym_col].astype(str).str.strip().str.upper()))
    logging.info("protein.info: %d proteins mapped to symbols", len(id2sym))

    links_df = pd.read_csv(links_gz, sep=" ", compression="gzip",
                           usecols=["protein1", "protein2", "combined_score"])
    edges = _links_to_edges(links_df, id2sym, confidence)
    del links_df  # 及时释放内存

    degree_series = pd.concat([edges["protein_a"], edges["protein_b"]]).value_counts()
    degree_df = pd.DataFrame({"gene": degree_series.index.astype(str),
                              "degree": degree_series.values.astype(int)})

    stem = os.path.join(config.INTERACTOME_CACHE_DIR,
                        _file_stem(species, confidence, string_version))
    report("writing", None)
    edges_path = _write_edges_auto(edges, stem)
    deg_path = stem + "_degrees.csv"
    # 原子写:度表与边表同属三件套,中断不能留半截文件顶替旧缓存
    deg_tmp = deg_path + ".tmp"
    degree_df.to_csv(deg_tmp, index=False, encoding="utf-8-sig")
    os.replace(deg_tmp, deg_path)

    meta = {
        "species": species,
        "taxon": taxon,
        "string_version": string_version,
        "confidence": confidence,
        "n_nodes": int(len(degree_df)),
        "n_edges": int(len(edges)),
        "built_at": datetime.now().isoformat(timespec="seconds"),
        "source_urls": [b + links_rel for b in config.INTERACTOME_DOWNLOAD_SOURCES],
        "files": {
            "edges": os.path.basename(edges_path),
            "degrees": os.path.basename(deg_path),
        },
    }
    with open(_meta_path(species, confidence, string_version), "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    if not keep_raw:
        shutil.rmtree(raw_dir, ignore_errors=True)

    clear_runtime_caches()
    logging.info("Interactome built: %d nodes / %d edges (%s, confidence>=%d)",
                 meta["n_nodes"], meta["n_edges"], string_version, confidence)
    return meta


def clean_old_version_caches(species: str = "Human") -> int:
    """删除该物种其它 STRING 版本/置信度组合的背景网络缓存(边表/度表/构建记录三件套)。
    供数据管理页在新版本构建成功后调用(确认弹窗已告知清理代价):当前生效组合一律保留,
    切回旧组合需重新构建。只删各 _meta.json 登记的文件,无构建记录的孤儿文件不动。返回删除文件数。"""
    from WY import config

    cur_conf = int(config.INTERACTOME_CONFIDENCE)
    cur_ver = config.STRING_VERSION
    removed = 0
    try:
        meta_names = [n for n in os.listdir(_cache_dir()) if n.endswith("_meta.json")]
    except OSError:
        return 0
    for name in meta_names:
        path = os.path.join(_cache_dir(), name)
        try:
            with open(path, "r", encoding="utf-8") as f:
                meta = json.load(f)
            if not isinstance(meta, dict):
                raise ValueError("meta is not a dict")
        except Exception:
            # 读不出的记录无法判断归属,保守保留
            logging.exception("Interactome meta unreadable, kept: %s", path)
            continue
        if meta.get("species") != species:
            continue
        try:
            is_current = (int(meta.get("confidence", -1)) == cur_conf
                          and meta.get("string_version") == cur_ver)
        except (TypeError, ValueError):
            is_current = False  # 置信度都解析不出的记录必不属于当前组合
        if is_current:
            continue  # 当前生效组合,保留
        # 非当前组合:尽力删除登记的边表/度表,再删除记录本身;清理绝不抛异常,
        # 避免传到页面把成功的构建误报为失败
        files = meta.get("files")
        if isinstance(files, dict):
            for key in ("edges", "degrees"):
                fname = files.get(key)
                if not fname:
                    continue
                try:
                    os.remove(os.path.join(_cache_dir(), fname))
                    removed += 1
                except OSError:
                    logging.exception("Failed to remove old interactome cache: %s", fname)
        try:
            os.remove(path)
            removed += 1
        except OSError:
            logging.exception("Failed to remove old interactome meta: %s", path)
    if removed:
        logging.info("Old interactome caches removed: %d file(s) (active: %s, confidence>=%d)",
                     removed, cur_ver, cur_conf)
    return removed


# ---------- 加载侧(页面/算法用) ----------


def require_meta(species: str = "Human") -> dict:
    """读取构建记录;缓存未构建时报错提示,绝不静默降级。"""
    from WY import config

    meta = _read_meta(species, config.INTERACTOME_CONFIDENCE, config.STRING_VERSION)
    if meta is None:
        raise InteractomeNotBuiltError(
            f"Interactome cache not built for {species}. "
            "Run: python scripts/build_interactome.py")
    return meta


@functools.lru_cache(maxsize=4)
def _load_graph_cached(edges_filename: str) -> "nx.Graph":
    """从磁盘加载背景网络图并设为只读;进程内只读一次盘,之后复用。"""
    import networkx as nx

    df = _read_edges_auto(os.path.join(_cache_dir(), edges_filename))
    graph = nx.from_pandas_edgelist(df, "protein_a", "protein_b", edge_attr="score")
    nx.freeze(graph)  # 设为只读,防止后续计算误改网络
    logging.info("Interactome graph loaded: %d nodes / %d edges (%s)",
                 graph.number_of_nodes(), graph.number_of_edges(), edges_filename)
    return graph


def load_interactome(species: str = "Human") -> "nx.Graph":
    """加载全蛋白背景互作网络(只读,自动复用缓存);未构建时报错。"""
    meta = require_meta(species)
    return _load_graph_cached(meta["files"]["edges"])


@functools.lru_cache(maxsize=4)
def _load_degrees_cached(degrees_filename: str) -> pd.DataFrame:
    return pd.read_csv(os.path.join(_cache_dir(), degrees_filename))


def load_degree_table(species: str = "Human") -> pd.DataFrame:
    """加载各基因的连接数(度)表,供随机对照使用;自动复用缓存。"""
    meta = require_meta(species)
    return _load_degrees_cached(meta["files"]["degrees"])


def interactome_age_days(species: str = "Human") -> int | None:
    """背景网络构建至今的天数;未构建返回 None。"""
    from WY import config
    from datetime import datetime

    meta = _read_meta(species, config.INTERACTOME_CONFIDENCE, config.STRING_VERSION)
    if meta is None:
        return None
    try:
        built_at = datetime.fromisoformat(meta["built_at"])
    except (KeyError, ValueError):
        return None
    return (datetime.now() - built_at).days


def is_stale(species: str = "Human", days: int | None = None) -> bool:
    """背景网络是否过于陈旧(超过 config 设定的天数阈值,默认365天);陈旧只是提示该检查 STRING 是否已发布新版本。
    未构建时返回 False。"""
    from WY import config

    days = days if days is not None else config.INTERACTOME_STALE_DAYS
    age = interactome_age_days(species)
    return age is not None and age > days


@functools.lru_cache(maxsize=4)
def _known_genes(species: str) -> frozenset:
    return frozenset(load_degree_table(species)["gene"])


def get_network_coverage(genes: list[str], species: str = "Human") -> tuple[list[str], list[str]]:
    """把基因列表分成"在背景网络中"和"不在"两组;基因名自动去空格并统一大写,结果去重且保持原顺序。"""
    known = _known_genes(species)
    seen: set[str] = set()
    inside: list[str] = []
    outside: list[str] = []
    for raw in genes:
        g = str(raw).strip().upper()
        if not g or g in seen:
            continue
        seen.add(g)
        (inside if g in known else outside).append(g)
    return inside, outside


def clear_runtime_caches() -> None:
    """清空内存中的加载缓存,重建背景网络后调用。"""
    _load_graph_cached.cache_clear()
    _load_degrees_cached.cache_clear()
    _known_genes.cache_clear()
    # 邻近度/RWR 的 LCC 上下文缓存也依赖背景网络:不清的话同进程重建后
    # 仍会用旧图算零分布(函数内导入,避免与 network_algorithms 循环依赖)
    from . import network_algorithms

    network_algorithms._lcc_context_cached.cache_clear()
