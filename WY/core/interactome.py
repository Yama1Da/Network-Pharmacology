"""背景互作组模块：STRING v12.0 全网络本地缓存（邻近度 Z-score / RWR / DIAMOnD 的数学地基）。

职责边界：
- 本模块只管"全量背景网络"的下载、缓存与加载；查询型 PPI（输入基因之间的边）
  仍归 core.ppi，两者互不替代；
- 下载与构建只由 scripts/build_interactome.py 调用（CLI 离线构建），Streamlit
  页面只读缓存——绝不把大文件下载放进同步 UI；
- 缓存缺失时 load_interactome / load_degree_table 抛 InteractomeNotBuiltError，
  由调用方翻译为用户提示。绝不允许静默降级为"拿查询子网络当背景"：那会让
  度匹配随机化与显著性计算产生系统性偏差。

缓存内容（INTERACTOME_CACHE_DIR，文件名含 taxon/置信度/版本，改参自然换文件）：
  string_{taxon}_edges_{conf}_{ver}.parquet   边表 protein_a / protein_b / score(0-1000)
  string_{taxon}_degrees_{conf}_{ver}.csv     度数表 gene / degree
  string_{taxon}_meta_{conf}_{ver}.json       构建元数据（含实际文件名，见 files 键）
"""

import functools
import json
import logging
import os

import pandas as pd


class InteractomeNotBuiltError(RuntimeError):
    """互作组缓存尚未构建（提示运行 scripts/build_interactome.py）。"""


def _taxon(species: str) -> str:
    """校验物种已启用并返回 NCBI taxon id（如 'Human' -> '9606'）。"""
    import config

    if species not in config.INTERACTOME_ENABLED_SPECIES:
        raise ValueError(
            f"species {species!r} not enabled for interactome; "
            f"enabled: {config.INTERACTOME_ENABLED_SPECIES}")
    return str(config.SPECIES_MAP[species]["ncbi_taxon"])


def _file_stem(species: str, confidence: int, string_version: str) -> str:
    return f"string_{_taxon(species)}_{confidence}_{string_version}"


def _meta_path(species: str, confidence: int, string_version: str) -> str:
    import config

    return os.path.join(
        config.INTERACTOME_CACHE_DIR,
        _file_stem(species, confidence, string_version) + "_meta.json")


def _read_meta(species: str, confidence: int, string_version: str) -> dict | None:
    """读构建元数据；缺失/损坏返回 None（缓存视为不存在）。"""
    path = _meta_path(species, confidence, string_version)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            meta = json.load(f)
    except Exception:
        logging.exception("Interactome meta file unreadable: %s", path)
        return None
    # 元数据里登记的实际文件必须都在（parquet 失败回退 pickle 时文件名不同）
    for key in ("edges", "degrees"):
        name = meta.get("files", {}).get(key)
        if not name or not os.path.exists(os.path.join(_cache_dir(), name)):
            logging.warning("Interactome cache file missing for key %s: %r", key, name)
            return None
    return meta


def _cache_dir() -> str:
    import config

    return config.INTERACTOME_CACHE_DIR


def is_built(species: str = "Human",
             confidence: int | None = None,
             string_version: str | None = None) -> bool:
    """互作组缓存是否可用（元数据存在且登记的边表/度数表文件齐全）。"""
    import config

    confidence = int(confidence if confidence is not None else config.INTERACTOME_CONFIDENCE)
    string_version = string_version or config.STRING_VERSION
    return _read_meta(species, confidence, string_version) is not None


# ---------- 构建侧（CLI 专用） ----------


def _download_gz(url_path: str, dest: str,
                 progress_cb=None, progress_range: tuple[float, float] = (0.0, 1.0),
                 progress_label: str = "download") -> str:
    """下载 STRING gz 文件到 dest，返回本地路径。

    - url_path 为不含域名的下载路径（基址已含 /download 前缀，见 build_interactome）；
    - 候选源按 config 配置依次尝试（中国镜像优先、官方兜底）；
    - dest 已存在且非空时直接跳过下载（支持手动预下载/中断后重跑）；
    - 流式写 .part 临时文件 + 原子替换，带重试；全部失败抛 RuntimeError；
    - progress_cb(stage, fraction, detail)：下载进度回调（Streamlit 进度条用）。
      fraction 落在 progress_range 区间内；stage 取 progress_label；
      detail 为 "x/y MB"（语言中立）或 "cached"（跳过下载）。
    """
    import time

    import requests

    import config

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    lo, hi = progress_range
    if os.path.exists(dest) and os.path.getsize(dest) > 1024:
        logging.info("Download skipped, file already present: %s", dest)
        if progress_cb:
            progress_cb(progress_label, hi, "cached")
        return dest

    bases = [getattr(config, "STRING_DOWNLOAD_MIRROR_BASE", None),
             config.STRING_DOWNLOAD_BASE]
    bases = [b.rstrip("/") for b in bases if b]
    last_err: Exception | None = None
    for base in bases:
        url = base + url_path
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
            except Exception as exc:
                last_err = exc
                logging.exception("Download failed (attempt %d/%d): %s",
                                  attempt, config.API_MAX_RETRIES, url)
                if attempt < config.API_MAX_RETRIES:
                    time.sleep(config.API_RETRY_DELAY)
        logging.warning("Download source exhausted, trying next: %s failed", base)
    raise RuntimeError(f"download failed from all sources: {url_path} ({last_err})")


def _info_column_names(info_df: pd.DataFrame) -> tuple[str, str]:
    """从 protein.info 表头解析（ENSP列, 基因符号列），找不到时抛 ValueError。

    按表头取列而非按列序号——STRING 各版本列名/列序会微调
    （实测 v12.0 info.gz 表头为：#string_protein_id / preferred_name / protein_size / annotation）。
    """
    id_candidates = ("protein_external_id", "#protein_external_id",
                     "string_protein_id", "#string_protein_id", "protein_id")
    sym_candidates = ("preferredName", "preferred_name", "symbol")
    id_col = next((c for c in id_candidates if c in info_df.columns), None)
    sym_col = next((c for c in sym_candidates if c in info_df.columns), None)
    if not id_col or not sym_col:
        raise ValueError(f"unexpected protein.info columns: {list(info_df.columns)}")
    return id_col, sym_col


def _links_to_edges(links_df: pd.DataFrame, id2sym: dict, confidence: int) -> pd.DataFrame:
    """
    把 STRING links 全表转为符号级边表（纯函数，便于自检）。

    步骤（顺序有讲究）：
      1. 先按 combined_score 过滤（0-1000 整数），把千万级行砍到百万级再做映射；
      2. ENSP -> 基因符号映射，映射不到的丢弃；
      3. 去自环（一个 ENSP 对的两个蛋白塌缩到同一符号）；
      4. 无向边统一 (小, 大) 方向后去重，保留最高分——多条 ENSP 对塌缩到同一
         符号对时必须只留一条，否则度数会被系统性高估。

    Returns
    -------
    pd.DataFrame
        列 protein_a / protein_b / score（score 为 STRING 0-1000 整数），按分数降序。
    """
    needed = {"protein1", "protein2", "combined_score"}
    if not needed.issubset(links_df.columns):
        raise ValueError(f"unexpected links columns: {list(links_df.columns)}")

    edges = links_df[links_df["combined_score"] >= confidence].copy()
    edges["protein_a"] = edges["protein1"].map(id2sym)
    edges["protein_b"] = edges["protein2"].map(id2sym)
    edges = edges.dropna(subset=["protein_a", "protein_b"])
    edges = edges[edges["protein_a"] != edges["protein_b"]]

    lo = pd.Series([min(a, b) for a, b in zip(edges["protein_a"], edges["protein_b"])],
                   index=edges.index)
    hi = pd.Series([max(a, b) for a, b in zip(edges["protein_a"], edges["protein_b"])],
                   index=edges.index)
    edges["_lo"], edges["_hi"] = lo, hi
    edges = (edges.sort_values("combined_score", ascending=False)
                  .drop_duplicates(subset=["_lo", "_hi"]))
    out = (edges[["_lo", "_hi", "combined_score"]]
           .rename(columns={"_lo": "protein_a", "_hi": "protein_b",
                            "combined_score": "score"})
           .reset_index(drop=True))
    logging.info("Interactome edges: %d (confidence>=%d)", len(out), confidence)
    return out


def _write_edges_auto(df: pd.DataFrame, stem: str) -> str:
    """边表落盘：优先 parquet（pyarrow），引擎不可用时回退 pickle。返回实际文件路径。"""
    try:
        path = stem + ".parquet"
        df.to_parquet(path, index=False)
        return path
    except Exception:
        logging.exception("Parquet write failed — falling back to pickle")
        path = stem + ".pkl"
        df.to_pickle(path)
        return path


def _read_edges_auto(path: str) -> pd.DataFrame:
    if path.endswith(".parquet"):
        return pd.read_parquet(path)
    return pd.read_pickle(path)


def build_interactome(species: str = "Human",
                      confidence: int | None = None,
                      string_version: str | None = None,
                      keep_raw: bool = False,
                      progress_cb=None) -> dict:
    """
    下载并构建全互作组缓存（CLI：scripts/build_interactome.py；或页面内一键构建）。

    幂等：同参数重跑直接覆盖同名文件。构建成功后清空运行期加载缓存。
    原始 gz 默认构建成功后删除（可 --keep-raw 保留调试），删了随时可重新下载。

    Parameters
    ----------
    progress_cb : callable | None
        进度回调 progress_cb(stage, fraction, detail)，供 UI 渲染进度条：
          stage ∈ {"download_links", "download_info", "processing", "writing"}（语言中立键），
          fraction ∈ [0,1]（None 表示不确定进度，如解析清洗阶段），
          detail 为 "x/y MB" / "cached" / ""。回调异常被吞掉（仅记日志），不影响构建。
    """
    import shutil
    from datetime import datetime

    import config

    confidence = int(confidence if confidence is not None else config.INTERACTOME_CONFIDENCE)
    string_version = string_version or config.STRING_VERSION
    taxon = _taxon(species)
    os.makedirs(config.INTERACTOME_CACHE_DIR, exist_ok=True)
    raw_dir = os.path.join(config.INTERACTOME_CACHE_DIR, "raw")

    # 基址（config.STRING_DOWNLOAD_*）已含 /download 前缀，相对路径不再重复携带
    links_rel = f"/protein.links.{string_version}/{taxon}.protein.links.{string_version}.txt.gz"
    info_rel = f"/protein.info.{string_version}/{taxon}.protein.info.{string_version}.txt.gz"

    def report(stage: str, fraction: float | None, detail: str = "") -> None:
        """进度回调安全封装：回调异常只记日志，绝不打断构建本身。"""
        if progress_cb is None:
            return
        try:
            progress_cb(stage, fraction, detail)
        except Exception:
            logging.exception("Interactome progress callback raised; ignoring")

    logging.info("Downloading STRING %s interactome (taxon %s, confidence>=%d)...",
                 string_version, taxon, confidence)
    # 进度区间分配：互作全网络占大头（0.05~0.95），映射表很小（0.95~1.0）
    links_gz = _download_gz(links_rel, os.path.join(raw_dir, os.path.basename(links_rel)),
                            progress_cb=report, progress_range=(0.05, 0.95),
                            progress_label="download_links")
    info_gz = _download_gz(info_rel, os.path.join(raw_dir, os.path.basename(info_rel)),
                           progress_cb=report, progress_range=(0.95, 1.0),
                           progress_label="download_info")

    report("processing", None)
    # ---- info: ENSP -> gene symbol ----
    info_df = pd.read_csv(info_gz, sep="\t", compression="gzip")
    id_col, sym_col = _info_column_names(info_df)
    id2sym = dict(zip(info_df[id_col].astype(str),
                      info_df[sym_col].astype(str).str.strip().str.upper()))
    logging.info("protein.info: %d proteins mapped to symbols", len(id2sym))

    # ---- links: 过滤 -> 映射 -> 去自环 -> 去重 ----
    links_df = pd.read_csv(links_gz, sep=" ", compression="gzip",
                           usecols=["protein1", "protein2", "combined_score"])
    edges = _links_to_edges(links_df, id2sym, confidence)
    del links_df  # 及时释放千万级行的内存

    degree_series = pd.concat([edges["protein_a"], edges["protein_b"]]).value_counts()
    degree_df = pd.DataFrame({"gene": degree_series.index.astype(str),
                              "degree": degree_series.values.astype(int)})

    stem = os.path.join(config.INTERACTOME_CACHE_DIR,
                        _file_stem(species, confidence, string_version))
    report("writing", None)
    edges_path = _write_edges_auto(edges, stem)
    deg_path = stem + "_degrees.csv"
    degree_df.to_csv(deg_path, index=False, encoding="utf-8-sig")

    meta = {
        "species": species,
        "taxon": taxon,
        "string_version": string_version,
        "confidence": confidence,
        "n_nodes": int(len(degree_df)),
        "n_edges": int(len(edges)),
        "built_at": datetime.now().isoformat(timespec="seconds"),
        "source_urls": [config.STRING_DOWNLOAD_MIRROR_BASE + links_rel,
                        config.STRING_DOWNLOAD_BASE + links_rel],
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


# ---------- 加载侧（页面 / 算法用） ----------


def require_meta(species: str = "Human") -> dict:
    """取构建元数据；缓存未构建时抛 InteractomeNotBuiltError（绝不静默降级）。"""
    import config

    meta = _read_meta(species, config.INTERACTOME_CONFIDENCE, config.STRING_VERSION)
    if meta is None:
        raise InteractomeNotBuiltError(
            f"Interactome cache not built for {species}. "
            "Run: python scripts/build_interactome.py")
    return meta


@functools.lru_cache(maxsize=4)
def _load_graph_cached(edges_filename: str) -> "nx.Graph":
    """按文件名加载并冻结图（lru_cache：Streamlit 每次 rerun 都会进页面，不能重复读盘）。"""
    import networkx as nx

    df = _read_edges_auto(os.path.join(_cache_dir(), edges_filename))
    graph = nx.from_pandas_edgelist(df, "protein_a", "protein_b", edge_attr="score")
    nx.freeze(graph)  # 防止后续代码意外加/删边（可变图是隐性偏差来源）
    logging.info("Interactome graph loaded: %d nodes / %d edges (%s)",
                 graph.number_of_nodes(), graph.number_of_edges(), edges_filename)
    return graph


def load_interactome(species: str = "Human") -> "nx.Graph":
    """加载全互作组图（只读、已 freeze、进程内缓存）。缺缓存抛 InteractomeNotBuiltError。"""
    meta = require_meta(species)
    return _load_graph_cached(meta["files"]["edges"])


@functools.lru_cache(maxsize=4)
def _load_degrees_cached(degrees_filename: str) -> pd.DataFrame:
    return pd.read_csv(os.path.join(_cache_dir(), degrees_filename))


def load_degree_table(species: str = "Human") -> pd.DataFrame:
    """加载度数表（gene / degree）。度匹配随机化（M3）每次都要用，独立缓存。"""
    meta = require_meta(species)
    return _load_degrees_cached(meta["files"]["degrees"])


def interactome_age_days(species: str = "Human") -> int | None:
    """背景网络构建至今的天数（meta.built_at 距今）；未构建返回 None。"""
    import config
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
    """背景网络是否超过陈旧阈值（config.INTERACTOME_STALE_DAYS，默认365天）。

    语义提醒：stale ≠ 需要重下同版本（同版本文件是冻结的），而是提示
    "构建至今较久，建议检查 STRING 是否已发布新版本（改 STRING_VERSION 后重建）"。
    未构建时返回 False（此时应报"未构建"而非"陈旧"，避免双重误导）。
    """
    import config

    days = days if days is not None else config.INTERACTOME_STALE_DAYS
    age = interactome_age_days(species)
    return age is not None and age > days


@functools.lru_cache(maxsize=4)
def _known_genes(species: str) -> frozenset:
    return frozenset(load_degree_table(species)["gene"])


def get_network_coverage(genes: list[str], species: str = "Human") -> tuple[list[str], list[str]]:
    """
    把基因列表按"是否存在于背景网络"二分（M3/M4 的覆盖率守卫用此函数）。

    输入会做 strip + 大写归一化；保持去重、保序。
    """
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
    """清空运行期加载缓存（重建缓存后调用；CLI 进程退出亦自然失效）。"""
    _load_graph_cached.cache_clear()
    _load_degrees_cached.cache_clear()
    _known_genes.cache_clear()
