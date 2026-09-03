"""靶点汇总：整合药物靶点（ChEMBL+DGIdb）与疾病靶点（Open Targets+DisGeNET），并取交集。"""

import logging
import os

import pandas as pd

from .caching import open_csv_cache, save_csv_cache
from .chembl import fetch_drug_targets_chembl
from .disgenet import fetch_disease_targets_disgenet
from .dgidb import fetch_drug_targets_dgidb
from .opentargets import get_disease_targets
from .text_utils import remember_name_mapping, safe_name


def _max_score_map(df: pd.DataFrame, score_col: str) -> dict[str, float]:
    """按基因符号分组取分数最大值；符号统一大写，分数统一转数值（非数值置空）。"""
    if df.empty or score_col not in df.columns:
        return {}
    tmp = df.dropna(subset=["approvedSymbol"]).copy()
    tmp["approvedSymbol"] = tmp["approvedSymbol"].str.upper()
    tmp[score_col] = pd.to_numeric(tmp[score_col], errors="coerce")
    return tmp.groupby("approvedSymbol")[score_col].max().to_dict()


def _merge_source_pair(
    name: str,
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    *,
    cache_path: str,
    kind: str,
    src_a: str,
    src_b: str,
    score_a: tuple[str, str],
    score_b: tuple[str, str],
    sources_desc: tuple[str, str],
    save: bool = True,
) -> pd.DataFrame:
    """两个数据源的靶点按基因符号取并集去重，并保留各源分数最大值。

    kind 为 "drug"/"disease"，用于日志；src_a/src_b 为来源标记列名；
    score_a/score_b 为 (输出列名, 源分数列名)；两源都没有数据时返回空表。
    save=False 时不写缓存(某路数据源失败时,避免把残缺结果永久缓存)。
    """
    genes_a = set(df_a["approvedSymbol"].dropna().str.upper()) if not df_a.empty else set()
    genes_b = set(df_b["approvedSymbol"].dropna().str.upper()) if not df_b.empty else set()
    all_genes = sorted(genes_a | genes_b)
    if not all_genes:
        logging.warning("No %s targets from either source for '%s'", kind, name)
        return pd.DataFrame()

    col_a, raw_a = score_a
    col_b, raw_b = score_b
    map_a = _max_score_map(df_a, raw_a)
    map_b = _max_score_map(df_b, raw_b)

    records: list[dict] = []
    for g in all_genes:
        in_a = g in genes_a
        in_b = g in genes_b
        records.append({
            "approvedSymbol": g,
            src_a: int(in_a),
            src_b: int(in_b),
            "source_count": int(in_a) + int(in_b),
            col_a: map_a.get(g),
            col_b: map_b.get(g),
        })

    merged = pd.DataFrame(records)
    if save:
        save_csv_cache(merged, cache_path, skip_empty=False)
    elif not merged.empty:
        logging.warning(
            "Merged %s '%s' NOT cached because some sources failed", kind, name)
    logging.info(
        "Merged %s '%s': %d targets (%s %d, %s %d)",
        kind, name, len(merged), sources_desc[0], len(genes_a), sources_desc[1], len(genes_b),
    )
    return merged


def get_drug_targets_merged(drug_name: str) -> pd.DataFrame:
    """整合 ChEMBL 与 DGIdb 两个数据库的药物靶点，按基因符号取并集去重。

    为每个基因标记来自哪些数据库（source_chmbl / source_dgidb），
    并保留各库的活性/交互分数。两个库都没有数据时返回空表；
    结果存为本地 CSV，已有缓存则直接复用。
    """
    import config

    cache_path = os.path.join(config.MERGED_CACHE_DIR, f"merged_drug_{safe_name(drug_name)}.csv")
    remember_name_mapping(drug_name)  # 记录安全文件名→原始名,供管理页还原
    cached = open_csv_cache(cache_path, label=f"merged drug targets for '{drug_name}'")
    if cached is not None:
        return cached

    chmbl_df, _ = fetch_drug_targets_chembl(drug_name)
    dgidb_df = fetch_drug_targets_dgidb(drug_name)
    return _merge_source_pair(
        drug_name, chmbl_df, dgidb_df,
        cache_path=cache_path, kind="drug",
        src_a="source_chmbl", src_b="source_dgidb",
        score_a=("score_pchembl", "pchembl_value"),
        score_b=("score_dgidb", "interaction_score"),
        sources_desc=("ChEMBL", "DGIdb"),
    )


def get_disease_targets_merged(disease_name: str) -> pd.DataFrame:
    """整合 Open Targets 与 DisGeNET 两个数据库的疾病靶点，按基因符号取并集去重。

    为每个基因标记来自哪些数据库（source_ot / source_disgenet），
    并保留各库的关联分数。两个库都没有数据时返回空表；
    结果存为本地 CSV，已有缓存则直接复用。
    """
    import config

    cache_path = os.path.join(config.MERGED_CACHE_DIR, f"merged_disease_{safe_name(disease_name)}.csv")
    remember_name_mapping(disease_name)  # 记录安全文件名→原始名,供管理页还原
    cached = open_csv_cache(cache_path, label=f"merged disease targets for '{disease_name}'")
    if cached is not None:
        return cached

    ot_df = get_disease_targets(disease_name)
    disgenet_df = fetch_disease_targets_disgenet(disease_name)
    return _merge_source_pair(
        disease_name, ot_df, disgenet_df,
        cache_path=cache_path, kind="disease",
        src_a="source_ot", src_b="source_disgenet",
        score_a=("score_ot", "association_score"),
        score_b=("score_disgenet", "disgenet_score"),
        sources_desc=("OpenTargets", "DisGeNET"),
    )


def get_common_targets(drug_name: str, disease_name: str, top_n: int = 1000) -> list[str]:
    """
    取药物靶点与疾病靶点的交集基因，作为后续 PPI 网络分析的输入。

    按字母序排序，最多返回 top_n 个（默认 1000，对齐 STRING 查询上限）；
    任一侧无数据时返回空列表。
    """
    drug_df, disease_df, _failed = get_targets_concurrent(drug_name, disease_name)
    # 提取药物和疾病两侧的基因符号
    drug_set = set(drug_df["approvedSymbol"].dropna()) if not drug_df.empty else set()
    disease_set = set(disease_df["approvedSymbol"].dropna()) if not disease_df.empty else set()
    # 取交集，按字母序排序并截断到 top_n
    return sorted(drug_set & disease_set)[:top_n]




# 各数据源单路等待期限(秒):超时先按空表继续,线程留在后台跑完并写缓存,
# 下次分析即可命中。ChEMBL 慢药物的批量分页留足余量;
# DisGeNET 偶发长时间不返回数据,期限卡住最坏等待。
_SOURCE_DEADLINE_SECONDS = {
    "chembl": 600,
    "dgidb": 120,
    "opentargets": 120,
    "disgenet": 240,
}


def get_targets_concurrent(
    drug_name: str, disease_name: str, progress_cb=None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """并发获取药物与疾病两侧靶点(四路数据源同时请求),供分析页一次拿全。

    语义与分别串行调用 get_drug_targets_merged + get_disease_targets_merged 一致:
    两侧 merged 缓存都命中时直接返回;否则只对缺失的一侧并发联网。

    progress_cb: 可选进度回调,每 0.5 秒或某路完成时被调用一次(始终在主线程,
    可安全更新 Streamlit 组件),签名 progress_cb(done: list[str], total: int),
    done 为已完成的数据源 key(chembl/dgidb/opentargets/disgenet)。

    返回 (药物靶点表, 疾病靶点表, 失败数据源 key 列表)。
    失败的数据源按空表参与展示,且该侧不写 merged 缓存——避免把残缺结果
    永久缓存;失败源下次分析自动重试,成功的源已缓存、秒回。

    并发写安全性:四路数据源的来源 CSV 文件名互不相同,各自单线程写自己的文件;
    merged 整合表只在主线程写一次;DisGeNET 限流记录文件内部有锁保护;
    共享 HTTP 缓存 sqlite 由 requests-cache 的线程本地连接保证安全。
    """
    from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor
    from concurrent.futures import wait as wait_futures

    import time

    import config

    drug_cache = os.path.join(config.MERGED_CACHE_DIR, f"merged_drug_{safe_name(drug_name)}.csv")
    disease_cache = os.path.join(
        config.MERGED_CACHE_DIR, f"merged_disease_{safe_name(disease_name)}.csv")

    drug_df = open_csv_cache(drug_cache, label=f"merged drug targets for '{drug_name}'")
    disease_df = open_csv_cache(disease_cache, label=f"merged disease targets for '{disease_name}'")
    if drug_df is not None and disease_df is not None:
        return drug_df, disease_df, []

    pool = ThreadPoolExecutor(max_workers=4, thread_name_prefix="target-fetch")
    futures: dict = {}
    if drug_df is None:
        futures["chembl"] = pool.submit(fetch_drug_targets_chembl, drug_name)
        futures["dgidb"] = pool.submit(fetch_drug_targets_dgidb, drug_name)
    if disease_df is None:
        futures["opentargets"] = pool.submit(get_disease_targets, disease_name)
        futures["disgenet"] = pool.submit(fetch_disease_targets_disgenet, disease_name)

    total = len(futures)
    results: dict = {}
    failed: list[str] = []
    done_keys: list[str] = []
    # 不等待仍在后台跑的线程:超时后让它们自己跑完并写缓存,供下次分析复用
    try:
        started = time.monotonic()
        pending = set(futures.values())
        key_of = {f: k for k, f in futures.items()}
        while pending:
            done, pending = wait_futures(pending, timeout=0.5, return_when=FIRST_COMPLETED)
            for fut in done:
                key = key_of[fut]
                try:
                    results[key] = fut.result()
                except Exception:
                    logging.exception("%s fetch thread raised; treating as failed", key)
                    failed.append(key)
                done_keys.append(key)
            # 单路超期:按失败处理,线程留在后台继续跑完并写来源缓存
            for fut in list(pending):
                key = key_of[fut]
                deadline = _SOURCE_DEADLINE_SECONDS[key]
                if time.monotonic() - started >= deadline:
                    pending.discard(fut)
                    failed.append(key)
                    done_keys.append(key)
                    logging.warning(
                        "%s fetch exceeded %ss deadline; continuing without it "
                        "(background thread keeps running and will cache the result)",
                        key, deadline)
            if progress_cb:
                progress_cb(list(done_keys), total)
    finally:
        pool.shutdown(wait=False)

    if drug_df is None:
        pair = results.get("chembl")
        chmbl_df = pair[0] if isinstance(pair, tuple) and pair[0] is not None else pd.DataFrame()
        # ChEMBL 请求中止时返回 (空表, None);查无此药返回 (空表, 路径)。前者算失败
        if "chembl" not in failed and not (isinstance(pair, tuple) and pair[1]):
            failed.append("chembl")
        dgidb_df = results.get("dgidb")
        if not isinstance(dgidb_df, pd.DataFrame):
            dgidb_df = pd.DataFrame()
        drug_df = _merge_source_pair(
            drug_name, chmbl_df, dgidb_df,
            cache_path=drug_cache, kind="drug",
            src_a="source_chmbl", src_b="source_dgidb",
            score_a=("score_pchembl", "pchembl_value"),
            score_b=("score_dgidb", "interaction_score"),
            sources_desc=("ChEMBL", "DGIdb"),
            save="chembl" not in failed and "dgidb" not in failed,
        )

    if disease_df is None:
        ot_df = results.get("opentargets")
        disgenet_df = results.get("disgenet")
        ot_df = ot_df if isinstance(ot_df, pd.DataFrame) else pd.DataFrame()
        disgenet_df = disgenet_df if isinstance(disgenet_df, pd.DataFrame) else pd.DataFrame()
        disease_df = _merge_source_pair(
            disease_name, ot_df, disgenet_df,
            cache_path=disease_cache, kind="disease",
            src_a="source_ot", src_b="source_disgenet",
            score_a=("score_ot", "association_score"),
            score_b=("score_disgenet", "disgenet_score"),
            sources_desc=("OpenTargets", "DisGeNET"),
            save="opentargets" not in failed and "disgenet" not in failed,
        )

    return (drug_df if drug_df is not None else pd.DataFrame(),
            disease_df if disease_df is not None else pd.DataFrame(),
            failed)
