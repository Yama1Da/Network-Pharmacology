"""靶点汇总：整合药物靶点（ChEMBL+DGIdb）与疾病靶点（Open Targets+DisGeNET），并取交集。"""

import logging

import pandas as pd

from .caching import open_named_csv_cache, save_csv_cache
from .chembl import fetch_drug_targets_chembl
from .disgenet import fetch_disease_targets_disgenet
from .dgidb import fetch_drug_targets_dgidb
from .opentargets import get_disease_targets
from .text_utils import remember_name_mapping


def _max_score_map(df: pd.DataFrame, score_col: str) -> dict[str, float]:
    """按基因符号分组取分数最大值；符号统一大写，分数统一转数值（非数值置空）。"""
    if df.empty or score_col not in df.columns:
        return {}
    tmp = df.dropna(subset=["approvedSymbol"]).copy()
    tmp["approvedSymbol"] = tmp["approvedSymbol"].str.upper()
    tmp[score_col] = pd.to_numeric(tmp[score_col], errors="coerce")
    return tmp.groupby("approvedSymbol")[score_col].max().to_dict()


# ======== 数据源降噪过滤(阈值集中配在 config 的"靶点降噪过滤"段) ========

# Open Targets 只保留这些生物类型:交集与后续分子对接都只关心蛋白靶点
_OT_KEEP_BIOTYPES = {"protein_coding"}


def _filter_dgidb_targets(df: pd.DataFrame) -> pd.DataFrame:
    """DGIdb 降噪：保留交互分数达到下限、且有至少 1 条 PMID 文献支撑的记录;
    单独卡分数或单独要求 PMID 都挡不住低质量预测,两条件需同时满足。"""
    import config

    if df.empty:
        return df
    out = df
    if "interaction_score" in out.columns:
        score = pd.to_numeric(out["interaction_score"], errors="coerce")
        out = out[score >= config.DGIDB_MIN_SCORE]
    if "interaction_pmids" in out.columns:
        has_pmid = out["interaction_pmids"].fillna("").astype(str).str.strip() != ""
        out = out[has_pmid]
    return out


def _filter_ot_targets(df: pd.DataFrame) -> pd.DataFrame:
    """Open Targets 降噪：总关联分数达到下限,且只保留蛋白编码基因。"""
    import config

    if df.empty:
        return df
    out = df
    if "association_score" in out.columns:
        score = pd.to_numeric(out["association_score"], errors="coerce")
        out = out[score >= config.OT_MIN_SCORE]
    if "biotype" in out.columns:
        out = out[out["biotype"].fillna("").isin(_OT_KEEP_BIOTYPES)]
    return out


def _filter_chembl_targets(df: pd.DataFrame, corroborated_genes: set[str],
                           *, blacklist_only: bool = False) -> pd.DataFrame:
    """ChEMBL 降噪：先剔黑名单,再按 pChEMBL 分层保留——保留：基因的最大 pChEMBL 达到阈值;或全部活性记录都
    没有 pChEMBL 定量值、但被 DGIdb 佐证(corroborated_genes 为通过 DGIdb 过滤的大写基因集合);剔除：仅有弱定量
    活性(0 < pChEMBL < 阈值)的基因。blacklist_only=True 时只做黑名单过滤,供整侧过滤后为空时的兜底。"""
    import config

    if df.empty:
        return df
    out = df.dropna(subset=["approvedSymbol"])
    symbols = out["approvedSymbol"].str.upper()
    out = out[~symbols.isin({g.upper() for g in config.CHEMBL_TARGET_BLACKLIST})]
    if blacklist_only or "pchembl_value" not in out.columns:
        return out

    symbols = out["approvedSymbol"].str.upper()
    pchem = pd.to_numeric(out["pchembl_value"], errors="coerce")
    max_by_gene = pchem.groupby(symbols).max()
    keep = set(max_by_gene[max_by_gene >= config.CHEMBL_PCHEMBL_THRESHOLD].index)
    keep |= {g for g in max_by_gene[max_by_gene.isna()].index if g in corroborated_genes}
    return out[symbols.isin(keep)]


def _prepare_drug_sources(drug_name: str, chmbl_df: pd.DataFrame,
                          dgidb_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """对药物两个数据源做降噪过滤,返回过滤后的 (ChEMBL, DGIdb) 两表。DGIdb 先过滤(分数+PMID),通过的基因集
    再用于佐证 ChEMBL 的无定量活性基因;两源过滤后都为空但原始数据非空时(DGIdb 失败、全部活性只有弱定量值等),
    回退为 ChEMBL 仅黑名单过滤,避免整侧结果为空导致交集无意义。"""
    dgidb_f = _filter_dgidb_targets(dgidb_df)
    corroborated = (
        set(dgidb_f["approvedSymbol"].dropna().str.upper()) if not dgidb_f.empty else set()
    )
    chmbl_f = _filter_chembl_targets(chmbl_df, corroborated)
    if chmbl_f.empty and dgidb_f.empty and not (chmbl_df.empty and dgidb_df.empty):
        logging.warning(
            "Drug '%s': noise filters emptied both sources; "
            "falling back to blacklist-only ChEMBL filtering", drug_name)
        chmbl_f = _filter_chembl_targets(chmbl_df, set(), blacklist_only=True)
    logging.info(
        "Drug noise filter '%s': ChEMBL %d->%d rows, DGIdb %d->%d rows",
        drug_name, len(chmbl_df), len(chmbl_f), len(dgidb_df), len(dgidb_f))
    return chmbl_f, dgidb_f


def _prepare_disease_sources(disease_name: str, ot_df: pd.DataFrame,
                             disgenet_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """对疾病两个数据源做降噪过滤,返回过滤后的 (OpenTargets, DisGeNET) 两表;
    DisGeNET 结果本身按分数降序且只取第一页(相当于 top100),不再额外过滤。"""
    ot_f = _filter_ot_targets(ot_df)
    logging.info(
        "Disease noise filter '%s': OpenTargets %d->%d rows",
        disease_name, len(ot_df), len(ot_f))
    return ot_f, disgenet_df


# 两个类别的合并参数(来源标记列/分数列/日志描述),串行与并发共四处调用点共用
_DRUG_MERGE = dict(
    kind="drug", src_a="source_chmbl", src_b="source_dgidb",
    score_a=("score_pchembl", "pchembl_value"),
    score_b=("score_dgidb", "interaction_score"),
    sources_desc=("ChEMBL", "DGIdb"),
)
_DISEASE_MERGE = dict(
    kind="disease", src_a="source_ot", src_b="source_disgenet",
    score_a=("score_ot", "association_score"),
    score_b=("score_disgenet", "disgenet_score"),
    sources_desc=("OpenTargets", "DisGeNET"),
)


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
    """两个数据源的靶点按基因符号取并集去重，并保留各源分数最大值。kind 为 "drug"/"disease"，用于日志;
    src_a/src_b 为来源标记列名；score_a/score_b 为 (输出列名, 源分数列名)；两源都没有数据时返回空表。
    save=False 时不写缓存(某路数据源失败时,避免把残缺结果永久缓存)。"""
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
    """整合 ChEMBL 与 DGIdb 两个数据库的药物靶点，按基因符号取并集去重;整合前先做降噪过滤(ChEMBL 黑名单+
    pChEMBL 阈值与 DGIdb 佐证、DGIdb 分数+PMID 双门槛,见 _prepare_drug_sources)。为每个基因标记来自哪些数据库
    (source_chmbl / source_dgidb)并保留各库的活性/交互分数;两个库都没有数据时返回空表;结果存为本地 CSV,已有缓存则直接复用。"""
    import config

    # merged 缓存文件名带降噪过滤规则版本号(config.MERGED_DRUG_PREFIX),规则或阈值变更后旧缓存自动失效
    cache_path, cached = open_named_csv_cache(
        config.MERGED_CACHE_DIR, config.MERGED_DRUG_PREFIX, drug_name,
        label=f"merged drug targets for '{drug_name}'",
        remember=remember_name_mapping,  # 记录安全文件名→原始名,供管理页还原
    )
    if cached is not None:
        return cached

    chmbl_df, _ = fetch_drug_targets_chembl(drug_name)
    dgidb_df = fetch_drug_targets_dgidb(drug_name)
    chmbl_df, dgidb_df = _prepare_drug_sources(drug_name, chmbl_df, dgidb_df)
    return _merge_source_pair(drug_name, chmbl_df, dgidb_df,
                              cache_path=cache_path, **_DRUG_MERGE)


def get_disease_targets_merged(disease_name: str) -> pd.DataFrame:
    """整合 Open Targets 与 DisGeNET 两个数据库的疾病靶点，按基因符号取并集去重;整合前先对 Open Targets 做
    降噪过滤(分数下限+只留蛋白编码;DisGeNET 本身按分数取第一页,不过滤,见 _prepare_disease_sources)。
    为每个基因标记来自哪些数据库(source_ot / source_disgenet)并保留各库的关联分数;两个库都没有数据时返回空表;结果存为本地 CSV,已有缓存则直接复用。"""
    import config

    cache_path, cached = open_named_csv_cache(
        config.MERGED_CACHE_DIR, config.MERGED_DISEASE_PREFIX, disease_name,
        label=f"merged disease targets for '{disease_name}'",
        remember=remember_name_mapping,  # 记录安全文件名→原始名,供管理页还原
    )
    if cached is not None:
        return cached

    ot_df = get_disease_targets(disease_name)
    disgenet_df = fetch_disease_targets_disgenet(disease_name)
    ot_df, disgenet_df = _prepare_disease_sources(disease_name, ot_df, disgenet_df)
    return _merge_source_pair(disease_name, ot_df, disgenet_df,
                              cache_path=cache_path, **_DISEASE_MERGE)


def get_common_targets(drug_name: str, disease_name: str, top_n: int = 1000) -> list[str]:
    """取药物靶点与疾病靶点的交集基因，作为后续 PPI 网络分析的输入;按字母序排序,最多返回 top_n 个
    (默认 1000,对齐 STRING 查询上限);任一侧无数据时返回空列表。"""
    drug_df, disease_df, _failed = get_targets_concurrent(drug_name, disease_name)
    drug_set = set(drug_df["approvedSymbol"].dropna()) if not drug_df.empty else set()
    disease_set = set(disease_df["approvedSymbol"].dropna()) if not disease_df.empty else set()
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
    """并发获取药物与疾病两侧靶点(四路数据源同时请求),供分析页一次拿全;语义与分别串行调用 get_drug_targets_merged + get_disease_targets_merged 一致:两侧 merged 缓存都命中时直接返回,否则只对缺失的一侧并发联网。progress_cb: 可选进度回调,每 0.5 秒或某路完成时被调用一次(始终在主线程,可安全更新 Streamlit 组件),签名 progress_cb(done: list[str], total: int),done 为已完成的数据源 key(chembl/dgidb/opentargets/disgenet)。
    返回 (药物靶点表, 疾病靶点表, 失败数据源 key 列表);失败的数据源按空表参与展示,且该侧不写 merged 缓存——避免把残缺结果永久缓存;失败源下次分析自动重试,成功的源已缓存、秒回。
    并发写安全性:四路数据源的来源 CSV 文件名互不相同,各自单线程写自己的文件;merged 整合表只在主线程写一次;DisGeNET 限流记录文件内部有锁保护;共享 HTTP 缓存 sqlite 由 requests-cache 的线程本地连接保证安全。"""
    from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor
    from concurrent.futures import wait as wait_futures

    import time

    import config

    drug_cache, drug_df = open_named_csv_cache(
        config.MERGED_CACHE_DIR, config.MERGED_DRUG_PREFIX, drug_name,
        label=f"merged drug targets for '{drug_name}'")
    disease_cache, disease_df = open_named_csv_cache(
        config.MERGED_CACHE_DIR, config.MERGED_DISEASE_PREFIX, disease_name,
        label=f"merged disease targets for '{disease_name}'")
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
        chmbl_df, dgidb_df = _prepare_drug_sources(drug_name, chmbl_df, dgidb_df)
        drug_df = _merge_source_pair(
            drug_name, chmbl_df, dgidb_df,
            cache_path=drug_cache, save="chembl" not in failed and "dgidb" not in failed,
            **_DRUG_MERGE,
        )

    if disease_df is None:
        ot_df = results.get("opentargets")
        disgenet_df = results.get("disgenet")
        ot_df = ot_df if isinstance(ot_df, pd.DataFrame) else pd.DataFrame()
        disgenet_df = disgenet_df if isinstance(disgenet_df, pd.DataFrame) else pd.DataFrame()
        ot_df, disgenet_df = _prepare_disease_sources(disease_name, ot_df, disgenet_df)
        disease_df = _merge_source_pair(
            disease_name, ot_df, disgenet_df,
            cache_path=disease_cache,
            save="opentargets" not in failed and "disgenet" not in failed,
            **_DISEASE_MERGE,
        )

    return (drug_df if drug_df is not None else pd.DataFrame(),
            disease_df if disease_df is not None else pd.DataFrame(),
            failed)
