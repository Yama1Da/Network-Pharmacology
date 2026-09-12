"""ChEMBL 数据源：查询药物的全部人类活性靶点（直接 REST 批量分页），并提供缓存有效期调整。"""

import logging

import pandas as pd
import streamlit as st

from .caching import open_named_csv_cache, save_csv_cache
from .http_client import request_with_retry
from .text_utils import remember_name_mapping

# ChEMBL REST 基地址:加 .json 后缀直接返回 JSON,
# 用统一的 request_with_retry 直连,不再依赖 chembl_webresource_client
# (该客户端每页固定 20 条、逐靶点单查,慢药物要发数百个请求)。
CHEMBL_API_URL = "https://www.ebi.ac.uk/chembl/api/data"

# 单页记录数:activity 记录小,官方上限 1000/页;target 记录带组件与同义词,体积大,25/批
# (批次越大单次响应越慢,ChEMBL 服务端慢时容易读超时,实测 25/批更稳)
_ACTIVITY_PAGE_SIZE = 1000
_TARGET_PAGE_SIZE = 25

# 管理页设置的缓存有效期(小时);None 表示沿用共享 HTTP 缓存的全局 TTL
_chembl_cache_hours: float | None = None


def set_chembl_cache_ttl(hours: int | None) -> None:
    """调整 ChEMBL 请求结果的缓存有效期(默认与全局 HTTP 缓存一致,24 小时)。
    hours 为缓存小时数:设得越小,数据越新、联网请求越频繁;超过有效期的旧数据会自动重新联网获取。
    """
    global _chembl_cache_hours
    _chembl_cache_hours = hours


def _sort_drug_df_by_pchembl(df: pd.DataFrame) -> pd.DataFrame:
    """按 pChEMBL 活性值从高到低排列药物靶点（无活性值的排最后）。"""
    if df.empty or "pchembl_value" not in df.columns:
        return df
    df = df.copy()
    df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
    return (
        df.sort_values("pchembl_value", ascending=False, na_position="last")
        .reset_index(drop=True)
    )


def _abort_fetch() -> tuple[pd.DataFrame, None]:
    """统一失败出口:界面提示后返回空表(不写缓存,避免把失败当有效数据)。"""
    try:
        st.toast("靶点获取错误，停止信息获取", icon="⚠️")
    except Exception:
        # 并发线程里没有 Streamlit 脚本上下文,toast 不可用时退化为日志
        logging.warning("靶点获取错误，停止信息获取")
    return pd.DataFrame(), None


def _chembl_get(path: str, params: dict, label: str) -> dict | None:
    """向 ChEMBL REST 发一次带缓存的重试 GET;全部失败返回 None。
    ChEMBL 服务端偶发很慢,超时与重试次数都比轻量接口放宽,避免个别慢请求导致整个靶点抓取中止。
    """
    import config

    resp = request_with_retry(
        f"{CHEMBL_API_URL}/{path}", method="get", params=params,
        max_retries=config.CHEMBL_API_MAX_RETRIES, timeout=config.CHEMBL_API_TIMEOUT,
        retry_delay=config.API_RETRY_DELAY * 2,
        label=label, use_cache=True, cache_ttl_hours=_chembl_cache_hours,
        throttle_key="chembl",
    )
    return resp.json() if resp is not None else None


def _find_chembl_id(drug_name: str) -> str | None:
    """按药物名查找化合物编号:先精确匹配,找不到时改用不区分大小写的模糊搜索。"""
    name = drug_name.strip()
    payload = _chembl_get("molecule.json", {"pref_name__iexact": name, "limit": 1},
                          f"ChEMBL molecule search for '{name}'")
    mols = (payload or {}).get("molecules", [])
    if not mols:
        payload = _chembl_get("molecule.json", {"pref_name__icontains": name, "limit": 1},
                              f"ChEMBL fuzzy molecule search for '{name}'")
        mols = (payload or {}).get("molecules", [])
    if not mols:
        return None
    return mols[0].get("molecule_chembl_id")


def _fetch_all_activities(chembl_id: str) -> list[dict] | None:
    """取该药物在人类中的全部活性记录,按 1000 条/页翻页;失败返回 None 表示中止。"""
    acts: list[dict] = []
    offset = 0
    while True:
        payload = _chembl_get(
            "activity.json",
            {"molecule_chembl_id": chembl_id, "organism": "Homo sapiens",
             "limit": _ACTIVITY_PAGE_SIZE, "offset": offset},
            f"ChEMBL activity page (offset {offset})",
        )
        if payload is None:
            return None
        batch = payload.get("activities", [])
        acts.extend(batch)
        # 按服务端实际生效的页大小翻页(防止服务端上限低于请求值时提前停止)
        page_size = (payload.get("page_meta") or {}).get("limit") or _ACTIVITY_PAGE_SIZE
        if len(batch) < page_size:
            return acts
        offset += page_size


def _fetch_targets_batched(tids: list[str]) -> dict[str, dict] | None:
    """按编号批量取靶点详情(25 个/批,一次请求顶过去 25 次);失败返回 None 表示中止。
    每批一次请求即可:批内编号数就是返回条数上限,不需要翻页;实测对 __in 过滤请求越界的 offset 会得到服务端 500,故绝不发空页请求。
    """
    targets: dict[str, dict] = {}
    for start in range(0, len(tids), _TARGET_PAGE_SIZE):
        batch = tids[start:start + _TARGET_PAGE_SIZE]
        payload = _chembl_get(
            "target.json",
            {"target_chembl_id__in": ",".join(batch), "limit": _TARGET_PAGE_SIZE},
            f"ChEMBL target batch ({len(batch)} ids)",
        )
        if payload is None:
            return None
        for t in payload.get("targets", []):
            tid = t.get("target_chembl_id")
            if tid:
                targets[tid] = t
    return targets


def _extract_gene_name(t: dict) -> str | None:
    """从靶点信息中提取基因符号，没有则退而用 HGNC 名称，都没有时返回 None"""
    if "target_components" in t:
        for comp in t["target_components"]:
            for syn in comp.get("target_component_synonyms", []):
                if syn.get("syn_type") == "GENE_SYMBOL":
                    return syn.get("component_synonym")
            for xref in comp.get("target_component_xrefs", []):
                if xref.get("xref_src_db") == "HGNC":
                    return xref.get("xref_name")
    return None


def fetch_drug_targets_chembl(drug_name: str) -> tuple[pd.DataFrame, str | None]:
    """从 ChEMBL 获取药物的全部活性靶点（全靶点模式）：不按活性值过滤，靶点更全，但会包含 ALB、CYP450 等弱关联/非特异性结合的蛋白；仅保留人类靶点，按 pChEMBL 活性值降序排列。
    请求全部批量分页:activity 1000 条/页、靶点详情 25 个/批,慢药物(如 metformin 数千条活性)也只需约 10 个请求;已有本地缓存时直接复用，请求失败时返回空表。
    """
    import config

    csv_path, cached = open_named_csv_cache(
        config.DRUG_CACHE_DIR, "ChEMBL_", drug_name,
        label=f"ChEMBL targets for '{drug_name}'",
        remember=remember_name_mapping,  # 记录安全文件名→原始名,供管理页还原
    )
    if cached is not None:
        return _sort_drug_df_by_pchembl(cached), csv_path

    try:
        chembl_id = _find_chembl_id(drug_name)
        if not chembl_id:
            # 查无此药属正常业务情况:返回空表并带上路径(不写文件),
            # 调用方据此区分"没有数据"与"请求失败"
            logging.warning("No ChEMBL molecule found for '%s'", drug_name)
            return pd.DataFrame(), csv_path

        # 第二步：取该药物在人类中的全部活性记录（不按 pChEMBL 过滤，靶点更全）。
        # 如需只保留有 pChEMBL 定量值的高置信度靶点，可在参数中
        # 加回 pchembl_value__isnull=False。
        acts = _fetch_all_activities(chembl_id)
        if acts is None:
            return _abort_fetch()  # 拉取失败即停止，且不保存缓存

        # 保序去重(dict.fromkeys 比 循环+list 成员判断 快,慢药物有数千条活性记录)
        tids = list(dict.fromkeys(
            act.get("target_chembl_id") for act in acts if act.get("target_chembl_id")))

        targets = _fetch_targets_batched(tids) if tids else {}
        if targets is None:
            return _abort_fetch()  # 拉取失败即停止，且不保存缓存

        rows: list[dict] = []
        seen: set[tuple] = set()
        for act in acts:
            tid = act.get("target_chembl_id")
            if not tid:
                continue

            t = targets.get(tid)
            if not t:
                # 正常不会走到这里；若发生说明数据异常，停止并提示
                logging.error("Target details missing for '%s', aborting fetch", tid)
                return _abort_fetch()

            gene_name = _extract_gene_name(t)
            if not gene_name:
                continue

            pref_name = t.get("pref_name", "") or ""
            organism = t.get("organism", "") or ""
            if organism != "Homo sapiens":
                # 双重保险：只保留人类靶点
                continue
            assay_type = act.get("assay_type", "") or ""
            pchembl_value = act.get("pchembl_value")

            key = (gene_name, assay_type, pchembl_value)
            if key in seen:
                continue
            seen.add(key)

            rows.append({
                "approvedSymbol": gene_name,
                "target_id": tid,
                "approvedName": pref_name,
                "organism": organism,
                "assay_type": assay_type,
                "pchembl_value": pchembl_value,
            })

        df = pd.DataFrame(rows).reset_index(drop=True)
        # 不设数量上限，按 pChEMBL 活性值降序排列（无值的排最后）
        df = _sort_drug_df_by_pchembl(df)

        # 保存到本地缓存（空结果也落盘，避免反复查询同一查不到的药物）
        save_csv_cache(df, csv_path, skip_empty=False)

        return df, csv_path

    except Exception:
        logging.exception("fetch_drug_targets_chembl failed for '%s'", drug_name)
        return pd.DataFrame(), None
