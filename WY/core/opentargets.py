"""Open Targets 数据源：按疾病名称查询关联靶点。"""

import logging
import os
import threading

import pandas as pd

from .caching import load_json_file, open_named_csv_cache, save_csv_cache, save_json_file
from .http_client import _opentargets_post
from .text_utils import _safe_round, remember_name_mapping

# 进程内成功结果 memo:只记成功抓取的结果,网络失败不缓存(下次自动重试)。
# 原实现用 @st.cache_data,失败返回的空表也被记忆化,直到重启前同一疾病永不重试。
_success_memo: dict[tuple, pd.DataFrame] = {}
_memo_lock = threading.Lock()


def find_disease_id(disease_name: str,
                    prefer_prefixes: tuple[str, ...] = ("MONDO", "UMLS")) -> str | None:
    """按疾病英文名搜索,返回疾病编号(EFO/MONDO 等);找不到返回 None。
    同名疾病常有多个词表编号(如 MONDO 与 EFO),命中多个时按 prefer_prefixes 优先返回(下游 DisGeNET 只认 MONDO/UMLS 等编号,不认 EFO)。
    """
    search_query = """
    query SearchDisease($queryString: String!) {
        search(queryString: $queryString) {
            hits {
                id
                entity
                name
            }
        }
    }
    """
    payload = _opentargets_post(search_query, {"queryString": disease_name.strip()},
                                use_cache=True)
    if payload is None:
        return None

    # 查询出错时，错误信息记录在 payload["errors"] 中，不会抛异常
    if "errors" in payload or "data" not in payload:
        logging.error("Open Targets search query failed: %s", payload.get("errors"))
        return None

    disease_hits = [h for h in payload["data"]["search"]["hits"]
                    if h.get("entity") == "disease"]
    for prefix in prefer_prefixes:
        for hit in disease_hits:
            hit_id = hit.get("id") or ""
            if hit_id.startswith(prefix + "_"):
                return hit_id
    return disease_hits[0]["id"] if disease_hits else None


def _get_disease_targets_impl(disease_name: str,
                              need_rows: int | None = None) -> tuple[pd.DataFrame, str]:
    """通过 Open Targets 查询疾病的关联靶点(不做进程内记忆化),返回 (表, 状态)。
    状态与管理页语义一致:"ok"(有数据)/"empty"(查无此病或 0 靶点)/"failed"(请求失败)。
    先按疾病英文名找到疾病编号，再按关联分数从高到低固定抓取
    config.OT_FETCH_SIZE 条(截断到分析所需的 Top-N 在整合侧完成,见 core.targets)。need_rows 为本次分析
    至少需要的条数(即 config.OT_TOP_N),原始缓存行数不足时自动重新抓取,否则直接复用缓存。
    查询失败或找不到该疾病时返回空表；结果存为本地 CSV，已有缓存则直接复用。
    """
    import config

    csv_path, cached = open_named_csv_cache(
        config.DISEASE_CACHE_DIR, "OpenTargets_", disease_name,
        label=f"disease targets for '{disease_name}'",
        remember=remember_name_mapping,  # 记录安全文件名→原始名,供管理页还原
    )
    # 隐藏哨兵文件记录上次是否已取完全部关联靶点:疾病关联总数不足抓取上限时,缓存行数永远到不了
    # need_rows,若只按行数判断命中会导致每次分析都重新联网
    meta_path = os.path.join(os.path.dirname(csv_path), "." + os.path.basename(csv_path) + ".meta")
    meta = load_json_file(meta_path, {}) or {}
    if cached is not None and (need_rows is None or meta.get("complete")
                               or len(cached) >= need_rows):
        return cached, "ok"

    try:
        disease_id = find_disease_id(disease_name)
        if disease_id is None:
            logging.warning("No disease found for '%s'", disease_name)
            return pd.DataFrame(), "empty"

        # 第二步：用疾病编号查询关联靶点列表（已按关联分数从高到低排列）。
        # score 为总关联分数（0~1，越大关联越强），novelty 为新颖性，
        # datatypeScores 为各证据类型（遗传关联、体细胞突变、文献等）的分数。
        disease_query = """
        query DiseaseTargets($diseaseId: String!, $size: Int!) {
            disease(efoId: $diseaseId) {
                id
                name
                associatedTargets(page: {index: 0, size: $size}) {
                    rows {
                        score
                        novelty
                        datatypeScores {
                            id
                            score
                        }
                        target {
                            id
                            approvedSymbol
                            approvedName
                            biotype
                            genomicLocation {
                                chromosome
                                start
                                end
                            }
                            targetClass {
                                label
                            }
                        }
                    }
                }
            }
        }
        """
        payload = _opentargets_post(disease_query,
                                    {"diseaseId": disease_id, "size": config.OT_FETCH_SIZE},
                                    use_cache=True)
        if payload is None:
            return pd.DataFrame(), "failed"

        if "errors" in payload or "data" not in payload:
            logging.error("Open Targets disease query failed: %s", payload.get("errors"))
            return pd.DataFrame(), "failed"

        rows = payload["data"]["disease"]["associatedTargets"]["rows"]
        rows_sorted = sorted(rows, key=lambda r: r.get("score", 0.0),
                             reverse=True)[:config.OT_FETCH_SIZE]

        records: list[dict] = []
        for r in rows_sorted:
            tgt = r.get("target", {})
            gloc = tgt.get("genomicLocation", {}) or {}
            # 靶点类别可能有多个，只取第一个
            tc = tgt.get("targetClass") or []
            # 各证据类型分数整理成"类型名: 分数"的形式，便于逐项取值
            dt_scores = {s.get("id"): s.get("score") for s in r.get("datatypeScores", []) if s.get("id")}
            records.append({
                "approvedSymbol": tgt.get("approvedSymbol"),      # 官方基因符号
                "ensembl_id": tgt.get("id"),                       # ENSG 编号
                "approvedName": tgt.get("approvedName"),           # 基因全名
                "biotype": tgt.get("biotype"),                     # 生物类型
                "chromosome": gloc.get("chromosome"),              # 所在染色体
                "genomic_start": gloc.get("start"),                # 基因组起始位置
                "genomic_end": gloc.get("end"),                    # 基因组结束位置
                "targetClass": tc[0].get("label") if tc else None, # 靶点类别（取第一个）
                "association_score": _safe_round(r.get("score", 0.0)),           # 总关联分数
                "novelty": _safe_round(r.get("novelty", 0.0)),                   # 新颖性
                "genetic_association": _safe_round(dt_scores.get("genetic_association")),  # 遗传关联分数
                "somatic_mutation": _safe_round(dt_scores.get("somatic_mutation")),        # 体细胞突变分数
                "literature": _safe_round(dt_scores.get("literature")),                    # 文献分数
                "animal_model": _safe_round(dt_scores.get("animal_model")),                # 动物模型分数
                "affected_pathway": _safe_round(dt_scores.get("affected_pathway")),        # 通路分数
                "clinical": _safe_round(dt_scores.get("clinical")),                        # 临床证据分数
            })

        df = pd.DataFrame(records)
        # 保存到本地缓存（空结果也缓存，避免反复查询同一无效疾病名）
        save_csv_cache(df, csv_path, skip_empty=False)
        # 抓取行数不足单次上限 => 已拿全该疾病的关联靶点(总数不足),记入哨兵避免日后反复重抓
        save_json_file(meta_path, {"complete": len(df) < config.OT_FETCH_SIZE})
        logging.info("Disease '%s' (%s): %d targets fetched", disease_name, disease_id, len(df))
        return df, "ok"

    except Exception:
        # 网络错误等异常统一在此兜底：记录日志并返回空表，不让页面崩溃
        logging.exception("get_disease_targets failed for '%s'", disease_name)
        return pd.DataFrame(), "failed"


def get_disease_targets(disease_name: str, need_rows: int | None = None) -> pd.DataFrame:
    """进程内 memo 化入口:只有成功抓取的结果进入 memo,失败(网络/HTTP/解析)返回空表
    但不缓存,下一次同疾病分析自动重试;跨重启持久化仍由本地 CSV 缓存负责。
    并发安全:读写均在锁内。
    """
    key = (str(disease_name), None if need_rows is None else int(need_rows))
    with _memo_lock:
        if key in _success_memo:
            return _success_memo[key]
    df, status = _get_disease_targets_impl(disease_name, need_rows)
    if status == "ok":
        with _memo_lock:
            _success_memo[key] = df
    return df


def invalidate_disease_targets_cache(disease_name: str) -> None:
    """失效指定疾病名的进程内 memo,供管理页强制更新前调用。
    只清该名称的条目,不影响其他疾病(旧实现 .clear() 全量清空,会在批量更新时
    干扰其他会话/疾病正在使用的缓存)。
    """
    name = str(disease_name)
    with _memo_lock:
        for k in [k for k in _success_memo if k[0] == name]:
            del _success_memo[k]
