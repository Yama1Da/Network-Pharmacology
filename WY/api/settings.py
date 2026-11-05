"""设置(App.py 齿轮职责的 TS 化):四栏侧边栏。"""

from fastapi import APIRouter
from fastapi import FastAPI, HTTPException, Request
from pydantic import BaseModel
from WY import config
from WY.core.management import (enrichment_display_name, save_api_env,
                                save_ui_settings, update_all_cache_files)
from .common import (_settings_snapshot)


router = APIRouter()


# ① 数据来源(抓取规模/STRING 版本) ② 分析算法(降噪阈值/黑名单)
# ③ 分子对接(口袋阈值/工具路径) ④ 通用设置(密钥/CLI 产物目录)

class ApiEnvBody(BaseModel):
    changes: dict




@router.get("/api/settings")
def settings_get() -> dict:
    return _settings_snapshot()


@router.post("/api/settings/api-env")
def settings_api_env(body: ApiEnvBody) -> dict:
    """密钥/工具路径写入:走 core.management.save_api_env(白名单+消毒+set_key+运行时回填)。"""
    return {"applied": save_api_env(body.changes)}


class StringVersionBody(BaseModel):
    version: str


@router.post("/api/settings/string-version")
def settings_string_version(body: StringVersionBody) -> dict:
    """STRING 版本切换:格式校验 → 下载站探活 → 激活(core.activate_string_version,
    版本常量/派生 URL/互作组运行时缓存一并处理)。"""
    from WY.core.interactome import string_version_exists
    from WY.core.management import activate_string_version

    target = config.normalized_string_version(body.version)
    if not target:
        raise HTTPException(status_code=400, detail=f"bad version format: {body.version!r}")
    if target == config.STRING_VERSION:
        return {"version": target, "changed": False}
    if not string_version_exists(target):
        raise HTTPException(status_code=400, detail=f"version {target} not found on download site")
    activate_string_version(target)
    return {"version": target, "changed": True}


class UiSettingsBody(BaseModel):
    changes: dict


# 放开全部运行时可调分析参数(core.UI_SETTINGS_KEYS 的 UI 子集;LOG_LEVEL/STRING_VERSION 走各自端点)。
# 值类型在此校验,防垃圾值经 setattr 污染运行时配置(save_ui_settings 只做持久化侧校验)。
_ALLOWED_UI_KEYS = (
    "POCKET_CUTOFF", "PPI_CONFIDENCE", "OT_TOP_N", "DISGENET_PAGES",
    "CHEMBL_PCHEMBL_THRESHOLD", "DGIDB_MIN_SCORE", "OT_MIN_SCORE",
    "CHEMBL_CORROBORATION_REQUIRED", "ENRICHMENT_PVALUE", "LOUVAIN_RESOLUTION",
    "UI_FONT_SCALE", "UI_THEME", "UI_LANG",
    # 统计稳健性/研究语义参数(缓存键均含这些值,改动后旧缓存自动失效重算)
    "RWR_ALPHA", "GSEA_PERMUTATION_NUM", "PROXIMITY_N_PERM", "RWR_NULL_PERMS",
    "GSEA_SET_MIN_SIZE", "GSEA_SET_MAX_SIZE", "GSEA_FDR_SCREEN",
    # 分子对接
    "EXHAUSTIVENESS", "DOCKING_CPU",
)
_BLACKLIST_KEY = "CHEMBL_TARGET_BLACKLIST"
_SOURCES_KEY = "ENRICHMENT_SOURCES"
# 数值范围(与前端 Stepper 的 min/max 一致,双保险)
_UI_RANGES: dict[str, tuple[float, float]] = {
    "POCKET_CUTOFF": (3, 6), "PPI_CONFIDENCE": (0.15, 0.9),
    "OT_TOP_N": (200, 1000), "DISGENET_PAGES": (1, 10),
    "CHEMBL_PCHEMBL_THRESHOLD": (0, 1), "DGIDB_MIN_SCORE": (0, 1),
    "OT_MIN_SCORE": (0, 1), "ENRICHMENT_PVALUE": (0.01, 0.2),
    "LOUVAIN_RESOLUTION": (0.5, 2),
    "UI_FONT_SCALE": (0.8, 1.5),
    # α→0 只剩种子、α→1 近乎不收敛,必须卡在退化区间之外
    "RWR_ALPHA": (0.1, 0.9), "GSEA_FDR_SCREEN": (0.05, 0.5),
    # 置换数/搜索深度:上限即耗时上限(满档单次分析十余分钟级),不开放更高档位
    "GSEA_PERMUTATION_NUM": (100, 5000), "PROXIMITY_N_PERM": (50, 2000),
    "RWR_NULL_PERMS": (50, 2000),
    "GSEA_SET_MIN_SIZE": (5, 100), "GSEA_SET_MAX_SIZE": (100, 1000),
    "EXHAUSTIVENESS": (1, 32), "DOCKING_CPU": (0, 64),
}
# 这些键影响降噪/抓取规模,值变化后旧参数算出的 merged 整合表不再可信 → 自动清缓存
_PURGE_MERGED_KEYS = {
    "CHEMBL_PCHEMBL_THRESHOLD", "DGIDB_MIN_SCORE", "OT_MIN_SCORE",
    "CHEMBL_CORROBORATION_REQUIRED", "CHEMBL_TARGET_BLACKLIST",
    "OT_TOP_N", "DISGENET_PAGES",
}


def settings_ui_save(changes: dict) -> dict:
    """应用 UI 设置的公共实现:类型/范围校验 → setattr → 持久化 → 值变化时同步派生值并清 merged。"""
    from WY.core.management import purge_merged_cache

    validated: dict = {}
    for k, v in changes.items():
        if k == _BLACKLIST_KEY:
            if not (isinstance(v, list) and all(isinstance(i, str) and i.strip() for i in v)):
                raise HTTPException(status_code=400, detail=f"{k} must be a list of gene symbols")
            validated[k] = [i.strip().upper() for i in v]
        elif k == _SOURCES_KEY:
            # 富集条目类型:取值表校验归一;空集会让 g:Profiler/Enrichr 请求非法,必须拒绝
            srcs = config.normalized_enrichment_sources(v)
            if not srcs:
                raise HTTPException(
                    status_code=400,
                    detail=f"{k} must be a non-empty subset of {config.ENRICHMENT_SOURCE_VALUES}")
            validated[k] = srcs
        elif k in _ALLOWED_UI_KEYS:
            if k in ("CHEMBL_CORROBORATION_REQUIRED",):
                if not isinstance(v, bool):
                    raise HTTPException(status_code=400, detail=f"{k} must be a boolean")
            elif k == "UI_THEME":
                if not (isinstance(v, str) and v.strip().lower() in ("light", "dark")):
                    raise HTTPException(status_code=400, detail=f"{k} must be 'light' or 'dark'")
                v = v.strip().lower()
            elif k == "UI_LANG":
                if not (isinstance(v, str) and v.strip().lower() in config.UI_LANG_VALUES):
                    raise HTTPException(status_code=400,
                                        detail=f"{k} must be one of {config.UI_LANG_VALUES}")
                v = v.strip().lower()
            else:
                if isinstance(v, bool) or not isinstance(v, (int, float)):
                    raise HTTPException(status_code=400, detail=f"{k} must be a number")
                lo, hi = _UI_RANGES[k]
                if not (lo <= v <= hi):
                    raise HTTPException(status_code=400,
                                        detail=f"{k} out of range [{lo}, {hi}]: {v}")
            validated[k] = v
        else:
            raise HTTPException(status_code=400,
                                detail=f"key not allowed: {k} (allowed: {_ALLOWED_UI_KEYS + (_BLACKLIST_KEY, _SOURCES_KEY)})")
    if not validated:
        raise HTTPException(status_code=400, detail="no changes")

    # GSEA 基因集规模窗成对校验:单改一边时与当前另一边比,夹逼区间不允许倒挂/相撞
    _lo = validated.get("GSEA_SET_MIN_SIZE", config.GSEA_SET_MIN_SIZE)
    _hi = validated.get("GSEA_SET_MAX_SIZE", config.GSEA_SET_MAX_SIZE)
    if _lo >= _hi:
        raise HTTPException(status_code=400,
                            detail=f"GSEA_SET_MIN_SIZE({_lo}) must be < GSEA_SET_MAX_SIZE({_hi})")

    def _differs(k: str, v) -> bool:
        old = getattr(config, k, None)
        if isinstance(old, tuple) and isinstance(v, list):
            return list(old) != v
        if isinstance(v, float) and isinstance(old, (int, float)) and not isinstance(old, bool):
            return abs(old - v) >= 1e-9
        return old != v

    changed = {k: v for k, v in validated.items() if _differs(k, v)}
    for k, v in validated.items():
        setattr(config, k, v)
    save_ui_settings(validated)
    if "PPI_CONFIDENCE" in changed:
        # 置信度被覆盖后重算背景互作组阈值,保证与 PPI 边筛选始终一致(config 末尾同款派生)
        config.INTERACTOME_CONFIDENCE = int(config.PPI_CONFIDENCE * 1000)
    purged = 0
    if any(k in _PURGE_MERGED_KEYS for k in changed):
        purged = purge_merged_cache()
    return {"snapshot": _settings_snapshot(), "changed": sorted(changed), "purged_merged": purged}


@router.post("/api/settings/ui")
def settings_ui(body: UiSettingsBody) -> dict:
    return settings_ui_save(body.changes)
