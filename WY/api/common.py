"""API 层共享 helper(被多个路由模块复用)。"""

import os

from fastapi import Response

from WY import config


_ALLOWED_LEVELS = ("INFO", "DEBUG")


_LOG_TAIL_MAX = 2 * 1024 * 1024  # 单次最多回传 2MB 尾部,防止历史大日志拖垮响应


_ZIP_NAMES = {
    "drug_sources": "Drug_sources.zip",
    "disease_sources": "Disease_sources.zip",
    "ppi": "PPI_files.zip",
    "enrichment": "Enrichment_files.zip",
    "stars": "Starred_results.zip",
    "downloads": "Data_files.zip",
}


def _resolve_under(base: str, rel: str, *, allow_root: bool = False) -> str | None:
    """realpath 包含校验(2026-11 P1-A 口径的公共件):把 rel 拼到 base 下并解析
    符号链接/../穿越,逃逸出 base 时返回 None。allow_root 决定恰好等于 base 本身
    (请求根目录)是否放行;Windows 盘符绝对路径会让 os.path.join 整段丢弃 base
    前缀,realpath 后同样逃逸,一并拦下。"""
    base_real = os.path.realpath(base)
    full = os.path.realpath(os.path.join(base_real, rel))
    contained = full.startswith(base_real + os.sep) or (allow_root and full == base_real)
    return full if contained else None


def _rel_downloads(path: str) -> str | None:
    """downloads 绝对路径 → 相对路径(统一正斜杠,前端拼 /api/artifacts?rel=…);
    空路径/越界/解析失败返回 None。"""
    if not path:
        return None
    try:
        real = os.path.realpath(path)
        base = os.path.realpath(config.DOWNLOADS_DIR)
        if real.startswith(base + os.sep):
            return os.path.relpath(real, base).replace(os.sep, "/")
    except OSError:
        pass
    return None


def _attachment(fname: str, disposition: str = "attachment") -> str:
    """Content-Disposition 头(自 retrieve.py 升为公共件)。文件名含非 ASCII(中文药名等)
    时按 RFC 5987 补 filename*:FastAPI 响应头按 latin-1 编码,直接放中文会 500。
    disposition 支持 attachment(下载)/ inline(内嵌展示)。"""
    if fname.isascii():
        return f'{disposition}; filename="{fname}"'
    from urllib.parse import quote

    fallback = fname.encode("ascii", "replace").decode().replace('"', "_")
    return f'{disposition}; filename="{fallback}"; filename*=UTF-8\'\'{quote(fname)}'


def _zip_response(data: bytes, filename: str) -> Response:
    """字节 zip 的统一下载响应:attachment 头走 _attachment(非 ASCII 文件名安全)。"""
    return Response(content=data, media_type="application/zip",
                    headers={"Content-Disposition": _attachment(filename)})


def _df_rows(df, limit: int | None = None) -> list[dict]:
    if df is None or len(df) == 0:
        return []
    if limit is not None:
        df = df.head(limit)
    return df.to_dict(orient="records")


def _settings_snapshot() -> dict:
    """设置快照:密钥只回布尔(不回值),路径明文,分析参数回当前生效值。"""
    from WY.core.docking import docking_timeout

    return {
        "keys": {"disgenet": bool(config.DisGeNET_API_KEY),
                 "ncbi": bool(config.NCBI_API_KEY)},
        "vina_dir": config.VINA_DIR,
        "obabel_dir": config.OBABEL_DIR,
        "cli_results_dir": config.CLI_RESULTS_DIR,
        "string_version": config.STRING_VERSION,
        "known_versions": list(config.STRING_KNOWN_VERSIONS),
        "ot_top_n": config.OT_TOP_N,
        "disgenet_pages": config.DISGENET_PAGES,
        "chembl_blacklist": list(config.CHEMBL_TARGET_BLACKLIST),
        "chembl_blacklist_options": list(config.CHEMBL_BLACKLIST_OPTIONS),
        "chembl_pchembl_threshold": config.CHEMBL_PCHEMBL_THRESHOLD,
        "dgidb_min_score": config.DGIDB_MIN_SCORE,
        "ot_min_score": config.OT_MIN_SCORE,
        "chembl_corroboration_required": config.CHEMBL_CORROBORATION_REQUIRED,
        "enrichment_pvalue": config.ENRICHMENT_PVALUE,
        "louvain_resolution": config.LOUVAIN_RESOLUTION,
        "ppi_confidence": config.PPI_CONFIDENCE,
        "pocket_cutoff": config.POCKET_CUTOFF,
        # 统计稳健性/研究语义参数(面板可调,随 ui_settings.json 迁移)
        "rwr_alpha": config.RWR_ALPHA,
        "gsea_permutation_num": config.GSEA_PERMUTATION_NUM,
        "proximity_n_perm": config.PROXIMITY_N_PERM,
        "rwr_null_perms": config.RWR_NULL_PERMS,
        "gsea_set_min_size": config.GSEA_SET_MIN_SIZE,
        "gsea_set_max_size": config.GSEA_SET_MAX_SIZE,
        "gsea_fdr_screen": config.GSEA_FDR_SCREEN,
        "enrichment_sources": list(config.ENRICHMENT_SOURCES),
        "exhaustiveness": config.EXHAUSTIVENESS,
        "docking_cpu": config.DOCKING_CPU,
        "num_modes": config.NUM_MODES,
        # 对接页参数行展示用:随搜索深度放宽的有效超时
        "docking_timeout": docking_timeout(),
        "ui_font_scale": config.UI_FONT_SCALE,
        "ui_theme": config.UI_THEME,
        "ui_lang": config.UI_LANG,
    }
