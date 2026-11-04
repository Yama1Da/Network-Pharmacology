"""API 层共享 helper(被多个路由模块复用)。"""

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
