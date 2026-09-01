"""Utils 兼容门面（facade）。

原单文件 Utils.py（2700+ 行）已按业务域拆分到 core/ 包：

    core/http_client.py     远端 API 公共层（重试 POST、端点常量）
    core/text_utils.py      safe_name / _safe_round 等小工具
    core/logging_setup.py   日志系统（控制台 + 按天轮转）
    core/caching.py         CSV 缓存读写与 TTL 判断
    core/opentargets.py     Open Targets（疾病/药物靶点）
    core/disgenet.py        DisGeNET（疾病靶点 + 主动限流）
    core/chembl.py          ChEMBL（药物靶点 + 缓存 TTL）
    core/dgidb.py           DGIdb（药物-基因交互）
    core/targets.py         多来源靶点整合（merged / 交集）
    core/files.py           缓存/下载/日志列表与 zip 打包
    core/ppi.py             STRING PPI：抓取、建图、指标、渲染、落盘
    core/enrichment.py      富集分析：g:Profiler + Enrichr 降级、绘图
    core/docking.py         分子对接：Vina 全流程与结果管理

本文件仅做 re-export，保持 `import Utils` + `Utils.xxx` 的既有调用方式，
App.py 与 pages/ 下所有页面无需任何改动。新代码建议直接从 core 子模块导入。
"""

from core.caching import _is_cache_fresh, _read_cache_csv
from core.chembl import (
    _sort_drug_df_by_pchembl,
    fetch_drug_targets_chembl,
    set_chembl_cache_ttl,
)
from core.dgidb import fetch_drug_targets_dgidb
from core.disgenet import (
    _acquire_disgenet_rate_limit_slot,
    fetch_disease_targets_disgenet,
)
from core.docking import (
    _docking_fingerprint,
    _pdb_coords,
    _split_models_python,
    _strip_receptor_pdb,
    _unique_path,
    clean_docking_cache,
    fetch_ligand_sdf_3d,
    fetch_receptor_pdb_candidates,
    infer_box_from_cofactor,
    parse_vina_scores,
    prepare_ligand_pdbqt,
    prepare_receptor_pdbqt,
    prune_docking_results,
    run_vina,
    split_vina,
    star_docking_result,
    write_docking_metadata,
)
from core.enrichment import (
    ENRICHR_LIB_MAP,
    _enrichment_cache_key,
    _normalize_enrichment_enrichr,
    _normalize_enrichment_gprofiler,
    fetch_enrichment_enrichr,
    fetch_enrichment_gprofiler,
    render_enrichment_plots,
    run_enrichment,
)
from core.files import (
    list_cache_files,
    list_download_dirs,
    list_log_files,
    make_dir_zip,
    make_ppi_zip,
)
from core.http_client import OPENTARGETS_API_URL, _opentargets_post
from core.logging_setup import _DailyRotatingFileHandler, setup_logging
from core.opentargets import (
    _search_drug,
    fetch_drug_targets_opentargets,
    get_disease_targets,
)
from core.ppi import (
    _string_post,
    build_ppi_graph,
    compute_ppi_metrics,
    fetch_string_ppi,
    mark_hub_genes,
    render_pyvis_html,
    render_static_plot,
    save_ppi_results,
)
from core.targets import (
    get_common_targets,
    get_disease_targets_merged,
    get_drug_targets_merged,
)
from core.text_utils import _safe_round, safe_name
