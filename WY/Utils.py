"""Utils 工具入口。

具体功能(靶点获取、PPI、富集分析、分子对接等)都在 core/ 包中实现,
这里统一转出,页面代码仍按 import Utils + Utils.xxx 的方式调用。
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
