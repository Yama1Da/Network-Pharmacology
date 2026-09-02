"""Utils 工具入口。

具体功能(靶点获取、PPI、富集分析、分子对接等)都在 core/ 包中实现,
这里统一转出,页面代码仍按 import Utils + Utils.xxx 的方式调用。
"""

from core.chembl import fetch_drug_targets_chembl, set_chembl_cache_ttl
from core.dgidb import fetch_drug_targets_dgidb
from core.disgenet import fetch_disease_targets_disgenet
from core.docking import (
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
from core.enrichment import render_enrichment_plots, run_enrichment
from core.files import (
    list_cache_files,
    list_download_dirs,
    list_log_files,
    make_dir_zip,
    make_ppi_zip,
)
from core.logging_setup import setup_logging
from core.opentargets import get_disease_targets
from core.ppi import (
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
from core.text_utils import safe_name
