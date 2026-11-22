# Third-Party Notices

This project incorporates the following third-party open-source components.

## Python Dependencies
<!-- Regenerate with: pip-licenses --format=markdown --with-urls --order=name
     --packages <names from WY/requirements.txt>
     Caveat: pydeseq2 metadata embeds the full MIT text — flatten that cell to "MIT". -->

| Name               | Version  | License                                                 | Homepage                                         |
|--------------------|----------|---------------------------------------------------------|--------------------------------------------------|
| PubChemPy          | 1.0.5    | MIT                                                     | https://github.com/mcs07/PubChemPy               |
| anndata            | 0.13.4   | BSD-3-Clause                                            | https://github.com/scverse/anndata               |
| fastapi            | 0.162.0  | MIT                                                     | https://github.com/fastapi/fastapi               |
| gemmi              | 0.7.5    | Mozilla Public License 2.0 (MPL 2.0)                    | https://github.com/project-gemmi/gemmi           |
| gprofiler-official | 1.0.0    | BSD License                                             | https://biit.cs.ut.ee/gprofiler/                 |
| gseapy             | 1.3.1    | BSD-3-Clause                                            | https://github.com/zqfang/gseapy                 |
| matplotlib         | 3.11.2   | Python Software Foundation License                      | https://matplotlib.org                           |
| matplotlib-venn    | 1.1.2    | MIT License                                             | https://github.com/konstantint/matplotlib-venn   |
| meeko              | 0.8.0    | GNU Lesser General Public License v2 or later (LGPLv2+) | https://github.com/forlilab/meeko                |
| networkx           | 3.6.1    | BSD-3-Clause                                            | https://networkx.org/                            |
| numpy              | 2.5.1    | BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0      | https://numpy.org                                |
| pandas             | 2.3.3    | BSD License                                             | https://pandas.pydata.org                        |
| py3Dmol            | 2.5.5    | MIT License                                             | https://3dmol.org                                |
| pyarrow            | 25.0.1   | Apache-2.0                                              | https://arrow.apache.org/                        |
| pydantic           | 2.13.4   | MIT                                                     | https://github.com/pydantic/pydantic             |
| pydeseq2           | 0.5.4    | MIT                                                     | https://github.com/owkin/PyDESeq2                |
| python-dotenv      | 1.2.3    | BSD-3-Clause                                            | https://github.com/theskumar/python-dotenv       |
| pyvis              | 0.3.2    | BSD                                                     | https://github.com/WestHealth/pyvis              |
| rdkit              | 2026.3.5 | BSD-3-Clause                                            | https://github.com/kuelumbus/rdkit-pypi          |
| requests           | 2.34.2   | Apache Software License                                 | https://github.com/psf/requests                  |
| requests-cache     | 1.3.3    | BSD-2-Clause                                            | https://github.com/requests-cache/requests-cache |
| scikit-learn       | 1.9.1    | BSD-3-Clause                                            | https://scikit-learn.org                         |
| scipy              | 1.18.0   | BSD License                                             | https://scipy.org/                               |
| uvicorn            | 0.52.4   | BSD-3-Clause                                            | https://uvicorn.dev/                             |
| xgboost            | 3.4.1    | Apache-2.0                                              | https://github.com/dmlc/xgboost                  |

## External Executables (Not Bundled)

### Open Babel
- **License**: GPL-2.0
- **Usage**: Invoked as a separate executable via `subprocess`; not linked, embedded, or distributed with this project.
- **License Text**: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html

### AutoDock Vina
- **License**: Apache-2.0
- **Usage**: Invoked as a separate executable via `subprocess`; not bundled with this project.
- **License Text**: https://www.apache.org/licenses/LICENSE-2.0

## External Data Sources (Phase 2)

Second-phase modules consume the following external datasets. None of them is
bundled or redistributed with this repository; data is fetched at runtime (or
referenced as a versioned local snapshot) and provenance fields are preserved.

### CCLE (Cancer Cell Line Encyclopedia, via DepMap Portal)
- **License**: CC BY 4.0
- **Usage**: cell-line to cancer-type lineage mapping for the DrugComb synergy
  benchmark (phase-2 §2.3, decision B5: CCLE single-source primary matching).
  Unmapped combinations are NOT dropped — they are reported under the
  pan-cancer control module stratum.
- **Source**: https://depmap.org/portal/

### DrugComb
- **License**: academic use with citation requirement (terms recorded during the
  phase-2 G1 preflight per plan §2.3)
- **Usage**: gold-standard drug-combination labels (ZIP ≥ +10 / ≤ −10, majority
  vote across cell lines) for the synergy-recall benchmark.
- **Citation**: Liu H, et al. DrugComb: an integrative cancer drug combination
  data portal. Nucleic Acids Res. 2020;48(W1):W174-W181. The exact citation
  format follows the publisher's current requirement as verified at preflight.
- **Source**: https://drugcomb.fimm.fi/

### CancerMine
- **License**: CC0
- **Usage**: literature-mined cancer driver/oncogene/tumor-suppressor genes,
  unioned with Open Targets cancer targets (CC BY 4.0) to build the pan-cancer
  control gene set (phase-2 §1.9). A filtered, dated snapshot is kept at
  `WY/resources/cancermine_cancer_genes.csv` for traceability.
- **Source**: Zenodo concept DOI 10.5281/zenodo.1156241

> Note: NCG (Network of Cancer Genes), originally planned for the pan-cancer
> gene set, was NOT adopted — its site has been unavailable since 2026-11-18
> (latest release NCGHD 7.0, 2022). CancerMine (CC0) is used instead. Cellosaurus
> is likewise NOT used anywhere in this project (CC BY-NC is incompatible with
> the Apache-2.0 distribution stance); cell-line metadata relies on CCLE only.
