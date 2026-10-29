[简体中文](README.md) | EN

# Network Pharmacology Analysis Platform (NPAP)

[![CI](https://github.com/Yama1Da/Network-Pharmacology/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/Yama1Da/Network-Pharmacology/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.13%2C%20%3C3.14-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

A one-stop network pharmacology platform: enter a drug and disease name and it
automatically runs the full pipeline of
"target acquisition → target integration → PPI network construction → enrichment analysis → molecular docking",
and can combine GEO expression profiles with multi-algorithm machine-learning consensus
to narrow the target list worth wet-lab validation. The project mainly serves network
pharmacology researchers and undergraduates. The web UI visualizes the whole workflow;
the CLI supports batch tabulation, but does not produce figures.

## Features

- **Disease/drug target acquisition**: integrates four data sources — DisGeNET, DGIdb,
  Open Targets, and ChEMBL — with automatic deduplication, and filters out
  non-specific binding proteins (albumin, CYP metabolic enzymes, etc.) to reduce noise;
  the filter conditions are user-adjustable.
- **PPI network**: builds protein–protein interaction networks based on STRING
  (v12.0 by default; the data version is switchable in Settings), with topology analysis
  identifying MCC core targets and multiple visualization styles.
- **Target prioritization**: caches the full STRING network locally as a background
  interactome, supporting proximity Z-score and network-wide RWR ranking; one-click build
  on the Data Management page, with mirror/official dual-source fallback and segmented
  concurrent downloads.
- **Enrichment analysis**: GO / KEGG pathway enrichment; automatically falls back to
  Enrichr when g:Profiler fails. Both sources have circuit-breaker degradation and
  automatic recovery, and when both are unavailable the page shows a
  "automatic retry in about X minutes" notice.
- **Molecular docking**: integrates AutoDock Vina, with RDKit + Meeko for automatic
  ligand preparation, built-in 3D result viewing, and focused-pocket display
  (surroundings dimmed).
- **Expression profiles and ML consensus screening**: searches and downloads GEO disease
  expression datasets, automatically performs differential expression analysis (DEG),
  then cross-screens with six algorithms — LASSO / Elastic Net / Ridge / SVM-RFE /
  Random Forest + Boruta / XGBoost — with consensus scoring + CV-ROC evaluation, and
  intersects the result with the platform-ranked targets to output a wet-lab priority
  shortlist.
- **Engineering**: request-level caching, local DisGeNET rate limiting (within the
  academic quota), enrichment dual-source circuit breaking, bilingual Chinese/English UI,
  daily rolling logs, and a target data management page.

## Requirements

Python dependencies need no manual installation: the launcher scripts (double-click
`start_npap.bat` on Windows, run `./start_npap.sh` on Linux/macOS) create the virtual
environment and install them on first run. Minimum requirements:

- Python 3.13 (exact version range in the Python badge above)
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina/releases) and
  [Open Babel](https://openbabel.org/): only the "Molecular Docking" feature needs them.
  They are discovered automatically once on PATH, or can be specified via
  `VINA_DIR` / `OBABEL_DIR` on the web Settings page or in `WY/API.env`
- A [DisGeNET](https://www.disgenet.org/) API key for DisGeNET, a disease–gene association
  database and one of the disease target sources; academic accounts can apply for free.
  Fill it in on the web Settings page and it is saved to `WY/API.env` automatically.
  The project also works without it — only the DisGeNET target table stays empty, all
  other data sources are unaffected

## Quick Start

Steps 1 and 2 are required; once they are done the main workflow runs in the browser.
Steps 3 and 4 are optional.

### 1. Install Python 3.13

On Windows, install from [python.org](https://www.python.org/downloads/) and tick
"Add python.exe to PATH"; on macOS / Linux:

```bash
brew install python@3.13                       # macOS
sudo apt install python3.13 python3.13-venv    # Debian/Ubuntu
sudo dnf install python3.13                    # Fedora
```

### 2. One-click launch

On Windows, double-click `start_npap.bat` (to stop, double-click `stop_npap.bat`);
on Linux/macOS run:

```bash
./start_npap.sh    # stop: ./stop_npap.sh
```

The first run creates the virtual environment and installs dependencies automatically
(a few minutes; the Tsinghua PyPI mirror is used by default — clear the `NPAP_PIP_INDEX`
environment variable to switch to the official source). When done it opens
`http://127.0.0.1:8000` in your browser automatically; open that address manually if
it does not.

### 3. Install the molecular docking tools (only the "Molecular Docking" feature needs them)

On Windows, download the installers from the links in "Requirements";
on Linux/macOS:

```bash
conda install -c conda-forge vina openbabel
```

### 4. (Optional) command-line batch analysis

The launcher has already installed the project in editable mode, so once the virtual
environment is activated, `npap` is available everywhere:

```bash
source .venv/bin/activate   # Windows: .venv\Scripts\activate
npap --help                 # see "Command-Line Batch Analysis (CLI)" below
```

## Usage Workflow

1. **Start an analysis on the home page**: enter a drug name and a disease name
   (the default example is metformin and Type 2 diabetes mellitus; only English disease
   names are supported for now), then click "Analyze". The program queries the four
   data sources concurrently.
2. **Review the target overview**: once the analysis finishes, the page shows the number
   of drug / disease / shared targets along with a Venn diagram. Expand the collapsible
   panels to view per-source target detail tables (sorted by source count and score).
3. **Data management**: open it from the "Data Management" menu at the top-left to view,
   add, remove, or correct the fetched disease and drug targets. The same page builds
   the background interactome in one click (a prerequisite for proximity / RWR analyses).
4. **Data network**: build a STRING interaction network from the integrated targets
   (multiple visualization styles available) and run GO / KEGG enrichment analysis.
5. **Molecular docking**: select core targets for AutoDock Vina docking and inspect the
   3D binding poses and focused-pocket view online.
6. **Expression profiles and ML consensus screening** (optional): search and select a
   GEO dataset (or upload a local expression matrix), automatically run DEG and
   six-algorithm consensus screening, and download the shortlist as a wet-lab priority
   reference.
7. **Logs and settings**: switch between Chinese/English and view run logs at the
   top-right; the gear button at the bottom-left adjusts analysis parameters (e.g.,
   target filter thresholds) and the STRING data version — when STRING releases a new
   version, just change the version number in Settings and rebuild the background
   network, with no program update needed.

All intermediate results are cached automatically (see `WY/config.py`): target pulls from
the four databases share a unified 60-day cache (upstream release cycles range from about
a quarter to two and a half years, so 60 days never misses a new release while keeping
quota consumption to a minimum), and integrated tables, PPI, enrichment and other analysis
caches default to 7 days. For an immediate refresh, click "Update" on a management page
to force a fresh download.

## Command-Line Batch Analysis (CLI)

`WY/cli.py` is a thin wrapper over `core/` and shares all caches and settings with the
web UI: analyses the web has already run return instantly in the CLI, and thresholds
changed in the web UI apply to the CLI automatically. The CLI targets batch tabulation;
it does not plot and does not dock — paper figures come from the web UI.

```bash
cd WY
source .venv/bin/activate   # Windows: .venv\Scripts\activate (the launcher created the env)

# Single-step commands: each writes CSV/JSON files into --out (default ./results;
# set CLI_RESULTS_DIR in API.env once for a permanent output directory)
python cli.py targets   --drug aspirin --disease "colorectal cancer"   # targets from 4 sources
python cli.py ppi       --genes genes.txt                              # STRING network + MCC core targets
python cli.py enrich    --genes common_targets.csv [--pvalue 0.05] [--sources GO:BP,KEGG]
python cli.py proximity --drug aspirin --disease "colorectal cancer"   # proximity Z-score
python cli.py rwr       --seeds core_targets.csv                       # network-wide RWR ranking

# Full pipeline for one pair: one result subdirectory + summary.json (parameter snapshot
# and per-step status, fully traceable)
python cli.py analyze --drug aspirin --disease "colorectal cancer"

# Batch: pairs.csv must contain 'drug' and 'disease' columns, one pair per row;
# one subdirectory per pair + batch_summary.csv
python cli.py batch --pairs pairs.csv --out results
python cli.py batch --pairs pairs.csv --out results --resume   # checkpoint resume
```

- **Exit codes**: 0 when everything succeeds, 1 when anything fails; failure reasons are
  printed directly (results are never silently incomplete). A single failed pair never
  aborts the batch, and `--resume` skips pairs that already succeeded and retries failed
  ones.
- **Gene lists** are interchangeable across `--genes` (`ppi`/`enrich`) and `--seeds`
  (`rwr`), in any of three forms: a .txt with one gene per line (`#` comments allowed),
  a .csv (auto-detects the `gene`/`approvedSymbol` column), or a comma-separated literal
  list. `targets`/`proximity`/`analyze`/`batch` take drug/disease names only; the two
  groups connect through the produced files (e.g.
  `enrich --genes common_targets.csv`).
- proximity / RWR require the background interactome: when it is missing they are skipped
  with a notice and never trigger the large download automatically.
- English disease names only; run batches sequentially — do not start multiple batch
  processes at once.
- `python cli.py serve` starts the web server directly (and opens the browser).

**Install as a global `npap` command (optional)**:

```bash
cd WY
pip install -e ..        # -e (editable) is required: caches/logs/API.env are all resolved
                         # relative to the WY/ source directory; a regular install would
                         # write them into site-packages
npap --help              # the entry script lands in the env's Scripts\ (Windows) or bin/
```

After installing, `npap` is available in any directory while the environment is active
(identical on all three platforms); to use it without activation, add the environment's
`Scripts\` (Windows) or `bin/` (Linux/macOS) directory to PATH. For global use, set
`CLI_RESULTS_DIR` in `API.env` once; otherwise, without `--out`, results land in the
`results` folder under the current directory.

## Configuration

| File           | Purpose                                                                                                                                               |
|----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| `WY/API.env`   | Secrets and local machine paths (`DisGeNET_API_KEY`, `VINA_DIR`, `OBABEL_DIR`, `CLI_RESULTS_DIR`); not committed to git                                  |
| `WY/config.py` | Data source URLs, rate limits/retries, circuit-breaker cooldown ladders, cache TTL, target filter thresholds, STRING version and download concurrency |

For behavior changes see `WY/config.py`; for secrets see `WY/API.env`.

## FAQ

- **DisGeNET returns 429?** The app already applies local rate limiting within the academic
  quota (50 requests/min); if it still triggers, just wait a moment. The quota is
  25,000 requests/day, and a normal single analysis consumes only 1 request. The project
  keeps a 10% margin, so the effective limits are 45 requests/min and 20,000 requests/day.
  If your account has a different tier's quota, adjust `DISGENET_RATE_LIMIT_PER_MINUTE`,
  `_PER_DAY`, and `_PER_MONTH` in `WY/config.py` yourself.
- **No enrichment results?** Repeated g:Profiler failures trip a circuit breaker and switch
  to Enrichr; when both sources are down the page shows an "automatic retry in about
  X minutes" notice. The cooldown ladders are configured in
  `GPROFILER_CIRCUIT_COOLDOWN_LADDER` / `ENRICHR_CIRCUIT_COOLDOWN_LADDER`.
- **Changed the filter thresholds but results look the same?** Increment
  `TARGET_FILTER_VERSION` accordingly; old integrated caches then expire automatically
  and are rebuilt.
- **Background interactome downloads slowly?** The CN mirror usually carries only the
  latest release; other versions automatically fall back to the official site. Raise
  `INTERACTOME_DOWNLOAD_CONCURRENCY` (default 6) in `WY/config.py` for more parallel
  segments (the official site throttles per connection), or use a version the mirror
  already carries.
- **Molecular docking can't find Vina / obabel?** Make sure the executable is on PATH, or
  point `VINA_DIR` / `OBABEL_DIR` in `WY/API.env` to the corresponding directory (a full
  path to the executable also works).
- **CLI results not where you expected?** The output location priority is
  `--out` > `CLI_RESULTS_DIR` (in `WY/API.env` or a system environment variable of the
  same name) > the `results` folder in the current directory; `--out` still overrides for
  a single run.

## Project Structure

```
WY/
├── api.py                  # FastAPI entry: /api endpoints + hosts the web/dist build output
├── cli.py                  # CLI entry: serve / targets / ppi / enrich / proximity / rwr / analyze / batch
├── config.py               # All tunable parameters (paths, rate limits/retries, cache TTL, filter thresholds, STRING version)
├── env_discovery.py        # Vina / Open Babel executable discovery
├── API.env                 # Secrets and local paths (not committed; auto-created on first save in Settings)
├── core/                   # Domain and orchestration layer, no UI framework dependencies
│   ├── targets / sources / chembl / dgidb / disgenet / opentargets
│   │                       #   four-source fetching and integrated noise reduction
│   ├── ppi / enrichment    #   STRING network, GO/KEGG enrichment (dual-source circuit breaking)
│   ├── interactome / network_algorithms    #   background interactome, proximity / RWR
│   ├── docking / docking_session           #   AutoDock Vina docking
│   ├── geo / expression / ml_screening / ml_pipeline / ml_results
│   │                       #   GEO expression profiles + ML consensus screening
│   ├── management / jobs / analysis / files  #   data-management orchestration, background jobs, result files
│   └── http_client / caching / logging_setup / text_utils  # infrastructure
├── web/                    # TS + Vue 3 frontend source (Vite build; dist/ hosted by api.py)
├── i18n/                   # Chinese/English dictionaries and lookup API
├── tools/                  # Standalone tools: build_interactome / geo_ml_screen / net_block / export_i18n_json
├── benchmark/              # Benchmarks: LIT-PCBA / cross-dataset robustness / ML scale / recovery
└── Static/                 # Runtime caches and analysis results
```

Dependencies are strictly one-way: `web(TS+Vue) / api.py → core/ → i18n, config`.
`core/` does not depend on any UI framework and can be imported, reused, and tested
headlessly; `cli.py` depends only on `core/` and `config`, sharing caches and settings
with the web UI.

## Development & Tests

```bash
pytest                            # full test suite, purely local, all external requests stubbed
pytest --cov --cov-report=html    # coverage report (writes htmlcov/; open it for a line-by-line view)
pytest -p WY.tools.net_block      # network-blocking audit: blocks non-loopback connections, proving zero external dependency
```

- Coverage policy lives in `pyproject.toml`: measures the `WY` package (excluding `tools/`,
  `benchmark/`, `web/`, `Static/`); **the suite fails below 75% coverage** (identical
  policy locally and in CI).
- ML-related tests run a slim configuration by default; set `NPAP_ML_TEST_FULL=1` for the
  full parameter set.
- CI (`.github/workflows/ci.yml`): every push / PR installs dependencies from scratch and
  runs the full suite; a full-parameter scheduled run also executes every Monday to catch
  statistical behavior drift from dependency upgrades.
- External data sources (DGIdb / DisGeNET / ChEMBL / Open Targets / STRING / g:Profiler /
  Enrichr / GEO) are all stubbed, and writes are sandboxed by `tests/conftest.py`: no API
  quota is consumed and the real `Static/` caches are untouched.
- Background interactome build: prefer the one-click build on the Data Management page,
  or run `WY/tools/build_interactome.py`.
- Docking benchmark validation: `WY/benchmark/lit_pcba_benchmark.py` (other benchmark
  scripts live in the same directory).

## Disclaimer

This project is for research purposes only; analysis results should be verified against
the official databases and experimental validation. When using the data, please follow
the citation requirements of each data source (STRING, DisGeNET, ChEMBL, DGIdb,
Open Targets, GEO). The third-party component license list is in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). If you have any questions about this
project, please contact the project author.
`NPAP` is this project's own abbreviation and is not affiliated with any third-party
organisation, project or trademark.
