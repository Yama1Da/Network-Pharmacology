[简体中文](README.md) | EN

# Network Pharmacology Analysis Platform

A one-stop network pharmacology analysis tool built on Streamlit: enter a drug and a disease
name, and it automatically completes the full workflow of
"target acquisition → target integration → PPI network construction → enrichment analysis → molecular docking".
Batch processing via the CLI is also supported (figures are not generated in CLI mode).

## Features

- **Disease/drug target acquisition**: integrates four data sources — DisGeNET, DGIdb,
  Open Targets, and ChEMBL — with automatic deduplication, and filters out
  non-specific binding proteins (albumin, CYP metabolic enzymes, etc.) to reduce noise;
  the noise-filtering conditions are user-adjustable
- **PPI network**: builds protein–protein interaction networks based on STRING
  (v12.0 by default; the data version is switchable in Settings), with multiple visualization styles
- **Background interactome**: a local cache of the full STRING network powering proximity
  Z-score / RWR target prioritization; one-click build on the Data Management page,
  with mirror/official multi-source fallback and segmented concurrent downloads
- **Enrichment analysis**: GO / KEGG pathway enrichment; automatically falls back to
  Enrichr when g:Profiler is unreachable. Both sources have circuit breakers (escalating
  cooldown ladders with automatic probing and recovery), and when both are down the page
  shows an "automatic retry in about X minutes" notice
- **Molecular docking**: integrates AutoDock Vina, with RDKit + Meeko for automatic
  ligand preparation, built-in 3D result viewing, and focused-pocket display (surroundings dimmed)
- **Engineering**: request-level caching, local DisGeNET rate limiting (within the academic
  quota), circuit-breaker failover and automatic recovery for enrichment data sources,
  bilingual Chinese/English UI, daily rotating logs, and data management pages
  for drug/disease targets

## Requirements

The requirements below apply to the local conda / pip routes; with Docker deployment, the host
only needs Docker and a DisGeNET key (see Docker Deployment below).

- Python >= 3.13, < 3.14 (development baseline 3.13.x; the pinned dependencies do not
  ship prebuilt wheels for 3.14 yet)
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) executable (per-platform
  install steps in Quick Start step 3; set `VINA_DIR` in `WY/API.env` — rename
  `API.env.example` to `API.env` first — or add its folder to PATH)
- The `obabel` command line from [Open Babel](https://openbabel.org/) (receptor conversion
  in molecular docking runs it as a subprocess; per-platform install steps in Quick Start
  step 3; set `OBABEL_DIR` or add its folder to PATH)
- A [DisGeNET](https://www.disgenet.org/) API key (academic accounts are free to apply;
  keys from other account tiers also work)

## Quick Start

```bash
# 1. Install dependencies (choose one)

# Option A (recommended): one-shot conda environment for dependencies (installs obabel;
# works on Windows/Linux/macOS). Keys and local paths are NOT part of the env —
# steps 2 and 3 below are still required.
cd WY
conda env create -f environment.yml
conda activate npy

# Option B: pip + venv
cd WY
python -m venv .venv
.venv\Scripts\activate        # Windows (Linux/macOS: source .venv/bin/activate)
pip install -r requirements.txt
# Receptor conversion additionally needs the obabel command line (invoked as a
# subprocess, not the Python bindings): install the official Open Babel installer on
# Windows, or openbabel via apt/dnf on Linux, or `conda install -c conda-forge openbabel`;
# make sure obabel is on PATH, or point OBABEL_DIR in API.env to its folder

# 2. Configure API.env (identical for conda and pip: keys are runtime config, not part
#    of any environment)
Copy API.env.example to API.env (keep the template so future updates still arrive),
edit it and fill in DisGeNET_API_KEY

# 3. Install AutoDock Vina and Open Babel (conda-forge has no win-64 vina build; obabel
#    is needed for receptor conversion)
Windows: download vina_1.2.7_win.exe (optionally together with vina_split_1.2.7_win.exe) 
#          from https://github.com/ccsb-scripps/AutoDock-Vina/releases into any folder,
#          then point VINA_DIR in API.env to that folder (or add it to PATH);
#          vina_split is optional — the app falls back to a Python splitter.
#          For Open Babel, grab the Windows installer (includes obabel.exe) from
#          https://openbabel.org/ and let it add itself to PATH; otherwise point
#          OBABEL_DIR in API.env to its folder.
Linux/macOS: conda install -c conda-forge vina   # installs vina and vina_split    
#          obabel: conda install -c conda-forge openbabel (or install openbabel via apt/dnf)

# 4. Launch
start_npap.bat  (or npap serve)

# 5. (Optional) install the command-line batch analysis as a global `npap` command
#    (see "Command-Line Batch Analysis (CLI)")
pip install -e ..        # run inside the WY directory; -e (editable) is required
npap --help              # usable from any directory while the environment is active
```

## Usage Workflow

1. **Start an analysis on the home page**: enter a drug name and a disease name
   (the default example is metformin and Type 2 diabetes mellitus; Chinese disease
   names are not supported yet), then click "Analyze". The program queries
   ChEMBL / DGIdb / Open Targets / DisGeNET concurrently.
2. **Review the target overview**: once the analysis finishes, the home page shows the number of
   drug targets, disease targets, and shared targets, along with a Venn diagram.
   Expand the collapsible panels to view per-source target detail tables
   (sorted by source count and score).
3. **Data management**: use the "Data Management" menu at the top-left of the home page
   to open the information pages, where you can view, add, remove, or correct
   the fetched disease and drug targets. The same page builds the background interactome
   in one click (a prerequisite for proximity / RWR analyses).
4. **PPI network and enrichment**: open the "Data Network" page from the menu to build a STRING
   interaction network from the integrated targets (multiple visualization styles available)
   and run GO / KEGG enrichment analysis.
5. **Molecular docking**: open the "Molecular Docking" page from the menu, select core targets,
   run AutoDock Vina docking, and inspect the 3D binding poses and focused-pocket view online.
6. **Logs and settings**: switch between Chinese/English and view run logs at the top-right;
   the gear button at the bottom-left opens global settings (e.g., target filter thresholds)
   and the STRING data version (when STRING releases a new version, just enter the tag in
   Settings and rebuild the background network — no program update needed).
   Changing settings automatically invalidates the integrated caches and rebuilds them.

All intermediate results are cached automatically (see `WY/config.py`), so repeated analyses
return in seconds. Target pulls from the four databases (ChEMBL / DGIdb / Open Targets /
DisGeNET) share a unified 60-day cache: the upstream release cycles range from about a quarter
to ~2.5 years (Open Targets releases quarterly, ChEMBL roughly every 7 months, DISGENET every
1~2 quarters, DGIdb major versions ~2.5 years apart), so 60 days never misses a new release
while keeping repeated requests — notably DisGeNET's metered API quota — to a minimum.
Integrated tables, PPI, enrichment and other analysis caches default to 7 days. For an
immediate refresh, hit "Update" on a management page to force a fresh download.

## Command-Line Batch Analysis (CLI)

`WY/cli.py` is a thin wrapper over `core/`: it shares all caches and settings with the web UI
(pairs the web has already analyzed return instantly; thresholds changed in the web gear apply
to the CLI automatically). It targets wet-lab + literature-comparison batch tabulation.
No plots, no docking — paper figures come from the web UI (all downloadable); the CLI
produces pure data output.

```bash
cd WY
conda activate npy        # or activate .venv

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

# Batch: pairs.csv must contain 'drug' and 'disease' columns (case-insensitive; extra
# columns like Excel row numbers are ignored; edit in Excel and save as CSV),
# one pair per row (no cross product); one subdirectory per pair + batch_summary.csv
python cli.py batch --pairs D:\cases\pairs.csv --out D:\cases\results
python cli.py batch --pairs D:\cases\pairs.csv --out D:\cases\results --resume   # checkpoint resume
```

- **Exit codes**: 0 when everything succeeds, 1 when anything fails (ready for CI/self-check
  scripts); failures print a concrete reason (both enrichment sources down, disease name with
  no returns, background interactome missing, etc.) — results are never silently incomplete.
- **One failed pair never aborts the batch**; `--resume` skips pairs that already succeeded
  and retries failed ones.
- Gene lists are **interchangeable** across `--genes` (`ppi`/`enrich`) and `--seeds` (`rwr`),
  in any of three forms: a .txt with one gene per line (`#` comments allowed), a .csv
  (auto-detects the `gene`/`approvedSymbol` column), or a comma-separated literal list.
  Name-based commands (`targets`/`proximity`/`analyze`/`batch`) take drug/disease names only;
  the two groups connect through the produced CSV files
  (e.g. `enrich --genes <pair dir>\common_targets.csv`).
- proximity / RWR require the background interactome: when it is missing they are skipped with
  a notice and never trigger the large download automatically (build it via the web UI or
  `scripts/build_interactome.py`).
- English disease names only, same as the web UI; run batches sequentially — do not start
  multiple batch processes at once.

**Install as a global `npap` command (optional; during development `python cli.py` works as is; Windows/Linux/macOS alike)**:

```bash
cd WY
pip install -e ..        # conda users: conda activate npy first; -e (editable) is required, see below
npap --help              # the entry script lands in the env's Scripts\ (Windows) or bin/ (Linux/macOS)
```

- **Editable install (`-e`) is mandatory**: the project is designed to run in place — caches,
  logs, and `API.env` are all resolved relative to the `WY/` source directory; a regular
  (non-editable) install would write them into site-packages.
- After installing, `npap` is available in any directory **while the environment is active**
  (identical on all three platforms); to use it in a terminal without activation, add the
  environment's `Scripts\` (Windows) or `bin/` (Linux/macOS) directory to PATH.
- Windows-only install-free transitional script: copy `npap.bat.example` as `npap.bat` and
  put its folder on PATH (machine-specific paths, kept out of git like `API.env`; the default
  works unmodified for the standard layout, conda and other layouts should replace the python
  path inside with an absolute one). On Linux/macOS, skip the local environment entirely with
  Docker (see Docker Deployment below).
- For global use, set `CLI_RESULTS_DIR` in `API.env` once: without `--out`, results land in
  the `results` folder under the **current directory**.

## Configuration

| File           | Purpose                                                                                                                                               |
|----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| `WY/API.env`   | Secrets and executable paths (`DisGeNET_API_KEY`, `VINA_DIR`, `OBABEL_DIR`, `CLI_RESULTS_DIR`); not committed to git                                  |
| `WY/config.py` | Data source URLs, rate limits/retries, circuit-breaker cooldown ladders, cache TTL, target filter thresholds, STRING version and download concurrency |

> Change behavior in `WY/config.py`, put secrets in `WY/API.env`: `WY/API.env` only holds
> sensitive information and local machine paths that stay out of git, while request control,
> cache validity periods, and target noise-filtering thresholds are centralized in `WY/config.py`.

## Project Structure

```
WY/
├── api.py                  # FastAPI: /api endpoints + built frontend hosting
├── cli.py                  # Command-line batch analysis entry (targets/ppi/enrich/proximity/rwr/analyze/batch)
├── config.py               # All tunable parameters (paths, rate limits/retries, cache TTL, filter thresholds)
├── API.env                 # Secrets and local paths (not committed; see API.env.example)
├── pages/                  # Streamlit subpages (target management, data network, enrichment, docking, logs, etc.)
├── core/                   # Domain & orchestration layer: data-source fetching, PPI, enrichment, docking, cache-file management
│   ├── management.py       #   Target data file refresh orchestration, settings persistence, cache-name display
│   ├── targets.py          #   Multi-source target fetching and integration
│   ├── ppi.py / enrichment.py / docking.py / interactome.py / network_algorithms.py
│   └── http_client.py / caching.py / text_utils.py / logging_setup.py  # infrastructure
├── ui/                     # UI component layer
│   ├── widgets.py          #   Reusable cross-page widgets (file tables, batch action bars, search box, etc.)
│   └── source_management.py  # Full-page rendering of the drug/disease target management pages
├── i18n/                   # Chinese/English translation dictionaries and lookup API
├── scripts/                # Offline self-check & benchmark scripts (see "Development & Self-Checks")
└── Static/                 # Runtime caches and analysis results
```

Dependencies are strictly one-way: `web (TS+Vue) / api.py → core/ → i18n, config`.
`core/` and `i18n/` do not depend on any UI framework (Streamlit) and can be imported,
reused, and tested without a UI environment; `cli.py` likewise depends only on `core/`
and `config`, sharing caches and settings with the web UI.

## Development & Tests

```bash
pytest tests/                          # full unit-test suite (120 tests, purely local, external requests stubbed)
python scripts/build_interactome.py    # download the background interactome if unsure it is already fetched
python scripts/lit_pcba_benchmark.py   # LIT-PCBA docking-benchmark validation script
```

The test suite runs purely locally (all external requests are stubbed) and consumes no
data-source quota. Run it as a regression after changing anything in `core/`.
For **build_interactome.py**, prefer clicking "One-click background interactome build" on the
Data Management page, which ensures the background interactome has been downloaded.

## FAQ

- **DisGeNET returns 429?** The app already applies local rate limiting within the academic quota
  (50 requests/min). If it still triggers, just wait a moment; the quota is 25,000 requests/day,
  and a normal single analysis consumes only 1 request. If your account has a different tier's
  quota, adjust `DISGENET_RATE_LIMIT_PER_MINUTE`, `DISGENET_RATE_LIMIT_PER_DAY`, and
  `DISGENET_RATE_LIMIT_PER_MONTH` in `WY/config.py` yourself.
- **No enrichment results?** Repeated g:Profiler failures trip a circuit breaker and switch to
  Enrichr; when both sources are down the page shows an "automatic retry in about X minutes"
  notice. The cooldown ladders (starting at 30 s, escalating to a 30-minute cap) are configured
  in `GPROFILER_CIRCUIT_COOLDOWN_LADDER` / `ENRICHR_CIRCUIT_COOLDOWN_LADDER`.
- **Changed the filter thresholds but results look the same?** Threshold changes require
  incrementing `TARGET_FILTER_VERSION` accordingly; old integrated caches then expire
  automatically and are rebuilt.
- **Background interactome downloads slowly?** The CN mirror usually carries only the latest
  release (currently v12.0 only, as measured); other versions automatically fall back to the
  official download site (slow per connection). Raise `INTERACTOME_DOWNLOAD_CONCURRENCY`
  (default 6) in `WY/config.py` to speed things up — the official site throttles per
  connection, so more segments scale nearly linearly; or stick to a version the mirror
  already carries.
- **Molecular docking can't find Vina / obabel?** Make sure the executable is on PATH, or
  point `VINA_DIR` / `OBABEL_DIR` in `WY/API.env` to the corresponding directory (a full
  path to the executable also works). The `obabel` used for receptor conversion comes from
  Open Babel — install steps in Quick Start step 3.
- **CLI results not where you expected?** Without `--out`, the output location priority is:
  the `--out` argument > `CLI_RESULTS_DIR` (in `WY/API.env` or an environment variable of the
  same name) > the `results` folder in the current directory; `--out` still overrides for a
  single run — no need to remove a configured `CLI_RESULTS_DIR`.
- **Docker container won't start or the page won't open?** Check `docker compose ps` and
  `docker logs network-pharmacy` first; on a remote server open `http://<server IP>:8000`
  (not localhost) and make sure port 8000 is allowed through the firewall/security group.

## Disclaimer

This project is for research purposes only; analysis results should be verified against the
official databases and experimental validation. When using the data, please follow the citation
requirements of each data source (STRING, DisGeNET, ChEMBL, DGIdb, Open Targets).
If you have any questions about this project, please contact the project author.
