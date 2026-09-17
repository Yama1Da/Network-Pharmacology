[简体中文](README.md) | EN

# Network Pharmacology Analysis Platform

A one-stop network pharmacology analysis tool built on Streamlit: enter a drug and a disease
name, and it automatically completes the full workflow of
"target acquisition → target integration → PPI network construction → enrichment analysis → molecular docking".

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

- Python >= 3.13, < 3.14 (development baseline 3.13.x; the pinned dependencies do not
  ship prebuilt wheels for 3.14 yet)
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) executable (per-platform
  install steps in Quick Start step 3; set `VINA_DIR` in `WY/API.env` — rename
  `API.env.example` to `API.env` first — or add its folder to PATH)
- A [DisGeNET](https://www.disgenet.org/) API key (academic accounts are free to apply;
  keys from other account tiers also work)

## Quick Start

```bash
# 1. Install dependencies (choose one)

# Option A (recommended): one-shot conda environment for dependencies (installs openbabel;
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
# openbabel is not in requirements.txt (its pip build often fails on platforms without
# a compiler toolchain): the app starts fine without it — only the receptor conversion
# in molecular docking needs it (lazy import). Install later if needed:
# conda install -c conda-forge openbabel    (or: pip install openbabel)

# 2. Configure API.env (identical for conda and pip: keys are runtime config, not part
#    of any environment)
Copy API.env.example to API.env (keep the template so future updates still arrive),
edit it and fill in DisGeNET_API_KEY

# 3. Install AutoDock Vina (conda-forge has no win-64 build; pick per platform)
# Windows: download vina_1.2.7_windows_x86_64.exe from
#          https://github.com/ccsb-scripps/AutoDock-Vina/releases into any folder,
#          then point VINA_DIR in API.env to that folder (or add it to PATH);
#          vina_split has no official Windows build — the app falls back to a
#          Python splitter, so it is optional.
# Linux/macOS: conda install -c conda-forge vina   # installs vina and vina_split

# 4. Launch
streamlit run App.py
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

All intermediate results are cached automatically (7 days by default, see `WY/config.py`),
so repeated analyses return in seconds. (The 7-day default exists because DGIdb updates
irregularly.)

## Configuration

| File | Purpose |
|---|---|
| `WY/API.env` | Secrets and executable paths (`DisGeNET_API_KEY`, `VINA_DIR`); not committed to git |
| `WY/config.py` | Data source URLs, rate limits/retries, circuit-breaker cooldown ladders, cache TTL, target filter thresholds, STRING version and download concurrency |

> Change behavior in `WY/config.py`, put secrets in `WY/API.env`: `WY/API.env` only holds
> sensitive information and local machine paths that stay out of git, while request control,
> cache validity periods, and target noise-filtering thresholds are centralized in `WY/config.py`.

## Project Structure

```
WY/
├── App.py                  # Home page: start an analysis, target overview, settings entry
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
├── scripts/                # Offline self-check scripts (see "Development & Self-Checks")
└── Static/                 # Runtime caches and analysis results
```

Dependencies are strictly one-way: `pages / App.py → ui/ → core/ → i18n, config`.
`core/` and `i18n/` do not depend on any UI framework (Streamlit) and can be imported,
reused, and tested without a UI environment.

## Development & Self-Checks

```bash
cd WY
python scripts/selfcheck_docking.py    # molecular docking module self-check
python scripts/selfcheck_network.py    # network algorithms & enrichment module self-check
python scripts/build_interactome.py    # download the background interactome if unsure it is already fetched
```

The first two self-check scripts run purely locally (all external requests are stubbed) and
consume no data-source quota. Run them as a regression after changing anything in `core/`.
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
- **Molecular docking can't find Vina?** Make sure `vina` is on PATH or the `VINA_DIR` environment
  variable exists, or point `VINA_DIR` in `WY/API.env` to the directory containing Vina
  (a full path to the executable also works).

## Disclaimer

This project is for research purposes only; analysis results should be verified against the
official databases and experimental validation. When using the data, please follow the citation
requirements of each data source (STRING, DisGeNET, ChEMBL, DGIdb, Open Targets).
