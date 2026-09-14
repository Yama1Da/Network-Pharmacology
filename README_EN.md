[简体中文](README.md) | EN

# Network Pharmacology Analysis Platform

A one-stop network pharmacology analysis tool built on Streamlit: enter a drug and a disease
name, and it automatically completes the full workflow of
"target acquisition → target integration → PPI network construction → enrichment analysis → molecular docking".

## Features

- **Disease/drug target acquisition**: integrates four data sources — DisGeNET, DGIdb,
  Open Targets, and ChEMBL — with automatic deduplication, and filters out
  non-specific binding proteins (albumin, CYP metabolic enzymes, etc.) to reduce noise
- **PPI network**: builds protein–protein interaction networks based on STRING v12.0,
  with multiple visualization styles
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

- Python ≥ 3.11
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) (add to PATH, or set `VINA_DIR`
  in `WY/API.env` — rename `API.env.example` to `API.env` first)
- A [DisGeNET](https://www.disgenet.org/) API key (academic accounts are free to apply;
  keys from other account tiers also work)

## Quick Start

```bash
# 1. Create a virtual environment
python -m venv .venv
.venv\Scripts\activate        # Windows (Linux/macOS: source .venv/bin/activate)

# 2. Install dependencies
pip install -r requirements.txt
# `openbabel>=3.2.0` can be installed with pip, but conda is easier:
# conda install -c conda-forge openbabel

# 3. Configure API.env
Edit API.env and fill in DisGeNET_API_KEY and the VINA_DIR path

# 4. Launch
streamlit run App.py  # run from the terminal
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
   the fetched disease and drug targets.
4. **PPI network and enrichment**: open the "Data Network" page from the menu to build a STRING
   interaction network from the integrated targets (multiple visualization styles available)
   and run GO / KEGG enrichment analysis.
5. **Molecular docking**: open the "Molecular Docking" page from the menu, select core targets,
   run AutoDock Vina docking, and inspect the 3D binding poses and focused-pocket view online.
6. **Logs and settings**: switch between Chinese/English and view run logs at the top-right;
   the gear button at the bottom-left opens global settings (e.g., target filter thresholds).
   Changing settings automatically invalidates the integrated caches and rebuilds them.

All intermediate results are cached automatically (7 days by default, see `WY/config.py`),
so repeated analyses return in seconds.

## Configuration

| File | Purpose |
|---|---|
| `WY/API.env` | Secrets and executable paths (`DisGeNET_API_KEY`, `VINA_DIR`); not committed to git |
| `WY/config.py` | Data source URLs, rate limits/retries, circuit-breaker cooldown ladders, cache TTL, target filter thresholds |

> Change behavior in `WY/config.py`, put secrets in `WY/API.env`: `WY/API.env` only holds
> sensitive information and local machine paths that stay out of git, while request control,
> cache validity periods, and target noise-filtering thresholds are centralized in `WY/config.py`.

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
- **Molecular docking can't find Vina?** Make sure `vina` is on PATH or the `VINA_DIR` environment
  variable exists, or point `VINA_DIR` in `WY/API.env` to the directory containing Vina
  (a full path to the executable also works).

## Disclaimer

This project is for research purposes only; analysis results should be verified against the
official databases and experimental validation. When using the data, please follow the citation
requirements of each data source (STRING, DisGeNET, ChEMBL, DGIdb, Open Targets).
