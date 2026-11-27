# NPAP Command-Line (CLI) Manual
[简体中文](CLI_Manual_zh.md) | [English](#english-version)

---

## English version

### 0. Prerequisites

1. Environment installed via the launch scripts (`start_npap.bat` / `start_npap.sh`);
   virtualenv lives at `WY/.venv`:
   ```bash
   cd WY
   source .venv/bin/activate        # Windows: .venv\Scripts\activate
   python cli.py --help             # or npap --help once installed globally
   ```
2. **Background interactome** (required by proximity / rwr / modules / synergy):
   click "Build" on the web Data page, or run `WY/tools/build_interactome.py`.
   Without it these commands warn and exit — they never trigger large downloads.
3. **DisGeNET API key** (optional): without it the disease target list loses one
   source; everything else works.
4. Disease names must be **in English** (e.g. `Type 2 diabetes mellitus`).

### 1. Command overview

| Command     | Purpose                                                                    | Typical runtime*                             |
|-------------|----------------------------------------------------------------------------|----------------------------------------------|
| `serve`     | Launch the web server (opens browser)                                      | seconds                                      |
| `targets`   | Target profiles from 4 sources (ChEMBL/DGIdb/Open Targets/DisGeNET)        | ~1 min (instant on cache)                    |
| `ppi`       | STRING network + topology metrics + MCC core targets                       | sec–min                                      |
| `enrich`    | GO/KEGG enrichment (g:Profiler, auto-fallback to Enrichr)                  | ~30 s                                        |
| `proximity` | Drug–disease module proximity Z-score                                      | ~2 min (200 perms)                           |
| `rwr`       | RWR whole-network ranking from seed genes                                  | ~1 min                                       |
| `analyze`   | Full pipeline for one pair (targets→intersection→PPI→enrich→proximity→RWR) | 3–5 min                                      |
| `batch`     | Many pairs + resume + batch summary                                        | pairs × analyze                              |
| `modules`   | Disease-module consensus (RWR/DIAMOnD/neighbor) + permutation significance | minutes on first run (try small perms first) |
| `synergy`   | Network-separation synergy prediction (single pair or batch)               | minutes/pair                                 |

\* First runs hit the network; cache hits (TTL 7 days, target sources 60 days)
are much faster.

### 2. Conventions

- **Output directory precedence**: `--out` > `CLI_RESULTS_DIR` (in `WY/API.env`
  or the same-named env var) > `./results` in the current directory.
- **Exit codes**: 0 = all succeeded; 1 = any failure (reason printed, results
  are never silently incomplete).
- **Gene lists, three forms** (`--genes`/`--seeds` interchangeable):
  ① `.txt` one gene per line (`#` comments allowed); ② `.csv` (auto-detects the
  `gene`/`approvedSymbol` column); ③ comma-separated string, e.g. `"TP53,EGFR"`.
  `targets`/`proximity`/`analyze`/`batch` take drug/disease names; the two
  command groups chain through output files.
- **Same origin as the web UI**: identical caches and Settings thresholds.
- **Batch discipline**: run sequentially; do not run multiple batch processes
  concurrently (cache write contention).

### 3. Command reference

#### 3.1 `targets` — four-source target profiles

```bash
python cli.py targets --drug aspirin --disease "colorectal cancer"
python cli.py targets --drug aspirin                # drug side only
```
Outputs: `<drug>_targets.csv` / `<disease>_targets.csv` /
`<drug>_<disease>_common_targets.csv`.
Notes: empty DisGeNET results — check the key in `WY/API.env` first; a failed
source is treated as empty with a warning and the **exit code is 1** (the
profile is incomplete for literature cross-checks).

#### 3.2 `ppi` — STRING network and core targets

```bash
python cli.py ppi --genes common_targets.csv
```
Outputs: edge list / metrics table (degree, betweenness, …) / MCC core-target
table. Gene counts above the STRING limit are truncated with a warning; shares
the edge cache with the web "manual genes" mode.

#### 3.3 `enrich` — GO/KEGG enrichment

```bash
python cli.py enrich --genes common_targets.csv --pvalue 0.05 --sources GO:BP,KEGG
```
Outputs: enrichment long-table CSV. g:Profiler failures trip a circuit breaker
that switches to Enrichr; if both are down the page reports the retry countdown.

#### 3.4 `proximity` — proximity Z-score

```bash
python cli.py proximity --drug aspirin --disease "colorectal cancer"
```
Outputs: `proximity.json` (Z-score, empirical p, degree-matched null parameters).
Exits 1 with a warning when the interactome is missing.

#### 3.5 `rwr` — RWR whole-network ranking

```bash
python cli.py rwr --seeds core_targets.csv
python cli.py rwr --seeds "TP53,EGFR,MTOR"
```
Outputs: `rwr.csv` (gene, score, is_seed; descending score = priority).

#### 3.6 `analyze` — full pipeline for one pair

```bash
python cli.py analyze --drug aspirin --disease "colorectal cancer"
```
Outputs: a result subdirectory + `summary.json` (**parameter snapshot + per-step
status/timing**, fully traceable). One failed step does not block the rest;
failed steps are listed at the end.

#### 3.7 `batch` — batch pipeline

```bash
python cli.py batch --pairs pairs.csv --out results
python cli.py batch --pairs pairs.csv --out results --resume   # resume
```
`pairs.csv` needs `drug` and `disease` columns (one pair per row). Outputs: one
subdirectory per pair + `batch_summary.csv` (status/duration/failure reason,
one row per pair). `--resume` skips succeeded pairs and retries failed ones;
any failure → exit code 1.

#### 3.8 `modules` — disease-module three-method consensus 

```bash
python cli.py modules --genes disease_targets.csv                     # defaults N=100, with significance
python cli.py modules --genes disease_targets.csv --top-n 20 --n-perms 10   # small trial run
python cli.py modules --genes disease_targets.csv --n-perms 0        # skip significance
```
Parameters: `--top-n` expansion size N per method (default 100, shared by all
three); `--k` consistency threshold (default 2: recalled by ≥2 methods);
`--n-perms` permutation rounds (default 200; **the DIAMOnD null is minutes-scale
on real graphs — try small round counts first**).
Outputs: `consensus.csv` (per-method ranks + Borda consensus + consistent flag) /
`modules.csv` (long-format members) / `significance.csv` (per-method permutation
p) / `summary.json` (parameter snapshot; seeds outside the LCC are recorded in
`seeds_out_of_lcc`, never silently dropped).

#### 3.9 `synergy` — network-separation synergy prediction 

```bash
# single pair
python cli.py synergy --a Paeonol --b Luteolin --disease "breast carcinoma"
# no disease anchor (pure separation) is also valid
python cli.py synergy --a Paeonol --b Luteolin
# batch: CSV must contain drug_a; drug_b/disease columns optional (fall back to --b/--disease)
python cli.py synergy --input pairs.csv --out results
python cli.py synergy --input pairs.csv --b melatonin --top-n 50 --n-perm 200
```
Behaviour: each pair archives a contract JSON into the shared `synergy-results/`
store (visible on the web results page); batches additionally produce
`synergy_batch_summary.csv` (**including the five E6 reachability columns**:
r_tier_a / r_tier_b / r_high_conf / r_sort_key / r_sort_key_note — sort in
Excel by r_sort_key descending for the E6 global ranking). `--top-n` target
slice (0 = full sets, default 100); `--n-perm` permutation rounds (≥10,
default 200). Per-target detail CSVs live in each run directory
(`reachability_<drug>.csv`). The CLI plots nothing; underpowered pairs carry a
`power_limited` flag.

#### 3.10 `serve` — launch the web server

```bash
python cli.py serve                    # opens the browser automatically
python cli.py serve --port 8000 --no-browser
```

### 4. End-to-end tutorials

**Tutorial A: complete network-pharmacology analysis for one drug**
```bash
python cli.py analyze --drug aspirin --disease "colorectal cancer"
# split into targets → ppi/enrich → proximity/rwr to inspect intermediate steps
```

**Tutorial B: batch screening with resume**
```bash
# 1) prepare pairs.csv (columns: drug, disease)
python cli.py batch --pairs pairs.csv --out results --resume
# 2) inspect failed rows in batch_summary.csv
```

**Tutorial C: synergy batch + reachability ranking**
```bash
# 1) prepare pairs.csv (columns: drug_a, drug_b, disease; disease may be empty)
python cli.py synergy --input pairs.csv --out results
# 2) open synergy_batch_summary.csv, sort by r_sort_key descending = E6 ranking
#    (r_high_conf=True rows = both drugs evidence-backed, highest trust)
```

### 5. Common errors

| Symptom                                  | Cause & fix                                                                                               |
|------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| `background interactome not built`       | prerequisite of proximity/rwr/modules/synergy: build on the Data page or run `tools/build_interactome.py` |
| `cannot write 'xxx.csv': file is locked` | Excel holds the output file — close it and rerun                                                          |
| DisGeNET empty / 429                     | check the key in `WY/API.env`; the client rate-limits (45 req/min), wait and retry                        |
| `no genes parsed from --genes`           | encoding/column issue: csv needs a `gene` or `approvedSymbol` column                                      |
| results not in the expected directory    | precedence `--out` > `CLI_RESULTS_DIR` > `./results`; set `CLI_RESULTS_DIR` once in `API.env`             |
| `synergy:no_targets_*` for one pair      | no targets in any source (likely a typo); other pairs continue, exit code 1                               |
| results unchanged after threshold edit   | bump the matching `*_VERSION` in `WY/config.py` so old caches invalidate                                  |

### 6. Global `npap` command (optional)

```bash
cd WY && pip install -e ..    # -e is required: caches/logs resolve relative to WY/
npap --help                   # usable anywhere inside the activated env
```
Set `CLI_RESULTS_DIR` once in `API.env`, otherwise results land in `./results`
of the current directory when `--out` is omitted.
