# QUICKSTART — running the pipeline locally

A copy-paste guide. Works on Linux and macOS. For Windows, use WSL2.

---

## 0. What you should have

After extracting the release, the folder `taxon_aware_amr/` should contain:

```
taxon_aware_amr/
├── README.md
├── QUICKSTART.md                  ← this file
├── CHANGELOG.md
├── LICENSE
├── pyproject.toml                 ← Python package metadata
├── envs/
│   └── environment.yml            ← conda environment (real-mode dependencies)
├── examples/
│   └── kenya_biodiversity_input.tsv   ← demo input
├── src/taxon_aware_amr/           ← all the Python code (14 modules)
│   ├── __init__.py
│   ├── __main__.py
│   ├── amr.py
│   ├── assembly.py
│   ├── cache.py
│   ├── cli.py
│   ├── final_output.py
│   ├── fixtures.py
│   ├── gene_count.py
│   ├── input.py
│   ├── ncbi.py
│   ├── orchestrator.py
│   ├── reporting.py
│   ├── taxonomy.py
│   └── virulence.py
└── tests/                         ← test suite
    ├── run_all.py
    ├── test_taxonomy.py
    ├── test_step2_pipeline.py
    ├── test_step3_pipeline.py
    ├── test_step4_pipeline.py
    └── test_pipeline_full.py
```

Sanity check from inside the extracted folder:

```bash
ls envs/environment.yml src/taxon_aware_amr/__init__.py
```

Both files must exist. If they don't, the archive was extracted incompletely
— re-extract.

---

## 1. Choose your install path

There are two ways to run the pipeline. **Start with Path A** even if you
plan to run for real — it confirms the package is wired correctly in 30
seconds with no external tools required.

| Path | What you need | What you can do |
|---|---|---|
| **A — Mock mode** | Python ≥ 3.10 + `pip` | Runs the entire pipeline using bundled fixtures. No NCBI, no AMRFinderPlus, no ABRicate needed. Confirms the install works. |
| **B — Real mode** | conda/mamba + the env file | Runs against real NCBI, real AMRFinderPlus, real ABRicate. This is production. |

---

## 2. Path A — Mock mode (smoke test, ~1 minute)

This needs **only Python** — nothing from conda, no NCBI, no biology tools.

```bash
# Go into the extracted folder
cd taxon_aware_amr

# Make a Python virtual environment (keeps your system Python clean)
python3 -m venv .venv
source .venv/bin/activate         # macOS / Linux
# On Windows PowerShell: .venv\Scripts\Activate.ps1

# Install the package in editable mode (no external dependencies needed)
pip install -e .

# Verify the CLI is installed
taxon_aware_amr --help

# Run the smoke test
taxon_aware_amr examples/kenya_biodiversity_input.tsv \
    --out demo_results \
    --cache demo_cache \
    --use-mock \
    --auto-replace-suppressed
```

You should see a streaming table of decisions, ending with:

```
Done. 16 rows processed.
  Per-row TSV:    demo_results/decisions.tsv
  Markdown report: demo_results/run_report.md
  SQLite cache:    demo_cache/cache.sqlite
```

**Inspect the outputs:**

```bash
# The final dataframe (matches the original task spec exactly)
cat demo_results/final_output.tsv

# Same data with all the audit columns
cat demo_results/final_output_full.tsv

# Human-readable report — open in any markdown viewer
less demo_results/run_report.md
```

If this works, the install is correct. You're ready for Path B.

---

## 3. Path B — Real mode (production)

Real mode calls NCBI's `datasets` CLI, AMRFinderPlus, and ABRicate. These are
all installed via conda from the bioconda channel.

### 3a. Install Miniforge (one-time, if you don't have conda)

If you have any of `conda`, `mamba`, or `micromamba` already, skip this.

```bash
# Linux
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# Follow the prompts; restart your shell after install.

# macOS (Apple Silicon)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
```

### 3b. Create the conda environment (one-time, ~10 minutes)

```bash
cd taxon_aware_amr

# Create the env from the bundled YAML. Use mamba if you have it, it's faster.
conda env create -f envs/environment.yml
# (or:  mamba env create -f envs/environment.yml)

# Activate it
conda activate taxon_aware_amr

# Install the package itself
pip install -e .

# Verify all four CLIs are on PATH
which datasets amrfinder abricate taxon_aware_amr
```

### 3c. Download tool databases (one-time, ~500 MB)

AMRFinderPlus and ABRicate ship without their reference databases — you fetch
them once after install:

```bash
# AMRFinderPlus DB (~250 MB)
amrfinder --update

# ABRicate's VFDB (and AMR DBs, for cross-checking)
# Modern ABRicate (v1.0+) ships with VFDB built-in. Verify:
abricate --list
# You should see vfdb, card, ncbi, resfinder, plasmidfinder, ecoli_vf, ...

# If VFDB isn't listed (older ABRicate), update:
abricate-get_db --db vfdb
```

### 3d. Run on real data

```bash
# Confirm everything works end-to-end (uses real NCBI now — drop --use-mock)
taxon_aware_amr examples/kenya_biodiversity_input.tsv \
    --out results_real_demo \
    --auto-replace-suppressed \
    --threads 8 \
    -v

# Run on your own input table
taxon_aware_amr /path/to/your_surveillance_input.tsv \
    --out results/run_2026-05 \
    --auto-replace-suppressed \
    --threads 8
```

---

## 4. Input file format

A tab-separated file with a header row. Required columns (any order):

| Column | Example | Notes |
|---|---|---|
| `taxon_id` | `1392` | Verbatim from your source table. Suffixes like `1392_ENA` allowed — they're audit-only. Alias `goat_taxon_id` accepted. |
| `genome_id` | `GCA_964341285.1` | NCBI Assembly accession. Use `none` for rows without an assembly. |
| `species_name` | `Bacillus_anthracis` | Free-form. |
| `category` | `Non_eukaryotic_microbe` | Audit-only. NCBI lineage is authoritative; if `category` says "microbe" but NCBI says it's a fungus, NCBI wins and the disagreement is logged. |

Optional columns are preserved and propagated to the audit output:
`ecological_role`, `endemic_status`, `iucn_status`, `genome_assembly_level`.

See `examples/kenya_biodiversity_input.tsv` for a working example.

---

## 5. What you get back

In your `--out` directory:

| File | Use |
|---|---|
| `final_output.tsv` | **The primary deliverable.** Minimal schema matching the task spec: `genome_id`, `total_gene_count`, `virulence_genes`, `virulence_gene_count`, `amr_genes`, `amr_gene_count`. |
| `final_output_full.tsv` | Same rows, 22 columns including all audit context (input verbatim, NCBI canonical, drug classes, virulence sources, status, reasons, warnings). |
| `decisions.tsv` | Per-row decision log streamed during the run. Joinable back to your input. |
| `run_report.md` | Executive summary on top, per-step diagnostics below. Open in a markdown viewer or paste into a report. |

In your `--cache` directory:

| File | Use |
|---|---|
| `cache.sqlite` | NCBI summaries, file inventory, AMR/virulence results. Re-runs use the cache; truly free. Inspect with `sqlite3`. |

---

## 6. Useful flags

```
--use-mock                  Use bundled fixtures instead of calling NCBI/AMRFinder/ABRicate.
                            For smoke tests and CI. No network or tool installs required.
--auto-replace-suppressed   When NCBI says "Replaced by GCA_xxx", chase the replacement
                            (up to 3 hops). Both accessions are logged.
--keep-files                Don't delete downloaded files after each row. Useful for
                            debugging a single problem row.
--threads N                 Threads for AMRFinderPlus + ABRicate (default: 4).
--refresh-cache             Re-query NCBI / re-run tools even for cached rows.
-o, --out PATH              Output directory (default: ./results).
-c, --cache PATH            Cache directory (default: ~/.cache/taxon_aware_amr).
--datasets-binary PATH      Path to datasets binary if not on PATH.
--amrfinder-binary PATH     Path to amrfinder binary if not on PATH.
--abricate-binary PATH      Path to abricate binary if not on PATH.
-v / -vv                    INFO / DEBUG logging.
```

---

## 7. Verifying the install with the test suite

```bash
cd taxon_aware_amr
source .venv/bin/activate              # or: conda activate taxon_aware_amr
python tests/run_all.py
```

Expected output:

```
PASSED: 5/5
  ✓ test_taxonomy.py
  ✓ test_step2_pipeline.py
  ✓ test_step3_pipeline.py
  ✓ test_step4_pipeline.py
  ✓ test_pipeline_full.py
```

If any test fails, the install isn't right — don't proceed to real data.

---

## 8. Common pitfalls

**"taxon_aware_amr: command not found"** — your virtual env or conda env
isn't activated. Run `source .venv/bin/activate` or
`conda activate taxon_aware_amr`.

**"datasets: command not found"** in real mode — the conda env wasn't
activated, or the env build is incomplete. Run `conda activate taxon_aware_amr`
then `which datasets`.

**"amrfinder: command not found"** — same fix. AMRFinderPlus is on the
bioconda channel; if `conda env create` skipped it (e.g. network blip), run
`mamba install -n taxon_aware_amr -c bioconda ncbi-amrfinderplus` to add it.

**"AMRFinderPlus reports no hits and no errors"** — you probably forgot
`amrfinder --update` to download the database. Run it once.

**Conda environment build fails on M1/M2 Mac** — bioconda has limited
arm64 coverage. Either run via Rosetta:
`CONDA_SUBDIR=osx-64 conda env create -f envs/environment.yml`,
or run the pipeline in a Linux Docker container.

**Pipeline hangs on the first row** — NCBI rate-limited you. Reduce
parallelism, wait 5 minutes, retry. The cache will skip rows already done.

**A row keeps failing with the same error** — failed lookups are cached
for 24 h so transient outages don't spam NCBI. To retry now, delete that
row from the cache:

```bash
sqlite3 demo_cache/cache.sqlite \
  "DELETE FROM assembly_summary WHERE accession = 'GCA_XXXXXX.1';"
```

---

## 9. Inspecting the SQLite cache

The cache is plain SQLite — query it with anything that speaks SQL:

```bash
sqlite3 demo_cache/cache.sqlite

sqlite> .tables
assembly_summary  file_inventory  run_log  schema_meta  tool_results

sqlite> SELECT accession, status, error_message
   ...> FROM assembly_summary WHERE status != 'ok';

sqlite> SELECT accession, tool, status, db_version
   ...> FROM tool_results;

sqlite> SELECT accession, file_type, size_bytes,
   ...>   CASE WHEN deleted_at IS NULL THEN 'on-disk' ELSE 'deleted' END as state
   ...> FROM file_inventory;
```

To reset everything:

```bash
rm -rf demo_cache/   # or ~/.cache/taxon_aware_amr/ for the default
```

---

## 10. Next steps

- Read `README.md` for the project's scientific rationale and design choices.
- Read `CHANGELOG.md` for the version history.
- Inspect `src/taxon_aware_amr/` — each module has a docstring at the top
  explaining what it does. The orchestrator (`orchestrator.py`) is the
  best place to start reading the code.
