# taxon_aware_amr

**Taxon-aware AMR and virulence detection for heterogeneous surveillance inputs.**

Version 1.0.0 · Author: Julien A. Nguinkal · BNITM (Hamburg)

---

## What problem this solves

Surveillance input tables in One Health programmes mix wild plants, livestock,
wildlife, fungi, protozoa, viruses and bacteria in a single sheet. The standard
AMR / virulence tools — AMRFinderPlus, CARD, ResFinder, VFDB — are curated
**for bacterial pathogens only**. Running them indiscriminately across a
heterogeneous table produces spurious homology hits to host transporters,
P450s, ABC families, or contamination — not biology.

This pipeline routes each row through a **taxon-aware decision tree** that
runs bacterial AMR/virulence detection only where it's biologically defensible,
and produces a clean dataframe matching the original task specification while
preserving a full audit trail.

## What's in the box

```
taxon_aware_amr/
├── pyproject.toml              installable package (`pip install -e .`)
├── README.md                   this file
├── QUICKSTART.md               copy-paste commands for laptop install + run
├── CHANGELOG.md                version history
├── LICENSE                     MIT
├── envs/environment.yml        conda dependencies (datasets, AMRFinderPlus,
│                                ABRicate, bakta, …)
├── examples/                   demo input TSV from a Kenya surveillance table
├── src/taxon_aware_amr/        focused modules + CLI
│   ├── taxonomy.py              Step 1 — routing decisions
│   ├── reporting.py             Step 1 — stdout + TSV streaming
│   ├── assembly.py              Step 2 — AssemblyRecord
│   ├── cache.py                 Step 2 — SQLite cache
│   ├── ncbi.py                  Step 2 — datasets-cli wrapper
│   ├── input.py                 Step 2 — input TSV reader
│   ├── gene_count.py            Step 3 — NCBI summary + GFF gene counting
│   ├── amr.py                   Step 4 — AMRFinderPlus combined-mode wrapper
│   ├── virulence.py             Step 5 — ABRicate-VFDB + cached AMRFinder --plus
│   ├── final_output.py          Step 6 — final dataframe assembly
│   ├── orchestrator.py          ties Steps 1-7 together
│   ├── cli.py                   argparse entry point
│   ├── fixtures.py              demo data for `--use-mock`
│   └── __main__.py              `python -m taxon_aware_amr`
└── tests/
    ├── test_taxonomy.py
    ├── test_step2_pipeline.py
    ├── test_step3_pipeline.py
    ├── test_step4_pipeline.py
    └── test_pipeline_full.py    master end-to-end (Steps 1-7)
```

## What the pipeline produces

For an input TSV with columns `taxon_id, genome_id, species_name, category`
(and optional audit columns), the pipeline writes to the output directory:

| File | Contents |
|---|---|
| `final_output.tsv` | **Minimal schema matching the original task spec:** `genome_id, total_gene_count, virulence_genes, virulence_gene_count, amr_genes, amr_gene_count`. Non-bacterial rows show `not_applicable`; rows without a usable assembly show NA. |
| `final_output_full.tsv` | Extended schema with all audit columns (input verbatim, NCBI canonical, drug classes, virulence sources, analysis status, reason, hint warnings, gene-count source). |
| `decisions.tsv` | Per-row decision log with every field on `TaxonDecision`. Joinable back to input. |
| `run_report.md` | Human-readable executive summary on top, followed by per-step diagnostics, audit warnings, and a method-notes appendix. Suitable for attaching to a surveillance report. |
| `cache.sqlite` | SQLite cache of assembly summaries, file inventory, and tool results (AMRFinderPlus, ABRicate). Re-runs are free. |

## The seven steps

| Step | Module | What it does | Tools |
|---|---|---|---|
| 1 | `taxonomy.py` | Routes each row to bacterial / non-bacterial / skipped based on **NCBI Assembly lineage** (input taxon_id is audit-only) | — |
| 2 | `ncbi.py` + `cache.py` | Fetches NCBI summaries via `datasets`, caches in SQLite, downloads only what's needed, deletes after success | `datasets` |
| 3 | `gene_count.py` | Reads `annotation_info.stats.gene_counts.total` from cached summary (zero download); falls back to GFF3 parsing | — |
| 4 | `amr.py` | Runs AMRFinderPlus in **combined mode** (`-n -p -g --plus`) on bacterial rows with `--organism` mapping for point mutations | `amrfinder` |
| 5 | `virulence.py` | Merges VFDB hits (via `abricate`) with cached AMRFinderPlus `--plus` virulence rows; dedup by gene symbol | `abricate` |
| 6 | `final_output.py` | Writes the minimal + extended final TSVs | — |
| 7 | `orchestrator.py` | Prepends an executive summary to `run_report.md` with the deliverable preview | — |

## Key design choices

- **GCA-first.** The NCBI Assembly accession is the only thing used for routing. The input table's `taxon_id` field can carry breed/strain suffixes (`9915_1`, `1392_ENA`) and is preserved verbatim as an audit field.
- **Lineage is authoritative.** Even when the input `category` column says `Non_eukaryotic_microbe`, the NCBI lineage decides bacterial vs viral routing. Disagreements are logged but never override.
- **Selective downloads.** A bacterial row pulls genome + protein + GFF3 (~10 MB); a non-bacterial row with NCBI annotation pulls nothing. Plant and mammal genomes in the input table contribute **zero bytes** to the run.
- **AMRFinderPlus combined mode.** Genome + protein + GFF together give NCBI's recommended highest-sensitivity configuration: protein for confirmed gene calls, nucleotide for fragments and unannotated ORFs, GFF to map nucleotide hits back to gene names. Protein-only is faster but misses gene fragments and partial point mutations.
- **Cleanup-after-success.** Each row's downloaded files are deleted as soon as it completes. The SQLite inventory remembers what was downloaded for full traceability.
- **Suppressed-assembly chase.** When NCBI returns "Replaced by GCA_xxx", the pipeline can chase the replacement (`--auto-replace-suppressed`). Both accessions are recorded.
- **Mock backends.** `--use-mock` runs the entire pipeline against bundled fixtures with no external tools required — useful for smoke-testing the install and for CI.

## Install + run

See `QUICKSTART.md` for full commands. TL;DR:

```bash
conda env create -f envs/environment.yml
conda activate taxon_aware_amr
pip install -e .

# Smoke test with bundled mock data
taxon_aware_amr examples/kenya_biodiversity_input.tsv --use-mock --auto-replace-suppressed

# Real run
taxon_aware_amr my_surveillance_input.tsv --auto-replace-suppressed --threads 8 -v
```

## What this pipeline does NOT do (current scope)

- Antifungal resistance markers (CYP51 / efflux). Fungal assemblies receive `not_applicable`.
- Antiparasitic resistance markers (*Plasmodium kelch13*, etc.). Protozoan assemblies receive `not_applicable`.
- Antiviral mutation screening. Viral assemblies receive `not_applicable`.
- De novo annotation for bacterial assemblies without NCBI annotation. (Bakta is in the conda env and ready to wire in if needed.)

Each of these is a defensible Step-8 extension; the current pipeline declines them by design rather than producing biologically meaningless output.

## Authoring

Julien A. Nguinkal — BNITM (One Health Disease Control Group), Hamburg.

## License

MIT. See `LICENSE`.
