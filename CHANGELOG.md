# Changelog

All notable changes to this project.

## 1.0.1 — 2026-05-13

### Fixed
- **Bacteria, archaea, and viruses were incorrectly classified as
  `skipped_taxon_unresolved` against real NCBI responses.** Two related
  bugs in `ncbi._parse_summary_json` and `taxonomy.lineage_to_group`:
  - The NCBI Taxonomy renamed `superkingdom` → `domain` in 2023. The
    parser only enumerated the older rank set, so the discriminative
    "Bacteria" / "Archaea" string at the top rank was never read.
  - Modern viral classifications omit the literal string "Viruses" and
    instead expose `realm` (Riboviria, Duplodnaviria, ...) and the
    canonical `-viricota` / `-viricetes` phylum/class suffix family.
    `lineage_to_group` only recognised the legacy "Viruses" string.

  `_parse_summary_json` now enumerates 20 ranks (incl. `domain`, `realm`,
  `acellular_root`, subranks) and falls through to any extra ranks the
  CLI might add. `lineage_to_group` now recognises six viral realms and
  the `-viricota` / `-viricetes` suffix family.

- The mock pipeline path was unaffected by these bugs (mock fixtures
  used the legacy rank names), so the test suite passed even with the
  production bug. v1.0.1 keeps the mock path identical.

### Added
- DEBUG-level diagnostic log of the raw `lineage_obj` when extraction
  returns empty. Run with `-vv` to see what NCBI actually returned for
  a row that fails to route — helps catch any further schema drift.

## 1.0.0 — 2026-05-13

Initial production release.

### Added
- **Step 5** — virulence detection. Merges ABRicate-VFDB hits with the
  `Element type == VIRULENCE` rows already in cache from AMRFinderPlus's
  `--plus` run. Deduplicated by gene symbol (case-insensitive).
  Bacterial rows only; non-bacterial assemblies receive `not_applicable`.
- **Step 6** — final output assembly. Writes two TSVs per run:
  - `final_output.tsv` with the minimal schema matching the original task
    specification (`genome_id`, `total_gene_count`, `virulence_genes`,
    `virulence_gene_count`, `amr_genes`, `amr_gene_count`).
  - `final_output_full.tsv` with all audit columns.
- **Step 7** — executive summary prepended to `run_report.md`, with a final
  output preview table and deliverable list. Suitable for direct inclusion
  in a surveillance report.
- **Master test** — `tests/test_pipeline_full.py` exercises Steps 1–7 in a
  single end-to-end run.

### Changed
- **AMRFinderPlus upgraded to combined mode** (`-n -p -g --plus`) for
  highest sensitivity. NCBI's recommended configuration. The Step 2
  bacterial download plan now retrieves genome + protein + GFF3 instead
  of protein-only.
- Cache schema gained a `tool_results` table reused by both AMR (Step 4)
  and virulence (Step 5).

### Fixed
- Removed dead code from earlier orchestrator restructurings.

## 0.4.0-step3 — Step 3 complete
- Gene counting from NCBI summary (primary) or downloaded GFF3 (fallback).
- Auto-replace of NCBI-suppressed assemblies via `--auto-replace-suppressed`,
  chasing "Replaced by GCA_xxx" messages.
- CLI exposed as the `taxon_aware_amr` console script.
- `--use-mock` flag and bundled demo fixtures for offline trials.
- `pyproject.toml` for `pip install -e .`.
- `QUICKSTART.md` documentation.

## 0.3.0-step2 — Step 2 complete
- NCBI `datasets` CLI wrapper (`DatasetsCLIClient`) with `MockNCBIClient`
  for offline testing.
- SQLite-backed cache (`AssemblyCache`) with three tables:
  `assembly_summary`, `file_inventory`, `run_log`.
- Selective downloads via `datasets --include`. Per-row cleanup after success.
- Failure taxonomy: `not_found` / `suppressed` / `network_error` / `parse_error`,
  with NCBI's error message surfaced in the audit warnings.
- Input TSV reader (`input.py`) accepting both `taxon_id` and `goat_taxon_id`.

## 0.2.0-step1-refactor — GCA-centric routing
- Refactored to use the NCBI Assembly accession (GCA/GCF) as the primary key,
  with input `taxon_id` and `category` demoted to audit-only fields.
- Introduced `AssemblyRecord` as the single source of truth.

## 0.1.0-step1 — initial taxon classifier
- `TaxonGroup`, `AnalysisStatus`, `TaxonDecision`.
- `lineage_to_group()` and `decide()` pure routing functions.
- `DecisionReporter` streaming to stdout, TSV log, and markdown summary.
