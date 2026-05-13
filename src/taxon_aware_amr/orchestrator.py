"""
taxon_aware_amr.orchestrator
============================

Wires the pieces together for an end-to-end run.

Flow (Steps 1 + 2 + 3)
----------------------
For each input row::

    1. normalise_genome_id       → cleaned accession or None
    2. ncbi.fetch_summary        → AssemblyFetch (uses cache)
       └─ if 'suppressed' and --auto-replace-suppressed:
              try fetching the replacement accession, log the substitution
    3. taxonomy.decide           → TaxonDecision (Step 1 routing)
    4. plan_downloads            → minimal set of files needed
    5. ncbi.download_files       → fetch them (if any)
    6. gene_count.count_genes    → Step 3: fill total_gene_count
    7. (Steps 4-5 plug in here: AMR, virulence)
    8. reporter.log              → stdout + decisions.tsv
    9. cleanup_files             → free disk, keep inventory
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .assembly      import normalise_genome_id
from .amr           import (
    AMRClient, AMRResult, amrfinder_organism_for,
)
from .cache         import AssemblyCache
from .final_output  import FinalRow, build_final_row, write_minimal_tsv, write_full_tsv
from .gene_count    import count_genes
from .input         import InputRow, read_input
from .ncbi          import AssemblyFetch, DownloadResult, NCBIClient, VALID_INCLUDES
from .reporting     import DecisionReporter
from .taxonomy      import AnalysisStatus, TaxonDecision, TaxonGroup, decide
from .virulence     import (
    ABRicateClient, VirulenceResult, detect_virulence,
)

logger = logging.getLogger(__name__)


# ===========================================================================
# Auto-replace helpers
# ===========================================================================

# NCBI suppression messages embed the replacement accession in free text.
# Real-world examples:
#   "Replaced by GCA_022413745.1"
#   "Replaced by accession GCA_022413745.1"
#   "This assembly has been replaced by GCF_000123456.2"
_REPLACEMENT_RE = re.compile(
    r"replaced\s+by\s+(?:accession\s+)?(GC[AF]_\d+\.\d+)",
    re.IGNORECASE,
)


def extract_replacement_accession(message: Optional[str]) -> Optional[str]:
    """Pull the replacement GCA/GCF out of an NCBI suppression message."""
    if not message:
        return None
    m = _REPLACEMENT_RE.search(message)
    return m.group(1) if m else None


# ===========================================================================
# Per-row execution plan
# ===========================================================================

@dataclass
class RowPlan:
    """What we will (or won't) download for a given decision."""
    needs_download: bool
    includes:       tuple[str, ...]
    reason:         str


def plan_downloads(
    decision: TaxonDecision, gene_count_from_summary: Optional[int],
) -> RowPlan:
    """Decide the minimal set of files to download.

    Rules
    -----
    1. Skipped rows                          → nothing.
    2. Bacterial rows                        → genome FASTA + protein FASTA + GFF3.
                                               AMRFinderPlus is invoked in combined
                                               mode (-n -p -g) for highest sensitivity:
                                               protein for annotated genes, nucleotide
                                               for fragments/unannotated ORFs, GFF to
                                               map nucleotide hits back to genes.
                                               GFF also serves as Step 3 gene-count
                                               fallback if NCBI lacks counts.
    3. Non-bacterial with NCBI gene count    → nothing — count from the summary.
    4. Non-bacterial WITHOUT NCBI gene count → gff3 only.
    """
    if decision.status not in (AnalysisStatus.OK_BACTERIAL,
                               AnalysisStatus.OK_GENECOUNT_ONLY):
        return RowPlan(False, (), "skipped row — nothing to retrieve")

    if decision.status is AnalysisStatus.OK_BACTERIAL:
        return RowPlan(True, ("genome", "protein", "gff3"),
                       "bacterial: genome+protein+gff3 for AMRFinderPlus combined mode")

    if gene_count_from_summary is not None:
        return RowPlan(False, (),
                       "non-bacterial: gene count in NCBI summary; no download")
    return RowPlan(True, ("gff3",),
                   "non-bacterial: gff3 only, for gene counting")


# ===========================================================================
# Per-row outcome bookkeeping
# ===========================================================================

@dataclass
class RowOutcome:
    decision:    TaxonDecision
    fetch:       AssemblyFetch
    plan:        RowPlan
    download:    Optional[DownloadResult]
    bytes_freed: int
    final_row:   Optional[FinalRow] = None


# ===========================================================================
# Fetch (with optional auto-replace of suppressed accessions)
# ===========================================================================

def fetch_with_replacement(
    accession: str,
    *,
    client: NCBIClient,
    auto_replace_suppressed: bool,
    max_hops: int = 3,
) -> tuple[AssemblyFetch, list[tuple[str, str]]]:
    """Fetch the summary; if suppressed and the message names a replacement,
    optionally chase it. Returns (final_fetch, hops).

    ``hops`` is a list of (from, to) pairs documenting every substitution
    made. Empty list when no auto-replace happened.
    """
    fetch = client.fetch_summary(accession)
    hops: list[tuple[str, str]] = []

    if not auto_replace_suppressed:
        return fetch, hops

    seen = {accession}
    while (fetch.status == "suppressed" and len(hops) < max_hops):
        replacement = extract_replacement_accession(fetch.error_message)
        if not replacement or replacement in seen:
            break
        hops.append((fetch.accession or accession, replacement))
        seen.add(replacement)
        fetch = client.fetch_summary(replacement)

    return fetch, hops


# ===========================================================================
# Run a single row
# ===========================================================================

def run_row(
    row: InputRow,
    *,
    client: NCBIClient,
    cache: AssemblyCache,
    download_dir: Path,
    keep_files: bool,
    auto_replace_suppressed: bool,
    amr_client: Optional[AMRClient] = None,
    abricate_client: Optional[ABRicateClient] = None,
) -> RowOutcome:
    """Run Steps 1 → 5 for a single input row."""
    clean_id = normalise_genome_id(row.genome_id)

    # ---------- 1. Fetch (optionally chasing replacements) ----------
    if clean_id is None:
        fetch = AssemblyFetch(accession="", status="no_genome_id",
                              assembly_record=None, gene_count_total=None,
                              error_message=None, from_cache=False)
        hops: list[tuple[str, str]] = []
    else:
        fetch, hops = fetch_with_replacement(
            clean_id, client=client,
            auto_replace_suppressed=auto_replace_suppressed,
        )

    # ---------- 2. Routing decision ----------
    decision = decide(
        assembly=fetch.assembly_record,
        input_genome_id=clean_id,
        input_taxon_id=row.taxon_id,
        input_species_name=row.species_name,
        input_category=row.category,
    )
    if fetch.status not in ("ok", "no_genome_id") and fetch.error_message:
        decision.hint_warnings.append(
            f"NCBI status: {fetch.status} — {fetch.error_message}"
        )
    if hops:
        chain_str = " → ".join([clean_id] + [h[1] for h in hops])
        decision.hint_warnings.append(
            f"Auto-replaced suppressed accession: {chain_str}"
        )

    # ---------- 3. Download plan ----------
    plan = plan_downloads(decision, fetch.gene_count_total)
    download: Optional[DownloadResult] = None
    if plan.needs_download and decision.ncbi_accession is not None:
        download = client.download_files(
            decision.ncbi_accession, plan.includes, download_dir,
        )
        if download.skipped:
            decision.hint_warnings.append(f"Download skipped: {download.reason}")

    # ---------- 4. Step 3: gene counting ----------
    gff_path: Optional[Path] = None
    if download and "gff3" in download.files:
        gff_path = download.files["gff3"]
    count, source = count_genes(
        summary_gene_count=fetch.gene_count_total,
        gff_path=gff_path,
    )
    # Only meaningful for rows where the assembly was actually fetched.
    if decision.status in (AnalysisStatus.OK_BACTERIAL,
                           AnalysisStatus.OK_GENECOUNT_ONLY):
        decision.total_gene_count = count
        decision.gene_count_source = source
        if count is None:
            decision.hint_warnings.append(
                f"Gene count unavailable: {source}"
            )

    # ---------- 5. Step 4: AMR detection (bacterial only, combined mode) ----------
    amr_result: Optional[AMRResult] = None
    if decision.status is AnalysisStatus.OK_BACTERIAL:
        protein_path = (download.files.get("protein")
                        if download and download.files else None)
        genome_path  = (download.files.get("genome")
                        if download and download.files else None)
        gff_path     = (download.files.get("gff3")
                        if download and download.files else None)
        if amr_client is None:
            decision.amr_status = "not_attempted"
            decision.amr_method = "no_client_configured"
            decision.amr_genes = ""
            decision.hint_warnings.append(
                "AMR client not configured — skipping AMRFinderPlus.")
        elif protein_path is None:
            decision.amr_status = "failed"
            decision.amr_method = "no_protein_fasta"
            decision.hint_warnings.append(
                "Protein FASTA missing for AMR detection — check Step 2 plan.")
        else:
            organism_flag = amrfinder_organism_for(decision.ncbi_organism_name)
            amr_result = amr_client.run(
                decision.ncbi_accession, protein_path,
                organism_flag=organism_flag,
                genome_path=genome_path,
                gff_path=gff_path,
            )
            decision.amr_status = amr_result.status
            decision.amr_method = f"amrfinderplus_{amr_result.db_version}"
            decision.amr_organism_flag = organism_flag or "auto:none"
            if amr_result.status == "ok":
                decision.amr_genes = "; ".join(amr_result.amr_gene_symbols())
                decision.amr_gene_count = len(amr_result.amr_hits)
                decision.amr_drug_classes = "; ".join(amr_result.amr_drug_classes())
            else:
                decision.amr_genes = ""
                decision.amr_gene_count = None
                decision.hint_warnings.append(
                    f"AMRFinderPlus failed: {amr_result.error_message}")
    elif decision.status is AnalysisStatus.OK_GENECOUNT_ONLY:
        decision.amr_status = "not_applicable"
        decision.amr_genes = "not_applicable"
        decision.amr_method = "not_applicable"
    # Skipped rows leave all amr_* fields blank/None (= NA in TSV).

    # ---------- 6. Step 5: virulence detection ----------
    virulence_genes_str = ""
    virulence_count: Optional[int] = None
    virulence_sources_str = ""
    virulence_methods_str = ""
    if decision.status is AnalysisStatus.OK_BACTERIAL:
        vresult = detect_virulence(
            decision.ncbi_accession,
            cache=cache,
            genome_path=(download.files.get("genome")
                         if download and download.files else None),
            abricate_client=abricate_client,
        )
        if vresult.status in ("ok", "failed") and vresult.hits:
            genes = vresult.gene_symbols(deduplicate=True)
            virulence_genes_str = "; ".join(genes)
            virulence_count = len(genes)
            virulence_sources_str = "; ".join(vresult.sources_used)
            virulence_methods_str = "; ".join(
                f"{src}@{ver}" for src, ver in vresult.db_versions.items() if ver
            )
        elif vresult.status == "ok" and not vresult.hits:
            virulence_genes_str = ""
            virulence_count = 0
            virulence_sources_str = "; ".join(vresult.sources_used)
        if vresult.error_message:
            decision.hint_warnings.append(vresult.error_message)
    elif decision.status is AnalysisStatus.OK_GENECOUNT_ONLY:
        virulence_genes_str = "not_applicable"

    outcome_final_row = build_final_row(
        decision,
        virulence_genes=virulence_genes_str,
        virulence_count=virulence_count,
        virulence_sources=virulence_sources_str,
        virulence_methods=virulence_methods_str,
    )

    # ---------- 7. Cleanup ----------
    bytes_freed = 0
    if download and download.files and not keep_files:
        bytes_freed = client.cleanup_files(decision.ncbi_accession)

    return RowOutcome(
        decision=decision, fetch=fetch, plan=plan,
        download=download, bytes_freed=bytes_freed,
        final_row=outcome_final_row,
    )


# ===========================================================================
# Run the whole pipeline
# ===========================================================================

def run_pipeline(
    input_path: Path,
    *,
    client: NCBIClient,
    cache: AssemblyCache,
    out_dir: Path,
    download_dir: Path,
    keep_files: bool = False,
    auto_replace_suppressed: bool = False,
    amr_client: Optional[AMRClient] = None,
    abricate_client: Optional[ABRicateClient] = None,
) -> dict:
    """Run Steps 1 → 7 over an input TSV. Returns summary statistics."""
    out_dir.mkdir(parents=True, exist_ok=True)
    download_dir.mkdir(parents=True, exist_ok=True)

    run_id = cache.start_run(notes=f"input={input_path}")
    reporter = DecisionReporter(
        tsv_path=out_dir / "decisions.tsv",
        summary_path=out_dir / "run_report.md",
    )

    total_downloaded_bytes = 0
    total_freed_bytes      = 0
    n_rows = 0
    plan_summary: dict[str, int] = {}
    gene_count_sources: dict[str, int] = {}
    replacements: list[tuple[str, str, str]] = []
    amr_stats: dict[str, int] = {"ok": 0, "failed": 0, "not_applicable": 0,
                                 "not_attempted": 0, "no_amr_found": 0}
    drug_class_tally: dict[str, int] = {}
    amr_rows: list[dict] = []
    virulence_rows: list[dict] = []
    virulence_stats = {"ok_with_hits": 0, "ok_clean": 0,
                       "not_applicable": 0, "failed": 0}
    final_rows: list[FinalRow] = []

    with reporter:
        for row in read_input(input_path):
            n_rows += 1
            outcome = run_row(
                row,
                client=client, cache=cache,
                download_dir=download_dir,
                keep_files=keep_files,
                auto_replace_suppressed=auto_replace_suppressed,
                amr_client=amr_client,
                abricate_client=abricate_client,
            )
            reporter.log(outcome.decision)
            if outcome.final_row:
                final_rows.append(outcome.final_row)

            plan_key = ("no_download" if not outcome.plan.needs_download
                        else "+".join(outcome.plan.includes))
            plan_summary[plan_key] = plan_summary.get(plan_key, 0) + 1
            if outcome.download and outcome.download.files:
                total_downloaded_bytes += outcome.download.bytes_
            total_freed_bytes += outcome.bytes_freed

            src = outcome.decision.gene_count_source or "—"
            gene_count_sources[src] = gene_count_sources.get(src, 0) + 1

            for w in outcome.decision.hint_warnings:
                if w.startswith("Auto-replaced suppressed accession:"):
                    final_acc = w.split("→")[-1].strip()
                    replacements.append((
                        outcome.decision.input_genome_id or "",
                        final_acc,
                        outcome.decision.input_species_name,
                    ))
                    break

            d = outcome.decision
            if d.amr_status:
                if d.amr_status == "ok" and d.amr_gene_count == 0:
                    amr_stats["no_amr_found"] += 1
                else:
                    amr_stats[d.amr_status] = amr_stats.get(d.amr_status, 0) + 1
                if d.amr_status == "ok" and d.amr_drug_classes:
                    for cl in d.amr_drug_classes.split("; "):
                        drug_class_tally[cl] = drug_class_tally.get(cl, 0) + 1
                if d.amr_status in ("ok", "failed"):
                    amr_rows.append({
                        "accession": d.ncbi_accession,
                        "organism": d.ncbi_organism_name,
                        "organism_flag": d.amr_organism_flag,
                        "gene_count": d.amr_gene_count,
                        "genes": d.amr_genes,
                        "classes": d.amr_drug_classes,
                        "status": d.amr_status,
                    })

            # Virulence tally
            fr = outcome.final_row
            if fr:
                if fr.virulence_genes == "not_applicable":
                    virulence_stats["not_applicable"] += 1
                elif fr.virulence_gene_count and fr.virulence_gene_count > 0:
                    virulence_stats["ok_with_hits"] += 1
                    virulence_rows.append({
                        "accession": fr.genome_id,
                        "organism": fr.ncbi_organism_name,
                        "genes": fr.virulence_genes,
                        "count": fr.virulence_gene_count,
                        "sources": fr.virulence_sources,
                    })
                elif d.status is AnalysisStatus.OK_BACTERIAL:
                    virulence_stats["ok_clean"] += 1

    # ---- Step 6: write the final output dataframes ----
    write_minimal_tsv(final_rows, out_dir / "final_output.tsv")
    write_full_tsv(final_rows,    out_dir / "final_output_full.tsv")

    cache.finish_run(run_id, n_rows=n_rows)
    cache_stats = cache.stats()

    _write_full_report(
        summary_path=out_dir / "run_report.md",
        n_rows=n_rows,
        plan_summary=plan_summary,
        gene_count_sources=gene_count_sources,
        replacements=replacements,
        downloaded=total_downloaded_bytes,
        freed=total_freed_bytes,
        cache_stats=cache_stats,
        keep_files=keep_files,
        auto_replace_suppressed=auto_replace_suppressed,
        amr_stats=amr_stats,
        drug_class_tally=drug_class_tally,
        amr_rows=amr_rows,
        virulence_stats=virulence_stats,
        virulence_rows=virulence_rows,
        final_rows=final_rows,
        out_dir=out_dir,
    )

    return {
        "n_rows": n_rows,
        "plan_summary": plan_summary,
        "gene_count_sources": gene_count_sources,
        "replacements": replacements,
        "total_downloaded_bytes": total_downloaded_bytes,
        "total_freed_bytes": total_freed_bytes,
        "cache_stats": cache_stats,
        "amr_stats": amr_stats,
        "drug_class_tally": drug_class_tally,
        "virulence_stats": virulence_stats,
        "final_rows": final_rows,
    }


def _write_full_report(
    summary_path: Path,
    n_rows: int,
    plan_summary: dict[str, int],
    gene_count_sources: dict[str, int],
    replacements: list[tuple[str, str, str]],
    downloaded: int,
    freed: int,
    cache_stats: dict,
    keep_files: bool,
    auto_replace_suppressed: bool,
    amr_stats: dict[str, int],
    drug_class_tally: dict[str, int],
    amr_rows: list[dict],
    virulence_stats: dict[str, int],
    virulence_rows: list[dict],
    final_rows: list[FinalRow],
    out_dir: Path,
) -> None:
    """Step 7: rewrite the markdown report with an executive summary on top.

    Note: this is called *after* DecisionReporter has already written the
    Steps 1-only body. We read that, prepend the executive summary, and
    append the Step 2-6 sections.
    """
    # Capture the body that DecisionReporter already wrote.
    body = summary_path.read_text(encoding="utf-8") if summary_path.exists() else ""

    with summary_path.open("w", encoding="utf-8") as fh:
        w = fh.write

        # ============== Executive summary (Step 7) ==============
        w("# Taxon-aware AMR / virulence pipeline — final report\n\n")
        w("## Executive summary\n\n")
        bact_ok = amr_stats.get("ok", 0)
        bact_none = amr_stats.get("no_amr_found", 0)
        bact_failed = amr_stats.get("failed", 0)
        n_bacterial = bact_ok + bact_none + bact_failed
        n_eukvir = amr_stats.get("not_applicable", 0)
        n_skipped = n_rows - n_bacterial - n_eukvir
        w(f"- **Input rows:** {n_rows}\n")
        w(f"- **Bacterial assemblies analysed:** {n_bacterial} "
          f"({bact_ok} with AMR hits, {bact_none} clean, {bact_failed} failed)\n")
        w(f"- **Non-bacterial assemblies (gene count only):** {n_eukvir}\n")
        w(f"- **Skipped rows** (no genome, suppressed, not found): {n_skipped}\n")
        if drug_class_tally:
            top = sorted(drug_class_tally.items(), key=lambda x: -x[1])[:5]
            w(f"- **Top drug classes detected:** "
              f"{', '.join(f'`{c}` ({n})' for c, n in top)}\n")
        if virulence_stats["ok_with_hits"]:
            w(f"- **Bacterial rows with virulence hits:** "
              f"{virulence_stats['ok_with_hits']}\n")
        if replacements:
            w(f"- **Auto-replaced suppressed accessions:** {len(replacements)}\n")
        w(f"- **Downloaded this run:** {_human_bytes(downloaded)}; "
          f"**freed by cleanup:** {_human_bytes(freed)}\n")
        w(f"- **Cache:** {cache_stats['summaries_ok']} successful summaries cached "
          f"({cache_stats['summaries_failed']} failed, will re-try after TTL)\n\n")
        w("**Deliverables:**\n")
        w(f"- `{out_dir.name}/final_output.tsv` — minimal schema "
          f"(genome_id, total_gene_count, virulence_genes, virulence_gene_count, "
          f"amr_genes, amr_gene_count)\n")
        w(f"- `{out_dir.name}/final_output_full.tsv` — extended schema with "
          f"all audit columns\n")
        w(f"- `{out_dir.name}/decisions.tsv` — per-row decision log\n")
        w(f"- `{out_dir.name}/run_report.md` — this document\n\n")

        # Final dataframe preview, as a markdown table.
        w("### Final output preview\n\n")
        w("| genome_id | species | total_gene_count | amr_gene_count | "
          "virulence_gene_count |\n")
        w("|---|---|---:|---:|---:|\n")
        for r in final_rows:
            gc = "" if r.total_gene_count is None else r.total_gene_count
            ac = "" if r.amr_gene_count is None else r.amr_gene_count
            vc = "" if r.virulence_gene_count is None else r.virulence_gene_count
            ag = "_n/a_" if r.amr_genes == "not_applicable" else ac
            vg = "_n/a_" if r.virulence_genes == "not_applicable" else vc
            w(f"| `{r.genome_id}` | {r.species_name} | {gc} | {ag} | {vg} |\n")
        w("\n---\n\n")

        # ============== The body that DecisionReporter wrote (Step 1) ==========
        w(body)

        # ============== Steps 2-5 sections ==============
        w("\n## Step 2 — retrieval, caching, cleanup\n\n")
        w("### Per-row download plan\n\n")
        w("Files downloaded depend on what each row needs. Bacterial rows need "
          "genome FASTA + protein FASTA + GFF3 for AMRFinderPlus combined mode; "
          "non-bacterial rows whose gene count is already in the NCBI summary "
          "download nothing.\n\n")
        w("| Files retrieved | Rows |\n|---|---:|\n")
        for plan, n in sorted(plan_summary.items(), key=lambda x: -x[1]):
            label = "_(none)_" if plan == "no_download" else f"`{plan}`"
            w(f"| {label} | {n} |\n")
        w("\n### Disk usage\n\n")
        w(f"- Downloaded this run: **{_human_bytes(downloaded)}**\n")
        w(f"- Freed by cleanup: **{_human_bytes(freed)}** "
          f"{'(cleanup disabled with --keep-files)' if keep_files else ''}\n")
        w(f"- Currently on disk: **{_human_bytes(cache_stats['bytes_on_disk'])}** "
          f"across {cache_stats['files_present']} files\n\n")
        w("### Cache health\n\n")
        w(f"- Summaries cached: **{cache_stats['summaries_total']}** "
          f"({cache_stats['summaries_ok']} ok, {cache_stats['summaries_failed']} failed)\n\n")

        if replacements:
            w("### Suppressed-assembly substitutions\n\n")
            if auto_replace_suppressed:
                w("| Input accession | Final accession used | Species |\n|---|---|---|\n")
                for inp, fin, sp in replacements:
                    w(f"| `{inp}` | `{fin}` | {sp} |\n")
                w("\n")

        w("## Step 3 — total gene counts\n\n")
        w("| Source | Rows |\n|---|---:|\n")
        for src, n in sorted(gene_count_sources.items(), key=lambda x: -x[1]):
            label = src if src else "_(not applicable)_"
            w(f"| `{label}` | {n} |\n")
        w("\n- `ncbi_summary`: from cached `annotation_info.stats.gene_counts.total`.\n")
        w("- `gff3_parsed`: GFF3 features of type `gene` or `pseudogene`.\n\n")

        w("## Step 4 — AMR gene detection (AMRFinderPlus)\n\n")
        w(f"- Bacterial rows analysed: **{n_bacterial}**\n")
        w(f"- Non-bacterial rows (AMR not applicable): **{n_eukvir}**\n")
        w("- AMRFinderPlus invoked in **combined mode** (`-n -p -g`) for highest "
          "sensitivity. `--plus` is on, so stress + virulence rows are also in "
          "cache for Step 5.\n\n")
        if amr_rows:
            w("### Per-row AMR profile\n\n")
            w("| Accession | Organism | --organism | Genes | Classes |\n")
            w("|---|---|---|---:|---|\n")
            for r in amr_rows:
                if r["status"] != "ok":
                    continue
                org_flag = r["organism_flag"] if r["organism_flag"] != "auto:none" else "_(none)_"
                classes = r["classes"] or "—"
                w(f"| `{r['accession']}` | *{r['organism']}* | `{org_flag}` | "
                  f"{r['gene_count']} | {classes} |\n")
                if r["genes"]:
                    w(f"|  |  |  |  | _Genes:_ {r['genes']} |\n")
            w("\n")
        if drug_class_tally:
            w("### Drug-class summary\n\n")
            w("| Drug class | Rows hit |\n|---|---:|\n")
            for cl, n in sorted(drug_class_tally.items(), key=lambda x: -x[1]):
                w(f"| `{cl}` | {n} |\n")
            w("\n")

        w("## Step 5 — virulence gene detection\n\n")
        w(f"- Bacterial rows with virulence hits: **{virulence_stats['ok_with_hits']}**\n")
        w(f"- Bacterial rows clean (no virulence): **{virulence_stats['ok_clean']}**\n")
        w(f"- Non-bacterial rows (virulence not applicable): "
          f"**{virulence_stats['not_applicable']}**\n")
        w("- Evidence merged from **two sources**: ABRicate-VFDB (curated virulence "
          "factor database) and AMRFinderPlus `--plus` virulence rows (already in "
          "cache from Step 4). Deduplicated by gene symbol.\n\n")
        if virulence_rows:
            w("### Per-row virulence profile\n\n")
            w("| Accession | Organism | Genes | Count | Sources |\n|---|---|---|---:|---|\n")
            for r in virulence_rows:
                w(f"| `{r['accession']}` | *{r['organism']}* | {r['genes']} | "
                  f"{r['count']} | {r['sources']} |\n")
            w("\n")

def _human_bytes(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"
