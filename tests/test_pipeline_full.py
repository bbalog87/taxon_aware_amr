"""
tests/test_pipeline_full.py
===========================

Master end-to-end test for the v1.0.0 release.

Exercises Steps 1 → 7 in a single run using only mock clients (no NCBI,
no AMRFinderPlus, no ABRicate required). Validates:

    * Routing decisions (Step 1)
    * Selective downloads with the new combined-mode plan (Step 2)
    * Gene counts from summary + GFF (Step 3)
    * AMR detection in combined mode with --organism (Step 4)
    * Virulence detection merged from VFDB + AMRFinderPlus --plus (Step 5)
    * Final output dataframes — minimal + full schemas (Step 6)
    * Executive summary atop run_report.md (Step 7)
"""

from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from taxon_aware_amr import (   # noqa: E402
    AssemblyCache,
    MockAMRClient,
    MockABRicateClient,
    MockNCBIClient,
    run_pipeline,
    MINIMAL_COLUMNS,
)
from taxon_aware_amr.fixtures import (    # noqa: E402
    DEMO_SUMMARIES, DEMO_LINEAGES,
    demo_amr_hits_as_records, demo_vfdb_hits_as_records,
)


def _patch_real_gff_and_genome(mock_client: MockNCBIClient) -> None:
    """Make the mock write real-looking files so Step 3 and Step 4 paths work."""
    original_dl = mock_client.download_files

    def dl(accession, includes, out_dir_):
        result = original_dl(accession, includes, out_dir_)
        if "gff3" in result.files:
            n = {"GCA_001049615.1": 2872}.get(accession, 5000)
            with open(result.files["gff3"], "w") as fh:
                fh.write("##gff-version 3\n")
                for i in range(n):
                    fh.write(f"chr\tRefSeq\tgene\t{i*100+1}\t{i*100+99}\t.\t+\t.\tID=g{i}\n")
            result.bytes_ = result.files["gff3"].stat().st_size
        if "genome" in result.files:
            result.files["genome"].write_text(">chr1\nACGTACGTACGT\n")
        if "protein" in result.files:
            result.files["protein"].write_text(">p1\nMACG\n")
        return result

    mock_client.download_files = dl  # type: ignore[assignment]


def main():
    project = Path(__file__).resolve().parent.parent
    input_path = project / "examples" / "kenya_biodiversity_input.tsv"
    out_dir    = project / "results_final"
    cache_dir  = project / "results_final" / "cache"
    download_dir = cache_dir / "genomes"

    if out_dir.exists():
        shutil.rmtree(out_dir)

    print("=" * 78)
    print("End-to-end v1.0.0 pipeline (Steps 1–7)")
    print("=" * 78)

    cache = AssemblyCache(cache_dir / "cache.sqlite")
    ncbi = MockNCBIClient(
        cache=cache, download_dir=download_dir,
        mock_summaries=DEMO_SUMMARIES, mock_lineages=DEMO_LINEAGES,
    )
    _patch_real_gff_and_genome(ncbi)
    amr  = MockAMRClient(cache=cache, mock_hits=demo_amr_hits_as_records())
    vfdb = MockABRicateClient(cache=cache, mock_hits=demo_vfdb_hits_as_records())

    stats = run_pipeline(
        input_path,
        client=ncbi, cache=cache,
        out_dir=out_dir, download_dir=download_dir,
        keep_files=False,
        auto_replace_suppressed=True,
        amr_client=amr,
        abricate_client=vfdb,
    )

    print()
    print("=" * 78)
    print("Final output dataframe (minimal schema)")
    print("=" * 78)

    # ---- Read the minimal final TSV ----
    with (out_dir / "final_output.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert list(reader.fieldnames) == list(MINIMAL_COLUMNS), \
            f"unexpected columns: {reader.fieldnames}"
        final_rows = list(reader)
    assert len(final_rows) == 16, f"expected 16 rows, got {len(final_rows)}"

    # Render a compact view
    print(f"{'genome_id':<22} {'genes':>6} {'amr':>4} {'vir':>4} species")
    print("-" * 78)
    by_acc = {}
    for r in final_rows:
        gid = r["genome_id"]
        by_acc[gid] = r
        sp = ""  # not in minimal schema; just for display
        gc = r["total_gene_count"] or "—"
        ac = r["amr_gene_count"] or ("n/a" if r["amr_genes"]=="not_applicable" else "—")
        vc = r["virulence_gene_count"] or ("n/a" if r["virulence_genes"]=="not_applicable" else "—")
        print(f"{gid:<22} {gc:>6} {ac:>4} {vc:>4}")

    print()
    print("=" * 78)
    print("Verification")
    print("=" * 78)

    # ---- Verify bacterial rows have BOTH AMR and virulence hits ----
    ban = by_acc["GCA_964341285.1"]   # B. anthracis
    assert ban["amr_gene_count"] == "2"
    assert ban["virulence_gene_count"] == "5"   # pagA, lef, cya, capB, capC
    assert "pagA" in ban["virulence_genes"]
    print("  ✓ B. anthracis: 2 AMR + 5 VFDB virulence genes (anthrax toxin + capsule)")

    mtb = by_acc["GCA_900654255.1"]   # M. tuberculosis
    assert mtb["amr_gene_count"] == "5"
    assert mtb["virulence_gene_count"] == "4"   # esxA, esxB, eccCb1, mmpL7
    print(f"  ✓ M. tuberculosis: 5 AMR + 4 VFDB virulence genes (ESX-1 secretion)")

    mrsa = by_acc["GCA_001049615.1"]  # S. aureus
    assert mrsa["amr_gene_count"] == "6"
    # MRSA has 8 VFDB hits (clfA, clfB, spa, hla, hlb, sak, sea, icaA)
    # plus AMRFinder --plus contributed `hlb` — dedup means hlb appears once.
    # So unique virulence count is 8.
    assert mrsa["virulence_gene_count"] == "8", \
        f"expected 8, got {mrsa['virulence_gene_count']}"
    print(f"  ✓ MRSA: 6 AMR + 8 unique virulence genes "
          f"(VFDB ∪ AMRFinder --plus, deduplicated on `hlb`)")

    # ---- Verify non-bacterial rows are "not_applicable" ----
    zea = by_acc["GCA_964199775.1"]   # Zea mays
    assert zea["amr_genes"] == "not_applicable"
    assert zea["virulence_genes"] == "not_applicable"
    assert zea["total_gene_count"] == "39591"
    print("  ✓ Zea mays: AMR + virulence not_applicable, gene count from NCBI summary")

    # ---- Verify auto-replaced row ----
    diceros = by_acc["GCA_022413745.1"]  # the replacement of GCA_015501595.1
    assert diceros["total_gene_count"] == "21845"
    assert diceros["amr_genes"] == "not_applicable"
    print("  ✓ Diceros bicornis: auto-replaced GCA_015501595.1 → GCA_022413745.1")

    # ---- Verify skipped rows ----
    skipped = [r for r in final_rows if r["genome_id"] == "none"
               or r["genome_id"] == "GCA_999999999.9"]
    for r in skipped:
        assert r["total_gene_count"] == ""
        assert r["amr_genes"] in ("NA", "")
    print(f"  ✓ {len(skipped)} skipped rows have empty gene counts and NA for AMR/vir")

    # ---- Verify the FULL TSV has the audit columns ----
    with (out_dir / "final_output_full.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames
        for must_have in ("input_genome_id", "input_taxon_id", "species_name",
                          "taxon_group", "ncbi_organism_name",
                          "amr_drug_classes", "virulence_sources",
                          "analysis_status", "reason", "hint_warnings"):
            assert must_have in fields, f"missing column {must_have}"
    print(f"  ✓ final_output_full.tsv has all {len(fields)} columns including audit")

    # ---- Verify executive summary in run_report.md ----
    report = (out_dir / "run_report.md").read_text()
    assert "## Executive summary" in report
    assert "### Final output preview" in report
    assert "## Step 5 — virulence" in report
    print("  ✓ run_report.md has executive summary + Step 5 section")

    # ---- Verify pipeline-level stats ----
    assert stats["amr_stats"]["ok"] == 3
    assert stats["virulence_stats"]["ok_with_hits"] == 3
    print(f"  ✓ pipeline: 3 AMR-positive bacterial rows, "
          f"3 virulence-positive bacterial rows")

    print()
    print(f"All v1.0.0 checks passed. Final deliverables in: {out_dir}")
    print(f"  - {out_dir.name}/final_output.tsv         (minimal schema)")
    print(f"  - {out_dir.name}/final_output_full.tsv    (extended)")
    print(f"  - {out_dir.name}/decisions.tsv            (per-row decisions)")
    print(f"  - {out_dir.name}/run_report.md            (executive summary)")


if __name__ == "__main__":
    main()
