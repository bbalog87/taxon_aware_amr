"""
tests/test_step3_pipeline.py
============================

Validates Step 3 features:

    1. ``total_gene_count`` populated from the NCBI summary when available
       (zero downloads — these rows hit `gene_count_source = "ncbi_summary"`).
    2. ``total_gene_count`` from a parsed GFF3 fallback when the summary
       lacks counts (Staphylococcus aureus in the fixtures).
    3. ``--auto-replace-suppressed`` chases "Replaced by GCA_xxx" messages
       and substitutes the replacement assembly, logging the substitution.
"""

from __future__ import annotations

import shutil
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from taxon_aware_amr import (        # noqa: E402
    AnalysisStatus,
    AssemblyCache,
    MockNCBIClient,
    count_genes_from_gff,
    extract_replacement_accession,
    run_pipeline,
)
from taxon_aware_amr.fixtures import DEMO_SUMMARIES, DEMO_LINEAGES   # noqa: E402


def _make_fake_gff(path: Path, n_genes: int, n_pseudogenes: int = 0) -> None:
    """Write a tiny but valid GFF3 with N gene rows."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["##gff-version 3"]
    for i in range(n_genes):
        lines.append(f"chr1\tRefSeq\tgene\t{i*100+1}\t{i*100+99}\t.\t+\t.\tID=gene{i};Name=g{i}")
    for i in range(n_pseudogenes):
        lines.append(f"chr1\tRefSeq\tpseudogene\t{i*100+5001}\t{i*100+5099}\t.\t+\t.\tID=pseudo{i}")
    # Throw in some non-gene rows to make sure they're skipped:
    for i in range(3):
        lines.append(f"chr1\tRefSeq\tCDS\t{i*100+1}\t{i*100+99}\t.\t+\t0\tID=cds{i}")
    path.write_text("\n".join(lines) + "\n")


def test_gff_parser():
    tmp = Path(__file__).resolve().parent / "_tmp_test_step3"
    tmp.mkdir(parents=True, exist_ok=True)
    gff = tmp / "demo.gff"
    _make_fake_gff(gff, n_genes=2700, n_pseudogenes=50)
    n = count_genes_from_gff(gff)
    assert n == 2750, f"expected 2750 gene+pseudogene, got {n}"
    shutil.rmtree(tmp)
    print(f"  ✓ GFF parser counts gene+pseudogene: {n}")


def test_replacement_regex():
    assert extract_replacement_accession(
        "Suppressed by NCBI: Replaced by GCA_022413745.1") == "GCA_022413745.1"
    assert extract_replacement_accession(
        "Replaced by accession GCF_000123456.2") == "GCF_000123456.2"
    assert extract_replacement_accession(
        "Some unrelated message") is None
    assert extract_replacement_accession(None) is None
    print("  ✓ replacement-accession regex handles all known message shapes")


def main():
    project = Path(__file__).resolve().parent.parent
    input_path = project / "examples" / "kenya_biodiversity_input.tsv"
    out_dir    = project / "results_step3"
    cache_dir  = project / "results_step3" / "cache"
    download_dir = cache_dir / "genomes"

    if out_dir.exists():
        shutil.rmtree(out_dir)

    print("=" * 78)
    print("Step 3 unit checks")
    print("=" * 78)
    test_gff_parser()
    test_replacement_regex()

    print()
    print("=" * 78)
    print("Step 3 end-to-end run (auto-replace ON, GFF parsed where needed)")
    print("=" * 78)

    cache  = AssemblyCache(cache_dir / "cache.sqlite")
    client = MockNCBIClient(
        cache=cache, download_dir=download_dir,
        mock_summaries=DEMO_SUMMARIES,
        mock_lineages=DEMO_LINEAGES,
    )

    # Inject realistic GFFs into the mock download flow:
    # the mock writes placeholder files, but we want the GFF to look real
    # so count_genes_from_gff can parse it. Patch in a hook.
    original_dl = client.download_files

    def download_with_real_gff(accession, includes, out_dir_):
        result = original_dl(accession, includes, out_dir_)
        if "gff3" in result.files:
            # Overwrite the placeholder with a small but valid GFF.
            # Use a per-accession synthetic count so the test can verify it.
            synthetic_count = {"GCA_001049615.1": 2872}.get(accession, 5000)
            _make_fake_gff(result.files["gff3"], synthetic_count)
            result.bytes_ = result.files["gff3"].stat().st_size
        return result

    client.download_files = download_with_real_gff  # type: ignore[assignment]

    stats = run_pipeline(
        input_path,
        client=client, cache=cache,
        out_dir=out_dir, download_dir=download_dir,
        keep_files=False,
        auto_replace_suppressed=True,
    )

    # Read the resulting TSV back to verify the gene counts
    import csv
    rows = {}
    with (out_dir / "decisions.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["input_species_name"]] = r

    print()
    print("Verification")
    print("-" * 78)

    # ---- Gene counts from NCBI summary (no download) ----
    assert rows["Bacillus_anthracis"]["gene_count_source"] == "ncbi_summary"
    assert rows["Bacillus_anthracis"]["total_gene_count"] == "5544"
    print(f"  ✓ Bacillus anthracis: count from summary  → "
          f"{rows['Bacillus_anthracis']['total_gene_count']}")

    assert rows["Zea_mays"]["gene_count_source"] == "ncbi_summary"
    assert rows["Zea_mays"]["total_gene_count"] == "39591"
    print(f"  ✓ Zea mays: count from summary            → "
          f"{rows['Zea_mays']['total_gene_count']}")

    # ---- Gene count from GFF3 (bacterial assembly without NCBI annotation) ----
    assert rows["Staphylococcus_aureus"]["gene_count_source"] == "gff3_parsed"
    assert rows["Staphylococcus_aureus"]["total_gene_count"] == "2872"
    print(f"  ✓ Staphylococcus aureus: count from GFF3   → "
          f"{rows['Staphylococcus_aureus']['total_gene_count']} (synthetic)")

    # ---- Skipped row: no gene count ----
    assert rows["Cercocebus_galeritus"]["total_gene_count"] == ""
    assert rows["Cercocebus_galeritus"]["gene_count_source"] == ""
    print(f"  ✓ Cercocebus galeritus (no genome): count = NA")

    # ---- Auto-replaced row ----
    diceros = rows["Diceros_bicornis"]
    # routing should now succeed via the replacement
    assert diceros["status"] == AnalysisStatus.OK_GENECOUNT_ONLY.value
    assert diceros["ncbi_accession"] == "GCA_022413745.1"   # the replacement
    assert diceros["input_genome_id"] == "GCA_015501595.1"  # original preserved
    assert "Auto-replaced" in diceros["hint_warnings"]
    print(f"  ✓ Diceros bicornis auto-replaced: "
          f"{diceros['input_genome_id']} → {diceros['ncbi_accession']} "
          f"(count {diceros['total_gene_count']})")

    # Pipeline-level stats
    assert len(stats["replacements"]) == 1
    print(f"  ✓ pipeline recorded {len(stats['replacements'])} substitution")

    print()
    print(f"All Step 3 checks passed. Reports in: {out_dir}")


if __name__ == "__main__":
    main()
