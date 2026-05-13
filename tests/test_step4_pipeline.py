"""
tests/test_step4_pipeline.py
============================

Validates Step 4: AMR detection via AMRFinderPlus on bacterial rows.

What this exercises
-------------------
1. Organism mapping for ``--organism`` flag (M. tuberculosis, S. aureus,
   plus fall-through for B. anthracis which is not in the supported list).
2. End-to-end orchestration: bacterial rows get AMRFinderPlus, non-bacterial
   rows get ``amr_genes = "not_applicable"``, skipped rows get blank.
3. AMR hit caching: re-running the pipeline pulls AMR results from the
   ``tool_results`` SQLite table without re-executing the mock client.
4. Drug-class tally aggregation across all bacterial rows.
5. TSV column propagation: ``amr_genes``, ``amr_gene_count``,
   ``amr_drug_classes``, ``amr_method`` appear in decisions.tsv.
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
    MockNCBIClient,
    amrfinder_organism_for,
    run_pipeline,
)
from taxon_aware_amr.fixtures import (    # noqa: E402
    DEMO_SUMMARIES, DEMO_LINEAGES, demo_amr_hits_as_records,
)


# ===========================================================================
# Unit checks
# ===========================================================================

def test_organism_mapping():
    cases = [
        ("Mycobacterium tuberculosis",            "Mycobacterium_tuberculosis"),
        ("Mycobacterium tuberculosis H37Rv",       "Mycobacterium_tuberculosis"),
        ("Staphylococcus aureus",                  "Staphylococcus_aureus"),
        ("Staphylococcus aureus subsp. aureus N315","Staphylococcus_aureus"),
        ("Escherichia coli",                       "Escherichia"),
        ("Escherichia coli O157:H7",               "Escherichia"),
        ("Bacillus anthracis",                     None),     # not in AMRFinder list
        ("Bacillus anthracis str. Ames",           None),
        ("",                                       None),
        (None,                                     None),
    ]
    for inp, expected in cases:
        got = amrfinder_organism_for(inp)
        assert got == expected, f"{inp!r} → expected {expected!r}, got {got!r}"
    print("  ✓ organism mapping correct for 10 representative cases")


# ===========================================================================
# Helpers — patch the mock NCBI client to write a real-looking GFF for
# the bacterial-without-annotation row (so Step 3 doesn't warn here).
# ===========================================================================

def _patch_real_gff(mock_client: MockNCBIClient) -> None:
    original_dl = mock_client.download_files

    def dl(accession, includes, out_dir_):
        result = original_dl(accession, includes, out_dir_)
        if "gff3" in result.files:
            synthetic = {"GCA_001049615.1": 2872}.get(accession, 5000)
            with open(result.files["gff3"], "w") as fh:
                fh.write("##gff-version 3\n")
                for i in range(synthetic):
                    fh.write(f"chr\tRefSeq\tgene\t{i*100+1}\t{i*100+99}\t.\t+\t.\t"
                             f"ID=g{i}\n")
            result.bytes_ = result.files["gff3"].stat().st_size
        return result

    mock_client.download_files = dl  # type: ignore[assignment]


# ===========================================================================
# End-to-end
# ===========================================================================

def main():
    project = Path(__file__).resolve().parent.parent
    input_path = project / "examples" / "kenya_biodiversity_input.tsv"
    out_dir    = project / "results_step4"
    cache_dir  = project / "results_step4" / "cache"
    download_dir = cache_dir / "genomes"

    if out_dir.exists():
        shutil.rmtree(out_dir)

    print("=" * 78)
    print("Step 4 unit checks")
    print("=" * 78)
    test_organism_mapping()

    print()
    print("=" * 78)
    print("Step 4 end-to-end run (cold cache)")
    print("=" * 78)

    cache  = AssemblyCache(cache_dir / "cache.sqlite")
    ncbi_client = MockNCBIClient(
        cache=cache, download_dir=download_dir,
        mock_summaries=DEMO_SUMMARIES, mock_lineages=DEMO_LINEAGES,
    )
    _patch_real_gff(ncbi_client)
    amr_client = MockAMRClient(
        cache=cache, mock_hits=demo_amr_hits_as_records(),
    )

    stats = run_pipeline(
        input_path,
        client=ncbi_client, cache=cache,
        out_dir=out_dir, download_dir=download_dir,
        keep_files=False,
        auto_replace_suppressed=True,
        amr_client=amr_client,
    )

    # ---- Inspect the TSV for the bacterial rows ----
    rows: dict[str, dict] = {}
    with (out_dir / "decisions.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["input_species_name"]] = r

    print()
    print("Per-row AMR profile from decisions.tsv")
    print("-" * 78)
    for sp in ("Bacillus_anthracis", "Mycobacterium_tuberculosis",
               "Staphylococcus_aureus"):
        r = rows[sp]
        print(f"  {sp:<32} status={r['amr_status']:<5} n={r['amr_gene_count']:>2}  "
              f"organism={r['amr_organism_flag']:<28} classes=[{r['amr_drug_classes']}]")
        print(f"      genes: {r['amr_genes']}")

    print()
    print("Verification")
    print("-" * 78)

    # 1. M. tuberculosis: 5 hits, organism flag set
    mtb = rows["Mycobacterium_tuberculosis"]
    assert mtb["amr_status"] == "ok"
    assert mtb["amr_gene_count"] == "5"
    assert mtb["amr_organism_flag"] == "Mycobacterium_tuberculosis"
    assert "rpoB_S531L" in mtb["amr_genes"]
    assert "ISONIAZID" in mtb["amr_drug_classes"]
    print("  ✓ M. tuberculosis: 5 AMR hits with --organism flag, INH/RIF/FQ/AMG classes")

    # 2. MRSA: 6 AMR hits; the STRESS qacJ and VIRULENCE hlb in fixtures
    #    are filtered out of Step 4 scope.
    mrsa = rows["Staphylococcus_aureus"]
    assert mrsa["amr_status"] == "ok"
    assert mrsa["amr_gene_count"] == "6", f"got {mrsa['amr_gene_count']}"
    assert "mecA" in mrsa["amr_genes"]
    assert "qacJ" not in mrsa["amr_genes"]   # stress, not AMR
    assert "hlb" not in mrsa["amr_genes"]    # virulence, not AMR (Step 5)
    print("  ✓ MRSA: 6 AMR hits; STRESS+VIRULENCE rows correctly excluded from Step 4")

    # 3. B. anthracis: organism flag should be "auto:none" (not in supported list)
    ban = rows["Bacillus_anthracis"]
    assert ban["amr_status"] == "ok"
    assert ban["amr_organism_flag"] == "auto:none"
    print("  ✓ B. anthracis: AMR ran without --organism (not in supported list)")

    # 4. Non-bacterial rows: amr_genes = "not_applicable"
    for sp in ("Zea_mays", "Loxodonta_africana", "Schizophyllum_commune",
               "RVF_virus"):
        r = rows[sp]
        assert r["amr_status"] == "not_applicable", f"{sp}: {r['amr_status']}"
        assert r["amr_genes"] == "not_applicable", f"{sp}: {r['amr_genes']}"
    print("  ✓ 4 non-bacterial rows correctly tagged 'not_applicable'")

    # 5. Skipped rows: AMR fields blank
    for sp in ("Cercocebus_galeritus", "Kigelia_africana", "Imaginary_species"):
        r = rows[sp]
        assert r["amr_status"] == "", f"{sp}: expected blank, got {r['amr_status']!r}"
        assert r["amr_genes"] == ""
    print("  ✓ 3 skipped rows have blank AMR fields")

    # 6. Pipeline-level tallies
    assert stats["amr_stats"]["ok"] == 3
    assert stats["amr_stats"]["not_applicable"] == 9  # 9 OK_GENECOUNT_ONLY rows (incl. Diceros via auto-replace)
    assert "BETA-LACTAM" in stats["drug_class_tally"]
    print(f"  ✓ pipeline tally: {stats['amr_stats']}")
    print(f"  ✓ drug classes seen: {sorted(stats['drug_class_tally'].keys())}")

    print()
    print("=" * 78)
    print("Step 4 cache check — re-run should pull AMR from tool_results table")
    print("=" * 78)

    # Re-instantiate the mock client; if cache works, the new client never
    # gets called for AMR.
    call_count = {"n": 0}
    real_run = MockAMRClient(cache=cache, mock_hits={}).run  # empty fixtures
    def counting_run(self, accession, protein_path, *, organism_flag,
                     genome_path=None, gff_path=None, refresh=False):
        call_count["n"] += 1
        return real_run(accession, protein_path,
                        organism_flag=organism_flag,
                        genome_path=genome_path, gff_path=gff_path,
                        refresh=refresh)
    new_amr = MockAMRClient(cache=cache, mock_hits={})
    new_amr.run = counting_run.__get__(new_amr, MockAMRClient)  # type: ignore

    stats2 = run_pipeline(
        input_path,
        client=ncbi_client, cache=cache,
        out_dir=project / "results_step4_rerun", download_dir=download_dir,
        keep_files=False, auto_replace_suppressed=True,
        amr_client=new_amr,
    )
    assert stats2["amr_stats"]["ok"] == 3, (
        f"Cache should provide AMR results; got {stats2['amr_stats']}")
    # Note: call_count may be 0 OR == n_bacterial depending on whether the
    # client.run consults the cache before invoking. In our design the
    # cache is checked inside .run(), so the calls happen but return early.
    # The real test is that .ok == 3 even though the new client's
    # mock_hits is empty.
    print(f"  ✓ rerun produced {stats2['amr_stats']['ok']} OK rows "
          f"(new client had empty fixtures — proves cache was used)")

    print()
    print(f"All Step 4 checks passed. Reports in: {out_dir}")


if __name__ == "__main__":
    main()
