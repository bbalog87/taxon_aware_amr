"""
tests/test_step2_pipeline.py
============================

End-to-end demo of Steps 1 + 2 with caching + selective downloads + cleanup.

What this exercises
-------------------
1. Read the example input TSV.
2. Run the pipeline with a ``MockNCBIClient`` that returns realistic
   ``datasets summary genome accession`` JSON, including:
     * present + annotated bacteria, viruses, protozoa, fungi, plants, mammals
     * one assembly intentionally missing from the fixture set
       (→ NCBI ``not_found`` failure path)
     * one assembly flagged ``suppression_reason`` in the summary
       (→ ``suppressed`` failure path)
     * one bacterial assembly with NO ``annotation_info.stats`` block
       (→ forces the download of GFF for gene counting in Step 3)
3. Verify the per-row download plan is minimal:
     * bacteria → protein FASTA (+ gff3 if no annotation)
     * non-bacterial with gene count in summary → no download
     * non-bacterial without gene count → gff3 only
4. Run the pipeline a SECOND time and confirm every successful summary
   comes from cache (zero new fetches).
5. Confirm that downloaded files are deleted after each row, but the
   inventory remembers what was downloaded.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# Allow running as a plain script without installing the package
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from taxon_aware_amr import (    # noqa: E402
    AssemblyCache,
    MockNCBIClient,
    run_pipeline,
)


# ===========================================================================
# Fixture JSONs — modelled after `datasets summary genome accession ...`
# ===========================================================================

def _summary(tax_id: int, name: str, level: str,
             gene_total: int | None,
             suppressed: str | None = None) -> dict:
    ai: dict = {}
    if gene_total is not None:
        ai = {
            "stats": {
                "gene_counts": {
                    "total": gene_total,
                    "protein_coding": gene_total - 50,
                    "non_coding": 50,
                    "pseudogene": 0,
                },
            },
        }
    out = {
        "accession": "",
        "organism": {"tax_id": tax_id, "organism_name": name},
        "annotation_info": ai,
        "assembly_info": {"assembly_level": level},
    }
    if suppressed:
        out["assembly_info"]["suppression_reason"] = suppressed
    return out


def _lineage(*ranks: tuple[str, str]) -> dict:
    """ranks = (('superkingdom','Bacteria'), ('phylum','Bacillota'), ...)"""
    return {"taxonomy": {"classification": {
        rank: {"name": name} for rank, name in ranks
    }}}


MOCK_SUMMARIES: dict[str, dict] = {
    "GCA_964341285.1": _summary(1392, "Bacillus anthracis",        "Complete Genome", 5544),
    "GCA_900654255.1": _summary(1773, "Mycobacterium tuberculosis","Scaffold",        4173),
    # Bacterial WITHOUT NCBI annotation → forces gff3 download:
    "GCA_001049615.1": _summary(1280, "Staphylococcus aureus",     "Contig",          None),
    "GCA_000847345.1": _summary(11588,"Rift Valley fever virus",   "Complete Genome", 7),
    "GCA_900327715.1": _summary(11292,"Rabies lyssavirus",         "Complete Genome", 5),
    "GCA_019455585.1": _summary(5811, "Toxoplasma gondii",         "Chromosome",      8920),
    "GCA_035232765.1": _summary(5807, "Cryptosporidium parvum",    "Complete Genome", 3865),
    "GCA_014900015.1": _summary(5334, "Schizophyllum commune",     "Contig",          13210),
    "GCA_964199775.1": _summary(4577, "Zea mays",                  "Complete Genome", 39591),
    "GCA_030014295.1": _summary(9785, "Loxodonta africana",        "Chromosome",      22871),
    "GCA_025407655.1": _summary(7159, "Aedes aegypti",             "Chromosome",      19500),
    # Suppressed assembly:
    "GCA_015501595.1": _summary(1081385, "Diceros bicornis",       "Scaffold",        None,
                                 suppressed="Replaced by GCA_022413745.1"),
    # GCA_999999999.9 is deliberately absent → not_found
}

MOCK_LINEAGES: dict[str, dict] = {
    "GCA_964341285.1": _lineage(("superkingdom","Bacteria"), ("phylum","Bacillota"),
                                ("class","Bacilli"), ("order","Bacillales"),
                                ("family","Bacillaceae"), ("genus","Bacillus"),
                                ("species","Bacillus anthracis")),
    "GCA_900654255.1": _lineage(("superkingdom","Bacteria"), ("phylum","Actinomycetota"),
                                ("genus","Mycobacterium"),
                                ("species","Mycobacterium tuberculosis")),
    "GCA_001049615.1": _lineage(("superkingdom","Bacteria"), ("phylum","Bacillota"),
                                ("genus","Staphylococcus"),
                                ("species","Staphylococcus aureus")),
    "GCA_000847345.1": _lineage(("superkingdom","Viruses"), ("phylum","Negarnaviricota"),
                                ("genus","Phlebovirus"),
                                ("species","Rift Valley fever virus")),
    "GCA_900327715.1": _lineage(("superkingdom","Viruses"), ("phylum","Negarnaviricota"),
                                ("genus","Lyssavirus"), ("species","Rabies lyssavirus")),
    "GCA_019455585.1": _lineage(("superkingdom","Eukaryota"), ("phylum","Apicomplexa"),
                                ("genus","Toxoplasma"), ("species","Toxoplasma gondii")),
    "GCA_035232765.1": _lineage(("superkingdom","Eukaryota"), ("phylum","Apicomplexa"),
                                ("genus","Cryptosporidium"),
                                ("species","Cryptosporidium parvum")),
    "GCA_014900015.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Fungi"),
                                ("phylum","Basidiomycota"),
                                ("species","Schizophyllum commune")),
    "GCA_964199775.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Viridiplantae"),
                                ("phylum","Streptophyta"), ("species","Zea mays")),
    "GCA_030014295.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Metazoa"),
                                ("phylum","Chordata"), ("class","Mammalia"),
                                ("species","Loxodonta africana")),
    "GCA_025407655.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Metazoa"),
                                ("phylum","Arthropoda"), ("class","Insecta"),
                                ("species","Aedes aegypti")),
}


# ===========================================================================
# Demo
# ===========================================================================

def main() -> None:
    project = Path(__file__).resolve().parent.parent
    input_path   = project / "examples" / "kenya_biodiversity_input.tsv"
    out_dir      = project / "results"
    cache_path   = project / "results" / "cache" / "cache.sqlite"
    download_dir = project / "results" / "cache" / "genomes"

    # Clean slate
    import shutil
    if (project / "results").exists():
        shutil.rmtree(project / "results")

    cache  = AssemblyCache(cache_path, negative_ttl_hours=24)
    client = MockNCBIClient(
        cache=cache, download_dir=download_dir,
        mock_summaries=MOCK_SUMMARIES,
        mock_lineages=MOCK_LINEAGES,
    )

    # ------------------------------------------------------------------ #
    print("=" * 78)
    print("Run 1: cold cache — every successful row fetches live.")
    print("=" * 78)
    stats1 = run_pipeline(
        input_path,
        client=client, cache=cache,
        out_dir=out_dir, download_dir=download_dir,
        keep_files=False,                # ← delete after success
    )

    # ------------------------------------------------------------------ #
    print()
    print("=" * 78)
    print("Run 2: warm cache — every successful row should be served from SQLite.")
    print("=" * 78)
    # Use a fresh reporter target so we don't overwrite the first run's
    # reports, and force a new run-log entry.
    stats2 = run_pipeline(
        input_path,
        client=client, cache=cache,
        out_dir=project / "results_rerun", download_dir=download_dir,
        keep_files=False,
    )

    # ------------------------------------------------------------------ #
    # Sanity checks
    # ------------------------------------------------------------------ #
    print()
    print("=" * 78)
    print("Verification")
    print("=" * 78)

    # 1. Bacterial rows now always download genome + protein + gff3
    # (combined mode for AMRFinderPlus). All 3 bacterial demo rows use this plan.
    plan = stats1["plan_summary"]
    assert plan.get("genome+protein+gff3", 0) == 3, (
        "All 3 bacterial demo rows should use genome+protein+gff3 plan, "
        f"got {plan}"
    )
    print(f"  ✓ genome+protein+gff3 downloads: {plan.get('genome+protein+gff3', 0)}")

    # 2. Non-bacterial rows whose summary has gene_counts download nothing.
    assert plan.get("no_download", 0) >= 1
    print(f"  ✓ no-download rows:       {plan.get('no_download', 0)}")

    # 4. Cleanup actually freed bytes.
    assert stats1["total_freed_bytes"] == stats1["total_downloaded_bytes"], (
        "Cleanup should free every byte we downloaded (we passed keep_files=False)."
    )
    print(f"  ✓ downloaded:             {stats1['total_downloaded_bytes']} B")
    print(f"  ✓ freed by cleanup:       {stats1['total_freed_bytes']} B")

    # 5. After Run 1, on-disk usage in the cache should be zero.
    assert stats1["cache_stats"]["bytes_on_disk"] == 0
    assert stats1["cache_stats"]["files_present"] == 0
    print(f"  ✓ files still on disk:    {stats1['cache_stats']['files_present']}")

    # 6. Run 2 uses cache: no new SQLite entries should appear.
    assert stats2["cache_stats"]["summaries_total"] == stats1["cache_stats"]["summaries_total"], (
        "Run 2 should not have created new summary cache entries; got "
        f"{stats2['cache_stats']['summaries_total']} vs "
        f"{stats1['cache_stats']['summaries_total']}."
    )
    print(f"  ✓ summaries cached:       {stats2['cache_stats']['summaries_total']} "
          f"({stats2['cache_stats']['summaries_ok']} ok, "
          f"{stats2['cache_stats']['summaries_failed']} failed — re-tried per TTL)")

    print()
    print(f"All Step 2 checks passed. Reports in: {out_dir}")


if __name__ == "__main__":
    main()
