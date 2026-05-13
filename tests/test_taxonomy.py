"""
tests/test_taxonomy.py
======================

End-to-end demo of Step 1's GCA-centric architecture.

For each input row we:

    1. Normalise the genome_id (cleaning 'none' → None).
    2. Look up the AssemblyRecord by GCA accession (mock in Step 1).
    3. Call ``decide(...)`` — pure, deterministic.
    4. Pass the decision through a ``DecisionReporter`` that
       - streams to stdout,
       - writes ``results/decisions.tsv``,
       - writes ``results/run_report.md``.

The input cases cover every routing branch:

    * Bacteria             → OK_BACTERIAL
    * Virus (same input category as bacteria — proves GCA-routing works)
    * Protozoa, Fungus, Plant, Mammal, Insect → OK_GENECOUNT_ONLY
    * genome_id == 'none'  → SKIPPED_NO_GENOME
    * GCA not in NCBI      → SKIPPED_ASSEMBLY_NOT_FOUND
    * Audit warnings on input rows where taxon_id has a '_ENA' suffix or
      where the input category disagrees with NCBI lineage.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running as a plain script without installing the package
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from taxon_aware_amr import (              # noqa: E402
    AnalysisStatus,
    DecisionReporter,
    TaxonGroup,
    decide,
    mock_assembly_record,
    normalise_genome_id,
)


# ---------------------------------------------------------------------------
# Input cases — lifted verbatim from the Kenya surveillance table.
# Format: (input_taxon_id, input_genome_id, input_species_name, input_category)
# ---------------------------------------------------------------------------

CASES: list[tuple[str, str, str, str]] = [
    # --- Bacteria (the only rows where AMR is biologically valid) ---
    ("1392_ENA",   "GCA_964341285.1", "Bacillus_anthracis",         "Non_eukaryotic_microbe"),
    ("1773_ENA",   "GCA_900654255.1", "Mycobacterium_tuberculosis", "Non_eukaryotic_microbe"),
    ("1280_ENA",   "GCA_001049615.1", "Staphylococcus_aureus",      "Non_eukaryotic_microbe"),

    # --- Viruses (same input category as bacteria — GCA routing disambiguates) ---
    ("11588_ENA",  "GCA_000847345.1", "RVF_virus",                  "Non_eukaryotic_microbe"),
    ("11292_ENA",  "GCA_900327715.1", "Rabies_lyssavirus",          "Non_eukaryotic_microbe"),

    # --- Protozoa ---
    ("5811",       "GCA_019455585.1", "Toxoplasma_gondii",          "Protozoa"),
    ("5807",       "GCA_035232765.1", "Cryptosporidium_parvum",     "Protozoa"),

    # --- Fungus ---
    ("5334",       "GCA_014900015.1", "Schizophyllum commune",      "Fungus"),

    # --- Crop plant ---
    ("4577",       "GCA_964199775.1", "Zea_mays",                   "Crop"),

    # --- Wildlife mammal ---
    ("9785",       "GCA_030014295.1", "Loxodonta africana",         "Wildlife"),

    # --- Disease vector (mosquito) ---
    ("7159",       "GCA_025407655.1", "Aedes_aegypti",              "Insect"),

    # --- No genome_id ---
    ("9532",       "none",            "Cercocebus galeritus",       "Wildlife"),
    ("70070",      "none",            "Kigelia africana",           "Wild-plant"),
    ("none_1",     "none",            "Clavipitaceous Fungi",       "Fungus"),

    # --- GCA given but not in NCBI (simulated by an unknown accession) ---
    ("9999",       "GCA_999999999.9", "Imaginary species",          "Wildlife"),

    # --- Audit-warning case: input category disagrees with NCBI ---
    # We label Aedes aegypti as 'Wildlife' here to demonstrate the
    # warning, even though Wildlife happens to map to METAZOA correctly;
    # to trigger a real disagreement we mislabel it as 'Crop'.
    ("7159",       "GCA_025407655.1", "Aedes_aegypti (mislabel)",   "Crop"),
]


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run(out_dir: Path = Path("results")) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    reporter = DecisionReporter(
        tsv_path=out_dir / "decisions.tsv",
        summary_path=out_dir / "run_report.md",
    )

    with reporter:
        for input_taxon_id, raw_genome_id, species, category in CASES:
            clean_id = normalise_genome_id(raw_genome_id)
            assembly = mock_assembly_record(clean_id) if clean_id else None
            decision = decide(
                assembly=assembly,
                input_genome_id=clean_id,
                input_taxon_id=input_taxon_id,
                input_species_name=species,
                input_category=category,
            )
            reporter.log(decision)

    # ---------------------------- assertions ---------------------------- #
    decisions = {d.input_species_name: d for d in reporter.decisions}

    # Bacteria → AMR applicable
    assert decisions["Bacillus_anthracis"].status is AnalysisStatus.OK_BACTERIAL
    assert decisions["Bacillus_anthracis"].amr_applicable is True
    assert decisions["Bacillus_anthracis"].ncbi_tax_id == 1392
    # Audit: '1392_ENA' parses to 1392, matches NCBI → no warning
    assert decisions["Bacillus_anthracis"].taxon_id_agrees is True
    assert decisions["Bacillus_anthracis"].hint_warnings == []

    # Viruses (same input category as bacteria) → NOT applicable
    assert decisions["RVF_virus"].status is AnalysisStatus.OK_GENECOUNT_ONLY
    assert decisions["RVF_virus"].taxon_group is TaxonGroup.VIRUS
    assert decisions["RVF_virus"].amr_applicable is False

    # Protozoa, Fungus, Plant, Metazoa → not applicable
    for sp in ("Toxoplasma_gondii", "Cryptosporidium_parvum",
               "Schizophyllum commune", "Zea_mays",
               "Loxodonta africana", "Aedes_aegypti"):
        assert decisions[sp].amr_applicable is False
        assert decisions[sp].status is AnalysisStatus.OK_GENECOUNT_ONLY

    # No genome → SKIPPED_NO_GENOME
    assert decisions["Cercocebus galeritus"].status is AnalysisStatus.SKIPPED_NO_GENOME
    assert decisions["Kigelia africana"].status   is AnalysisStatus.SKIPPED_NO_GENOME
    assert decisions["Clavipitaceous Fungi"].status is AnalysisStatus.SKIPPED_NO_GENOME

    # Unknown accession → SKIPPED_ASSEMBLY_NOT_FOUND, with warning
    imag = decisions["Imaginary species"]
    assert imag.status is AnalysisStatus.SKIPPED_ASSEMBLY_NOT_FOUND
    assert imag.hint_warnings, "expected a 'not resolvable' warning"

    # Audit-warning case: mosquito mislabelled as Crop
    mislabel = decisions["Aedes_aegypti (mislabel)"]
    assert mislabel.status is AnalysisStatus.OK_GENECOUNT_ONLY  # routing still correct
    assert mislabel.taxon_group is TaxonGroup.METAZOA           # from lineage
    assert mislabel.category_agrees is False                    # audit caught it
    assert any("expected plant" in w for w in mislabel.hint_warnings)

    print("\nAll branch assertions passed.")


if __name__ == "__main__":
    run(out_dir=Path(__file__).resolve().parent.parent / "results")
