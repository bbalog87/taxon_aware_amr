"""
taxon_aware_amr.assembly
========================

Assembly metadata records and retrieval.

Why this module exists
----------------------
The single source of truth for routing decisions is the **NCBI Assembly record**,
addressed by its GCA / GCF accession. This module wraps that record in an
immutable dataclass (``AssemblyRecord``) and exposes a fetcher abstraction so
the rest of the pipeline can stay pure and testable.

We deliberately do NOT route decisions from the input table's
``taxon_id`` column because that column, in the surveillance input used here,
contains values that NCBI will not accept:

    * Breed / strain suffixes: ``9915_1``, ``9925_2``
    * Database-source tags:    ``1392_ENA``
    * Literal absence:         ``none_1``
    * Repeated rows for the same assembly in different ecological roles
      (e.g. ``1392_ENA`` appears for *B. anthracis* both as pathogen and as
      a soil reservoir).

The GCA accession, by contrast, is unique, versioned, and uniquely resolvable
through ``datasets summary genome accession <GCA_*>``.

How it is used
--------------
::

    record = fetch_assembly_record("GCA_964341285.1")    # Step 2: NCBI call
    # or
    record = mock_assembly_record("GCA_964341285.1")     # Step 1 demo

    decision = taxonomy.decide(
        assembly=record,                  # may also be None
        input_genome_id="GCA_964341285.1",
        input_taxon_id="1392_ENA",
        input_species_name="Bacillus_anthracis",
        input_category="Non_eukaryotic_microbe",
    )
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# AssemblyRecord
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AssemblyRecord:
    """Authoritative NCBI Assembly metadata, addressed by GCA accession.

    Attributes
    ----------
    accession
        ``GCA_xxxxxxxxx.N`` (or RefSeq ``GCF_xxxxxxxxx.N``). The version
        suffix is preserved.
    tax_id
        NCBI canonical tax_id reported by the assembly summary. This is
        always integer-valued.
    organism_name
        NCBI canonical scientific name, e.g. ``"Bacillus anthracis"``.
        May include strain in parentheses for some assemblies.
    lineage
        Tuple of scientific names from root to species. Immutable.
    assembly_level
        ``contig`` | ``scaffold`` | ``chromosome`` | ``complete_genome``.
    source
        How this record was obtained: ``"ncbi_datasets"`` (live) or
        ``"mock"`` (Step-1 demo fixture).
    """
    accession: str
    tax_id: int
    organism_name: str
    lineage: tuple[str, ...]
    assembly_level: str
    source: str = "ncbi_datasets"


# ---------------------------------------------------------------------------
# Fetcher abstraction
# ---------------------------------------------------------------------------

class AssemblyNotFound(Exception):
    """Raised when an accession cannot be resolved (suppressed, withdrawn,
    or typo'd)."""


def fetch_assembly_record(accession: str) -> AssemblyRecord:
    """Live NCBI fetch — implemented in Step 2.

    Will shell out to::

        datasets summary genome accession <accession> --as-json-lines

    parse the JSON, and assemble an ``AssemblyRecord``. For Step 1 the
    function is a stub so the rest of the pipeline can be exercised
    against the mock fixtures below.
    """
    raise NotImplementedError(
        "Live NCBI retrieval is implemented in Step 2 (taxon_aware_amr.ncbi). "
        "Use mock_assembly_record() for Step 1."
    )


# ---------------------------------------------------------------------------
# Mock records for Step 1 testing
# ---------------------------------------------------------------------------
#
# These mirror what `datasets summary genome accession <GCA>` returns for
# representative rows of the input table. They let us validate the routing
# logic end-to-end without hitting the network.
#
# Note: Some `assembly_level` values are normalised (e.g. NCBI returns
# "Complete Genome" → we use "complete_genome" to match the input file's
# convention).

_MOCK_ASSEMBLIES: dict[str, AssemblyRecord] = {
    # --- Bacteria (AMR + virulence applicable) ---
    "GCA_964341285.1": AssemblyRecord(
        accession="GCA_964341285.1",
        tax_id=1392,
        organism_name="Bacillus anthracis",
        lineage=("root", "cellular root", "Bacteria", "Bacillota", "Bacilli",
                 "Bacillales", "Bacillaceae", "Bacillus", "Bacillus cereus group"),
        assembly_level="complete_genome",
        source="mock",
    ),
    "GCA_900654255.1": AssemblyRecord(
        accession="GCA_900654255.1",
        tax_id=1773,
        organism_name="Mycobacterium tuberculosis",
        lineage=("root", "cellular root", "Bacteria", "Actinomycetota",
                 "Actinomycetes", "Mycobacteriales", "Mycobacteriaceae",
                 "Mycobacterium tuberculosis complex"),
        assembly_level="scaffold",
        source="mock",
    ),
    "GCA_001049615.1": AssemblyRecord(
        accession="GCA_001049615.1",
        tax_id=1280,
        organism_name="Staphylococcus aureus",
        lineage=("root", "cellular root", "Bacteria", "Bacillota", "Bacilli",
                 "Bacillales", "Staphylococcaceae", "Staphylococcus"),
        assembly_level="contig",
        source="mock",
    ),

    # --- Viruses (same input category as bacteria — proves GCA-routing works) ---
    "GCA_000847345.1": AssemblyRecord(
        accession="GCA_000847345.1",
        tax_id=11588,
        organism_name="Rift Valley fever virus",
        lineage=("root", "Viruses", "Riboviria", "Orthornavirae",
                 "Negarnaviricota", "Ellioviricetes", "Bunyavirales",
                 "Phenuiviridae", "Phlebovirus"),
        assembly_level="complete_genome",
        source="mock",
    ),
    "GCA_900327715.1": AssemblyRecord(
        accession="GCA_900327715.1",
        tax_id=11292,
        organism_name="Rabies lyssavirus",
        lineage=("root", "Viruses", "Riboviria", "Orthornavirae",
                 "Negarnaviricota", "Monjiviricetes", "Mononegavirales",
                 "Rhabdoviridae", "Lyssavirus"),
        assembly_level="complete_genome",
        source="mock",
    ),

    # --- Protozoa ---
    "GCA_019455585.1": AssemblyRecord(
        accession="GCA_019455585.1",
        tax_id=5811,
        organism_name="Toxoplasma gondii",
        lineage=("root", "cellular root", "Eukaryota", "Sar", "Alveolata",
                 "Apicomplexa", "Conoidasida", "Coccidia", "Eucoccidiorida",
                 "Eimeriorina", "Sarcocystidae", "Toxoplasma"),
        assembly_level="chromosome",
        source="mock",
    ),
    "GCA_035232765.1": AssemblyRecord(
        accession="GCA_035232765.1",
        tax_id=5807,
        organism_name="Cryptosporidium parvum",
        lineage=("root", "cellular root", "Eukaryota", "Sar", "Alveolata",
                 "Apicomplexa", "Conoidasida", "Coccidia", "Eucoccidiorida",
                 "Eimeriorina", "Cryptosporidiidae", "Cryptosporidium"),
        assembly_level="complete_genome",
        source="mock",
    ),

    # --- Fungus ---
    "GCA_014900015.1": AssemblyRecord(
        accession="GCA_014900015.1",
        tax_id=5334,
        organism_name="Schizophyllum commune",
        lineage=("root", "cellular root", "Eukaryota", "Opisthokonta",
                 "Fungi", "Dikarya", "Basidiomycota", "Agaricomycotina",
                 "Agaricomycetes", "Agaricales", "Schizophyllaceae"),
        assembly_level="contig",
        source="mock",
    ),

    # --- Plant ---
    "GCA_964199775.1": AssemblyRecord(
        accession="GCA_964199775.1",
        tax_id=4577,
        organism_name="Zea mays",
        lineage=("root", "cellular root", "Eukaryota", "Viridiplantae",
                 "Streptophyta", "Embryophyta", "Tracheophyta", "Magnoliopsida",
                 "Poales", "Poaceae", "Zea"),
        assembly_level="complete_genome",
        source="mock",
    ),

    # --- Mammal (Metazoa) ---
    "GCA_030014295.1": AssemblyRecord(
        accession="GCA_030014295.1",
        tax_id=9785,
        organism_name="Loxodonta africana",
        lineage=("root", "cellular root", "Eukaryota", "Opisthokonta",
                 "Metazoa", "Eumetazoa", "Bilateria", "Deuterostomia",
                 "Chordata", "Mammalia", "Proboscidea", "Elephantidae",
                 "Loxodonta"),
        assembly_level="chromosome",
        source="mock",
    ),

    # --- Insect (Metazoa) ---
    "GCA_025407655.1": AssemblyRecord(
        accession="GCA_025407655.1",
        tax_id=7159,
        organism_name="Aedes aegypti",
        lineage=("root", "cellular root", "Eukaryota", "Opisthokonta",
                 "Metazoa", "Eumetazoa", "Bilateria", "Protostomia",
                 "Ecdysozoa", "Arthropoda", "Hexapoda", "Insecta",
                 "Diptera", "Culicidae", "Aedes"),
        assembly_level="chromosome",
        source="mock",
    ),
}


def mock_assembly_record(accession: str) -> Optional[AssemblyRecord]:
    """Return a mock ``AssemblyRecord`` for Step 1 testing.

    Returns ``None`` if the accession is not in the fixture set, simulating
    the ``AssemblyNotFound`` situation that Step 2 will raise.
    """
    return _MOCK_ASSEMBLIES.get(accession)


# ---------------------------------------------------------------------------
# Input cleaning helpers
# ---------------------------------------------------------------------------

_GENOME_SENTINELS: frozenset[str] = frozenset({"", "none", "None", "NA", "NaN", "nan"})


def normalise_genome_id(raw: Optional[str]) -> Optional[str]:
    """Return a clean accession or None if the input means 'no genome'.

    Examples
    --------
    >>> normalise_genome_id("GCA_964341285.1")
    'GCA_964341285.1'
    >>> normalise_genome_id("none") is None
    True
    >>> normalise_genome_id("  ") is None
    True
    """
    if raw is None:
        return None
    s = raw.strip()
    if s in _GENOME_SENTINELS:
        return None
    return s
