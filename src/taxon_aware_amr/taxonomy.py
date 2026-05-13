"""
taxon_aware_amr.taxonomy
========================

Routing decisions for the AMR / virulence pipeline.

Design (read this before editing)
---------------------------------
This module is **pure**: it makes a decision given inputs, with no I/O.

A decision is made from exactly one source of truth — the ``AssemblyRecord``
fetched from NCBI by GCA accession. The input table's ``taxon_id`` and
``category`` columns are recorded for **audit only**; they never affect
routing. This avoids the well-known input quirks documented in
``assembly.py``.

Five possible outcomes per input row:

    ┌────────────────────────────────┬───────────────────────────────────────┐
    │ outcome                        │ trigger                               │
    ├────────────────────────────────┼───────────────────────────────────────┤
    │ SKIPPED_NO_GENOME              │ genome_id is 'none' / empty           │
    │ SKIPPED_ASSEMBLY_NOT_FOUND     │ GCA given but NCBI cannot resolve it  │
    │ SKIPPED_TAXON_UNRESOLVED       │ lineage matched no known group        │
    │ OK_BACTERIAL                   │ lineage contains 'Bacteria'           │
    │ OK_GENECOUNT_ONLY              │ lineage is virus / fungus / protozoa  │
    │                                │ / plant / metazoa / archaea           │
    └────────────────────────────────┴───────────────────────────────────────┘

Why bacterial-only for AMR/virulence
------------------------------------
AMRFinderPlus, CARD, ResFinder and VFDB are bacterial pathogen databases.
Running them against eukaryotic or viral assemblies returns spurious
homology to conserved transporters, P450s, ABC families, or assembly
contamination — not biology. See README for full rationale.

Audit fields
------------
Each decision carries:
    * ``input_*``           — exactly what was in the input row, for traceability.
    * ``ncbi_*``            — what NCBI reported, from the assembly record.
    * ``taxon_id_agrees``   — does the *core* of input_taxon_id match NCBI's?
                              (Core = digits before any '_' suffix.)
    * ``category_agrees``   — does the input category map to the resolved group?
    * ``hint_warnings``     — short notes for the run report.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, asdict, field
from enum import Enum
from typing import Optional

from .assembly import AssemblyRecord

logger = logging.getLogger(__name__)


# ===========================================================================
# Enums
# ===========================================================================

class TaxonGroup(str, Enum):
    """High-level taxonomic group used for routing decisions."""
    BACTERIA = "bacteria"
    ARCHAEA  = "archaea"
    VIRUS    = "virus"
    FUNGUS   = "fungus"
    PROTOZOA = "protozoa"
    PLANT    = "plant"
    METAZOA  = "metazoa"     # all animals (vertebrates + invertebrates)
    UNKNOWN  = "unknown"


class AnalysisStatus(str, Enum):
    """Per-row outcome."""
    OK_BACTERIAL                = "ok_bacterial"
    OK_GENECOUNT_ONLY           = "ok_genecount_only"
    SKIPPED_NO_GENOME           = "skipped_no_genome"
    SKIPPED_ASSEMBLY_NOT_FOUND  = "skipped_assembly_not_found"
    SKIPPED_TAXON_UNRESOLVED    = "skipped_taxon_unresolved"


# ===========================================================================
# Decision record
# ===========================================================================

@dataclass
class TaxonDecision:
    """A complete routing decision plus its audit trail."""

    # ----- Inputs as supplied (verbatim, for traceability) -----
    input_genome_id:    Optional[str]
    input_taxon_id:     str
    input_species_name: str
    input_category:     str

    # ----- NCBI authoritative fields (None when no assembly) -----
    ncbi_accession:      Optional[str]
    ncbi_tax_id:         Optional[int]
    ncbi_organism_name:  Optional[str]
    ncbi_assembly_level: Optional[str]

    # ----- Routing -----
    taxon_group:          TaxonGroup
    amr_applicable:       bool
    virulence_applicable: bool
    status:               AnalysisStatus

    # ----- Audit -----
    taxon_id_agrees:  bool       # input_taxon_id core matches ncbi_tax_id
    category_agrees:  bool       # input_category maps to resolved taxon_group
    hint_warnings:    list[str] = field(default_factory=list)

    # ----- Provenance / reason -----
    method_used: str = ""        # "assembly_lineage" | "no_assembly"
    reason:      str = ""        # human-readable explanation surfaced in report

    # ----- Step 3: gene counting (filled by orchestrator) -----
    total_gene_count:   Optional[int] = None
    gene_count_source:  str           = ""    # ncbi_summary | gff3_parsed | no_source | gff_parse_error...

    # ----- Step 4: AMR detection (filled by orchestrator) -----
    amr_status:         str           = ""    # ok | failed | not_applicable | not_attempted
    amr_genes:          str           = ""    # "; "-joined gene symbols, or "not_applicable"
    amr_gene_count:     Optional[int] = None
    amr_drug_classes:   str           = ""    # "; "-joined drug classes
    amr_method:         str           = ""    # "amrfinderplus_<dbver>" or "not_applicable"
    amr_organism_flag:  str           = ""    # value passed to --organism, or "" / "auto:none"

    def to_row(self) -> dict:
        """Flatten to a dict suitable for a pandas DataFrame row or TSV."""
        d = asdict(self)
        d["taxon_group"] = self.taxon_group.value
        d["status"] = self.status.value
        d["hint_warnings"] = "; ".join(self.hint_warnings)
        return d


# ===========================================================================
# Category-hint mapping (used for AUDIT only, never for routing)
# ===========================================================================
#
# Maps the input table's `category` column to the TaxonGroup the user
# expected. `None` means the category is intrinsically ambiguous (e.g.
# 'Non_eukaryotic_microbe' covers both bacteria and viruses).

CATEGORY_TO_EXPECTED_GROUP: dict[str, Optional[TaxonGroup]] = {
    "Wild-plant":             TaxonGroup.PLANT,
    "Crop":                   TaxonGroup.PLANT,
    "Wildlife":               TaxonGroup.METAZOA,
    "Wild-bird":              TaxonGroup.METAZOA,
    "Domestic-bird":          TaxonGroup.METAZOA,
    "Livestock":              TaxonGroup.METAZOA,
    "Mammal":                 TaxonGroup.METAZOA,
    "Reptile":                TaxonGroup.METAZOA,
    "Amphibians":             TaxonGroup.METAZOA,
    "Fish":                   TaxonGroup.METAZOA,
    "Insect":                 TaxonGroup.METAZOA,
    "Arachnid":               TaxonGroup.METAZOA,
    "Fungus":                 TaxonGroup.FUNGUS,
    "Protozoa":               TaxonGroup.PROTOZOA,
    "Non_eukaryotic_microbe": None,   # ambiguous — bacteria OR virus
}


# ===========================================================================
# Lineage → TaxonGroup
# ===========================================================================
#
# Protozoa is paraphyletic; we recognise it by characteristic eukaryotic
# clades that are neither Fungi, Viridiplantae, nor Metazoa.

_PROTOZOAN_CLADES: frozenset[str] = frozenset({
    "Apicomplexa",     # Plasmodium, Toxoplasma, Cryptosporidium
    "Kinetoplastea",   # Trypanosoma, Leishmania
    "Euglenozoa",      # superphylum incl. Kinetoplastea
    "Fornicata",       # Giardia
    "Amoebozoa",       # Entamoeba
    "Heterolobosea",   # Naegleria
    "Parabasalia",     # Trichomonas
})


# Modern viral realms (ICTV 2020+ classification, recognised by NCBI Taxonomy).
# Viruses are not assigned a superkingdom/domain in the new ranked lineage,
# so realms (or the -viricota / -viricetes suffix family) are the only
# top-level signal that lands in the classification we receive.
_VIRAL_REALMS: frozenset[str] = frozenset({
    "Riboviria", "Duplodnaviria", "Monodnaviria",
    "Varidnaviria", "Adnaviria", "Ribozyviria",
})


def lineage_to_group(lineage: tuple[str, ...] | list[str]) -> TaxonGroup:
    """Resolve an NCBI taxonomic lineage to a high-level group.

    The lineage is searched as a set — order doesn't matter, only
    membership. Tests are applied in priority order (bacteria first
    because some sequencing kits include bacterial spike-ins that
    could otherwise be mis-classified).

    Returns ``TaxonGroup.UNKNOWN`` if no rule fires; this indicates a
    missing rule (e.g. a previously unseen clade) and triggers
    ``SKIPPED_TAXON_UNRESOLVED`` downstream.
    """
    names = set(lineage)
    if "Bacteria"      in names: return TaxonGroup.BACTERIA
    if "Archaea"       in names: return TaxonGroup.ARCHAEA
    if "Viruses"       in names: return TaxonGroup.VIRUS
    if "Fungi"         in names: return TaxonGroup.FUNGUS
    if "Viridiplantae" in names: return TaxonGroup.PLANT
    if "Metazoa"       in names: return TaxonGroup.METAZOA
    if names & _PROTOZOAN_CLADES: return TaxonGroup.PROTOZOA
    # Viral fallback: NCBI's "Viruses" superkingdom does not appear in the
    # modern ranked classification; we identify viruses by realm membership
    # (Riboviria, Duplodnaviria, ...) or by the canonical phylum/class
    # suffixes -viricota / -viricetes.
    if names & _VIRAL_REALMS:
        return TaxonGroup.VIRUS
    for n in names:
        if n.endswith("viricota") or n.endswith("viricetes"):
            return TaxonGroup.VIRUS
    return TaxonGroup.UNKNOWN


# ===========================================================================
# Input taxon_id cleaning (audit-only)
# ===========================================================================

_TAXON_ID_CORE_RE = re.compile(r"^\s*(\d+)")


def _taxon_id_core(raw: str) -> Optional[int]:
    """Extract the leading integer from a possibly-decorated input taxon_id.

    Examples
    --------
    >>> _taxon_id_core("1392_ENA")
    1392
    >>> _taxon_id_core("9915_1")
    9915
    >>> _taxon_id_core("none_1") is None
    True
    >>> _taxon_id_core("  658858  ")
    658858
    """
    if raw is None:
        return None
    m = _TAXON_ID_CORE_RE.match(raw)
    return int(m.group(1)) if m else None


# ===========================================================================
# decide() — the routing function
# ===========================================================================

def decide(
    *,                                          # keyword-only for readability
    assembly:           Optional[AssemblyRecord],
    input_genome_id:    Optional[str],
    input_taxon_id:     str,
    input_species_name: str,
    input_category:     str,
) -> TaxonDecision:
    """Make a routing decision for a single input row.

    The function is **pure** — it does no I/O. Callers are responsible for
    fetching the ``AssemblyRecord`` (or passing ``None`` if no genome_id
    was supplied / the accession could not be resolved).

    Distinguishing the two None-assembly cases:
        * ``input_genome_id is None``    → no genome was requested  → SKIPPED_NO_GENOME
        * ``input_genome_id is not None``
          and ``assembly is None``       → NCBI couldn't resolve it → SKIPPED_ASSEMBLY_NOT_FOUND

    Returns
    -------
    TaxonDecision
    """
    # ----------------------------------------------------------------- #
    # Case A: No genome_id at all                                        #
    # ----------------------------------------------------------------- #
    if input_genome_id is None:
        return TaxonDecision(
            input_genome_id=None,
            input_taxon_id=input_taxon_id,
            input_species_name=input_species_name,
            input_category=input_category,
            ncbi_accession=None,
            ncbi_tax_id=None,
            ncbi_organism_name=None,
            ncbi_assembly_level=None,
            taxon_group=TaxonGroup.UNKNOWN,
            amr_applicable=False,
            virulence_applicable=False,
            status=AnalysisStatus.SKIPPED_NO_GENOME,
            taxon_id_agrees=False,        # nothing to compare against
            category_agrees=False,
            hint_warnings=[],
            method_used="no_assembly",
            reason="No genome accession provided in input.",
        )

    # ----------------------------------------------------------------- #
    # Case B: genome_id given but NCBI could not resolve it              #
    # ----------------------------------------------------------------- #
    if assembly is None:
        return TaxonDecision(
            input_genome_id=input_genome_id,
            input_taxon_id=input_taxon_id,
            input_species_name=input_species_name,
            input_category=input_category,
            ncbi_accession=None,
            ncbi_tax_id=None,
            ncbi_organism_name=None,
            ncbi_assembly_level=None,
            taxon_group=TaxonGroup.UNKNOWN,
            amr_applicable=False,
            virulence_applicable=False,
            status=AnalysisStatus.SKIPPED_ASSEMBLY_NOT_FOUND,
            taxon_id_agrees=False,
            category_agrees=False,
            hint_warnings=[
                f"Accession {input_genome_id} not resolvable at NCBI "
                f"(suppressed / withdrawn / typo?)."
            ],
            method_used="no_assembly",
            reason=f"NCBI assembly summary unavailable for {input_genome_id}.",
        )

    # ----------------------------------------------------------------- #
    # Case C: We have an assembly. Resolve taxon group from its lineage. #
    # ----------------------------------------------------------------- #
    group = lineage_to_group(assembly.lineage)

    # ---- Audit cross-checks (don't affect routing — only inform the report) ----
    warnings: list[str] = []

    input_core = _taxon_id_core(input_taxon_id)
    taxon_id_agrees = (input_core is not None) and (input_core == assembly.tax_id)
    if input_core is None:
        warnings.append(
            f"Input taxon_id {input_taxon_id!r} has no leading integer; "
            "unparseable for audit."
        )
    elif not taxon_id_agrees:
        warnings.append(
            f"Input taxon_id core ({input_core}) ≠ NCBI tax_id "
            f"({assembly.tax_id}) for {assembly.organism_name}."
        )

    expected = CATEGORY_TO_EXPECTED_GROUP.get(input_category)
    # If the category is intrinsically ambiguous (None), we cannot
    # disagree — the audit is a no-op.
    category_agrees = (expected is None) or (expected == group)
    if not category_agrees:
        warnings.append(
            f"Input category {input_category!r} expected "
            f"{expected.value if expected else 'unknown'}, "
            f"but NCBI lineage resolves to {group.value}."
        )

    # ---- Applicability + status from group ----
    if group == TaxonGroup.UNKNOWN:
        status = AnalysisStatus.SKIPPED_TAXON_UNRESOLVED
        amr_applicable = False
        virulence_applicable = False
        reason = (
            f"Lineage {list(assembly.lineage)!r} matched no known group rule. "
            "Add a rule in lineage_to_group() if this clade should be supported."
        )
    elif group == TaxonGroup.BACTERIA:
        status = AnalysisStatus.OK_BACTERIAL
        amr_applicable = True
        virulence_applicable = True
        reason = (
            f"Bacterial assembly ({assembly.organism_name}) — "
            "AMRFinderPlus and VFDB applicable."
        )
    else:
        status = AnalysisStatus.OK_GENECOUNT_ONLY
        amr_applicable = False
        virulence_applicable = False
        reason = (
            f"{group.value.capitalize()} assembly ({assembly.organism_name}) — "
            "bacterial AMR/virulence databases are not biologically applicable. "
            "Gene count from canonical annotation only."
        )

    return TaxonDecision(
        input_genome_id=input_genome_id,
        input_taxon_id=input_taxon_id,
        input_species_name=input_species_name,
        input_category=input_category,
        ncbi_accession=assembly.accession,
        ncbi_tax_id=assembly.tax_id,
        ncbi_organism_name=assembly.organism_name,
        ncbi_assembly_level=assembly.assembly_level,
        taxon_group=group,
        amr_applicable=amr_applicable,
        virulence_applicable=virulence_applicable,
        status=status,
        taxon_id_agrees=taxon_id_agrees,
        category_agrees=category_agrees,
        hint_warnings=warnings,
        method_used="assembly_lineage",
        reason=reason,
    )
