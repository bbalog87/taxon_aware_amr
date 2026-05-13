"""
taxon_aware_amr.final_output
============================

Step 6: assemble the final deliverable dataframes.

Two outputs
-----------
1. ``final_output.tsv`` — **minimal schema** matching the original task spec:

       genome_id
       total_gene_count
       virulence_genes
       virulence_gene_count
       amr_genes
       amr_gene_count

   This is the clean dataframe that downstream consumers asked for. Non-bacterial
   rows show ``not_applicable`` in the AMR and virulence columns; rows without
   a usable genome show ``NA``.

2. ``final_output_full.tsv`` — **extended schema** for transparency. Adds the
   audit columns produced by Steps 1–5 (taxon_group, ncbi_accession,
   amr_drug_classes, virulence_sources, analysis_status, reason, hint_warnings).

The minimal output uses the **GCA accession actually analysed** as
``genome_id`` (i.e., after auto-replace). The original input genome_id is
preserved in the extended output as ``input_genome_id`` for audit.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .taxonomy import AnalysisStatus, TaxonDecision

logger = logging.getLogger(__name__)


# ===========================================================================
# Per-row final record
# ===========================================================================

@dataclass
class FinalRow:
    """The merged result of Steps 1–5 for one input row, ready to serialise."""

    # Minimal schema (matches original task spec)
    genome_id:               str           # NCBI accession actually analysed
    total_gene_count:        Optional[int]
    virulence_genes:         str           # "; "-joined, or "not_applicable", or "NA"
    virulence_gene_count:    Optional[int]
    amr_genes:               str
    amr_gene_count:          Optional[int]

    # Extended schema
    input_genome_id:         Optional[str]
    input_taxon_id:          str
    species_name:            str
    taxon_group:             str
    ncbi_organism_name:      Optional[str]
    ncbi_tax_id:             Optional[int]
    ncbi_assembly_level:     Optional[str]
    amr_drug_classes:        str
    amr_organism_flag:       str
    amr_method:              str
    virulence_sources:       str           # "; "-joined source labels
    virulence_methods:       str           # "; "-joined "tool@db_version"
    analysis_status:         str
    reason:                  str
    hint_warnings:           str           # "; "-joined
    gene_count_source:       str

    # ---- Serialisation helpers ----
    def minimal_row(self) -> dict:
        return {
            "genome_id":            self.genome_id,
            "total_gene_count":     "" if self.total_gene_count is None
                                    else self.total_gene_count,
            "virulence_genes":      self.virulence_genes,
            "virulence_gene_count": "" if self.virulence_gene_count is None
                                    else self.virulence_gene_count,
            "amr_genes":            self.amr_genes,
            "amr_gene_count":       "" if self.amr_gene_count is None
                                    else self.amr_gene_count,
        }

    def full_row(self) -> dict:
        d = self.minimal_row()
        d.update({
            "input_genome_id":     self.input_genome_id or "",
            "input_taxon_id":      self.input_taxon_id,
            "species_name":        self.species_name,
            "taxon_group":         self.taxon_group,
            "ncbi_organism_name":  self.ncbi_organism_name or "",
            "ncbi_tax_id":         "" if self.ncbi_tax_id is None else self.ncbi_tax_id,
            "ncbi_assembly_level": self.ncbi_assembly_level or "",
            "amr_drug_classes":    self.amr_drug_classes,
            "amr_organism_flag":   self.amr_organism_flag,
            "amr_method":          self.amr_method,
            "virulence_sources":   self.virulence_sources,
            "virulence_methods":   self.virulence_methods,
            "analysis_status":     self.analysis_status,
            "reason":              self.reason,
            "hint_warnings":       self.hint_warnings,
            "gene_count_source":   self.gene_count_source,
        })
        return d


# ===========================================================================
# Build a FinalRow from a TaxonDecision + virulence info
# ===========================================================================

def build_final_row(
    decision: TaxonDecision,
    *,
    virulence_genes:   str = "",
    virulence_count:   Optional[int] = None,
    virulence_sources: str = "",
    virulence_methods: str = "",
) -> FinalRow:
    """Merge Step 1–5 outputs into one FinalRow."""
    # genome_id: use the NCBI accession actually analysed (post auto-replace).
    # Fall back to the input value for skipped-no-genome rows.
    if decision.ncbi_accession:
        genome_id = decision.ncbi_accession
    elif decision.input_genome_id:
        genome_id = decision.input_genome_id
    else:
        genome_id = "none"

    return FinalRow(
        genome_id=genome_id,
        total_gene_count=decision.total_gene_count,
        virulence_genes=virulence_genes or _na_value_for(decision, "virulence"),
        virulence_gene_count=virulence_count,
        amr_genes=decision.amr_genes or _na_value_for(decision, "amr"),
        amr_gene_count=decision.amr_gene_count,
        input_genome_id=decision.input_genome_id,
        input_taxon_id=decision.input_taxon_id,
        species_name=decision.input_species_name,
        taxon_group=decision.taxon_group.value,
        ncbi_organism_name=decision.ncbi_organism_name,
        ncbi_tax_id=decision.ncbi_tax_id,
        ncbi_assembly_level=decision.ncbi_assembly_level,
        amr_drug_classes=decision.amr_drug_classes,
        amr_organism_flag=decision.amr_organism_flag,
        amr_method=decision.amr_method,
        virulence_sources=virulence_sources,
        virulence_methods=virulence_methods,
        analysis_status=decision.status.value,
        reason=decision.reason,
        hint_warnings="; ".join(decision.hint_warnings),
        gene_count_source=decision.gene_count_source,
    )


def _na_value_for(decision: TaxonDecision, kind: str) -> str:
    """Choose the right empty-cell label for AMR/virulence columns."""
    if decision.status is AnalysisStatus.OK_GENECOUNT_ONLY:
        return "not_applicable"
    if decision.status is AnalysisStatus.OK_BACTERIAL:
        # Bacterial but no hits or analysis failed.
        return ""
    return "NA"


# ===========================================================================
# TSV writers
# ===========================================================================

MINIMAL_COLUMNS = (
    "genome_id", "total_gene_count",
    "virulence_genes", "virulence_gene_count",
    "amr_genes", "amr_gene_count",
)

FULL_COLUMNS = MINIMAL_COLUMNS + (
    "input_genome_id", "input_taxon_id", "species_name",
    "taxon_group", "ncbi_organism_name", "ncbi_tax_id", "ncbi_assembly_level",
    "amr_drug_classes", "amr_organism_flag", "amr_method",
    "virulence_sources", "virulence_methods",
    "analysis_status", "reason", "hint_warnings", "gene_count_source",
)


def write_minimal_tsv(rows: list[FinalRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=MINIMAL_COLUMNS,
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r.minimal_row())


def write_full_tsv(rows: list[FinalRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FULL_COLUMNS,
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r.full_row())
