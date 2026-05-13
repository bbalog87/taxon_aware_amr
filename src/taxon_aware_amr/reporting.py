"""
taxon_aware_amr.reporting
=========================

Decision logging and run reporting.

This module exists because routing decisions are useless if they're invisible.
For every input row, the user must be able to see:

    * **In the command** (stdout, while the pipeline runs)
        - what input was processed
        - what NCBI returned
        - what the routing decision was, and why

    * **In a report file** (persistent, machine + human readable)
        - a TSV with one row per decision (decisions.tsv)
        - a Markdown summary suitable for emailing or attaching to a report
          (run_report.md), grouping rows by status, listing audit warnings,
          and ending with a counts table

Both outputs are produced by a single ``DecisionReporter`` so the streaming
log and the persisted log can never drift.

Typical use
-----------
::

    reporter = DecisionReporter(
        tsv_path=Path("results/decisions.tsv"),
        summary_path=Path("results/run_report.md"),
    )
    with reporter:
        for row in input_rows:
            decision = decide(...)
            reporter.log(decision)
    # Markdown summary is written when the context exits.
"""

from __future__ import annotations

import csv
import logging
import sys
from collections import Counter
from dataclasses import fields
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO

from .taxonomy import AnalysisStatus, TaxonDecision, TaxonGroup

logger = logging.getLogger(__name__)


# ===========================================================================
# Human-readable status descriptions for the summary report
# ===========================================================================

_STATUS_DESCRIPTION: dict[AnalysisStatus, str] = {
    AnalysisStatus.OK_BACTERIAL: (
        "Bacterial assembly — AMRFinderPlus and VFDB will run."
    ),
    AnalysisStatus.OK_GENECOUNT_ONLY: (
        "Eukaryotic or viral assembly — gene count only; "
        "bacterial AMR/VFDB databases are not biologically applicable."
    ),
    AnalysisStatus.SKIPPED_NO_GENOME: (
        "No genome accession in input — nothing to analyse."
    ),
    AnalysisStatus.SKIPPED_ASSEMBLY_NOT_FOUND: (
        "Accession given but NCBI could not resolve it "
        "(suppressed / withdrawn / typo)."
    ),
    AnalysisStatus.SKIPPED_TAXON_UNRESOLVED: (
        "Assembly fetched but its NCBI lineage matched no known taxon-group "
        "rule. Pipeline rule needs extending."
    ),
}


# ===========================================================================
# DecisionReporter
# ===========================================================================

class DecisionReporter:
    """Streams routing decisions to stdout, a TSV log, and a Markdown summary.

    Parameters
    ----------
    tsv_path
        Where to write the per-row TSV log. One header row, then one row per
        input. Columns are taken from ``TaxonDecision.to_row()``.
    summary_path
        Where to write the Markdown summary report (produced on ``close()``
        or when the context manager exits).
    stream
        File-like object for the live log. Default: ``sys.stdout``. Pass a
        StringIO for testing or ``open(os.devnull, "w")`` for quiet mode.
    show_progress
        If ``True``, writes one line per decision to ``stream``.

    Attributes
    ----------
    decisions
        All decisions logged so far. Reset on ``close()``.
    """

    def __init__(
        self,
        tsv_path: Path,
        summary_path: Path,
        stream: TextIO = sys.stdout,
        show_progress: bool = True,
    ) -> None:
        self.tsv_path = Path(tsv_path)
        self.summary_path = Path(summary_path)
        self.stream = stream
        self.show_progress = show_progress
        self.decisions: list[TaxonDecision] = []
        self._tsv_handle: TextIO | None = None
        self._tsv_writer: csv.DictWriter | None = None
        self._started_at = datetime.now(timezone.utc)

        self.tsv_path.parent.mkdir(parents=True, exist_ok=True)
        self.summary_path.parent.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Context manager interface                                           #
    # ------------------------------------------------------------------ #
    def __enter__(self) -> "DecisionReporter":
        self._open_tsv()
        if self.show_progress:
            self._emit_header()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    # ------------------------------------------------------------------ #
    # Public API                                                          #
    # ------------------------------------------------------------------ #
    def log(self, decision: TaxonDecision) -> None:
        """Record a decision: stream it, write it to the TSV, and remember it."""
        self.decisions.append(decision)
        if self._tsv_writer is None:
            self._open_tsv()
        assert self._tsv_writer is not None
        self._tsv_writer.writerow(decision.to_row())
        if self.show_progress:
            self._emit_row(decision)

    def close(self) -> None:
        """Close the TSV and write the Markdown summary."""
        if self._tsv_handle is not None:
            self._tsv_handle.close()
            self._tsv_handle = None
            self._tsv_writer = None
        self._write_summary()
        if self.show_progress:
            self.stream.write(
                f"\nWrote {len(self.decisions)} decisions to {self.tsv_path}\n"
                f"Wrote run summary to {self.summary_path}\n"
            )

    # ------------------------------------------------------------------ #
    # Internals — TSV                                                     #
    # ------------------------------------------------------------------ #
    def _open_tsv(self) -> None:
        # Field order = order of TaxonDecision dataclass fields, so the
        # TSV column order tracks the source of truth automatically.
        fieldnames = [f.name for f in fields(TaxonDecision)]
        self._tsv_handle = self.tsv_path.open("w", newline="", encoding="utf-8")
        self._tsv_writer = csv.DictWriter(
            self._tsv_handle, fieldnames=fieldnames, delimiter="\t",
            extrasaction="ignore",
        )
        self._tsv_writer.writeheader()

    # ------------------------------------------------------------------ #
    # Internals — stdout streaming                                        #
    # ------------------------------------------------------------------ #
    _COLS = [
        ("status",       28),
        ("group",        10),
        ("input_id",     20),
        ("ncbi_accession", 18),
        ("ncbi_organism", 32),
        ("amr",            4),
        ("vir",            4),
        ("audit",          8),
    ]

    def _emit_header(self) -> None:
        header = " | ".join(name.ljust(w) for name, w in self._COLS)
        rule = "-+-".join("-" * w for _, w in self._COLS)
        self.stream.write(header + "\n" + rule + "\n")

    def _emit_row(self, d: TaxonDecision) -> None:
        audit = "ok"
        if d.hint_warnings:
            audit = "warn"
        cells = [
            d.status.value,
            d.taxon_group.value,
            (d.input_genome_id or "—"),
            (d.ncbi_accession or "—"),
            (d.ncbi_organism_name or d.input_species_name)[:32],
            "yes" if d.amr_applicable else "no",
            "yes" if d.virulence_applicable else "no",
            audit,
        ]
        line = " | ".join(str(c).ljust(w) for (_, w), c in zip(self._COLS, cells))
        self.stream.write(line + "\n")

    # ------------------------------------------------------------------ #
    # Internals — Markdown summary                                        #
    # ------------------------------------------------------------------ #
    def _write_summary(self) -> None:
        finished_at = datetime.now(timezone.utc)
        duration = finished_at - self._started_at
        total = len(self.decisions)
        by_status = Counter(d.status for d in self.decisions)
        by_group  = Counter(d.taxon_group for d in self.decisions)
        warned    = [d for d in self.decisions if d.hint_warnings]
        bacterial = [d for d in self.decisions if d.status is AnalysisStatus.OK_BACTERIAL]
        gene_only = [d for d in self.decisions if d.status is AnalysisStatus.OK_GENECOUNT_ONLY]
        skipped   = [d for d in self.decisions if d.status.value.startswith("skipped_")]

        with self.summary_path.open("w", encoding="utf-8") as fh:
            w = fh.write

            w("# Taxon-aware AMR / virulence pipeline — run summary\n\n")
            w(f"- **Started**:  `{self._started_at.isoformat(timespec='seconds')}`\n")
            w(f"- **Finished**: `{finished_at.isoformat(timespec='seconds')}`\n")
            w(f"- **Duration**: `{duration}`\n")
            w(f"- **Rows processed**: {total}\n")
            w(f"- **Per-row TSV log**: `{self.tsv_path}`\n\n")

            # ---- Status counts ---------------------------------------- #
            w("## Outcome counts\n\n")
            w("| Status | Count | Meaning |\n")
            w("|---|---:|---|\n")
            for status in AnalysisStatus:
                cnt = by_status.get(status, 0)
                if cnt == 0:
                    continue
                w(f"| `{status.value}` | {cnt} | {_STATUS_DESCRIPTION[status]} |\n")
            w("\n")

            # ---- Group counts ----------------------------------------- #
            w("## Resolved taxon groups\n\n")
            w("| Taxon group | Count |\n|---|---:|\n")
            for grp in TaxonGroup:
                cnt = by_group.get(grp, 0)
                if cnt:
                    w(f"| `{grp.value}` | {cnt} |\n")
            w("\n")

            # ---- Rows that will run AMR -------------------------------- #
            w(f"## Bacterial rows (AMR + virulence will run): {len(bacterial)}\n\n")
            if bacterial:
                w("| Input genome_id | NCBI accession | Organism | Assembly level |\n")
                w("|---|---|---|---|\n")
                for d in bacterial:
                    w(f"| `{d.input_genome_id}` | `{d.ncbi_accession}` | "
                      f"*{d.ncbi_organism_name}* | {d.ncbi_assembly_level} |\n")
                w("\n")

            # ---- Gene-count-only rows ---------------------------------- #
            w(f"## Gene-count-only rows: {len(gene_only)}\n\n")
            if gene_only:
                w("Bacterial AMR/VFDB **not run** for these — they are not bacterial. "
                  "Gene counts will be derived from the canonical annotation.\n\n")
                w("| Genome | Organism | Group |\n|---|---|---|\n")
                for d in gene_only:
                    w(f"| `{d.ncbi_accession}` | *{d.ncbi_organism_name}* | "
                      f"{d.taxon_group.value} |\n")
                w("\n")

            # ---- Skipped rows ----------------------------------------- #
            w(f"## Skipped rows: {len(skipped)}\n\n")
            if skipped:
                w("| Input genome_id | Input species | Reason |\n|---|---|---|\n")
                for d in skipped:
                    w(f"| `{d.input_genome_id or 'none'}` | "
                      f"{d.input_species_name} | {d.reason} |\n")
                w("\n")

            # ---- Audit warnings --------------------------------------- #
            w(f"## Audit warnings: {len(warned)}\n\n")
            if warned:
                w("These rows were routed correctly from the NCBI lineage, but the "
                  "input columns (`taxon_id` or `category`) disagreed with NCBI. "
                  "Routing is unaffected; this is a data-quality signal for the "
                  "input table.\n\n")
                for d in warned:
                    w(f"- **{d.input_species_name}** "
                      f"(`{d.input_genome_id}` → `{d.ncbi_accession}`)\n")
                    for ww in d.hint_warnings:
                        w(f"    - {ww}\n")
                w("\n")
            else:
                w("None.\n\n")

            # ---- Method footer ---------------------------------------- #
            w("## Method notes\n\n")
            w("- Routing key: **NCBI Assembly accession (GCA/GCF)**. "
              "Input `taxon_id` and `category` are audit fields only.\n")
            w("- AMR / virulence tools (AMRFinderPlus, VFDB) are run **only** "
              "on assemblies whose NCBI lineage contains `Bacteria`.\n")
            w("- Eukaryotic and viral assemblies receive `not_applicable` in the "
              "AMR / virulence columns, with the reason recorded per row.\n")
