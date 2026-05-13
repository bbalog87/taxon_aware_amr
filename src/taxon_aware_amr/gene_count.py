"""
taxon_aware_amr.gene_count
==========================

Compute the ``total_gene_count`` field of the output dataframe.

Two sources, in priority order
------------------------------
1. **NCBI assembly summary** — ``annotation_info.stats.gene_counts.total``.
   This is the authoritative value for assemblies NCBI has annotated. It
   matches the count you'd get by parsing the official GFF, and is already
   in the cached summary JSON, so no download is needed.

2. **GFF3 parsing** — fallback when the summary lacks counts (older or
   custom assemblies). We count features whose ``type`` column equals
   ``gene`` or, optionally, ``pseudogene``. NCBI's annotation puts every
   gene record (protein-coding + ncRNA + pseudogene) as a top-level
   ``gene`` feature, so a single pass through the GFF gives the right total.

When neither source is available
--------------------------------
Returns ``(None, "no_source")``. Callers should leave ``total_gene_count``
as NA in the output and record the reason in the row's audit warnings.

Bacteria without NCBI annotation
--------------------------------
For bacterial assemblies that lack annotation in NCBI, the right answer in
Step 3 is still NA from this module — but Step 4 will run Bakta on the
genome FASTA which produces its own GFF, and that GFF will be counted here
on the second pass. The function is designed to be reusable for that.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Optional, Iterable

logger = logging.getLogger(__name__)

# Feature types NCBI considers "genes" at the top level of a GFF3.
# `pseudogene` is included by NCBI in the total in `gene_counts.total`, so
# we include it too for consistency.
DEFAULT_GENE_TYPES: tuple[str, ...] = ("gene", "pseudogene")


def count_genes_from_gff(
    gff_path: Path,
    *,
    gene_types: Iterable[str] = DEFAULT_GENE_TYPES,
) -> int:
    """Count top-level gene features in a GFF3 file.

    Parameters
    ----------
    gff_path
        Path to a ``.gff``, ``.gff3``, or ``.gff.gz`` file.
    gene_types
        GFF3 type column values to count. Defaults to ``("gene", "pseudogene")``
        to match NCBI's gene_counts.total convention.

    Returns
    -------
    int
        Number of matching features.

    Raises
    ------
    FileNotFoundError, OSError
        On unreadable input.
    ValueError
        If the file is not a GFF (no '##gff' magic and no tab-separated rows).
    """
    types = frozenset(gene_types)
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    seen_data_line = False
    count = 0
    with opener(str(gff_path), "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seen_data_line = True
            if parts[2] in types:
                count += 1
    if not seen_data_line:
        raise ValueError(f"No GFF data rows found in {gff_path}")
    return count


def count_genes(
    *,
    summary_gene_count: Optional[int],
    gff_path: Optional[Path] = None,
) -> tuple[Optional[int], str]:
    """Resolve a row's total_gene_count from the best available source.

    Returns
    -------
    (count, source)
        ``source`` is one of:
            ``"ncbi_summary"``     — used the cached summary JSON
            ``"gff3_parsed"``      — parsed a downloaded GFF
            ``"gff_parse_error"``  — GFF was present but unparseable
            ``"no_source"``        — neither summary count nor GFF available
    """
    if summary_gene_count is not None:
        return summary_gene_count, "ncbi_summary"
    if gff_path is not None and gff_path.exists():
        try:
            return count_genes_from_gff(gff_path), "gff3_parsed"
        except (OSError, ValueError) as e:
            logger.warning("Failed to parse GFF %s: %s", gff_path, e)
            return None, f"gff_parse_error: {e}"
    return None, "no_source"
