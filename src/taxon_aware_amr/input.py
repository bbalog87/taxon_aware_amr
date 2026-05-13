"""
taxon_aware_amr.input
=====================

Reader for the surveillance input table.

Expected columns (header row, tab- or whitespace-separated):

    goat_taxon_id | taxon_id    (either name accepted)
    genome_id
    species_name
    category
    ecological_role
    endemic_status
    iucn_status
    genome_assembly_level

Only the first four are required for routing; the rest are kept verbatim as
audit / context fields and propagated to the output dataframe.
"""

from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

logger = logging.getLogger(__name__)


REQUIRED_COLUMNS = ("genome_id", "species_name", "category")
TAXON_ID_ALIASES = ("taxon_id", "goat_taxon_id", "tax_id")


@dataclass
class InputRow:
    taxon_id:               str
    genome_id:              str
    species_name:           str
    category:               str
    ecological_role:        Optional[str] = None
    endemic_status:         Optional[str] = None
    iucn_status:            Optional[str] = None
    genome_assembly_level:  Optional[str] = None
    extra:                  dict[str, str] = field(default_factory=dict)


def read_input(path: Path) -> Iterator[InputRow]:
    """Yield ``InputRow`` objects from a TSV / TXT input table.

    The reader auto-detects:
      * tab vs multi-space separators (the PDF table used spaces);
      * either ``taxon_id`` or ``goat_taxon_id`` as the first column name.

    Raises ``ValueError`` if required columns are missing.
    """
    with Path(path).open("r", encoding="utf-8") as fh:
        # Sniff the dialect to handle both TSV and whitespace-separated input.
        sample = fh.read(4096)
        fh.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;|") \
            if "\t" in sample or "," in sample else None

        if dialect is not None:
            reader = csv.DictReader(fh, dialect=dialect)
        else:
            # Fall back to splitting on any run of whitespace.
            header = fh.readline().split()
            reader = (dict(zip(header, line.split())) for line in fh if line.strip())

        # Resolve taxon_id column name
        first_row = next(iter(reader), None)
        if first_row is None:
            return
        fieldnames = list(first_row.keys())
        tax_col = next((c for c in TAXON_ID_ALIASES if c in fieldnames), None)
        if tax_col is None:
            raise ValueError(
                f"Input must contain one of {TAXON_ID_ALIASES} columns; "
                f"found {fieldnames}."
            )
        missing = [c for c in REQUIRED_COLUMNS if c not in fieldnames]
        if missing:
            raise ValueError(f"Input missing required columns: {missing}")

        # Re-emit the first row plus the rest
        def _all_rows():
            yield first_row
            yield from reader

        for raw in _all_rows():
            yield InputRow(
                taxon_id=raw.get(tax_col, "").strip(),
                genome_id=raw.get("genome_id", "").strip(),
                species_name=raw.get("species_name", "").strip(),
                category=raw.get("category", "").strip(),
                ecological_role=raw.get("ecological_role", "").strip() or None,
                endemic_status=raw.get("endemic_status", "").strip() or None,
                iucn_status=raw.get("iucn_status", "").strip() or None,
                genome_assembly_level=raw.get("genome_assembly_level", "").strip() or None,
                extra={k: v for k, v in raw.items()
                       if k not in (tax_col, "genome_id", "species_name", "category",
                                    "ecological_role", "endemic_status", "iucn_status",
                                    "genome_assembly_level")},
            )
