"""
taxon_aware_amr.virulence
=========================

Step 5: virulence gene detection on bacterial assemblies.

Two evidence sources, merged
----------------------------
1. **VFDB via ABRicate** — ``abricate --db vfdb <genome.fna>``. VFDB is the
   most comprehensive curated virulence factor database for bacterial
   pathogens (Chen et al., MGC.AC.CN). ABRicate is a lightweight BLAST
   wrapper that runs the screen quickly and outputs a stable TSV.

2. **AMRFinderPlus ``--plus`` cache rows** — AMRFinderPlus also reports
   selected virulence genes when ``--plus`` is on. Step 4 already invoked
   AMRFinderPlus with ``--plus`` and cached the full output. Step 5 pulls
   the ``Element type == 'VIRULENCE'`` rows out of the cache **without
   re-executing AMRFinderPlus**.

The two sources are deduplicated by gene symbol (case-insensitive). Both
methods are recorded per row so downstream readers see whether a gene was
confirmed by both DBs or only one.

Bacterial only
--------------
Same rule as AMR: VFDB is a bacterial-pathogen database. Eukaryotic and
viral assemblies receive ``virulence_genes = "not_applicable"``.
"""

from __future__ import annotations

import csv
import io
import json
import logging
import subprocess
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Protocol

from .amr   import AMRHit, _cached_result as cached_amr_result
from .cache import AssemblyCache

logger = logging.getLogger(__name__)


# ===========================================================================
# Types
# ===========================================================================

@dataclass
class VirulenceHit:
    gene_symbol:   str
    source:        str     # 'vfdb_abricate' | 'amrfinderplus_plus'
    product:       str
    db:            str
    pct_identity:  Optional[float]
    pct_coverage:  Optional[float]
    contig:        str
    start:         Optional[int]
    stop:          Optional[int]


@dataclass
class VirulenceResult:
    accession:     str
    status:        str            # ok | failed | not_applicable | not_attempted
    sources_used:  list[str] = field(default_factory=list)   # which DBs we consulted
    db_versions:   dict[str, str] = field(default_factory=dict)
    hits:          list[VirulenceHit] = field(default_factory=list)
    error_message: Optional[str] = None
    from_cache:    bool = False

    def gene_symbols(self, *, deduplicate: bool = True) -> list[str]:
        if not deduplicate:
            return [h.gene_symbol for h in self.hits]
        seen: list[str] = []
        seen_lower: set[str] = set()
        for h in self.hits:
            key = h.gene_symbol.lower()
            if key not in seen_lower:
                seen.append(h.gene_symbol)
                seen_lower.add(key)
        return seen


# ===========================================================================
# ABRicate output parser
# ===========================================================================
#
# ABRicate TSV columns (v1.x):
#   #FILE, SEQUENCE, START, END, STRAND, GENE, COVERAGE, COVERAGE_MAP, GAPS,
#   %COVERAGE, %IDENTITY, DATABASE, ACCESSION, PRODUCT, RESISTANCE

def parse_abricate_tsv(text: str, *, db_label: str = "vfdb") -> list[VirulenceHit]:
    """Parse ABRicate TSV output into VirulenceHit records."""
    out: list[VirulenceHit] = []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        def g(k: str) -> str:
            return row.get(k, "") or row.get(k.lstrip("#"), "")
        def gf(k: str) -> Optional[float]:
            v = g(k)
            try:
                return float(v) if v not in ("", "NA") else None
            except ValueError:
                return None
        def gi(k: str) -> Optional[int]:
            v = g(k)
            try:
                return int(v) if v not in ("", "NA") else None
            except ValueError:
                return None

        out.append(VirulenceHit(
            gene_symbol=g("GENE"),
            source=f"{db_label}_abricate",
            product=g("PRODUCT"),
            db=g("DATABASE") or db_label,
            pct_identity=gf("%IDENTITY"),
            pct_coverage=gf("%COVERAGE"),
            contig=g("SEQUENCE"),
            start=gi("START"),
            stop=gi("END"),
        ))
    return out


def _hits_to_json(hits: list[VirulenceHit]) -> str:
    return json.dumps([asdict(h) for h in hits])


def _hits_from_json(text: str) -> list[VirulenceHit]:
    return [VirulenceHit(**d) for d in json.loads(text)]


# ===========================================================================
# AMRFinderPlus --plus → VirulenceHit bridge
# ===========================================================================

def virulence_hits_from_amrfinder_cache(
    cache: AssemblyCache, accession: str,
) -> tuple[list[VirulenceHit], str]:
    """Pull ``Element type == 'VIRULENCE'`` rows out of Step 4's cache.

    Returns
    -------
    (hits, db_version)
    """
    cached = cached_amr_result(cache, accession)
    if cached is None or cached.status != "ok":
        return [], ""
    out: list[VirulenceHit] = []
    for h in cached.virulence_hits:
        out.append(VirulenceHit(
            gene_symbol=h.gene_symbol,
            source="amrfinderplus_plus",
            product=h.drug_class or "",   # AMRFinder repurposes 'Class' for virulence subtype
            db="amrfinderplus_plus",
            pct_identity=h.pct_identity,
            pct_coverage=h.pct_coverage,
            contig=h.contig,
            start=h.start,
            stop=h.stop,
        ))
    return out, cached.db_version


# ===========================================================================
# ABRicate client (real + mock)
# ===========================================================================

class ABRicateClient(Protocol):
    def run(
        self, accession: str, genome_path: Path, *,
        db: str = "vfdb", refresh: bool = False,
    ) -> tuple[list[VirulenceHit], str, Optional[str]]: ...
    """Returns (hits, db_version, error_message)."""


class ABRicateCLIClient:
    """Run ABRicate via the ``abricate`` binary."""

    def __init__(
        self, cache: AssemblyCache, *,
        binary: str = "abricate",
        threads: int = 4,
        timeout_seconds: int = 1200,
    ) -> None:
        self.cache = cache
        self.binary = binary
        self.threads = threads
        self.timeout_seconds = timeout_seconds

    def run(
        self, accession: str, genome_path: Path, *,
        db: str = "vfdb", refresh: bool = False,
    ) -> tuple[list[VirulenceHit], str, Optional[str]]:
        # Check cache (key: accession + tool name 'abricate_<db>')
        tool_key = f"abricate_{db}"
        if not refresh:
            cached = self._cached(accession, tool_key)
            if cached is not None:
                return cached

        cmd = [self.binary, "--db", db, "--threads", str(self.threads),
               "--quiet", str(genome_path)]
        try:
            proc = subprocess.run(
                cmd, check=True, capture_output=True, text=True,
                timeout=self.timeout_seconds,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
                FileNotFoundError) as e:
            msg = str(e)
            if isinstance(e, subprocess.CalledProcessError):
                msg = (e.stderr or e.stdout or "").strip().splitlines()[-1:] or [str(e)]
                msg = msg[0]
            self._persist(accession, tool_key, status="failed",
                          db_version="", hits=[], error_message=msg)
            return [], "", msg

        hits = parse_abricate_tsv(proc.stdout, db_label=db)
        db_version = self._detect_db_version(db)
        self._persist(accession, tool_key, status="ok",
                      db_version=db_version, hits=hits)
        return hits, db_version, None

    def _detect_db_version(self, db: str) -> str:
        try:
            proc = subprocess.run(
                [self.binary, "--list"], capture_output=True, text=True, timeout=10,
            )
            for line in proc.stdout.splitlines():
                if line.startswith(db) or db in line.split()[:2]:
                    # `abricate --list` shows: db, sequences, datadir, format, date
                    return line.strip()
            return "unknown"
        except (subprocess.SubprocessError, FileNotFoundError):
            return "unknown"

    def _cached(self, accession: str, tool: str):
        row = self.cache._conn.execute(
            "SELECT * FROM tool_results WHERE accession=? AND tool=?",
            (accession, tool),
        ).fetchone()
        if row is None:
            return None
        hits = _hits_from_json(row["parsed_json"]) if row["parsed_json"] else []
        return hits, row["db_version"] or "", row["error_message"]

    def _persist(self, accession: str, tool: str, *,
                 status: str, db_version: str,
                 hits: list[VirulenceHit],
                 error_message: Optional[str] = None) -> None:
        self.cache._conn.execute(
            """
            INSERT OR REPLACE INTO tool_results
                (accession, tool, completed_at, status, db_version, organism_flag,
                 raw_output, parsed_json, error_message)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (accession, tool,
             datetime.now(timezone.utc).isoformat(timespec="seconds"),
             status, db_version, "",
             None, _hits_to_json(hits), error_message),
        )


class MockABRicateClient:
    """Returns canned VFDB hits per accession; for offline testing."""

    def __init__(
        self, cache: AssemblyCache, *,
        mock_hits: dict[str, list[VirulenceHit]],
        db_version: str = "vfdb_2026-01 (mock)",
    ) -> None:
        self.cache = cache
        self.mock_hits = mock_hits
        self.db_version = db_version

    def run(
        self, accession: str, genome_path: Path, *,
        db: str = "vfdb", refresh: bool = False,
    ) -> tuple[list[VirulenceHit], str, Optional[str]]:
        tool_key = f"abricate_{db}"
        if not refresh:
            row = self.cache._conn.execute(
                "SELECT * FROM tool_results WHERE accession=? AND tool=?",
                (accession, tool_key),
            ).fetchone()
            if row is not None:
                hits = _hits_from_json(row["parsed_json"]) if row["parsed_json"] else []
                return hits, row["db_version"] or "", row["error_message"]
        hits = list(self.mock_hits.get(accession, []))
        self.cache._conn.execute(
            """
            INSERT OR REPLACE INTO tool_results
                (accession, tool, completed_at, status, db_version, organism_flag,
                 raw_output, parsed_json, error_message)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (accession, tool_key,
             datetime.now(timezone.utc).isoformat(timespec="seconds"),
             "ok", self.db_version, "", None, _hits_to_json(hits), None),
        )
        return hits, self.db_version, None


# ===========================================================================
# Orchestration helper — merge VFDB + AMRFinderPlus cache for one row
# ===========================================================================

def detect_virulence(
    accession: str,
    *,
    cache: AssemblyCache,
    genome_path: Optional[Path],
    abricate_client: Optional[ABRicateClient],
) -> VirulenceResult:
    """Compute the merged virulence profile for a bacterial accession.

    Always pulls AMRFinderPlus --plus virulence rows out of the Step 4 cache.
    Additionally runs ABRicate-VFDB when an abricate_client and a genome
    FASTA are available.
    """
    sources_used: list[str] = []
    db_versions: dict[str, str] = {}
    all_hits: list[VirulenceHit] = []
    error: Optional[str] = None

    # 1. AMRFinderPlus --plus virulence (already cached by Step 4)
    af_hits, af_dbver = virulence_hits_from_amrfinder_cache(cache, accession)
    if af_hits or af_dbver:
        sources_used.append("amrfinderplus_plus")
        db_versions["amrfinderplus_plus"] = af_dbver
        all_hits.extend(af_hits)

    # 2. ABRicate-VFDB
    if abricate_client is not None and genome_path is not None:
        vfdb_hits, vfdb_ver, vfdb_err = abricate_client.run(
            accession, genome_path, db="vfdb",
        )
        sources_used.append("vfdb_abricate")
        db_versions["vfdb_abricate"] = vfdb_ver
        all_hits.extend(vfdb_hits)
        if vfdb_err:
            error = f"ABRicate-VFDB failed: {vfdb_err}"
    elif abricate_client is not None and genome_path is None:
        error = ("ABRicate not run: no genome FASTA available "
                 "(check Step 2 download plan)")

    if not sources_used:
        return VirulenceResult(
            accession=accession, status="not_attempted",
            sources_used=[], db_versions={},
            error_message="No virulence source available",
        )

    return VirulenceResult(
        accession=accession,
        status="failed" if error else "ok",
        sources_used=sources_used,
        db_versions=db_versions,
        hits=all_hits,
        error_message=error,
    )
