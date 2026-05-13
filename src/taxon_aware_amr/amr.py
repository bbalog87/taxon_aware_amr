"""
taxon_aware_amr.amr
===================

Step 4: AMR gene detection on bacterial assemblies.

Tool
----
`AMRFinderPlus <https://github.com/ncbi/amr>`_ — NCBI's curated AMR
detection pipeline. Detects acquired AMR genes by BLAST/HMM against
the NCBI Pathogen Detection Reference Gene Catalog, and point
mutations conferring resistance for supported organisms.

We always invoke it with ``--plus`` to also capture stress (biocide,
metal) and virulence rows. Step 4 filters the cached output to
``Element type == "AMR"``; Step 5 reuses the same cache row for
virulence without re-executing AMRFinderPlus.

Input
-----
The protein FASTA already downloaded in Step 2 (``--include protein``).
AMRFinderPlus is more accurate against translated sequence than against
raw genomic DNA, and the protein file is a fraction of the size.

--organism flag
---------------
When the species is on AMRFinderPlus's supported organism list (M.
tuberculosis, S. aureus, Salmonella, Klebsiella, Vibrio cholerae, ...),
passing ``--organism`` enables higher-confidence point-mutation calling
for that lineage. We map NCBI organism_name → AMRFinderPlus organism
string via :data:`AMRFINDER_ORGANISM_MAP`, with prefix-matching for
strain-annotated names like "Bacillus anthracis str. Ames".

Output
------
A list of :class:`AMRHit` records and an :class:`AMRResult` summary
that the orchestrator merges into ``TaxonDecision`` and the run report.
"""

from __future__ import annotations

import csv
import io
import json
import logging
import re
import subprocess
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional, Protocol

from .cache import AssemblyCache

logger = logging.getLogger(__name__)


# ===========================================================================
# Types
# ===========================================================================

@dataclass
class AMRHit:
    """One row from AMRFinderPlus output."""
    gene_symbol:     str
    element_type:    str    # AMR | STRESS | VIRULENCE
    element_subtype: str    # AMR | POINT | etc.
    drug_class:      str
    subclass:        str
    method:          str    # BLASTX / HMM / POINTX / ALLELEX / ...
    pct_identity:    Optional[float]
    pct_coverage:    Optional[float]
    contig:          str
    start:           Optional[int]
    stop:            Optional[int]
    reference_acc:   str


@dataclass
class AMRResult:
    """Pipeline-level AMR result for a single accession."""
    accession:      str
    status:         str               # ok | failed | not_applicable | not_attempted
    tool:           str               # 'amrfinderplus'
    db_version:     str               # AMRFinderPlus DB version
    organism_flag:  str               # value passed to --organism, or ''
    hits_all:       list[AMRHit] = field(default_factory=list)
    error_message:  Optional[str] = None
    from_cache:     bool = False

    # ---------------- helpers ----------------
    @property
    def amr_hits(self) -> list[AMRHit]:
        """Subset of hits where Element type == 'AMR' (Step 4 scope)."""
        return [h for h in self.hits_all if h.element_type == "AMR"]

    @property
    def virulence_hits(self) -> list[AMRHit]:
        """Subset of hits where Element type == 'VIRULENCE' (used in Step 5)."""
        return [h for h in self.hits_all if h.element_type == "VIRULENCE"]

    @property
    def stress_hits(self) -> list[AMRHit]:
        return [h for h in self.hits_all if h.element_type == "STRESS"]

    def amr_gene_symbols(self) -> list[str]:
        return [h.gene_symbol for h in self.amr_hits]

    def amr_drug_classes(self) -> list[str]:
        seen: list[str] = []
        for h in self.amr_hits:
            if h.drug_class and h.drug_class not in seen:
                seen.append(h.drug_class)
        return seen


# ===========================================================================
# Organism mapping
# ===========================================================================
#
# Generated from `amrfinder -l` output and the NCBI Pathogen Detection
# organism table. Keys are exact NCBI organism names or prefix-matchable
# genus/species; values are the strings AMRFinderPlus expects on its
# `--organism` flag. Maintainers: keep this synced with `amrfinder -l`.

AMRFINDER_ORGANISM_MAP: dict[str, str] = {
    "Acinetobacter baumannii":          "Acinetobacter_baumannii",
    "Burkholderia cepacia":             "Burkholderia_cepacia",
    "Burkholderia pseudomallei":        "Burkholderia_pseudomallei",
    "Campylobacter":                    "Campylobacter",
    "Campylobacter jejuni":             "Campylobacter",
    "Campylobacter coli":               "Campylobacter",
    "Citrobacter freundii":             "Citrobacter_freundii",
    "Clostridioides difficile":         "Clostridioides_difficile",
    "Enterobacter asburiae":            "Enterobacter_asburiae",
    "Enterobacter cloacae":             "Enterobacter_cloacae",
    "Enterococcus faecalis":            "Enterococcus_faecalis",
    "Enterococcus faecium":             "Enterococcus_faecium",
    "Enterococcus hirae":               "Enterococcus_hirae",
    "Escherichia coli":                 "Escherichia",
    "Escherichia":                      "Escherichia",
    "Klebsiella oxytoca":               "Klebsiella_oxytoca",
    "Klebsiella pneumoniae":            "Klebsiella_pneumoniae",
    "Mycobacterium tuberculosis":       "Mycobacterium_tuberculosis",
    "Neisseria gonorrhoeae":            "Neisseria_gonorrhoeae",
    "Neisseria meningitidis":           "Neisseria_meningitidis",
    "Pseudomonas aeruginosa":           "Pseudomonas_aeruginosa",
    "Salmonella enterica":              "Salmonella",
    "Salmonella":                       "Salmonella",
    "Serratia marcescens":              "Serratia_marcescens",
    "Staphylococcus aureus":            "Staphylococcus_aureus",
    "Staphylococcus pseudintermedius":  "Staphylococcus_pseudintermedius",
    "Streptococcus agalactiae":         "Streptococcus_agalactiae",
    "Streptococcus pneumoniae":         "Streptococcus_pneumoniae",
    "Streptococcus pyogenes":           "Streptococcus_pyogenes",
    "Vibrio cholerae":                  "Vibrio_cholerae",
    "Vibrio parahaemolyticus":          "Vibrio_parahaemolyticus",
    "Vibrio vulnificus":                "Vibrio_vulnificus",
}


def amrfinder_organism_for(ncbi_organism_name: Optional[str]) -> Optional[str]:
    """Resolve an NCBI organism name to AMRFinderPlus's ``--organism`` value.

    Tries exact match, then genus-species (strips strain info), then genus.
    Returns ``None`` if no rule fires; the caller then runs AMRFinderPlus
    without ``--organism`` (still detects acquired genes, just not point
    mutations).

    Examples
    --------
    >>> amrfinder_organism_for("Mycobacterium tuberculosis")
    'Mycobacterium_tuberculosis'
    >>> amrfinder_organism_for("Bacillus anthracis str. Ames")  # no rule
    >>> amrfinder_organism_for("Staphylococcus aureus subsp. aureus N315")
    'Staphylococcus_aureus'
    """
    if not ncbi_organism_name:
        return None
    name = ncbi_organism_name.strip()
    if name in AMRFINDER_ORGANISM_MAP:
        return AMRFINDER_ORGANISM_MAP[name]
    parts = name.split()
    if len(parts) >= 2:
        gs = " ".join(parts[:2])
        if gs in AMRFINDER_ORGANISM_MAP:
            return AMRFINDER_ORGANISM_MAP[gs]
    if parts and parts[0] in AMRFINDER_ORGANISM_MAP:
        return AMRFINDER_ORGANISM_MAP[parts[0]]
    return None


# ===========================================================================
# Parsing AMRFinderPlus TSV output
# ===========================================================================
#
# AMRFinderPlus writes a header row with these columns (as of v4):
#   Protein identifier, Contig id, Start, Stop, Strand, Gene symbol,
#   Sequence name, Scope, Element type, Element subtype, Class, Subclass,
#   Method, Target length, Reference sequence length,
#   % Coverage of reference sequence, % Identity to reference sequence,
#   Alignment length, Accession of closest sequence, Name of closest sequence,
#   HMM id, HMM description

def parse_amrfinder_tsv(text: str) -> list[AMRHit]:
    """Parse AMRFinderPlus TSV output into a list of AMRHit records."""
    out: list[AMRHit] = []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    for row in reader:
        # Defensive: AMRFinderPlus column names can shift; canonicalise.
        def g(*aliases: str) -> str:
            for a in aliases:
                if a in row and row[a]:
                    return row[a]
            return ""

        def gf(*aliases: str) -> Optional[float]:
            v = g(*aliases)
            try:
                return float(v) if v not in ("", "NA") else None
            except ValueError:
                return None

        def gi(*aliases: str) -> Optional[int]:
            v = g(*aliases)
            try:
                return int(v) if v not in ("", "NA") else None
            except ValueError:
                return None

        out.append(AMRHit(
            gene_symbol=g("Gene symbol", "Element symbol"),
            element_type=g("Element type"),
            element_subtype=g("Element subtype"),
            drug_class=g("Class"),
            subclass=g("Subclass"),
            method=g("Method"),
            pct_identity=gf("% Identity to reference sequence",
                            "% Identity to reference"),
            pct_coverage=gf("% Coverage of reference sequence",
                            "% Coverage of reference"),
            contig=g("Contig id"),
            start=gi("Start"),
            stop=gi("Stop"),
            reference_acc=g("Accession of closest sequence",
                            "Closest reference accession"),
        ))
    return out


def _hits_to_json(hits: list[AMRHit]) -> str:
    return json.dumps([asdict(h) for h in hits])


def _hits_from_json(text: str) -> list[AMRHit]:
    return [AMRHit(**d) for d in json.loads(text)]


# ===========================================================================
# Client protocol + cache helpers
# ===========================================================================

class AMRClient(Protocol):
    def run(
        self, accession: str, protein_path: Path, *,
        organism_flag: Optional[str],
        refresh: bool = False,
    ) -> AMRResult: ...


def _cached_result(cache: AssemblyCache, accession: str) -> Optional[AMRResult]:
    """Look up a cached AMRFinderPlus run for this accession."""
    row = cache._conn.execute(   # cache exposes _conn intentionally
        "SELECT * FROM tool_results WHERE accession = ? AND tool = ?",
        (accession, "amrfinderplus"),
    ).fetchone()
    if row is None:
        return None
    hits = _hits_from_json(row["parsed_json"]) if row["parsed_json"] else []
    return AMRResult(
        accession=accession,
        status=row["status"],
        tool="amrfinderplus",
        db_version=row["db_version"] or "",
        organism_flag=row["organism_flag"] or "",
        hits_all=hits,
        error_message=row["error_message"],
        from_cache=True,
    )


def _persist_result(cache: AssemblyCache, r: AMRResult) -> None:
    cache._conn.execute(
        """
        INSERT OR REPLACE INTO tool_results
            (accession, tool, completed_at, status, db_version, organism_flag,
             raw_output, parsed_json, error_message)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            r.accession, r.tool,
            datetime.now(timezone.utc).isoformat(timespec="seconds"),
            r.status, r.db_version, r.organism_flag,
            None,                            # raw_output reserved for Step 5
            _hits_to_json(r.hits_all),
            r.error_message,
        ),
    )


# ===========================================================================
# Production client — shells out to amrfinder
# ===========================================================================

_VERSION_RE = re.compile(r"AMRFinderPlus database version[:\s]+(\S+)", re.IGNORECASE)


class AMRFinderPlusClient:
    """Run AMRFinderPlus via the ``amrfinder`` binary on bacterial assemblies.

    Invoked in **combined mode** (``-n genome -p protein -g annotation.gff``)
    for highest sensitivity: protein input catches confirmed gene calls fast,
    nucleotide input catches fragmented or unannotated ORFs, and the GFF maps
    nucleotide hits back to gene names. This is NCBI's recommended mode.

    Combined mode also enables the full POINTX (point mutation) catalogue for
    supported organisms (Mycobacterium tuberculosis, Salmonella, etc.) — some
    point mutations cannot be confidently called in protein-only mode.
    """

    def __init__(
        self, cache: AssemblyCache, *,
        binary: str = "amrfinder",
        threads: int = 4,
        timeout_seconds: int = 1800,
        plus: bool = True,
    ) -> None:
        self.cache = cache
        self.binary = binary
        self.threads = threads
        self.timeout_seconds = timeout_seconds
        self.plus = plus

    def run(
        self, accession: str,
        protein_path: Path,
        *,
        organism_flag: Optional[str],
        genome_path: Optional[Path] = None,
        gff_path: Optional[Path] = None,
        refresh: bool = False,
    ) -> AMRResult:
        if not refresh:
            cached = _cached_result(self.cache, accession)
            if cached is not None:
                return cached

        cmd = [self.binary, "--threads", str(self.threads), "-p", str(protein_path)]
        # Combined mode when we have all three inputs (the new default).
        if genome_path is not None:
            cmd += ["-n", str(genome_path)]
        if gff_path is not None:
            cmd += ["-g", str(gff_path), "--annotation_format", "ncbi"]
        if self.plus:
            cmd.append("--plus")
        if organism_flag:
            cmd += ["--organism", organism_flag]

        db_version = self._detect_db_version()
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
            r = AMRResult(
                accession=accession, status="failed",
                tool="amrfinderplus", db_version=db_version,
                organism_flag=organism_flag or "",
                error_message=msg,
            )
            _persist_result(self.cache, r)
            return r

        hits = parse_amrfinder_tsv(proc.stdout)
        r = AMRResult(
            accession=accession, status="ok",
            tool="amrfinderplus", db_version=db_version,
            organism_flag=organism_flag or "",
            hits_all=hits,
        )
        _persist_result(self.cache, r)
        return r

    def _detect_db_version(self) -> str:
        """Best-effort: parse `amrfinder --version` for the DB version."""
        try:
            proc = subprocess.run(
                [self.binary, "--version"], capture_output=True,
                text=True, timeout=15,
            )
            m = _VERSION_RE.search(proc.stdout + proc.stderr)
            if m:
                return m.group(1)
            return proc.stdout.strip().splitlines()[0] if proc.stdout else "unknown"
        except (subprocess.SubprocessError, FileNotFoundError):
            return "unknown"


# ===========================================================================
# Mock client — for offline testing
# ===========================================================================

class MockAMRClient:
    """In-memory AMR client that returns canned hits per accession.

    Same interface as ``AMRFinderPlusClient``. The bundled fixtures simulate
    realistic AMR profiles for the three bacterial demo accessions.
    """

    def __init__(
        self, cache: AssemblyCache, *,
        mock_hits: dict[str, list[AMRHit]],
        db_version: str = "2026-01-15.1 (mock)",
    ) -> None:
        self.cache = cache
        self.mock_hits = mock_hits
        self.db_version = db_version

    def run(
        self, accession: str,
        protein_path: Path,
        *,
        organism_flag: Optional[str],
        genome_path: Optional[Path] = None,
        gff_path: Optional[Path] = None,
        refresh: bool = False,
    ) -> AMRResult:
        if not refresh:
            cached = _cached_result(self.cache, accession)
            if cached is not None:
                return cached
        hits = list(self.mock_hits.get(accession, []))
        r = AMRResult(
            accession=accession, status="ok",
            tool="amrfinderplus", db_version=self.db_version,
            organism_flag=organism_flag or "",
            hits_all=hits,
        )
        _persist_result(self.cache, r)
        return r
