"""
taxon_aware_amr.ncbi
====================

NCBI client.

This module is the only place the pipeline touches the network. It exposes a
protocol (``NCBIClient``) with two production-grade implementations:

    * ``DatasetsCLIClient`` — wraps the official ``datasets`` CLI from NCBI.
    * ``MockNCBIClient``    — emits realistic JSON for offline testing.

Operations
----------
    fetch_summary(accession)              → AssemblyFetch
    download_files(accession, includes)   → DownloadResult
    cleanup_files(accession)              → int  (bytes freed on disk)

All operations consult and update the ``AssemblyCache``. Successful summaries
are cached indefinitely; failures are cached briefly so transient NCBI outages
don't poison the cache.

Failure taxonomy
----------------
We map raw NCBI failures to four categories so the run report can explain
*why* a row was skipped:

    not_found      — NCBI returned no record (typo, retracted accession)
    suppressed     — NCBI returned a record but flagged it suppressed/withdrawn
    network_error  — transient: timeout, connectivity, rate-limit
    parse_error    — datasets returned JSON we couldn't interpret

The first two are durable failures; the latter two are typically transient
and re-tried on the next run after the cache TTL expires.
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional, Protocol

from .assembly import AssemblyRecord
from .cache import AssemblyCache, SummaryEntry

logger = logging.getLogger(__name__)


# ===========================================================================
# Types
# ===========================================================================

VALID_INCLUDES: frozenset[str] = frozenset({
    "genome",       # FASTA of the assembled sequences
    "gff3",         # gene-level annotation
    "protein",      # predicted protein FASTA
    "cds",          # CDS nucleotide FASTA
    "rna",          # RNA FASTA
    "seq-report",   # per-sequence chromosome/scaffold map
})


@dataclass
class AssemblyFetch:
    """Result of a summary lookup, success or failure."""
    accession:        str
    status:           str                       # ok | not_found | suppressed | network_error | parse_error
    assembly_record:  Optional[AssemblyRecord]
    gene_count_total: Optional[int]             # from NCBI annotation_info.stats if present
    error_message:    Optional[str]
    from_cache:       bool


@dataclass
class DownloadResult:
    accession:  str
    files:      dict[str, Path] = field(default_factory=dict)   # file_type → path on disk
    bytes_:     int = 0
    skipped:    bool = False
    reason:     Optional[str] = None


# ===========================================================================
# Protocol
# ===========================================================================

class NCBIClient(Protocol):
    """Interface shared by DatasetsCLIClient and MockNCBIClient."""

    def fetch_summary(self, accession: str, *, refresh: bool = False) -> AssemblyFetch: ...
    def download_files(
        self, accession: str, includes: Iterable[str], out_dir: Path,
    ) -> DownloadResult: ...
    def cleanup_files(self, accession: str) -> int: ...


# ===========================================================================
# Parsing helpers — used by both real and mock clients
# ===========================================================================

def _parse_summary_json(
    accession: str, summary_obj: dict, lineage_obj: Optional[dict],
) -> tuple[Optional[AssemblyRecord], Optional[int], Optional[str]]:
    """Turn the raw datasets JSON into (AssemblyRecord, gene_count, error_msg).

    Designed to be lenient — datasets has changed its output shape across
    major versions. We try a handful of well-known paths and fall back to
    a parse_error if nothing works.
    """
    try:
        # Suppressed / withdrawn flag (newer datasets versions include this).
        assembly_info = summary_obj.get("assembly_info", {}) or {}
        if assembly_info.get("suppression_reason"):
            return (None, None,
                    f"Suppressed by NCBI: {assembly_info['suppression_reason']}")

        organism = summary_obj.get("organism", {}) or {}
        tax_id = organism.get("tax_id")
        organism_name = organism.get("organism_name") or organism.get("name")
        if not tax_id or not organism_name:
            return (None, None,
                    "Summary JSON lacks organism.tax_id or organism_name.")

        # Lineage — try the dedicated lineage record first, then a fallback
        # within the genome summary itself.
        lineage_names: list[str] = []
        if lineage_obj:
            # Modern `datasets summary taxonomy taxon X --as-json-lines` puts
            # ranked lineage under `taxonomy.classification`. Some shapes wrap
            # under `reports[0].taxonomy`. Cover both.
            tax_blob = (
                lineage_obj.get("taxonomy")
                or (lineage_obj.get("reports") or [{}])[0].get("taxonomy")
                or {}
            )
            classification = (
                tax_blob.get("classification")
                or lineage_obj.get("classification")
                or {}
            )
            # NCBI renamed `superkingdom` → `domain` in 2023, and uses `realm`
            # for viruses. Enumerate the known ranks in hierarchical order so
            # the lineage tuple reads top-down, then catch any extras the
            # datasets CLI might add.
            KNOWN_RANKS = (
                "domain", "superkingdom", "realm", "acellular_root",
                "kingdom", "subkingdom",
                "phylum", "subphylum",
                "class", "subclass",
                "order", "suborder",
                "family", "subfamily",
                "genus", "subgenus",
                "species", "subspecies", "strain",
            )
            seen_keys: set[str] = set()
            for rank in KNOWN_RANKS:
                node = classification.get(rank) if isinstance(classification, dict) else None
                if isinstance(node, dict) and node.get("name"):
                    lineage_names.append(node["name"])
                    seen_keys.add(rank)
            if isinstance(classification, dict):
                for rank, node in classification.items():
                    if rank in seen_keys:
                        continue
                    if isinstance(node, dict) and node.get("name"):
                        lineage_names.append(node["name"])
            # As a last-ditch fallback, some shapes expose a flat
            # `ranked_lineage` list of {rank, name} dicts.
            if not lineage_names:
                for item in tax_blob.get("ranked_lineage", []) or []:
                    if isinstance(item, dict) and item.get("name"):
                        lineage_names.append(item["name"])
        if not lineage_names:
            # Older summary shape: organism.lineage is a list of dicts or strings.
            for item in organism.get("lineage", []) or []:
                if isinstance(item, dict):
                    name = item.get("name")
                    if name:
                        lineage_names.append(name)
                elif isinstance(item, str):
                    lineage_names.append(item)

        if not lineage_names:
            # Diagnostic: emit the raw lineage_obj at DEBUG so users running
            # with -vv can see what NCBI actually returned and report a
            # response shape we don't yet handle.
            logger.debug(
                "Lineage extraction returned empty for %s. "
                "lineage_obj keys: %s; classification keys: %s. "
                "Raw lineage_obj (truncated): %r",
                accession,
                list(lineage_obj.keys()) if isinstance(lineage_obj, dict) else type(lineage_obj),
                list(classification.keys()) if isinstance(classification, dict) else type(classification),
                (json.dumps(lineage_obj)[:1000] + "...") if lineage_obj else None,
            )
            return (None, None,
                    "Could not extract a lineage from datasets summary "
                    "(neither taxonomy classification nor organism.lineage present).")

        # Assembly level — normalise to our snake_case convention.
        raw_level = assembly_info.get("assembly_level", "") or ""
        level = raw_level.lower().replace(" ", "_")

        record = AssemblyRecord(
            accession=accession,
            tax_id=int(tax_id),
            organism_name=organism_name,
            lineage=tuple(lineage_names),
            assembly_level=level or "unknown",
            source="ncbi_datasets",
        )

        # Gene count — when NCBI has annotated the assembly, we can grab the
        # total directly and skip the GFF download.
        gene_count = None
        ai = summary_obj.get("annotation_info") or {}
        stats = ai.get("stats") or {}
        gc = stats.get("gene_counts") or {}
        if "total" in gc:
            try:
                gene_count = int(gc["total"])
            except (TypeError, ValueError):
                gene_count = None

        return record, gene_count, None

    except (KeyError, TypeError, ValueError) as exc:
        return (None, None, f"parse_error: {exc!r}")


def _fetch_from_cache_entry(entry: SummaryEntry) -> AssemblyFetch:
    """Reconstruct an AssemblyFetch from a cached SummaryEntry."""
    if entry.status != "ok":
        return AssemblyFetch(
            accession=entry.accession,
            status=entry.status,
            assembly_record=None,
            gene_count_total=None,
            error_message=entry.error_message,
            from_cache=True,
        )
    summary_obj = json.loads(entry.summary_json) if entry.summary_json else {}
    lineage_obj = json.loads(entry.lineage_json) if entry.lineage_json else None
    record, gene_count, err = _parse_summary_json(
        entry.accession, summary_obj, lineage_obj,
    )
    if record is None:
        return AssemblyFetch(
            accession=entry.accession, status="parse_error",
            assembly_record=None, gene_count_total=None,
            error_message=err, from_cache=True,
        )
    return AssemblyFetch(
        accession=entry.accession, status="ok",
        assembly_record=record, gene_count_total=gene_count,
        error_message=None, from_cache=True,
    )


# ===========================================================================
# Production client — shells out to `datasets`
# ===========================================================================

class DatasetsCLIClient:
    """NCBI client built on the ``datasets`` CLI.

    External dependency: ``datasets`` (https://github.com/ncbi/datasets).

    Production use::

        cache  = AssemblyCache(Path("cache/cache.sqlite"))
        client = DatasetsCLIClient(cache=cache, download_dir=Path("cache/genomes"))
        fetch  = client.fetch_summary("GCA_964341285.1")
    """

    def __init__(
        self,
        cache: AssemblyCache,
        download_dir: Path,
        *,
        binary: str = "datasets",
        timeout_seconds: int = 180,
    ) -> None:
        self.cache = cache
        self.download_dir = Path(download_dir)
        self.binary = binary
        self.timeout_seconds = timeout_seconds
        self.download_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # fetch_summary                                                       #
    # ------------------------------------------------------------------ #
    def fetch_summary(self, accession: str, *, refresh: bool = False) -> AssemblyFetch:
        if not refresh:
            cached = self.cache.get_summary(accession)
            if cached:
                return _fetch_from_cache_entry(cached)

        # Live call: datasets summary genome accession <acc> --as-json-lines
        try:
            summary_text = self._run([
                self.binary, "summary", "genome", "accession", accession,
                "--as-json-lines",
            ])
        except _SubprocessFailure as e:
            self._persist_failure(accession, "network_error", str(e))
            return AssemblyFetch(accession, "network_error", None, None, str(e), False)

        if not summary_text.strip():
            msg = "NCBI returned no record for this accession."
            self._persist_failure(accession, "not_found", msg)
            return AssemblyFetch(accession, "not_found", None, None, msg, False)

        try:
            summary_obj = json.loads(summary_text.splitlines()[0])
        except json.JSONDecodeError as e:
            self._persist_failure(accession, "parse_error", str(e))
            return AssemblyFetch(accession, "parse_error", None, None, str(e), False)

        # Lineage call: datasets summary taxonomy taxon <tax_id>
        tax_id = (summary_obj.get("organism") or {}).get("tax_id")
        lineage_obj: Optional[dict] = None
        if tax_id:
            try:
                lineage_text = self._run([
                    self.binary, "summary", "taxonomy", "taxon", str(tax_id),
                    "--as-json-lines",
                ])
                if lineage_text.strip():
                    lineage_obj = json.loads(lineage_text.splitlines()[0])
            except (_SubprocessFailure, json.JSONDecodeError) as e:
                # Non-fatal: we have the summary, just no lineage. The parser
                # will try to fall back to organism.lineage.
                logger.warning("Lineage fetch failed for tax_id %s: %s", tax_id, e)

        record, gene_count, err = _parse_summary_json(
            accession, summary_obj, lineage_obj,
        )
        if record is None:
            # 'suppressed' is detected inside _parse_summary_json
            status = "suppressed" if err and err.startswith("Suppressed") else "parse_error"
            self._persist_failure(accession, status, err or "Unknown parse error")
            return AssemblyFetch(accession, status, None, None, err, False)

        # Success: cache it
        self.cache.put_summary(SummaryEntry(
            accession=accession,
            fetched_at=datetime.now(timezone.utc),
            status="ok",
            summary_json=json.dumps(summary_obj),
            lineage_json=json.dumps(lineage_obj) if lineage_obj else None,
            error_message=None,
        ))
        return AssemblyFetch(accession, "ok", record, gene_count, None, False)

    # ------------------------------------------------------------------ #
    # download_files                                                      #
    # ------------------------------------------------------------------ #
    def download_files(
        self, accession: str, includes: Iterable[str], out_dir: Path,
    ) -> DownloadResult:
        includes = list(includes)
        bad = [i for i in includes if i not in VALID_INCLUDES]
        if bad:
            raise ValueError(f"Unsupported include(s): {bad}. "
                             f"Allowed: {sorted(VALID_INCLUDES)}")

        target = Path(out_dir) / accession
        target.mkdir(parents=True, exist_ok=True)
        zip_path = target / "ncbi_dataset.zip"

        cmd = [
            self.binary, "download", "genome", "accession", accession,
            "--include", ",".join(includes),
            "--filename", str(zip_path),
            "--no-progressbar",
        ]
        try:
            self._run(cmd)
        except _SubprocessFailure as e:
            return DownloadResult(accession=accession, skipped=True, reason=str(e))

        # Unzip
        try:
            shutil.unpack_archive(str(zip_path), str(target))
            zip_path.unlink(missing_ok=True)
        except (shutil.ReadError, FileNotFoundError) as e:
            return DownloadResult(accession=accession, skipped=True,
                                  reason=f"Unzip failed: {e}")

        # Inventory what we actually got
        files: dict[str, Path] = {}
        total_bytes = 0
        data_dir = target / "ncbi_dataset" / "data" / accession
        for path in data_dir.rglob("*") if data_dir.exists() else []:
            if not path.is_file():
                continue
            ft = _classify_downloaded_file(path)
            if ft:
                files[ft] = path
                size = path.stat().st_size
                total_bytes += size
                self.cache.record_file(accession, ft, path, size)

        return DownloadResult(accession=accession, files=files, bytes_=total_bytes)

    # ------------------------------------------------------------------ #
    # cleanup_files                                                       #
    # ------------------------------------------------------------------ #
    def cleanup_files(self, accession: str) -> int:
        """Delete all on-disk files for this accession; return bytes freed."""
        rows = self.cache.list_files(accession, include_deleted=False)
        freed = 0
        for r in rows:
            try:
                # Delete the whole accession directory once, then break.
                acc_dir = self.download_dir / accession
                if acc_dir.exists():
                    freed += sum(p.stat().st_size for p in acc_dir.rglob("*") if p.is_file())
                    shutil.rmtree(acc_dir)
                break
            except OSError as e:
                logger.warning("Cleanup error for %s: %s", accession, e)
        self.cache.mark_files_deleted(accession)
        return freed

    # ------------------------------------------------------------------ #
    # Internals                                                           #
    # ------------------------------------------------------------------ #
    def _run(self, cmd: list[str]) -> str:
        logger.debug("EXEC: %s", " ".join(cmd))
        try:
            result = subprocess.run(
                cmd, check=True, capture_output=True, text=True,
                timeout=self.timeout_seconds,
            )
            return result.stdout
        except subprocess.TimeoutExpired as e:
            raise _SubprocessFailure(
                f"Timed out after {self.timeout_seconds}s: {e.cmd}") from e
        except subprocess.CalledProcessError as e:
            msg = (e.stderr or e.stdout or "").strip().splitlines()[-1:] or [""]
            raise _SubprocessFailure(
                f"datasets exited {e.returncode}: {msg[0]}") from e
        except FileNotFoundError as e:
            raise _SubprocessFailure(
                f"'{self.binary}' binary not found on PATH "
                "(install ncbi-datasets-cli via conda).") from e

    def _persist_failure(self, accession: str, status: str, msg: str) -> None:
        self.cache.put_summary(SummaryEntry(
            accession=accession,
            fetched_at=datetime.now(timezone.utc),
            status=status,
            summary_json=None,
            lineage_json=None,
            error_message=msg,
        ))


class _SubprocessFailure(RuntimeError):
    """Internal — raised by ``DatasetsCLIClient._run`` on any process failure."""


def _classify_downloaded_file(path: Path) -> Optional[str]:
    """Map a filename to one of the VALID_INCLUDES keys."""
    name = path.name.lower()
    if name.endswith(".gff") or name.endswith(".gff.gz") or "genomic.gff" in name:
        return "gff3"
    if "_protein.faa" in name or name.endswith(".faa"):
        return "protein"
    if "_cds_from_genomic" in name or "cds.fna" in name:
        return "cds"
    if "_rna_from_genomic" in name or "rna.fna" in name:
        return "rna"
    if name.endswith(".fna") or name.endswith(".fa") or name.endswith(".fasta"):
        return "genome"
    if "sequence_report" in name:
        return "seq-report"
    return None


# ===========================================================================
# Mock client — for offline testing and pipeline development
# ===========================================================================

class MockNCBIClient:
    """In-memory NCBI client that uses fixtures from ``assembly.mock_*``.

    Same interface as ``DatasetsCLIClient``. Useful for:
      * CI tests with no network.
      * Local pipeline development on a laptop without ``datasets`` installed.
      * Reproducing a past run from cached fixtures.
    """

    def __init__(
        self,
        cache: AssemblyCache,
        download_dir: Path,
        *,
        mock_summaries: dict[str, dict],
        mock_lineages: Optional[dict[str, dict]] = None,
        mock_gene_counts: Optional[dict[str, int]] = None,
        mock_gff_gene_counts: Optional[dict[str, int]] = None,
        unknown_accession_behaviour: str = "not_found",  # or "suppressed"
    ) -> None:
        self.cache = cache
        self.download_dir = Path(download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)
        self.summaries = mock_summaries
        self.lineages = mock_lineages or {}
        self.gene_counts = mock_gene_counts or {}
        # Gene counts to bake into the synthetic GFF3 produced by
        # download_files. Used by --use-mock to make the demo output show
        # realistic Staphylococcus aureus / Bacillus / Mycobacterium counts.
        self._mock_gff_gene_count = mock_gff_gene_counts or {}
        self.unknown_accession_behaviour = unknown_accession_behaviour

    def fetch_summary(self, accession: str, *, refresh: bool = False) -> AssemblyFetch:
        if not refresh:
            cached = self.cache.get_summary(accession)
            if cached:
                return _fetch_from_cache_entry(cached)

        if accession not in self.summaries:
            if self.unknown_accession_behaviour == "suppressed":
                self._persist_failure(accession, "suppressed",
                                      "Mock: accession marked suppressed.")
                return AssemblyFetch(accession, "suppressed", None, None,
                                     "Mock: accession marked suppressed.", False)
            msg = "Mock: accession not present in fixture set."
            self._persist_failure(accession, "not_found", msg)
            return AssemblyFetch(accession, "not_found", None, None, msg, False)

        summary_obj = self.summaries[accession]
        lineage_obj = self.lineages.get(accession)
        record, gene_count, err = _parse_summary_json(
            accession, summary_obj, lineage_obj,
        )
        # Allow tests to inject gene_count override (handy when summary lacks it)
        if gene_count is None and accession in self.gene_counts:
            gene_count = self.gene_counts[accession]

        if record is None:
            # Same categorisation logic as DatasetsCLIClient — keeps the mock
            # honest so report labels match what production would produce.
            status = "suppressed" if err and err.startswith("Suppressed") else "parse_error"
            self._persist_failure(accession, status, err or "unknown")
            return AssemblyFetch(accession, status, None, None, err, False)

        self.cache.put_summary(SummaryEntry(
            accession=accession,
            fetched_at=datetime.now(timezone.utc),
            status="ok",
            summary_json=json.dumps(summary_obj),
            lineage_json=json.dumps(lineage_obj) if lineage_obj else None,
            error_message=None,
        ))
        return AssemblyFetch(accession, "ok", record, gene_count, None, False)

    def download_files(
        self, accession: str, includes: Iterable[str], out_dir: Path,
    ) -> DownloadResult:
        """Simulate a selective download by writing tiny but parseable files.

        Each file type gets realistic-shaped content so downstream parsers
        (GFF, FASTA) work without further patching:

            gff3    → valid GFF3 with a small synthetic gene set
            genome  → 1-record FASTA placeholder
            protein → 1-record FASTA placeholder
            cds/rna → similar small FASTA
        """
        target = Path(out_dir) / accession
        target.mkdir(parents=True, exist_ok=True)
        files: dict[str, Path] = {}
        total_bytes = 0
        for ft in includes:
            if ft not in VALID_INCLUDES:
                raise ValueError(f"Unsupported include {ft!r}")
            path = target / f"{accession}_{ft}"
            if ft == "gff3":
                # Write a small valid GFF3 — n_genes drawn from the mock
                # gene_counts table if available, else default to 5000.
                n = self._mock_gff_gene_count.get(accession, 5000)
                lines = ["##gff-version 3"]
                for i in range(n):
                    lines.append(
                        f"chr1\tRefSeq\tgene\t{i*100+1}\t{i*100+99}\t.\t+\t."
                        f"\tID=gene{i};Name=g{i}"
                    )
                payload = ("\n".join(lines) + "\n").encode()
                path = path.with_suffix(".gff")
            elif ft in ("genome", "cds", "rna"):
                payload = f">chr1 mock {ft} {accession}\nACGTACGTACGTACGT\n".encode()
                path = path.with_suffix(".fna")
            elif ft == "protein":
                payload = f">protein_1 mock {accession}\nMACGTACGTACGT\n".encode()
                path = path.with_suffix(".faa")
            else:
                payload = f"# Mock {ft} for {accession}\n".encode()
                path = path.with_suffix(".txt")
            path.write_bytes(payload)
            sz = len(payload)
            files[ft] = path
            total_bytes += sz
            self.cache.record_file(accession, ft, path, sz)
        return DownloadResult(accession=accession, files=files, bytes_=total_bytes)

    def cleanup_files(self, accession: str) -> int:
        rows = self.cache.list_files(accession, include_deleted=False)
        freed = 0
        acc_dir = self.download_dir / accession
        if acc_dir.exists():
            for p in acc_dir.rglob("*"):
                if p.is_file():
                    freed += p.stat().st_size
            shutil.rmtree(acc_dir)
        # Mark files deleted whether or not the directory existed (idempotent).
        for _ in rows:
            pass
        self.cache.mark_files_deleted(accession)
        return freed

    def _persist_failure(self, accession: str, status: str, msg: str) -> None:
        self.cache.put_summary(SummaryEntry(
            accession=accession,
            fetched_at=datetime.now(timezone.utc),
            status=status,
            summary_json=None,
            lineage_json=None,
            error_message=msg,
        ))
