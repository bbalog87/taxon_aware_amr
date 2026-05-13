"""
taxon_aware_amr.cache
=====================

SQLite-backed cache for NCBI assembly metadata and downloaded-file bookkeeping.

Why a cache
-----------
GCA accessions are versioned and immutable (``GCA_964341285.1`` will never
return a different summary than it did yesterday), so we can cache summaries
indefinitely. NCBI's Datasets API also rate-limits aggressive callers; for an
80-row input table that we'll iterate on while developing the pipeline, an
on-disk cache is essential to avoid hammering the service.

Failed lookups are also cached, but with a shorter TTL (default 24 h), so a
transient outage doesn't poison the cache for typos and so genuine typos
don't keep firing live calls.

Schema
------
``assembly_summary``  one row per accession we have ever queried
    accession          TEXT PRIMARY KEY
    fetched_at         TEXT NOT NULL        -- ISO timestamp UTC
    status             TEXT NOT NULL        -- ok | not_found | suppressed | error
    summary_json       TEXT                 -- raw `datasets summary genome ...` JSON
    lineage_json       TEXT                 -- raw `datasets summary taxonomy ...` JSON
    error_message      TEXT                 -- human-readable failure reason

``file_inventory``  one row per downloaded artefact
    accession          TEXT NOT NULL
    file_type          TEXT NOT NULL        -- genome | gff3 | protein | cds | rna
    path               TEXT NOT NULL
    size_bytes         INTEGER NOT NULL
    downloaded_at      TEXT NOT NULL
    deleted_at         TEXT                 -- NULL while file is present on disk
    PRIMARY KEY (accession, file_type)

``run_log``  one row per pipeline invocation
    id                 INTEGER PRIMARY KEY AUTOINCREMENT
    started_at         TEXT NOT NULL
    finished_at        TEXT
    n_rows             INTEGER
    notes              TEXT
"""

from __future__ import annotations

import json
import logging
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


# ===========================================================================
# DTOs
# ===========================================================================

@dataclass
class SummaryEntry:
    """One cached assembly summary."""
    accession:     str
    fetched_at:    datetime
    status:        str                 # ok | not_found | suppressed | error
    summary_json:  Optional[str]
    lineage_json:  Optional[str]
    error_message: Optional[str]

    def is_stale(self, *, negative_ttl: timedelta) -> bool:
        """Failures expire after negative_ttl; successes never expire."""
        if self.status == "ok":
            return False
        age = datetime.now(timezone.utc) - self.fetched_at
        return age > negative_ttl


@dataclass
class FileInventoryRow:
    accession:     str
    file_type:     str
    path:          Path
    size_bytes:    int
    downloaded_at: datetime
    deleted_at:    Optional[datetime]


# ===========================================================================
# Cache
# ===========================================================================

class AssemblyCache:
    """SQLite-backed cache. Thread-safe enough for our single-process pipeline.

    Parameters
    ----------
    db_path
        Path to the SQLite file. Created on first use.
    negative_ttl_hours
        How long to remember failed lookups before re-trying NCBI. Default 24 h.
    """

    SCHEMA_VERSION = 1

    def __init__(self, db_path: Path, *, negative_ttl_hours: int = 24) -> None:
        self.db_path = Path(db_path)
        self.negative_ttl = timedelta(hours=negative_ttl_hours)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self.db_path), isolation_level=None)
        self._conn.row_factory = sqlite3.Row
        self._init_schema()

    # ------------------------------------------------------------------ #
    # Schema                                                              #
    # ------------------------------------------------------------------ #
    def _init_schema(self) -> None:
        with self._conn:
            self._conn.executescript("""
                CREATE TABLE IF NOT EXISTS schema_meta (
                    key   TEXT PRIMARY KEY,
                    value TEXT
                );

                CREATE TABLE IF NOT EXISTS assembly_summary (
                    accession      TEXT PRIMARY KEY,
                    fetched_at     TEXT NOT NULL,
                    status         TEXT NOT NULL,
                    summary_json   TEXT,
                    lineage_json   TEXT,
                    error_message  TEXT
                );

                CREATE TABLE IF NOT EXISTS file_inventory (
                    accession      TEXT NOT NULL,
                    file_type      TEXT NOT NULL,
                    path           TEXT NOT NULL,
                    size_bytes     INTEGER NOT NULL,
                    downloaded_at  TEXT NOT NULL,
                    deleted_at     TEXT,
                    PRIMARY KEY (accession, file_type)
                );

                -- Results from per-row tools (AMRFinderPlus, ABRicate VFDB).
                -- One row per (accession, tool). Lets re-runs skip re-execution
                -- and lets Step 5 reuse the AMRFinderPlus --plus output that
                -- Step 4 already produced.
                CREATE TABLE IF NOT EXISTS tool_results (
                    accession      TEXT NOT NULL,
                    tool           TEXT NOT NULL,
                    completed_at   TEXT NOT NULL,
                    status         TEXT NOT NULL,   -- ok | failed
                    db_version     TEXT,
                    organism_flag  TEXT,
                    raw_output     TEXT,
                    parsed_json    TEXT,
                    error_message  TEXT,
                    PRIMARY KEY (accession, tool)
                );

                CREATE TABLE IF NOT EXISTS run_log (
                    id           INTEGER PRIMARY KEY AUTOINCREMENT,
                    started_at   TEXT NOT NULL,
                    finished_at  TEXT,
                    n_rows       INTEGER,
                    notes        TEXT
                );
            """)
            self._conn.execute(
                "INSERT OR REPLACE INTO schema_meta(key, value) VALUES (?, ?)",
                ("version", str(self.SCHEMA_VERSION)),
            )

    # ------------------------------------------------------------------ #
    # Assembly summaries                                                  #
    # ------------------------------------------------------------------ #
    def get_summary(self, accession: str) -> Optional[SummaryEntry]:
        row = self._conn.execute(
            "SELECT * FROM assembly_summary WHERE accession = ?",
            (accession,),
        ).fetchone()
        if row is None:
            return None
        entry = SummaryEntry(
            accession=row["accession"],
            fetched_at=datetime.fromisoformat(row["fetched_at"]),
            status=row["status"],
            summary_json=row["summary_json"],
            lineage_json=row["lineage_json"],
            error_message=row["error_message"],
        )
        if entry.is_stale(negative_ttl=self.negative_ttl):
            logger.debug("Cache entry for %s is stale; will refetch.", accession)
            return None
        return entry

    def put_summary(self, entry: SummaryEntry) -> None:
        self._conn.execute(
            """
            INSERT OR REPLACE INTO assembly_summary
                (accession, fetched_at, status, summary_json, lineage_json, error_message)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                entry.accession,
                entry.fetched_at.isoformat(timespec="seconds"),
                entry.status,
                entry.summary_json,
                entry.lineage_json,
                entry.error_message,
            ),
        )

    # ------------------------------------------------------------------ #
    # File inventory                                                      #
    # ------------------------------------------------------------------ #
    def record_file(
        self, accession: str, file_type: str, path: Path, size_bytes: int,
    ) -> None:
        self._conn.execute(
            """
            INSERT OR REPLACE INTO file_inventory
                (accession, file_type, path, size_bytes, downloaded_at, deleted_at)
            VALUES (?, ?, ?, ?, ?, NULL)
            """,
            (accession, file_type, str(path), int(size_bytes),
             datetime.now(timezone.utc).isoformat(timespec="seconds")),
        )

    def mark_files_deleted(self, accession: str) -> list[FileInventoryRow]:
        """Mark all files for an accession as deleted; return what was deleted."""
        rows = self.list_files(accession, include_deleted=False)
        self._conn.execute(
            """
            UPDATE file_inventory
            SET    deleted_at = ?
            WHERE  accession = ? AND deleted_at IS NULL
            """,
            (datetime.now(timezone.utc).isoformat(timespec="seconds"), accession),
        )
        return rows

    def list_files(
        self, accession: str, *, include_deleted: bool = True,
    ) -> list[FileInventoryRow]:
        sql = "SELECT * FROM file_inventory WHERE accession = ?"
        if not include_deleted:
            sql += " AND deleted_at IS NULL"
        sql += " ORDER BY file_type"
        out = []
        for row in self._conn.execute(sql, (accession,)).fetchall():
            out.append(FileInventoryRow(
                accession=row["accession"],
                file_type=row["file_type"],
                path=Path(row["path"]),
                size_bytes=row["size_bytes"],
                downloaded_at=datetime.fromisoformat(row["downloaded_at"]),
                deleted_at=(datetime.fromisoformat(row["deleted_at"])
                            if row["deleted_at"] else None),
            ))
        return out

    # ------------------------------------------------------------------ #
    # Run log                                                             #
    # ------------------------------------------------------------------ #
    def start_run(self, notes: str = "") -> int:
        cur = self._conn.execute(
            "INSERT INTO run_log(started_at, notes) VALUES (?, ?)",
            (datetime.now(timezone.utc).isoformat(timespec="seconds"), notes),
        )
        return int(cur.lastrowid)

    def finish_run(self, run_id: int, n_rows: int, notes: str = "") -> None:
        self._conn.execute(
            """
            UPDATE run_log
            SET    finished_at = ?,
                   n_rows      = ?,
                   notes       = COALESCE(notes, '') || ?
            WHERE  id = ?
            """,
            (datetime.now(timezone.utc).isoformat(timespec="seconds"),
             n_rows, ("\n" + notes if notes else ""), run_id),
        )

    # ------------------------------------------------------------------ #
    # Stats / introspection                                               #
    # ------------------------------------------------------------------ #
    def stats(self) -> dict:
        c = self._conn.cursor()
        return {
            "summaries_total": c.execute(
                "SELECT COUNT(*) FROM assembly_summary").fetchone()[0],
            "summaries_ok": c.execute(
                "SELECT COUNT(*) FROM assembly_summary WHERE status = 'ok'").fetchone()[0],
            "summaries_failed": c.execute(
                "SELECT COUNT(*) FROM assembly_summary WHERE status != 'ok'").fetchone()[0],
            "files_present": c.execute(
                "SELECT COUNT(*) FROM file_inventory WHERE deleted_at IS NULL").fetchone()[0],
            "files_deleted": c.execute(
                "SELECT COUNT(*) FROM file_inventory WHERE deleted_at IS NOT NULL").fetchone()[0],
            "bytes_on_disk": c.execute(
                "SELECT COALESCE(SUM(size_bytes),0) FROM file_inventory "
                "WHERE deleted_at IS NULL").fetchone()[0],
            "bytes_ever_downloaded": c.execute(
                "SELECT COALESCE(SUM(size_bytes),0) FROM file_inventory").fetchone()[0],
        }

    def close(self) -> None:
        self._conn.close()
