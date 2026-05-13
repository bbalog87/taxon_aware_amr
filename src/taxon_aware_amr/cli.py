"""
taxon_aware_amr.cli
===================

Command-line entry point.

Usage::

    taxon_aware_amr INPUT [-o OUT] [-c CACHE]
                          [--auto-replace-suppressed]
                          [--keep-files]
                          [--use-mock]
                          [--datasets-binary PATH]
                          [--refresh-cache]
                          [-v | -vv]

See ``--help`` for the full list. Run::

    taxon_aware_amr examples/kenya_biodiversity_input.tsv --use-mock

for a complete end-to-end demo without any external tools.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Sequence

from .cache        import AssemblyCache
from .ncbi         import DatasetsCLIClient, MockNCBIClient
from .orchestrator import run_pipeline

DEFAULT_CACHE_DIR = Path.home() / ".cache" / "taxon_aware_amr"


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="taxon_aware_amr",
        description=(
            "Taxon-aware AMR / virulence detection pipeline. "
            "Routes each row to the right tool based on the NCBI assembly's "
            "lineage; bacterial AMR / virulence databases are run ONLY on "
            "bacterial assemblies."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Outputs in OUT/:\n"
            "  decisions.tsv   one row per input, machine-readable\n"
            "  run_report.md   human-readable summary\n\n"
            "Cache in CACHE/:\n"
            "  cache.sqlite    assembly summaries + file inventory + run log\n"
            "  genomes/        downloaded FASTA / GFF (deleted after use\n"
            "                  unless --keep-files)\n"
        ),
    )
    p.add_argument("input", type=Path,
                   help="Input TSV with at least taxon_id, genome_id, "
                        "species_name, category columns.")
    p.add_argument("-o", "--out", type=Path, default=Path("results"),
                   help="Output directory (default: ./results)")
    p.add_argument("-c", "--cache", type=Path, default=DEFAULT_CACHE_DIR,
                   help=f"Cache directory (default: {DEFAULT_CACHE_DIR})")

    p.add_argument("--auto-replace-suppressed", action="store_true",
                   help="When NCBI returns 'Suppressed → Replaced by GCA_xxx', "
                        "automatically fetch the replacement (logged as audit "
                        "warning). Default: off; suppressed rows are skipped.")
    p.add_argument("--keep-files", action="store_true",
                   help="Do not delete downloaded files after each row finishes. "
                        "Useful for debugging; expensive on a large input table.")
    p.add_argument("--use-mock", action="store_true",
                   help="Use the bundled mock NCBI client and demo fixtures "
                        "(works only with examples/kenya_biodiversity_input.tsv). "
                        "Lets you try the pipeline without installing datasets-cli.")
    p.add_argument("--datasets-binary", default="datasets",
                   metavar="PATH",
                   help="Path to the `datasets` CLI binary (default: 'datasets' "
                        "on PATH; ignored with --use-mock).")
    p.add_argument("--amrfinder-binary", default="amrfinder",
                   metavar="PATH",
                   help="Path to the `amrfinder` (AMRFinderPlus) CLI binary "
                        "(default: 'amrfinder' on PATH; ignored with --use-mock).")
    p.add_argument("--abricate-binary", default="abricate",
                   metavar="PATH",
                   help="Path to the `abricate` CLI binary "
                        "(default: 'abricate' on PATH; ignored with --use-mock).")
    p.add_argument("--threads", type=int, default=4, metavar="N",
                   help="Threads to pass to AMRFinderPlus (default: 4).")
    p.add_argument("--refresh-cache", action="store_true",
                   help="Ignore cached summaries and re-query NCBI. "
                        "(Not implemented in mock mode.)")
    p.add_argument("-v", "--verbose", action="count", default=0,
                   help="Increase log verbosity. -v: INFO, -vv: DEBUG.")
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)

    # ---- Logging ----
    level = [logging.WARNING, logging.INFO, logging.DEBUG][min(args.verbose, 2)]
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-7s %(name)s :: %(message)s",
        datefmt="%H:%M:%S",
    )

    # ---- Sanity checks ----
    if not args.input.exists():
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        return 2

    cache_path   = args.cache / "cache.sqlite"
    download_dir = args.cache / "genomes"

    print(f"Input:          {args.input}")
    print(f"Output:         {args.out}")
    print(f"Cache:          {args.cache}")
    print(f"NCBI backend:   {'mock fixtures' if args.use_mock else f'datasets CLI ({args.datasets_binary})'}")
    print(f"Auto-replace:   {'on' if args.auto_replace_suppressed else 'off'}")
    print(f"Keep files:     {'yes' if args.keep_files else 'no (cleanup after each row)'}")
    print()

    cache = AssemblyCache(cache_path)

    if args.use_mock:
        from .fixtures import (
            DEMO_SUMMARIES, DEMO_LINEAGES, DEMO_GFF_GENE_COUNTS,
            demo_amr_hits_as_records, demo_vfdb_hits_as_records,
        )
        from .amr import MockAMRClient
        from .virulence import MockABRicateClient
        client = MockNCBIClient(
            cache=cache, download_dir=download_dir,
            mock_summaries=DEMO_SUMMARIES,
            mock_lineages=DEMO_LINEAGES,
            mock_gff_gene_counts=DEMO_GFF_GENE_COUNTS,
        )
        amr_client = MockAMRClient(
            cache=cache, mock_hits=demo_amr_hits_as_records(),
        )
        abricate_client = MockABRicateClient(
            cache=cache, mock_hits=demo_vfdb_hits_as_records(),
        )
    else:
        from .amr import AMRFinderPlusClient
        from .virulence import ABRicateCLIClient
        client = DatasetsCLIClient(
            cache=cache, download_dir=download_dir,
            binary=args.datasets_binary,
        )
        amr_client = AMRFinderPlusClient(
            cache=cache, binary=args.amrfinder_binary,
            threads=args.threads,
        )
        abricate_client = ABRicateCLIClient(
            cache=cache, binary=args.abricate_binary,
            threads=args.threads,
        )

    try:
        summary = run_pipeline(
            input_path=args.input,
            client=client, cache=cache,
            out_dir=args.out, download_dir=download_dir,
            keep_files=args.keep_files,
            auto_replace_suppressed=args.auto_replace_suppressed,
            amr_client=amr_client,
            abricate_client=abricate_client,
        )
    finally:
        cache.close()

    print()
    print(f"Done. {summary['n_rows']} rows processed.")
    print(f"  Per-row TSV:    {args.out / 'decisions.tsv'}")
    print(f"  Markdown report: {args.out / 'run_report.md'}")
    print(f"  SQLite cache:    {cache_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
