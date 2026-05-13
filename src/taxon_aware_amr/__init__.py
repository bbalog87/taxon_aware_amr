"""taxon_aware_amr — taxon-aware AMR / virulence detection pipeline.

Steps:
    1. Taxon classifier + reporting.
    2. NCBI retrieval + SQLite cache + selective downloads + cleanup.
    3. Gene counting (NCBI summary preferred; GFF3 fallback) + auto-replace
       of suppressed assemblies.
    4. AMR detection — AMRFinderPlus in combined mode (-n -p -g) on bacterial
       assemblies, `--plus` on, supported-organism mapping for point mutations.
    5. Virulence detection — ABRicate-VFDB + reuse of AMRFinderPlus `--plus`
       virulence rows from cache, merged + deduplicated.
    6. Final output dataframe matching the original task spec.
    7. Executive summary + run report.
"""

from .assembly     import AssemblyRecord, mock_assembly_record, normalise_genome_id
from .taxonomy     import (
    TaxonGroup, AnalysisStatus, TaxonDecision,
    CATEGORY_TO_EXPECTED_GROUP, lineage_to_group, decide,
)
from .reporting    import DecisionReporter
from .cache        import AssemblyCache, SummaryEntry, FileInventoryRow
from .ncbi         import (
    NCBIClient, DatasetsCLIClient, MockNCBIClient,
    AssemblyFetch, DownloadResult, VALID_INCLUDES,
)
from .input        import InputRow, read_input
from .gene_count   import count_genes, count_genes_from_gff, DEFAULT_GENE_TYPES
from .amr          import (
    AMRClient, AMRFinderPlusClient, MockAMRClient,
    AMRHit, AMRResult,
    AMRFINDER_ORGANISM_MAP, amrfinder_organism_for,
    parse_amrfinder_tsv,
)
from .virulence    import (
    ABRicateClient, ABRicateCLIClient, MockABRicateClient,
    VirulenceHit, VirulenceResult,
    detect_virulence, parse_abricate_tsv,
    virulence_hits_from_amrfinder_cache,
)
from .final_output import (
    FinalRow, build_final_row, write_minimal_tsv, write_full_tsv,
    MINIMAL_COLUMNS, FULL_COLUMNS,
)
from .orchestrator import (
    RowPlan, RowOutcome,
    plan_downloads, run_row, run_pipeline,
    extract_replacement_accession,
)

__all__ = [
    # assembly
    "AssemblyRecord", "mock_assembly_record", "normalise_genome_id",
    # taxonomy
    "TaxonGroup", "AnalysisStatus", "TaxonDecision",
    "CATEGORY_TO_EXPECTED_GROUP", "lineage_to_group", "decide",
    # reporting
    "DecisionReporter",
    # cache
    "AssemblyCache", "SummaryEntry", "FileInventoryRow",
    # ncbi
    "NCBIClient", "DatasetsCLIClient", "MockNCBIClient",
    "AssemblyFetch", "DownloadResult", "VALID_INCLUDES",
    # input
    "InputRow", "read_input",
    # gene_count
    "count_genes", "count_genes_from_gff", "DEFAULT_GENE_TYPES",
    # amr
    "AMRClient", "AMRFinderPlusClient", "MockAMRClient",
    "AMRHit", "AMRResult",
    "AMRFINDER_ORGANISM_MAP", "amrfinder_organism_for",
    "parse_amrfinder_tsv",
    # virulence
    "ABRicateClient", "ABRicateCLIClient", "MockABRicateClient",
    "VirulenceHit", "VirulenceResult",
    "detect_virulence", "parse_abricate_tsv",
    "virulence_hits_from_amrfinder_cache",
    # final_output
    "FinalRow", "build_final_row", "write_minimal_tsv", "write_full_tsv",
    "MINIMAL_COLUMNS", "FULL_COLUMNS",
    # orchestrator
    "RowPlan", "RowOutcome",
    "plan_downloads", "run_row", "run_pipeline",
    "extract_replacement_accession",
]

__version__ = "1.0.1"

