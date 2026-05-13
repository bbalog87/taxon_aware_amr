"""
taxon_aware_amr.fixtures
========================

In-package mock NCBI fixtures for the bundled example input
(``examples/kenya_biodiversity_input.tsv``).

These let a user try the pipeline end-to-end with ``--use-mock`` before
installing the real ``ncbi-datasets-cli`` and ``amrfinderplus`` tooling::

    taxon_aware_amr examples/kenya_biodiversity_input.tsv --use-mock

The JSON shapes match what ``datasets summary genome accession`` returns
in real life, so the parsing code in ``ncbi.py`` is exercised against
production-shaped data.
"""

from __future__ import annotations


def _summary(tax_id: int, name: str, level: str,
             gene_total: int | None,
             suppressed: str | None = None) -> dict:
    """Build a `datasets summary genome accession ...` shaped record."""
    ai: dict = {}
    if gene_total is not None:
        ai = {"stats": {"gene_counts": {
            "total": gene_total,
            "protein_coding": max(0, gene_total - 50),
            "non_coding": 50,
            "pseudogene": 0,
        }}}
    out = {
        "accession": "",
        "organism": {"tax_id": tax_id, "organism_name": name},
        "annotation_info": ai,
        "assembly_info": {"assembly_level": level},
    }
    if suppressed:
        out["assembly_info"]["suppression_reason"] = suppressed
    return out


def _lineage(*ranks: tuple[str, str]) -> dict:
    """Build a `datasets summary taxonomy taxon ...` shaped record."""
    return {"taxonomy": {"classification": {
        rank: {"name": name} for rank, name in ranks
    }}}


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

DEMO_SUMMARIES: dict[str, dict] = {
    # ---- Bacteria ----
    "GCA_964341285.1": _summary(1392, "Bacillus anthracis",         "Complete Genome", 5544),
    "GCA_900654255.1": _summary(1773, "Mycobacterium tuberculosis", "Scaffold",        4173),
    # Bacterial WITHOUT NCBI annotation → forces gff3 download for gene count:
    "GCA_001049615.1": _summary(1280, "Staphylococcus aureus",      "Contig",          None),

    # ---- Viruses ----
    "GCA_000847345.1": _summary(11588, "Rift Valley fever virus", "Complete Genome", 7),
    "GCA_900327715.1": _summary(11292, "Rabies lyssavirus",       "Complete Genome", 5),

    # ---- Protozoa ----
    "GCA_019455585.1": _summary(5811, "Toxoplasma gondii",      "Chromosome",      8920),
    "GCA_035232765.1": _summary(5807, "Cryptosporidium parvum", "Complete Genome", 3865),

    # ---- Fungus ----
    "GCA_014900015.1": _summary(5334, "Schizophyllum commune", "Contig", 13210),

    # ---- Plant ----
    "GCA_964199775.1": _summary(4577, "Zea mays", "Complete Genome", 39591),

    # ---- Metazoa ----
    "GCA_030014295.1": _summary(9785, "Loxodonta africana", "Chromosome", 22871),
    "GCA_025407655.1": _summary(7159, "Aedes aegypti",      "Chromosome", 19500),

    # ---- Suppressed → has a replacement that resolves successfully ----
    "GCA_015501595.1": _summary(1081385, "Diceros bicornis", "Scaffold", None,
                                 suppressed="Replaced by GCA_022413745.1"),
    "GCA_022413745.1": _summary(1081385, "Diceros bicornis", "Chromosome", 21845),

    # GCA_999999999.9 is deliberately absent → not_found
}


# ---------------------------------------------------------------------------
# Lineages
# ---------------------------------------------------------------------------

DEMO_LINEAGES: dict[str, dict] = {
    "GCA_964341285.1": _lineage(("superkingdom","Bacteria"), ("phylum","Bacillota"),
                                ("class","Bacilli"), ("genus","Bacillus"),
                                ("species","Bacillus anthracis")),
    "GCA_900654255.1": _lineage(("superkingdom","Bacteria"), ("phylum","Actinomycetota"),
                                ("species","Mycobacterium tuberculosis")),
    "GCA_001049615.1": _lineage(("superkingdom","Bacteria"), ("phylum","Bacillota"),
                                ("species","Staphylococcus aureus")),
    "GCA_000847345.1": _lineage(("superkingdom","Viruses"),
                                ("phylum","Negarnaviricota"),
                                ("species","Rift Valley fever virus")),
    "GCA_900327715.1": _lineage(("superkingdom","Viruses"),
                                ("phylum","Negarnaviricota"),
                                ("species","Rabies lyssavirus")),
    "GCA_019455585.1": _lineage(("superkingdom","Eukaryota"),
                                ("phylum","Apicomplexa"),
                                ("species","Toxoplasma gondii")),
    "GCA_035232765.1": _lineage(("superkingdom","Eukaryota"),
                                ("phylum","Apicomplexa"),
                                ("species","Cryptosporidium parvum")),
    "GCA_014900015.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Fungi"),
                                ("species","Schizophyllum commune")),
    "GCA_964199775.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Viridiplantae"),
                                ("species","Zea mays")),
    "GCA_030014295.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Metazoa"),
                                ("species","Loxodonta africana")),
    "GCA_025407655.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Metazoa"),
                                ("species","Aedes aegypti")),
    "GCA_022413745.1": _lineage(("superkingdom","Eukaryota"), ("kingdom","Metazoa"),
                                ("species","Diceros bicornis")),
}


# ---------------------------------------------------------------------------
# Mock AMR hits — Step 4 fixtures for the bacterial demo accessions
# ---------------------------------------------------------------------------
# These mirror what AMRFinderPlus typically reports for the wild-type or
# common reference isolate of each species. Gene symbols and class strings
# match the NCBI Reference Gene Catalog conventions.

def _hit(gene: str, etype: str, drug_class: str, *, subclass: str = "",
         subtype: str = "", method: str = "BLASTX", pid: float = 99.5,
         cov: float = 100.0) -> dict:
    return {
        "gene_symbol": gene, "element_type": etype,
        "element_subtype": subtype or etype,
        "drug_class": drug_class, "subclass": subclass or drug_class,
        "method": method, "pct_identity": pid, "pct_coverage": cov,
        "contig": "chr1", "start": 1, "stop": 1000, "reference_acc": "WP_000000000.1",
    }


DEMO_AMR_HITS: dict[str, list[dict]] = {
    # Bacillus anthracis — typical wild-type carries intrinsic blaA + a
    # macrolide efflux gene; not really an "MDR pathogen" but Step 4 should
    # still surface what's there.
    "GCA_964341285.1": [
        _hit("blaA",     "AMR", "BETA-LACTAM",  subclass="CEPHALOSPORIN"),
        _hit("msrA",     "AMR", "MACROLIDE",    subclass="MACROLIDE"),
    ],
    # MDR M. tuberculosis — INH/RIF/FQ/aminoglycoside resistance with
    # point mutations (--organism flag enables POINTX detection).
    "GCA_900654255.1": [
        _hit("katG_S315T", "AMR", "ISONIAZID",   method="POINTX",
             subclass="ISONIAZID", subtype="POINT"),
        _hit("inhA",       "AMR", "ISONIAZID",   subclass="ISONIAZID"),
        _hit("rpoB_S531L", "AMR", "RIFAMYCIN",   method="POINTX",
             subclass="RIFAMPIN", subtype="POINT"),
        _hit("gyrA_D94G",  "AMR", "QUINOLONE",   method="POINTX",
             subclass="FLUOROQUINOLONE", subtype="POINT"),
        _hit("rrs_A1401G", "AMR", "AMINOGLYCOSIDE", method="POINTX",
             subclass="STREPTOMYCIN", subtype="POINT"),
    ],
    # MRSA — mecA (the classic resistance determinant) + the rest of the
    # SCCmec cassette neighbours, plus a couple of acquired genes.
    "GCA_001049615.1": [
        _hit("mecA",            "AMR", "BETA-LACTAM", subclass="METHICILLIN"),
        _hit("mecR1",           "AMR", "BETA-LACTAM", subclass="METHICILLIN"),
        _hit("blaZ",            "AMR", "BETA-LACTAM", subclass="PENICILLIN"),
        _hit("ermC",            "AMR", "MACROLIDE",   subclass="MLS"),
        _hit("tet(K)",          "AMR", "TETRACYCLINE",subclass="TETRACYCLINE"),
        _hit("aac(6')-aph(2'')","AMR", "AMINOGLYCOSIDE",
             subclass="GENTAMICIN"),
        # One stress and one virulence row → confirms --plus is on; these
        # are filtered out in Step 4 but visible to Step 5.
        _hit("qacJ",            "STRESS",    "QUATERNARY AMMONIUM"),
        _hit("hlb",              "VIRULENCE", "VIRULENCE"),
    ],
}


def demo_amr_hits_as_records() -> dict:
    """Convert demo hit dicts to AMRHit instances for MockAMRClient."""
    from .amr import AMRHit
    return {acc: [AMRHit(**d) for d in hits]
            for acc, hits in DEMO_AMR_HITS.items()}


# ---------------------------------------------------------------------------
# Mock VFDB virulence hits — Step 5 fixtures
# ---------------------------------------------------------------------------
# These mirror what ABRicate-VFDB typically reports for the reference isolate
# of each species. Realistic virulence factor names from VFDB.

def _vhit(gene: str, product: str, *, pid: float = 99.0, cov: float = 100.0) -> dict:
    return {
        "gene_symbol": gene, "source": "vfdb_abricate",
        "product": product, "db": "vfdb",
        "pct_identity": pid, "pct_coverage": cov,
        "contig": "chr1", "start": 1, "stop": 1000,
    }


DEMO_VFDB_HITS: dict[str, list[dict]] = {
    # Anthrax toxin: pagA (protective antigen), lef (lethal factor),
    # cya (edema factor) + capsule biosynthesis cap genes.
    "GCA_964341285.1": [
        _vhit("pagA", "Anthrax toxin moiety: protective antigen"),
        _vhit("lef",  "Anthrax toxin moiety: lethal factor"),
        _vhit("cya",  "Anthrax toxin moiety: edema factor"),
        _vhit("capB", "Poly-gamma-D-glutamate capsule biosynthesis"),
        _vhit("capC", "Poly-gamma-D-glutamate capsule biosynthesis"),
    ],
    # M. tuberculosis: ESX-1 type VII secretion system + mycolic acid synth.
    "GCA_900654255.1": [
        _vhit("esxA", "ESAT-6 (early secreted antigen target 6 kDa)"),
        _vhit("esxB", "CFP-10 (culture filtrate protein 10 kDa)"),
        _vhit("eccCb1","ESX-1 type VII secretion ATPase"),
        _vhit("mmpL7","Mycobacterial membrane protein large 7"),
    ],
    # MRSA: adhesins (clf), protein A (spa), haemolysins (hla, hlb),
    # staphylokinase (sak), enterotoxin A.
    "GCA_001049615.1": [
        _vhit("clfA", "Clumping factor A (fibrinogen-binding adhesin)"),
        _vhit("clfB", "Clumping factor B"),
        _vhit("spa",  "Immunoglobulin G-binding protein A"),
        _vhit("hla",  "Alpha-haemolysin"),
        _vhit("hlb",  "Beta-haemolysin"),
        _vhit("sak",  "Staphylokinase"),
        _vhit("sea",  "Staphylococcal enterotoxin A"),
        _vhit("icaA", "Intercellular adhesin / biofilm formation"),
    ],
}


def demo_vfdb_hits_as_records() -> dict:
    """Convert demo dicts to VirulenceHit instances for MockABRicateClient."""
    from .virulence import VirulenceHit
    return {acc: [VirulenceHit(**d) for d in hits]
            for acc, hits in DEMO_VFDB_HITS.items()}


# ---------------------------------------------------------------------------
# Mock GFF gene counts — used by MockNCBIClient to bake realistic gene
# counts into the synthetic GFF3 files it writes during --use-mock demos.
# Only matters for bacterial accessions whose NCBI summary lacks gene_counts
# (so Step 3 has to fall back to GFF parsing). The values below are realistic
# protein-coding gene counts for the reference isolate of each species.
# ---------------------------------------------------------------------------

DEMO_GFF_GENE_COUNTS: dict[str, int] = {
    "GCA_964341285.1": 5544,    # Bacillus anthracis ~5500 genes (Ames-like ref)
    "GCA_900654255.1": 4173,    # Mycobacterium tuberculosis H37Rv ~4000-4200
    "GCA_001049615.1": 2872,    # Staphylococcus aureus N315-like ~2870
}
