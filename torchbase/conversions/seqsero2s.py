"""SeqSero2S converter — local database files to torch format.

SeqSero2S is a curated variant of SeqSero2 (github.com/LSTUGA/SeqSero2S) that:
- Uses a corrected allele database (H_and_O_and_specific_genes.fasta)
- Adds 7-gene MLST for sequence typing (aroC, dnaN, hemD, hisD, purE, sucA, thrA)
- Follows the simplified Kauffmann-White-Le Minor (KWS) nomenclature

Database inputs expected:
  - Antigen FASTA: H_and_O_and_specific_genes.fasta (O-antigen, fliC, fljB alleles)
  - MLST FASTAs:  aroC.tfa  dnaN.tfa  hemD.tfa  hisD.tfa  purE.tfa  sucA.tfa  thrA.tfa
                  (from seqsero2s_db/kmer/)
  - MLST profiles: salmonella_profile.txt (TSV: ST, aroC, dnaN, hemD, hisD, purE, sucA, thrA)
  - Serotype profiles: serotypes TSV (Serotype, O, H1, H2)

Usage:
    torchtools convert seqsero2s \\
        --antigen-db H_and_O_and_specific_genes.fasta \\
        --mlst-profiles salmonella_profile.txt \\
        --serotype-profiles serotypes.tsv \\
        aroC.tfa dnaN.tfa hemD.tfa hisD.tfa purE.tfa sucA.tfa thrA.tfa
"""

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import List, IO, Optional

import toml

from torchbase.quality.kmer_analysis import analyze_locus
from torchbase.conversions.log import get_logger, TRACE

_log = get_logger("seqsero2s")


TAXA = ["Salmonella enterica"]
MLST_LOCI = ["aroC", "dnaN", "hemD", "hisD", "purE", "sucA", "thrA"]

_SS2S_RAW = "https://raw.githubusercontent.com/LSTUGA/SeqSero2S/master"

DOWNLOAD_SOURCES = {
    "repo": "https://github.com/LSTUGA/SeqSero2S",
    "antigen_db": (
        "H_and_O_and_specific_genes.fasta",
        f"{_SS2S_RAW}/seqsero2s_db/H_and_O_and_specific_genes.fasta",
    ),
    "mlst_profiles": (
        "salmonella_profile.txt",
        f"{_SS2S_RAW}/seqsero2s_db/salmonella_profile.txt",
    ),
    "sequences": [
        (f"{locus}.fasta", f"{_SS2S_RAW}/seqsero2s_db/kmer/{locus}.tfa")
        for locus in MLST_LOCI
    ],
    # Serotype profiles (Serotype/O/H1/H2/ST) are not in the repo; provide
    # --serotype-profiles manually.
}


def download_sources(dest_dir: Path) -> dict:
    """Download canonical SeqSero2S files to dest_dir.

    Returns dict with keys matching convert_local kwargs:
      antigen_db, mlst_profiles, sequences (list of open file handles).
    Serotype profiles must be provided separately via --serotype-profiles.
    """
    from torchbase.conversions import fetch_file

    dest_dir = Path(dest_dir)

    adb_name, adb_url = DOWNLOAD_SOURCES["antigen_db"]
    antigen_db_path = fetch_file(adb_url, dest_dir / adb_name)

    mlst_name, mlst_url = DOWNLOAD_SOURCES["mlst_profiles"]
    mlst_profiles_path = fetch_file(mlst_url, dest_dir / mlst_name)

    sequence_files = []
    for filename, url in DOWNLOAD_SOURCES["sequences"]:
        path = fetch_file(url, dest_dir / filename)
        sequence_files.append(open(path))

    return {
        "antigen_db": open(antigen_db_path),
        "mlst_profiles": open(mlst_profiles_path),
        "sequences": sequence_files,
    }


def convert_local(
    sequence_files: List[IO],
    antigen_db: Optional[IO] = None,
    mlst_profiles: Optional[IO] = None,
    serotype_profiles: Optional[IO] = None,
    output_path: str = ".",
    namespace: str = "seqsero2s",
    name: str = "seqsero2s",
    version: str = "1.0.0",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> str:
    """Convert local SeqSero2S database files to torch format.

    Args:
        sequence_files: Open file handles for MLST locus FASTA files
            (aroC.tfa, dnaN.tfa, etc. — extension does not matter).
        antigen_db: Open file handle for the combined antigen FASTA
            (H_and_O_and_specific_genes.fasta). Copied as antigen_db.fasta.
        mlst_profiles: Open file handle for salmonella_profile.txt (TSV:
            ST, aroC, dnaN, hemD, hisD, purE, sucA, thrA).
        serotype_profiles: Open file handle for serotype definitions TSV
            (Serotype, O, H1, H2). Used as the primary profiles.tsv.
        output_path: Directory in which to create the torch.
        namespace: Torch namespace (default: "seqsero2s").
        name: Torch name (default: "seqsero2s").
        version: Torch version string.
        kmer_size: K-mer size for quality analysis.
        overlap_threshold: Overlap similarity threshold.
        duplicate_threshold: Duplicate similarity threshold.

    Returns:
        Path to the created torch directory.
    """
    output_path = Path(output_path)
    _log.info("Starting seqsero2s conversion → %s/%s %s", namespace, name, version)
    torch_dir = output_path / namespace / name / f"{version}.torch"
    torch_dir.mkdir(parents=True, exist_ok=True)

    resources_dir = torch_dir / "_resources"
    resources_dir.mkdir(exist_ok=True)

    mlst_dir = resources_dir / "mlst"
    mlst_dir.mkdir(exist_ok=True)

    # Copy MLST locus files
    locus_names = []
    for fasta_fh in sequence_files:
        src = Path(fasta_fh.name)
        # Normalise extension to .fasta
        dest_name = src.stem + ".fasta"
        dest = mlst_dir / dest_name
        shutil.copy2(src, dest)
        locus_names.append(src.stem)
        _log.debug("  copying %s → %s (%d bytes)", src.name, dest, dest.stat().st_size)

    _log.info("  %d locus file(s): %s", len(locus_names), ", ".join(locus_names))

    # Copy antigen database FASTA
    if antigen_db is not None:
        dest = resources_dir / "antigen_db.fasta"
        shutil.copy2(antigen_db.name, dest)
        _log.info("  antigen db: %s", Path(antigen_db.name).name if antigen_db else "none")
    else:
        _log.info("  antigen db: none")

    # Write MLST profiles into resources
    if mlst_profiles is not None:
        dest = resources_dir / "mlst_profiles.tsv"
        shutil.copy2(mlst_profiles.name, dest)
        mlst_profiles.seek(0)
        mlst_rows = list(csv.reader(mlst_profiles, delimiter="\t"))
        mlst_profile_count = max(0, len(mlst_rows) - 1)
    else:
        mlst_profile_count = 0

    _log.info("  MLST profiles: %d rows", mlst_profile_count)

    # Write serotype profiles as the primary profiles.tsv
    profiles_dest = torch_dir / "profiles.tsv"
    if serotype_profiles is not None:
        shutil.copy2(serotype_profiles.name, profiles_dest)
        serotype_profiles.seek(0)
        rows = list(csv.reader(serotype_profiles, delimiter="\t"))
        serotype_count = max(0, len(rows) - 1)
    else:
        _write_stub_serotype_profiles(profiles_dest)
        serotype_count = 0

    _log.info("  serotype profiles: %d rows", serotype_count)
    _log.info("  profiles: %d rows", serotype_count)
    _log.debug("  profiles written to %s (%d rows)", profiles_dest, serotype_count)

    _log.info("  running k-mer quality analysis on %d loci (k=%d)", len(locus_names), kmer_size)
    quality_results = _run_quality_analysis(
        mlst_dir, kmer_size, overlap_threshold, duplicate_threshold
    )
    _log.info("  quality summary: %d suspect loci, %d duplicate pairs",
              quality_results.get("suspect_loci", 0), len(quality_results.get("duplicate_pairs", [])))

    now = datetime.now(timezone.utc).isoformat()
    metadata = {
        "namespace": namespace,
        "name": name,
        "version": version,
        "version_info": {"strategy": "snapshot", "timestamp": now},
        "typing": {
            "method": "serotyping",
            "scheme": "SeqSero2S simplified KWS + 7-gene MLST",
            "antigen_loci": ["O-antigen", "fliC", "fljB"],
            "mlst_loci": MLST_LOCI,
            "serotype_count": serotype_count,
            "mlst_profile_count": mlst_profile_count,
        },
        "description": {
            "short": "SeqSero2S Salmonella serotyping + MLST torch",
            "long": (
                "Salmonella enterica serotyping (simplified KWS) combined with "
                "7-gene MLST sequence typing (aroC, dnaN, hemD, hisD, purE, sucA, thrA). "
                "Corrected allele database from LSTUGA/SeqSero2S."
            ),
            "taxa": TAXA,
        },
        "data_quality": {
            "kmer_size": kmer_size,
            "overlap_threshold": overlap_threshold,
            "duplicate_threshold": duplicate_threshold,
        },
        "manifest": {
            "profiles": "profiles.tsv",
            "resources": [
                f.name for f in sorted(resources_dir.iterdir())
                if f.is_file() and not f.name.startswith(".")
            ],
        },
    }
    with open(torch_dir / "metadata.toml", "w") as f:
        toml.dump(metadata, f)
    _log.debug("  metadata.toml written")

    quality_report = _build_quality_report(locus_names, quality_results, kmer_size, overlap_threshold)
    with open(torch_dir / "quality.json", "w") as f:
        json.dump(quality_report, f, indent=2)
    _log.debug("  quality.json written")

    _log.info("Torch written: %s", torch_dir)
    return str(torch_dir)


def _write_stub_serotype_profiles(profiles_path: Path) -> None:
    with open(profiles_path, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerow(["Serotype", "O", "H1", "H2", "ST"])


def _build_quality_report(locus_names, quality_results, kmer_size, overlap_threshold):
    return {
        "loci": {
            locus: {
                "allele_count": quality_results["loci_results"].get(locus, {}).get("allele_count", 0),
                "kmer_size": kmer_size,
                "similarity_stats": quality_results["loci_results"].get(locus, {}).get("similarity_stats", {}),
                "threshold": quality_results["loci_results"].get(locus, {}).get("threshold", overlap_threshold),
                "threshold_method": "gap_detection",
                "suspect_pairs": quality_results["loci_results"].get(locus, {}).get("suspect_pairs", []),
            }
            for locus in locus_names
        },
        "summary": {
            "total_loci": len(locus_names),
            "suspect_loci": quality_results.get("suspect_loci", 0),
            "suspect_alleles": len(quality_results.get("duplicate_pairs", [])),
        },
    }


def _run_quality_analysis(mlst_dir, kmer_size, overlap_threshold, duplicate_threshold):
    results = {
        "total_loci": 0,
        "suspect_loci": 0,
        "similar_pairs": [],
        "duplicate_pairs": [],
        "loci_results": {},
    }
    for fasta_file in sorted(mlst_dir.glob("*.fasta")):
        locus_name = fasta_file.stem
        report = analyze_locus(
            fasta_file,
            k_size=kmer_size,
            overlap_threshold=overlap_threshold * 100,
            duplicate_threshold=duplicate_threshold * 100,
        )
        results["total_loci"] += 1
        entry = {
            "suspect_pairs": report.suspect_pairs,
            "threshold": report.threshold,
            "allele_count": len(report.similarities) + 1 if report.similarities else 0,
            "similarity_stats": {},
        }
        if report.similarities:
            sims = sorted(report.similarities)
            n = len(sims)
            entry["similarity_stats"] = {
                "min": sims[0],
                "median": sims[n // 2],
                "percentile_99": sims[min(int(n * 0.99), n - 1)],
            }
        results["loci_results"][locus_name] = entry
        _log.debug("    %s: %d alleles, %d suspect pairs (threshold=%.3f)",
                   locus_name, entry["allele_count"], len(report.suspect_pairs), report.threshold)
        for pair in report.suspect_pairs:
            _log.log(TRACE, "      suspect: %s ↔ %s  sim=%.4f  type=%s",
                     pair.get("allele1"), pair.get("allele2"),
                     pair.get("similarity", 0), pair.get("issue_type", "?"))
        if report.suspect_pairs:
            results["suspect_loci"] += 1
            for pair in report.suspect_pairs:
                bucket = "duplicate_pairs" if pair.get("issue_type") == "duplicate" else "similar_pairs"
                results[bucket].append({"locus": locus_name, **pair})
    return results
