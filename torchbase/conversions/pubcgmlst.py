"""PubMLST cgMLST converter — local files to torch format.

Converts a cgMLST scheme provided as local files (profiles TSV + per-locus
FASTA files) into a torch directory structure.

Usage:
    torchtools convert pubcgmlst profiles.tsv locus1.fasta locus2.fasta ...
"""

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, IO, Optional
import shutil

import toml

from torchbase.quality.kmer_analysis import analyze_locus
from torchbase.conversions.log import get_logger, TRACE

_log = get_logger("pubcgmlst")


def convert_local(
    scheme_file: IO,
    sequence_files: List[IO],
    output_path: str = ".",
    namespace: str = "pubcgmlst",
    name: Optional[str] = None,
    version: str = "1.0.0",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> str:
    """Convert local cgMLST scheme files to torch format.

    Args:
        scheme_file: Open file handle for profiles TSV (cgST/ID + locus allele columns)
        sequence_files: Open file handles for per-locus FASTA files
        output_path: Directory in which to create the torch
        namespace: Torch namespace (default: "pubcgmlst")
        name: Torch name; inferred from scheme filename if None
        version: Torch version string
        kmer_size: K-mer size for quality analysis
        overlap_threshold: Overlap threshold for quality analysis
        duplicate_threshold: Duplicate threshold for quality analysis

    Returns:
        Path to the created torch directory
    """
    output_path = Path(output_path)

    # Infer torch name from scheme filename if not provided
    if name is None:
        name = Path(scheme_file.name).stem.replace(" ", "_").lower()

    _log.info("Starting pubcgmlst conversion → %s/%s %s", namespace, name, version)

    torch_dir = output_path / namespace / name / f"{version}.torch"
    torch_dir.mkdir(parents=True, exist_ok=True)

    resources_dir = torch_dir / "_resources"
    resources_dir.mkdir(exist_ok=True)

    # Copy FASTA files into _resources/
    locus_names = []
    for fasta_fh in sequence_files:
        src = Path(fasta_fh.name)
        dest = resources_dir / src.name
        shutil.copy2(src, dest)
        locus_names.append(src.stem)
        _log.debug("  copying %s → %s (%d bytes)", src.name, dest, dest.stat().st_size)

    _log.info("  %d locus file(s): %s", len(locus_names), ", ".join(locus_names))

    # Copy profiles TSV
    profiles_dest = torch_dir / "profiles.tsv"
    shutil.copy2(scheme_file.name, profiles_dest)

    # Count profiles
    scheme_file.seek(0)
    reader = csv.reader(scheme_file, delimiter="\t")
    rows = list(reader)
    profile_count = max(0, len(rows) - 1)  # subtract header

    _log.info("  profiles: %d rows", profile_count)
    _log.debug("  profiles written to %s (%d rows)", profiles_dest, profile_count)

    # Run k-mer quality analysis
    _log.info("  running k-mer quality analysis on %d loci (k=%d)", len(locus_names), kmer_size)
    quality_results = _run_quality_analysis(
        resources_dir, kmer_size, overlap_threshold, duplicate_threshold
    )
    _log.info("  quality summary: %d suspect loci, %d duplicate pairs",
              quality_results.get("suspect_loci", 0), len(quality_results.get("duplicate_pairs", [])))

    now = datetime.now(timezone.utc).isoformat()

    metadata = {
        "namespace": namespace,
        "name": name,
        "version": version,
        "version_info": {
            "strategy": "snapshot",
            "timestamp": now,
        },
        "typing": {
            "method": "cgmlst",
            "loci_count": len(locus_names),
            "profiles_count": profile_count,
        },
        "data_quality": {
            "kmer_size": kmer_size,
            "overlap_threshold": overlap_threshold,
            "duplicate_threshold": duplicate_threshold,
        },
        "manifest": {
            "profiles": "profiles.tsv",
            "resources": [f.name for f in sorted(resources_dir.iterdir())
                          if f.is_file() and not f.name.startswith(".")],
        },
    }

    with open(torch_dir / "metadata.toml", "w", encoding="utf-8") as f:
        toml.dump(metadata, f)
    _log.debug("  metadata.toml written")

    quality_report = {
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

    with open(torch_dir / "quality.json", "w", encoding="utf-8") as f:
        json.dump(quality_report, f, indent=2)
    _log.debug("  quality.json written")

    _log.info("Torch written: %s", torch_dir)
    return str(torch_dir)


def _run_quality_analysis(
    resources_dir: Path,
    kmer_size: int,
    overlap_threshold: float,
    duplicate_threshold: float,
) -> dict:
    results = {
        "total_loci": 0,
        "suspect_loci": 0,
        "similar_pairs": [],
        "duplicate_pairs": [],
        "loci_results": {},
    }

    for fasta_file in sorted(resources_dir.glob("*.fasta")):
        locus_name = fasta_file.stem
        report = analyze_locus(
            fasta_file,
            k_size=kmer_size,
            overlap_threshold=overlap_threshold * 100,
            duplicate_threshold=duplicate_threshold * 100,
        )

        results["total_loci"] += 1
        locus_entry = {
            "suspect_pairs": report.suspect_pairs,
            "threshold": report.threshold,
            "allele_count": len(report.similarities) + 1 if report.similarities else 0,
            "similarity_stats": {},
        }
        if report.similarities:
            sims = sorted(report.similarities)
            n = len(sims)
            locus_entry["similarity_stats"] = {
                "min": sims[0],
                "median": sims[n // 2],
                "percentile_99": sims[min(int(n * 0.99), n - 1)],
            }

        results["loci_results"][locus_name] = locus_entry

        _log.debug("    %s: %d alleles, %d suspect pairs (threshold=%.3f)",
                   locus_name, locus_entry["allele_count"], len(report.suspect_pairs), report.threshold)
        for pair in report.suspect_pairs:
            _log.log(TRACE, "      suspect: %s ↔ %s  sim=%.4f  type=%s",
                     pair.get("allele1"), pair.get("allele2"),
                     pair.get("similarity", 0), pair.get("issue_type", "?"))

        if report.suspect_pairs:
            results["suspect_loci"] += 1
            for pair in report.suspect_pairs:
                if pair.get("issue_type") == "duplicate":
                    results["duplicate_pairs"].append({"locus": locus_name, **pair})
                else:
                    results["similar_pairs"].append({"locus": locus_name, **pair})

    return results
