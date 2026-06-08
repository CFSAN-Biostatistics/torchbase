"""SeqSero2 converter — local database files to torch format.

SeqSero2 serotypes Salmonella by typing wzx/wzy (O-antigen), fliC (H1-antigen),
and fljB (H2-antigen) loci.

Expected serotype table (TSV) columns: Serotype, O, H1, H2
Each column after Serotype is an antigen formula component; missing or ambiguous
antigens should use "-" or "?" per the White-Kauffmann-Le Minor scheme.

Usage:
    torchtools convert seqsero2 --profiles serotypes.tsv wzx_wzy.fasta fliC.fasta fljB.fasta
    torchtools convert seqsero2 --profiles serotypes.tsv --db Salmonella_antigen_db.fasta
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

_log = get_logger("seqsero2")


TAXA = ["Salmonella"]

DOWNLOAD_SOURCES = {
    "repo": "https://github.com/CDCgov/SeqSero2",
    "sequences": [
        ("wzx_wzy.fasta",  "https://raw.githubusercontent.com/CDCgov/SeqSero2/master/SeqSero2_package/data/Otype_variable_region.fasta"),
        ("fliC.fasta",     "https://raw.githubusercontent.com/CDCgov/SeqSero2/master/SeqSero2_package/data/H1_fliC_variant.fasta"),
        ("fljB.fasta",     "https://raw.githubusercontent.com/CDCgov/SeqSero2/master/SeqSero2_package/data/H2_fljB_variant.fasta"),
    ],
    # Serotype profiles are compiled from NCBI / literature data and are not
    # available as a simple download from the repo. Provide --profiles manually.
}


def download_sources(dest_dir: Path) -> dict:
    """Download canonical SeqSero2 FASTA files to dest_dir.

    Returns {'sequences': [list of open file handles]}.
    Serotype profiles must be provided separately via --profiles.
    """
    from torchbase.conversions import fetch_file

    dest_dir = Path(dest_dir)
    sequence_files = []
    for filename, url in DOWNLOAD_SOURCES["sequences"]:
        path = fetch_file(url, dest_dir / filename)
        sequence_files.append(open(path))
    return {"sequences": sequence_files}


def convert_local(
    sequence_files: List[IO],
    profiles_file: Optional[IO] = None,
    output_path: str = ".",
    namespace: str = "seqsero2",
    name: str = "seqsero2",
    version: str = "1.0.0",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> str:
    """Convert local SeqSero2 database files to torch format.

    Args:
        sequence_files: Open file handles for antigen gene FASTA files.
            Typically: wzx_wzy.fasta (or separate wzx.fasta + wzy.fasta),
            fliC.fasta, fljB.fasta.
        profiles_file: Open file handle for serotype definitions TSV.
            Columns: Serotype, O, H1, H2 (at minimum). If absent a stub
            header-only table is written so the torch is loadable.
        output_path: Directory in which to create the torch.
        namespace: Torch namespace (default: "seqsero2").
        name: Torch name (default: "seqsero2").
        version: Torch version string.
        kmer_size: K-mer size for quality analysis.
        overlap_threshold: Overlap similarity threshold for quality analysis.
        duplicate_threshold: Duplicate similarity threshold for quality analysis.

    Returns:
        Path to the created torch directory.
    """
    output_path = Path(output_path)
    _log.info("Starting seqsero2 conversion → %s/%s %s", namespace, name, version)
    torch_dir = output_path / namespace / name / f"{version}.torch"
    torch_dir.mkdir(parents=True, exist_ok=True)

    resources_dir = torch_dir / "_resources"
    resources_dir.mkdir(exist_ok=True)

    locus_names = []
    for fasta_fh in sequence_files:
        src = Path(fasta_fh.name)
        dest = resources_dir / src.name
        shutil.copy2(src, dest)
        locus_names.append(src.stem)
        _log.debug("  copying %s → %s (%d bytes)", src.name, dest, dest.stat().st_size)

    _log.info("  %d locus file(s): %s", len(locus_names), ", ".join(locus_names))

    profiles_dest = torch_dir / "profiles.tsv"
    if profiles_file is not None:
        shutil.copy2(profiles_file.name, profiles_dest)
        profiles_file.seek(0)
        rows = list(csv.reader(profiles_file, delimiter="\t"))
        profile_count = max(0, len(rows) - 1)
    else:
        _write_stub_profiles(profiles_dest, locus_names)
        profile_count = 0

    _log.info("  profiles: %d rows", profile_count)
    _log.debug("  profiles written to %s (%d rows)", profiles_dest, profile_count)

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
        "version_info": {"strategy": "snapshot", "timestamp": now},
        "typing": {
            "method": "serotyping",
            "scheme": "White-Kauffmann-Le Minor",
            "loci_count": len(locus_names),
            "profiles_count": profile_count,
        },
        "description": {
            "short": "SeqSero2 Salmonella serotyping torch",
            "long": (
                "Salmonella serotyping based on wzx/wzy (O-antigen), "
                "fliC (H1-antigen), and fljB (H2-antigen) loci"
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


def _write_stub_profiles(profiles_path: Path, locus_names: List[str]) -> None:
    with open(profiles_path, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerow(["Serotype", "O", "H1", "H2"])


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


def _run_quality_analysis(resources_dir, kmer_size, overlap_threshold, duplicate_threshold):
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
