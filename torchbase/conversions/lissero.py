"""LisSero converter — local database files to torch format.

LisSero serotypes Listeria monocytogenes based on the presence or absence of
eight PCR-target loci: prs, LMOSA, LMOSB, ORF2110, ORF2819, ldh, lin0764,
lin1118. Serogroup is determined by the binary pattern of detected loci.

Expected serogroup table (TSV) columns:
    Serogroup, prs, LMOSA, LMOSB, ORF2110, ORF2819, ldh, lin0764, lin1118
Values in each locus column should be 1 (present) or 0 (absent).

Usage:
    torchtools convert lissero --profiles serogroups.tsv \\
        prs.fasta LMOSA.fasta LMOSB.fasta ORF2110.fasta ORF2819.fasta \\
        ldh.fasta lin0764.fasta lin1118.fasta
"""

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import List, IO, Optional

import toml

from torchbase.quality.kmer_analysis import analyze_locus


TAXA = ["Listeria monocytogenes"]

# Canonical LisSero loci in expected order
CANONICAL_LOCI = ["prs", "LMOSA", "LMOSB", "ORF2110", "ORF2819", "ldh", "lin0764", "lin1118"]

_LISSERO_RAW = "https://raw.githubusercontent.com/MDU-PHL/LisSero/master"

DOWNLOAD_SOURCES = {
    "repo": "https://github.com/MDU-PHL/LisSero",
    "sequences": [
        (f"{locus}.fasta", f"{_LISSERO_RAW}/lissero/db/{locus}.fasta")
        for locus in CANONICAL_LOCI
    ],
    "profiles": (
        "Serogroups.tsv",
        f"{_LISSERO_RAW}/lissero/db/Serogroups.tsv",
    ),
}


def download_sources(dest_dir: Path) -> dict:
    """Download canonical LisSero database files to dest_dir.

    Returns {'sequences': [list of open file handles], 'profiles': open file handle}.
    """
    from torchbase.conversions import fetch_file

    dest_dir = Path(dest_dir)

    sequence_files = []
    for filename, url in DOWNLOAD_SOURCES["sequences"]:
        path = fetch_file(url, dest_dir / filename)
        sequence_files.append(open(path))

    prof_name, prof_url = DOWNLOAD_SOURCES["profiles"]
    prof_path = fetch_file(prof_url, dest_dir / prof_name)

    return {
        "sequences": sequence_files,
        "profiles": open(prof_path),
    }


def convert_local(
    sequence_files: List[IO],
    profiles_file: Optional[IO] = None,
    output_path: str = ".",
    namespace: str = "lissero",
    name: str = "lissero",
    version: str = "1.0.0",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> str:
    """Convert local LisSero database files to torch format.

    Args:
        sequence_files: Open file handles for per-locus FASTA files.
            Canonical loci: prs, LMOSA, LMOSB, ORF2110, ORF2819, ldh,
            lin0764, lin1118. Filenames should match locus names.
        profiles_file: Open file handle for serogroup definitions TSV.
            Columns: Serogroup, prs, LMOSA, LMOSB, ORF2110, ORF2819,
            ldh, lin0764, lin1118 (binary 1/0 per locus). If absent, a
            stub header-only table is written so the torch is loadable.
        output_path: Directory in which to create the torch.
        namespace: Torch namespace (default: "lissero").
        name: Torch name (default: "lissero").
        version: Torch version string.
        kmer_size: K-mer size for quality analysis.
        overlap_threshold: Overlap similarity threshold.
        duplicate_threshold: Duplicate similarity threshold.

    Returns:
        Path to the created torch directory.
    """
    output_path = Path(output_path)
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

    profiles_dest = torch_dir / "profiles.tsv"
    if profiles_file is not None:
        shutil.copy2(profiles_file.name, profiles_dest)
        profiles_file.seek(0)
        rows = list(csv.reader(profiles_file, delimiter="\t"))
        profile_count = max(0, len(rows) - 1)
    else:
        _write_stub_profiles(profiles_dest, locus_names or CANONICAL_LOCI)
        profile_count = 0

    quality_results = _run_quality_analysis(
        resources_dir, kmer_size, overlap_threshold, duplicate_threshold
    )

    now = datetime.now(timezone.utc).isoformat()
    metadata = {
        "namespace": namespace,
        "name": name,
        "version": version,
        "version_info": {"strategy": "snapshot", "timestamp": now},
        "typing": {
            "method": "serotyping",
            "scheme": "LisSero PCR-target serogroup",
            "loci_count": len(locus_names),
            "profiles_count": profile_count,
        },
        "description": {
            "short": "LisSero Listeria monocytogenes serogroup torch",
            "long": (
                "Listeria monocytogenes serogroup determination by "
                "presence/absence of eight PCR-target loci: "
                + ", ".join(CANONICAL_LOCI)
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

    quality_report = _build_quality_report(locus_names, quality_results, kmer_size, overlap_threshold)
    with open(torch_dir / "quality.json", "w") as f:
        json.dump(quality_report, f, indent=2)

    return str(torch_dir)


def _write_stub_profiles(profiles_path: Path, locus_names: List[str]) -> None:
    with open(profiles_path, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerow(["Serogroup"] + locus_names)


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
        if report.suspect_pairs:
            results["suspect_loci"] += 1
            for pair in report.suspect_pairs:
                bucket = "duplicate_pairs" if pair.get("issue_type") == "duplicate" else "similar_pairs"
                results[bucket].append({"locus": locus_name, **pair})
    return results
