"""ECTyper converter — local database files to torch format.

ECTyper serotypes E. coli (and Shigella) by typing O-antigen and H-antigen loci.
The ECTyper database FASTA uses headers prefixed with the antigen class:
    >O1-1-wzm_1   (O-antigen gene)
    >H2-1-fliC_1  (H-antigen gene)

The converter splits the combined database into O-antigen and H-antigen FASTA
files when a single combined FASTA is provided, or accepts pre-split files.

Expected serotype table (TSV) columns: Serotype, O, H
Missing antigens should use "-" or "ND" per ECTyper convention.

Usage:
    # Combined ECTyper database FASTA
    torchtools convert ectyper --profiles ectyper.tsv --db ECTyperDB.fasta

    # Pre-split per-antigen FASTA files
    torchtools convert ectyper --profiles ectyper.tsv O_antigen.fasta H_antigen.fasta
"""

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import List, IO, Optional

import toml

from torchbase.quality.kmer_analysis import analyze_locus


TAXA = ["Escherichia coli", "Shigella"]


def convert_local(
    sequence_files: List[IO],
    profiles_file: Optional[IO] = None,
    db_fasta: Optional[IO] = None,
    output_path: str = ".",
    namespace: str = "ectyper",
    name: str = "ectyper",
    version: str = "1.0.0",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> str:
    """Convert local ECTyper database files to torch format.

    Args:
        sequence_files: Open file handles for per-antigen FASTA files (O, H).
        profiles_file: Open file handle for serotype definitions TSV.
            Columns: Serotype, O, H (at minimum).
        db_fasta: Open file handle for a combined ECTyper DB FASTA. When
            provided, sequences are split into O-antigen and H-antigen files
            based on header prefixes (O* / H*). Takes precedence over
            sequence_files if both are given.
        output_path: Directory in which to create the torch.
        namespace: Torch namespace (default: "ectyper").
        name: Torch name (default: "ectyper").
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

    if db_fasta is not None:
        locus_names = _split_combined_fasta(db_fasta, resources_dir)
    else:
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
        _write_stub_profiles(profiles_dest)
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
            "scheme": "O:H antigen",
            "loci_count": len(locus_names),
            "profiles_count": profile_count,
        },
        "description": {
            "short": "ECTyper E. coli / Shigella serotyping torch",
            "long": "E. coli and Shigella serotyping based on O-antigen and H-antigen loci",
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


def _split_combined_fasta(db_fasta_fh: IO, resources_dir: Path) -> List[str]:
    """Split an ECTyper combined FASTA into O-antigen and H-antigen files.

    Headers starting with 'O' go to O_antigen.fasta; those starting with 'H'
    go to H_antigen.fasta. Unrecognised prefixes go to other_antigen.fasta.
    """
    buckets: dict = {"O": [], "H": [], "other": []}
    current_bucket = None
    current_lines: List[str] = []

    def flush(bucket, lines):
        if lines:
            buckets[bucket].extend(lines)

    for raw_line in db_fasta_fh:
        line = raw_line.rstrip("\n")
        if line.startswith(">"):
            if current_bucket is not None:
                flush(current_bucket, current_lines)
            header = line[1:]
            if header.upper().startswith("O"):
                current_bucket = "O"
            elif header.upper().startswith("H"):
                current_bucket = "H"
            else:
                current_bucket = "other"
            current_lines = [line]
        else:
            if current_bucket is not None:
                current_lines.append(line)

    if current_bucket is not None:
        flush(current_bucket, current_lines)

    locus_names = []
    bucket_map = {"O": "O_antigen", "H": "H_antigen", "other": "other_antigen"}
    for key, fname in bucket_map.items():
        if buckets[key]:
            dest = resources_dir / f"{fname}.fasta"
            dest.write_text("\n".join(buckets[key]) + "\n")
            locus_names.append(fname)

    return locus_names


def _write_stub_profiles(profiles_path: Path) -> None:
    with open(profiles_path, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerow(["Serotype", "O", "H"])


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
