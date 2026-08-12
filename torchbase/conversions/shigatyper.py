"""ShigaTyper converter — local FASTA files to torch format.

ShigaTyper identifies Shigella serotypes from presence/absence of marker
genes (species markers `ipaH`/`ipaB`/`EclacY`/`cadA`, plus per-serotype
O-antigen `wzx`/`wzy` and accessory loci) mapped by minimap2, evaluated
through a checkpoint decision cascade — not a lookup table, and not (yet)
something any built-in typing model in torchbase runs. See "Known gap" below.

If a profiles TSV is provided it is used as-is; otherwise a stub header-only
table is written so the torch is immediately loadable.

Usage:
    torchtools convert shigatyper --download
    torchtools convert shigatyper --profiles serotypes.tsv markers/*.fasta

Known gap: this converter produces the reference data — one locus per marker,
matching ShigaTyper's presence/absence model — but no typing model in
torchbase implements ShigaTyper's actual decision cascade (checkpoints,
mpileup depth-ratio tie-breaks, hard-coded exceptions like the S. boydii
9/15 vs EIEC check). The allelic model's Profile matching does not apply: a
serotype call here is not "which allele at each locus", it is a rule cascade
over which markers are present at all, evaluated in a fixed order with
several special cases. Encoding that cascade as profiles.tsv rows would be
fabricating logic that does not exist yet, and — per docs/operon-strategy-plan.md
§9 risk 1 — a plausible-looking wrong answer at full confidence is worse than
no answer. Treat this torch as reference data pending a real "shigatyper"
typing model (transcribed from shigatyper.py's `predict` cascade, the way
`torchbase/operon.py` transcribes StxTyper's), not as ready to type from.
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

_log = get_logger("shigatyper")


TAXA = ["Shigella"]

_SHIGATYPER_RAW = "https://raw.githubusercontent.com/CFSAN-Biostatistics/shigatyper/master"

DOWNLOAD_SOURCES = {
    "repo": "https://github.com/CFSAN-Biostatistics/shigatyper",
    # Upstream consolidated its per-locus FASTA files (wzx.fasta, wzy.fasta,
    # fliC.fasta -- none of which exist any more) into one multi-marker
    # reference. download_sources() below splits it back into one file per
    # marker so each becomes its own presence/absence locus, matching
    # convert_local's one-file-per-locus convention.
    "reference": (
        "ShigellaRef5.fasta",
        f"{_SHIGATYPER_RAW}/shigatyper/resources/ShigellaRef5.fasta",
    ),
    # Serotype profiles are not available as a standalone TSV in the repo,
    # and could not be transcribed as one even if they were -- see the
    # "Known gap" note above. Provide --profiles manually if you have one.
}


def _parse_fasta_records(text: str):
    """Yield (header, sequence) pairs from FASTA text."""
    header, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:].split()[0], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def download_sources(dest_dir: Path) -> dict:
    """Download ShigaTyper's consolidated reference and split it into
    one-marker-per-file loci.

    Returns {'sequences': [list of open file handles]}.
    Serotype profiles must be provided separately via --profiles, if at all
    (see the module docstring's "Known gap").
    """
    from torchbase.conversions import fetch_file

    dest_dir = Path(dest_dir)
    filename, url = DOWNLOAD_SOURCES["reference"]
    combined = fetch_file(url, dest_dir / filename)

    markers_dir = dest_dir / "markers"
    markers_dir.mkdir(exist_ok=True)
    sequence_files = []
    for marker, sequence in _parse_fasta_records(
        combined.read_text(encoding="utf-8")
    ):
        marker_path = markers_dir / f"{marker}.fasta"
        marker_path.write_text(f">{marker}\n{sequence}\n", encoding="utf-8")
        sequence_files.append(open(marker_path, encoding="utf-8"))

    _log.info("split %s into %d marker loci", filename, len(sequence_files))
    return {"sequences": sequence_files}


def convert_local(
    sequence_files: List[IO],
    profiles_file: Optional[IO] = None,
    output_path: str = ".",
    namespace: str = "hfp",
    name: str = "shigatyper",
    version: str = "1.0.0",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> str:
    """Convert local ShigaTyper FASTA files to torch format.

    Args:
        sequence_files: Open file handles for antigen gene FASTA files.
        profiles_file: Optional open file handle for serotype profiles TSV.
            Columns: Serotype, O, H (at minimum). If absent, a stub
            header-only table is written so the torch is loadable.
        output_path: Directory in which to create the torch.
        namespace: Torch namespace. Defaults to "hfp": a torch's namespace
            names the authority for the data, and the ShigaTyper database is
            an FDA Human Foods Program product
            (github.com/CFSAN-Biostatistics/ShigaTyper, from the center's
            former name). Compare "ncbi" for the StxTyper conversion and
            "pubmlst" for PubMLST schemes.
        name: Torch name (default: "shigatyper").
        version: Torch version string.
        kmer_size: K-mer size for quality analysis.
        overlap_threshold: Overlap similarity threshold.
        duplicate_threshold: Duplicate similarity threshold.

    Returns:
        Path to the created torch directory.
    """
    output_path = Path(output_path)
    _log.info("Starting shigatyper conversion → %s/%s %s", namespace, name, version)
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
            "scheme": "Shigella O:H antigen (presence/absence marker set)",
            "loci_count": len(locus_names),
            "profiles_count": profile_count,
        },
        "description": {
            "short": "ShigaTyper Shigella serotyping reference (no built-in typing model yet)",
            "long": (
                "Species and O-antigen/accessory marker sequences from "
                "CFSAN-Biostatistics/shigatyper's ShigellaRef5.fasta, one "
                "locus per marker. ShigaTyper calls a serotype via a "
                "checkpoint decision cascade over marker presence plus "
                "mpileup depth-ratio tie-breaks -- not a per-locus allele "
                "lookup -- and no typing model in torchbase implements that "
                "cascade yet. This torch is reference data pending one; see "
                "the module docstring in torchbase/conversions/shigatyper.py."
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
    with open(torch_dir / "metadata.toml", "w", encoding="utf-8") as f:
        toml.dump(metadata, f)
    _log.debug("  metadata.toml written")

    quality_report = _build_quality_report(locus_names, quality_results, kmer_size, overlap_threshold)
    with open(torch_dir / "quality.json", "w", encoding="utf-8") as f:
        json.dump(quality_report, f, indent=2)
    _log.debug("  quality.json written")

    _log.info("Torch written: %s", torch_dir)
    return str(torch_dir)


def _write_stub_profiles(profiles_path: Path, locus_names: List[str]) -> None:
    with open(profiles_path, "w", newline="", encoding="utf-8") as f:
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
