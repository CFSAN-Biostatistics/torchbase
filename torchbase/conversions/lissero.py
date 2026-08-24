"""LisSero converter — local database files to torch format.

LisSero serotypes Listeria monocytogenes from presence/absence of five
PCR-target genes -- Prs (species gate), lmo0737, lmo1118, ORF2110, ORF2819 --
against upstream's single combined `lissero/db/sequences.fasta`. If a
profiles TSV is provided it is used as-is; otherwise this converter writes
the canonical table below, a mechanical transcription of upstream's
`Serotype.report_maker()` decision tree (see `_build_canonical_profiles`),
so a torch built with `--download` is loadable *and* typeable.

torchtools convert lissero --download --output torches/

See docs/adr/0003-profile-matching-value-domain-agnostic.md.
"""

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, IO, Optional

import toml

from torchbase.conversions import parse_fasta_records
from torchbase.quality.kmer_analysis import analyze_locus
from torchbase.conversions.log import get_logger, TRACE

_log = get_logger("lissero")


TAXA = ["Listeria monocytogenes"]

# The five genes LisSero's Serotype.report_maker() actually inspects, in its
# own declared dict order. (This corrects an earlier, wrong locus list --
# prs/LMOSA/LMOSB/ORF2110/ORF2819/ldh/lin0764/lin1118 -- that named genes not
# present in LisSero's database at all.)
CANONICAL_LOCI = ["Prs", "lmo0737", "lmo1118", "ORF2110", "ORF2819"]

# A single-allele locus's one FASTA record is always named "{locus}_1" (see
# download_sources), so calls_from_presence reports allele_id="1" on a hit --
# matching that here means "Prs: 1" in profiles.tsv, the same allele-ID
# convention every other profiles.tsv in torchbase uses; ABSENT has no real
# allele to name and stays a true sentinel.
PRESENT, ABSENT, WILD = "1", "absent", "?"

_LISSERO_RAW = "https://raw.githubusercontent.com/MDU-PHL/LisSero/master"

DOWNLOAD_SOURCES = {
    "repo": "https://github.com/MDU-PHL/LisSero",
    # Upstream's database is ONE combined FASTA (`lissero/db/sequences.fasta`),
    # headers "{gene}~~{accession}:{coords} {description}"; there is no
    # per-locus file and no downloadable profiles table (the serotype rule
    # is Python source, transcribed into _build_canonical_profiles below).
    "reference": (
        "sequences.fasta",
        f"{_LISSERO_RAW}/lissero/db/sequences.fasta",
    ),
}


def download_sources(dest_dir: Path) -> dict:
    """Download LisSero's combined reference and split it into one file per gene.

    Returns {'sequences': [list of open file handles]}.
    """
    from torchbase.conversions import fetch_file

    dest_dir = Path(dest_dir)
    filename, url = DOWNLOAD_SOURCES["reference"]
    combined = fetch_file(url, dest_dir / filename)

    genes_dir = dest_dir / "genes"
    genes_dir.mkdir(exist_ok=True)
    sequence_files = []
    for header, sequence in parse_fasta_records(combined.read_text(encoding="utf-8")):
        # Upstream truncates a hit name at "~~" too (Serotype._blast_parse).
        gene = header.split("~~")[0]
        gene_path = genes_dir / f"{gene}.fasta"
        # "_1": an explicit allele id, matching extract_locus_and_allele's
        # "split at the last underscore" convention -- and the real token
        # calls_from_presence reports on a hit (see PRESENT above).
        gene_path.write_text(f">{gene}_1\n{sequence}\n", encoding="utf-8")
        sequence_files.append(open(gene_path, encoding="utf-8"))

    _log.info("split %s into %d gene loci", filename, len(sequence_files))
    return {"sequences": sequence_files}


# --------------------------------------------------------------------------
# Canonical serogroup table: a mechanical transcription of upstream's
# Serotype.report_maker() decision tree (lissero/scripts/Serotype.py) into
# profiles.tsv rows. torchbase.allele_calls.calls_from_presence's
# present/absent tokens plus torchbase.profile_match's first-matching-row
# semantics reproduce it exactly, with the same row priority as upstream's
# if/elif chain: Prs- is checked (and wins) unconditionally first.
# --------------------------------------------------------------------------


def _row(profile_id: str, **locus_values: str) -> Dict[str, str]:
    row = {"Serogroup": profile_id}
    for locus in CANONICAL_LOCI:
        row[locus] = locus_values.pop(locus, WILD)
    if locus_values:
        raise ValueError(f"unknown locus column(s): {sorted(locus_values)}")
    return row


def _build_canonical_profiles() -> List[Dict[str, str]]:
    """Build the canonical profiles.tsv rows, in upstream's own priority order."""
    return [
        _row("Nontypeable", Prs=ABSENT),
        _row("1/2c, 3c", Prs=PRESENT, lmo0737=PRESENT, ORF2819=ABSENT, ORF2110=ABSENT, lmo1118=PRESENT),
        _row("1/2a, 3a", Prs=PRESENT, lmo0737=PRESENT, ORF2819=ABSENT, ORF2110=ABSENT, lmo1118=ABSENT),
        _row("4b, 4d, 4e", Prs=PRESENT, ORF2819=PRESENT, lmo0737=ABSENT, lmo1118=ABSENT, ORF2110=PRESENT),
        _row("1/2b, 3b, 7", Prs=PRESENT, ORF2819=PRESENT, lmo0737=ABSENT, lmo1118=ABSENT, ORF2110=ABSENT),
        _row("Nontypeable", Prs=PRESENT, lmo0737=PRESENT, ORF2819=PRESENT, ORF2110=PRESENT, lmo1118=PRESENT),
        _row("4b, 4d, 4e*", Prs=PRESENT, lmo0737=PRESENT, ORF2819=PRESENT, ORF2110=PRESENT, lmo1118=ABSENT),
        # Everything else (Prs+, no combination above matched) is upstream's
        # own final `else: Nontypeable` -- left as the natural novel_profile
        # fallthrough rather than an explicit row.
    ]


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
        sequence_files: Open file handles for per-gene FASTA files. Canonical
            loci: Prs, lmo0737, lmo1118, ORF2110, ORF2819 (see
            `download_sources`). Filenames should match locus names.
        profiles_file: Optional open file handle for a serogroup definitions
            TSV, overriding the canonical table this converter would
            otherwise write (see module docstring). Columns: Serogroup, then
            one column per locus in `CANONICAL_LOCI`.
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
    _log.info("Starting lissero conversion → %s/%s %s", namespace, name, version)
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
    elif set(locus_names) == set(CANONICAL_LOCI):
        _write_canonical_profiles(profiles_dest)
        profile_count = len(_build_canonical_profiles())
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
            "scheme": "LisSero PCR-target serogroup",
            "loci_count": len(locus_names),
            "profiles_count": profile_count,
            "calling_mode": "presence_absence",
            "id_column": "Serogroup",
        },
        "description": {
            "short": "LisSero Listeria monocytogenes serogroup torch",
            "long": (
                "Listeria monocytogenes serogroup determination by "
                "presence/absence of five PCR-target loci: "
                + ", ".join(CANONICAL_LOCI) + ". profiles.tsv is a "
                "mechanical transcription of LisSero's "
                "Serotype.report_maker() decision tree into rows, matched "
                "by presence/absence (torchbase.allele_calls."
                "calls_from_presence); see "
                "docs/adr/0003-profile-matching-value-domain-agnostic.md."
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
        csv.writer(f, delimiter="\t").writerow(["Serogroup"] + locus_names)


def _write_canonical_profiles(profiles_path: Path) -> None:
    with open(profiles_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Serogroup"] + CANONICAL_LOCI, delimiter="\t")
        writer.writeheader()
        writer.writerows(_build_canonical_profiles())


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
