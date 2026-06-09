"""PubMLST Converter module for converting PubMLST schemes to torch format.

PubMLST changed their data licensing on 2025-01-01: allele sequences submitted
from that date onward require authentication and carry a non-commercial-only
license.  Alleles entered on or before 2024-12-31 remain freely accessible and
redistributable.  The default cutoff_date here reflects that boundary.

This module can build single-scheme or multi-scheme torches.  Pass a list of
scheme_ids to bundle multiple schemes (e.g., MLST + wgMLST for one organism)
into a single torch under the multi-scheme format.

Version string used: YYYY.MM.DD derived from the cutoff date, making the data
snapshot explicit (e.g., 2024.12.31 for the last freely-licensed snapshot,
matching tseemann/mlst v2.33.0).
"""

import json
import re
import tempfile  # noqa: F401  kept for future use
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Optional, Dict, List, Any, Sequence
import csv
import toml

from torchbase.conversions.bigsdb_client import BIGSdbClient, SchemeData, DatabaseInfo
from torchbase.conversions.log import get_logger, TRACE

_log = get_logger("pubmlst")

# PubMLST changed licensing for data submitted from this date onward.
PUBMLST_LICENSE_CUTOFF = date(2024, 12, 31)


def convert_all(
    base_url: str,
    output_path: str,
    namespace: str = "pubmlst",
    torch_name: str = "pubmlst",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
    cutoff_date: date = PUBMLST_LICENSE_CUTOFF,
    skip_errors: bool = True,
) -> str:
    """Convert every scheme in every PubMLST database into one multi-scheme torch.

    Enumerates all seqdef databases exposed by the API, then all schemes within
    each database, and writes them all into a single torch under the multi-scheme
    format.  Each scheme becomes one entry under schemes/ named
    {database}_{scheme_key} to avoid collisions across organisms.

    Args:
        base_url: PubMLST/BIGSdb REST API root (e.g., https://pubmlst.org/api)
        output_path: Directory where the torch tree will be written
        namespace: Torch namespace (default: "pubmlst")
        torch_name: Torch name (default: "pubmlst")
        kmer_size: K-mer size for quality analysis
        overlap_threshold: Overlap threshold for quality analysis
        duplicate_threshold: Duplicate threshold for quality analysis
        cutoff_date: Exclude alleles entered after this date
        skip_errors: If True, log failures per scheme and continue rather than
            aborting the whole run

    Returns:
        Path to the created torch directory

    Raises:
        BIGSdbError: If the database list cannot be fetched (always fatal)
    """
    import sys

    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    client = BIGSdbClient(base_url)
    torch_version = cutoff_date.strftime("%Y.%-m.%-d")

    torch_dir = output_path / namespace / torch_name / f"{torch_version}.torch"
    torch_dir.mkdir(parents=True, exist_ok=True)
    schemes_dir = torch_dir / "schemes"
    schemes_dir.mkdir(exist_ok=True)

    _log.info("Starting full PubMLST conversion (cutoff: %s)", cutoff_date)
    databases = client.list_databases()

    all_quality_results: Dict[str, Any] = {
        "total_loci": 0,
        "similar_pairs": [],
        "duplicate_pairs": [],
        "loci_results": {},
    }
    scheme_registry: Dict[str, Dict[str, Any]] = {}
    all_scheme_data: List[SchemeData] = []
    fetch_timestamp = datetime.now(timezone.utc).isoformat()
    provenance_entries: List[Dict[str, Any]] = []

    for db_info in databases:
        database = db_info.name
        try:
            scheme_list = client.list_schemes(database)
        except Exception as exc:
            if skip_errors:
                import logging as _stdlib_logging
                _log.warning("skipping database %s: %s", database, exc)
                print(f"[warn] skipping database {database}: {exc}", file=sys.stderr)
                continue
            raise

        _log.info("Processing database: %s (%d schemes)", database, len(scheme_list))
        for scheme_info in scheme_list:
            _log.debug("  scheme %d: %s", scheme_info.scheme_id, scheme_info.description)
            try:
                scheme_data = client.fetch_scheme(
                    database, scheme_info.scheme_id, cutoff_date=cutoff_date
                )
            except Exception as exc:
                if skip_errors:
                    print(
                        f"[warn] skipping {database}/scheme {scheme_info.scheme_id}: {exc}",
                        file=sys.stderr,
                    )
                    continue
                raise

            _log.debug("  fetched scheme %d: %s (%d loci)", scheme_info.scheme_id, scheme_data.metadata.name, len(scheme_data.loci))

            # Prefix with database name to avoid collisions across organisms.
            # Strip the common "pubmlst_" prefix and "_seqdef" suffix for readability.
            db_short = re.sub(r"^pubmlst_|_seqdef$", "", database)
            scheme_key = f"{db_short}_{_sanitize_name(scheme_data.metadata.name)}"

            organism_dir = schemes_dir / scheme_key
            organism_dir.mkdir(exist_ok=True)
            alleles_dir = organism_dir / "alleles"
            alleles_dir.mkdir(exist_ok=True)

            for locus in scheme_data.loci:
                fasta_text = client._fetch_alleles_fasta(
                    database, locus.locus_id, cutoff_date=cutoff_date
                )
                locus_fasta_path = alleles_dir / f"{locus.locus_id}.fasta"
                locus_fasta_path.write_text(fasta_text)
                _log.debug("    locus %s: wrote %s", locus.locus_id, locus_fasta_path)

            _write_profiles_tsv(
                organism_dir / "profiles.tsv", scheme_data.profiles.profiles
            )

            q = _run_quality_analysis(
                alleles_dir,
                kmer_size=kmer_size,
                overlap_threshold=overlap_threshold,
                duplicate_threshold=duplicate_threshold,
            )
            _log_quality_results(_log, scheme_key, q)
            all_quality_results["total_loci"] += q["total_loci"]
            all_quality_results["similar_pairs"].extend(q["similar_pairs"])
            all_quality_results["duplicate_pairs"].extend(q["duplicate_pairs"])
            all_quality_results["loci_results"].update(
                {f"{scheme_key}/{k}": v for k, v in q["loci_results"].items()}
            )

            scheme_registry[scheme_key] = {
                "organism": scheme_data.metadata.name,
                "database": database,
                "loci": [locus.locus_id for locus in scheme_data.loci],
            }
            all_scheme_data.append(scheme_data)
            provenance_entries.append(
                {"database": database, "scheme_id": scheme_info.scheme_id}
            )

    metadata = _generate_metadata(
        namespace=namespace,
        name=torch_name,
        version=torch_version,
        schemes=all_scheme_data,
        database_url=base_url,
        scheme_ids=[e["scheme_id"] for e in provenance_entries],
        scheme_registry=scheme_registry,
        fetch_timestamp=fetch_timestamp,
        kmer_size=kmer_size,
        overlap_threshold=overlap_threshold,
        duplicate_threshold=duplicate_threshold,
        quality_results=all_quality_results,
        cutoff_date=cutoff_date,
    )
    # Also record the per-database breakdown in provenance
    metadata["provenance"]["databases"] = provenance_entries

    metadata_path = torch_dir / "metadata.toml"
    with open(metadata_path, "w") as f:
        toml.dump(metadata, f)

    quality_report = _generate_quality_report(
        kmer_size=kmer_size,
        overlap_threshold=overlap_threshold,
        duplicate_threshold=duplicate_threshold,
        quality_results=all_quality_results,
    )
    quality_path = torch_dir / "quality.json"
    with open(quality_path, "w") as f:
        json.dump(quality_report, f, indent=2)

    _log.info("Torch written: %d schemes, %d total loci → %s",
              len(scheme_registry), all_quality_results['total_loci'], torch_dir)
    return str(torch_dir)


def convert_scheme(
    database_url: str,
    scheme_id: int,
    output_path: str,
    namespace: str = "pubmlst",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
    cutoff_date: date = PUBMLST_LICENSE_CUTOFF,
) -> str:
    """Convert a single PubMLST scheme to torch format.

    Convenience wrapper around convert_schemes() for the common case of a
    single scheme ID.
    """
    return convert_schemes(
        database_url=database_url,
        scheme_ids=[scheme_id],
        output_path=output_path,
        namespace=namespace,
        kmer_size=kmer_size,
        overlap_threshold=overlap_threshold,
        duplicate_threshold=duplicate_threshold,
        cutoff_date=cutoff_date,
    )


def convert_schemes(
    database_url: str,
    scheme_ids: Sequence[int],
    output_path: str,
    namespace: str = "pubmlst",
    torch_name: Optional[str] = None,
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
    cutoff_date: date = PUBMLST_LICENSE_CUTOFF,
) -> str:
    """Convert one or more PubMLST schemes into a multi-scheme torch.

    Each scheme_id becomes one subdirectory under schemes/ in the resulting
    torch.  When a single scheme_id is given the result is a valid single-entry
    multi-scheme torch (compatible with both single- and multi-scheme loaders).

    Args:
        database_url: Base URL of the PubMLST/BIGSdb REST API
        scheme_ids: One or more numeric scheme IDs to include
        output_path: Output directory where the torch tree will be written
        namespace: Namespace for the torch (default: "pubmlst")
        torch_name: Override the torch name; defaults to the name of the first
            scheme (sanitized)
        kmer_size: K-mer size for quality analysis
        overlap_threshold: Overlap threshold for quality analysis
        duplicate_threshold: Duplicate threshold for quality analysis
        cutoff_date: Exclude alleles entered after this date.  Defaults to
            2024-12-31, the last day of freely-redistributable PubMLST data.

    Returns:
        Path to the created torch directory

    Raises:
        ValueError: If scheme_ids is empty
        Exception: If scheme fetch or conversion fails
    """
    if not scheme_ids:
        raise ValueError("At least one scheme_id must be provided")

    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # If the user passed a full database URL (e.g. https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef)
    # split into base URL and database name; otherwise fall back to simple extraction.
    m = re.match(r'^(https?://[^/]+(?:/[^/]+)*)/db/([^/]+)/?$', database_url.rstrip('/'))
    if m:
        base_url = m.group(1)
        database_name = m.group(2)
    else:
        base_url = database_url
        database_name = _extract_database_name(database_url)

    _log.info("Converting %d scheme(s) from %s/%s (cutoff: %s)",
              len(scheme_ids), base_url, database_name, cutoff_date)
    client = BIGSdbClient(base_url)

    # Version is derived from the cutoff date to make the data snapshot explicit.
    # e.g. cutoff_date=2024-12-31 → version "2024.12.31"
    torch_version = cutoff_date.strftime("%Y.%-m.%-d")

    # Fetch all schemes first so we can derive a name before creating dirs.
    schemes: List[SchemeData] = []
    for sid in scheme_ids:
        sd = client.fetch_scheme(database_name, sid, cutoff_date=cutoff_date)
        _log.info("  scheme %d: %s (%d loci, %d profiles)", sd.metadata.scheme_id, sd.metadata.name,
                  len(sd.loci), sd.profiles.row_count if sd.profiles else 0)
        schemes.append(sd)

    derived_name = _sanitize_name(schemes[0].metadata.name)
    if torch_name is None:
        torch_name = derived_name

    torch_dir = output_path / namespace / torch_name / f"{torch_version}.torch"
    torch_dir.mkdir(parents=True, exist_ok=True)

    schemes_dir = torch_dir / "schemes"
    schemes_dir.mkdir(exist_ok=True)

    all_quality_results: Dict[str, Any] = {
        "total_loci": 0,
        "similar_pairs": [],
        "duplicate_pairs": [],
        "loci_results": {},
    }

    scheme_registry: Dict[str, Dict[str, Any]] = {}
    fetch_timestamp = datetime.now(timezone.utc).isoformat()

    for scheme_data in schemes:
        scheme_key = _sanitize_name(scheme_data.metadata.name)
        organism_dir = schemes_dir / scheme_key
        organism_dir.mkdir(exist_ok=True)
        alleles_dir = organism_dir / "alleles"
        alleles_dir.mkdir(exist_ok=True)

        allele_counts = {}
        for locus in scheme_data.loci:
            fasta_text = client._fetch_alleles_fasta(
                database_name, locus.locus_id, cutoff_date=cutoff_date
            )
            fasta_path = alleles_dir / f"{locus.locus_id}.fasta"
            fasta_path.write_text(fasta_text)
            allele_counts[locus.locus_id] = locus.alleles_count
            _log.debug("    locus %s: %d alleles", locus.locus_id, locus.alleles_count)

        profiles_path = organism_dir / "profiles.tsv"
        _write_profiles_tsv(profiles_path, scheme_data.profiles.profiles)

        q = _run_quality_analysis(
            alleles_dir,
            kmer_size=kmer_size,
            overlap_threshold=overlap_threshold,
            duplicate_threshold=duplicate_threshold,
        )
        _log_quality_results(_log, scheme_key, q)
        all_quality_results["total_loci"] += q["total_loci"]
        all_quality_results["similar_pairs"].extend(q["similar_pairs"])
        all_quality_results["duplicate_pairs"].extend(q["duplicate_pairs"])
        all_quality_results["loci_results"].update(q["loci_results"])

        scheme_registry[scheme_key] = {
            "organism": scheme_data.metadata.name,
            "loci": [locus.locus_id for locus in scheme_data.loci],
        }

    metadata = _generate_metadata(
        namespace=namespace,
        name=torch_name,
        version=torch_version,
        schemes=schemes,
        database_url=database_url,
        scheme_ids=list(scheme_ids),
        scheme_registry=scheme_registry,
        fetch_timestamp=fetch_timestamp,
        kmer_size=kmer_size,
        overlap_threshold=overlap_threshold,
        duplicate_threshold=duplicate_threshold,
        quality_results=all_quality_results,
        cutoff_date=cutoff_date,
    )

    metadata_path = torch_dir / "metadata.toml"
    with open(metadata_path, "w") as f:
        toml.dump(metadata, f)

    quality_report = _generate_quality_report(
        kmer_size=kmer_size,
        overlap_threshold=overlap_threshold,
        duplicate_threshold=duplicate_threshold,
        quality_results=all_quality_results,
    )
    quality_path = torch_dir / "quality.json"
    with open(quality_path, "w") as f:
        json.dump(quality_report, f, indent=2)

    _log.info("Torch written: %s", torch_dir)
    return str(torch_dir)


def _extract_database_name(url: str) -> str:
    """Derive a BIGSdb database name from a PubMLST API URL."""
    if "pubmlst" in url.lower():
        return "pubmlst"
    return "bigsdb"


def _sanitize_name(name: str) -> str:
    """Sanitize a scheme name for use in file paths."""
    return name.replace(" ", "_").replace("/", "_").lower()


def _write_profiles_tsv(profiles_path: Path, profiles: List[Dict[str, str]]) -> None:
    """Write profiles.tsv file."""
    if not profiles:
        with open(profiles_path, "w") as f:
            f.write("ST\n")
        return

    fieldnames = list(profiles[0].keys())
    with open(profiles_path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(profiles)


def _log_quality_results(log, scheme_key: str, q: Dict[str, Any]) -> None:
    """Emit quality summary at DEBUG and suspect pairs at TRACE."""
    suspect_loci = [k for k, v in q.get("loci_results", {}).items() if v.get("suspect_pairs")]
    duplicate_pairs = q.get("duplicate_pairs", [])
    log.debug("  quality: %s — %d loci, %d suspect loci, %d duplicate pairs",
              scheme_key, q.get("total_loci", 0), len(suspect_loci), len(duplicate_pairs))
    for locus_name, locus_data in q.get("loci_results", {}).items():
        for pair in locus_data.get("suspect_pairs", []):
            log.log(TRACE, "    suspect: %s allele %s ↔ %s (sim=%.3f, type=%s)",
                    locus_name,
                    pair.get("allele1"), pair.get("allele2"),
                    pair.get("similarity", 0), pair.get("issue_type", "?"))


def _run_quality_analysis(
    alleles_dir: Path,
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> Dict[str, Any]:
    """Run k-mer quality analysis on all loci in alleles_dir."""
    from torchbase.quality.kmer_analysis import analyze_locus

    results: Dict[str, Any] = {
        "total_loci": 0,
        "similar_pairs": [],
        "duplicate_pairs": [],
        "loci_results": {},
    }

    for fasta_file in sorted(alleles_dir.glob("*.fasta")):
        locus_name = fasta_file.stem
        report = analyze_locus(
            fasta_file,
            k_size=kmer_size,
            overlap_threshold=overlap_threshold * 100,
            duplicate_threshold=duplicate_threshold * 100,
        )

        results["total_loci"] += 1
        results["loci_results"][locus_name] = {
            "suspect_pairs": report.suspect_pairs,
            "threshold": report.threshold,
        }

        _log.debug("    %s: %d suspect pairs (threshold=%.3f)",
                   locus_name, len(report.suspect_pairs), report.threshold)
        for pair in report.suspect_pairs:
            _log.log(TRACE, "      suspect: %s ↔ %s  sim=%.4f  type=%s",
                     pair.get("allele1"), pair.get("allele2"),
                     pair.get("similarity", 0), pair.get("issue_type", "?"))
            entry = {
                "locus": locus_name,
                "allele1": pair["allele1"],
                "allele2": pair["allele2"],
                "similarity": pair["similarity"],
            }
            if pair["issue_type"] == "duplicate":
                results["duplicate_pairs"].append(entry)
            elif pair["issue_type"] == "overlap":
                results["similar_pairs"].append(entry)

    return results


def _generate_metadata(
    namespace: str,
    name: str,
    version: str,
    schemes: List[SchemeData],
    database_url: str,
    scheme_ids: List[int],
    scheme_registry: Dict[str, Any],
    fetch_timestamp: str,
    kmer_size: int,
    overlap_threshold: float,
    duplicate_threshold: float,
    quality_results: Dict[str, Any],
    cutoff_date: date,
) -> Dict[str, Any]:
    """Generate metadata.toml content."""
    total_profiles = sum(
        (s.profiles.row_count if s.profiles else 0) for s in schemes
    )

    return {
        "namespace": namespace,
        "name": name,
        "version": version,
        "version_info": {
            "strategy": "pubmlst-snapshot",
            "cutoff_date": cutoff_date.isoformat(),
            "fetch_timestamp": fetch_timestamp,
            # Alleles submitted to PubMLST after 2024-12-31 require
            # authentication and carry non-commercial-only terms.  This
            # snapshot contains only freely-redistributable pre-2025 data,
            # matching the dataset frozen in tseemann/mlst v2.33.0.
            "license_note": (
                "Alleles from PubMLST with date_entered <= "
                f"{cutoff_date.isoformat()} (freely redistributable; "
                "pre-2025 licensing terms apply)"
            ),
        },
        "provenance": {
            "source": "PubMLST",
            "database_url": database_url,
            "scheme_ids": scheme_ids,
            "fetch_date": fetch_timestamp,
        },
        "data_quality": {
            "kmer_analysis_performed": True,
            "kmer_size": kmer_size,
            "duplicate_threshold": duplicate_threshold,
            "overlap_threshold": overlap_threshold,
            "similar_alleles": quality_results.get("similar_pairs", []),
            "duplicate_alleles": quality_results.get("duplicate_pairs", []),
        },
        "typing": {
            "scheme_count": len(schemes),
            "loci_count": sum(len(s.loci) for s in schemes),
            "profiles_count": total_profiles,
        },
        "schemes": scheme_registry,
    }


def _generate_quality_report(
    kmer_size: int,
    overlap_threshold: float,
    duplicate_threshold: float,
    quality_results: Dict[str, Any],
) -> Dict[str, Any]:
    """Generate quality.json content."""
    return {
        "kmer_analysis": {
            "performed": True,
            "kmer_size": kmer_size,
            "parameters": {
                "overlap_threshold": overlap_threshold,
                "duplicate_threshold": duplicate_threshold,
            },
            "results": {
                "total_loci": quality_results.get("total_loci", 0),
                "similar_pairs": quality_results.get("similar_pairs", []),
                "duplicate_pairs": quality_results.get("duplicate_pairs", []),
            },
        },
    }
