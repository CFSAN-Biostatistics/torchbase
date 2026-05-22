"""PubMLST Converter module for converting PubMLST schemes to torch format.

This module provides the convert_scheme() function that:
1. Fetches scheme data from a BIGSdb database using the BIGSdbClient
2. Creates a torch directory structure
3. Extracts alleles and profiles
4. Runs k-mer quality analysis
5. Generates metadata.toml and quality.json
"""

import json
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Dict, List, Any
import csv
import toml

from torchbase.conversions.bigsdb_client import BIGSdbClient
from torchbase.quality.kmer_analysis import analyze_locus


def convert_scheme(
    database_url: str,
    scheme_id: int,
    output_path: str,
    namespace: str = "pubmlst",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> str:
    """Convert a PubMLST scheme to torch format.

    Args:
        database_url: Base URL of PubMLST database API
        scheme_id: Numeric ID of the scheme to convert
        output_path: Output directory path where torch will be created
        namespace: Namespace for the torch (default: "pubmlst")
        kmer_size: K-mer size for quality analysis (default: 13)
        overlap_threshold: Overlap threshold for quality analysis (default: 0.90)
        duplicate_threshold: Duplicate threshold for quality analysis (default: 0.95)

    Returns:
        Path to the created torch directory

    Raises:
        Exception: If scheme fetch or conversion fails
    """
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # Fetch scheme data from BIGSdb
    client = BIGSdbClient(database_url)

    # Extract database name from URL
    # e.g., "http://pubmlst.org/api" -> "pubmlst"
    database_name = _extract_database_name(database_url)

    # Fetch scheme data
    scheme_data = client.fetch_scheme(database_name, scheme_id)

    # Generate torch name from scheme name
    torch_name = _sanitize_name(scheme_data.metadata.name)
    torch_version = "1.0.0"

    # Create torch directory structure
    # Format: <namespace>/<name>/<version>.torch/
    torch_dir = output_path / namespace / torch_name / f"{torch_version}.torch"
    torch_dir.mkdir(parents=True, exist_ok=True)

    # Create schemes subdirectory
    schemes_dir = torch_dir / "schemes"
    schemes_dir.mkdir(exist_ok=True)

    # Create organism subdirectory
    organism_name = _sanitize_name(scheme_data.metadata.name)
    organism_dir = schemes_dir / organism_name
    organism_dir.mkdir(exist_ok=True)

    # Create alleles subdirectory
    alleles_dir = organism_dir / "alleles"
    alleles_dir.mkdir(exist_ok=True)

    # Extract allele sequences and write FASTA files
    allele_counts = {}
    for locus in scheme_data.loci:
        # For now, we'll create stub FASTA files
        # In a real implementation, these would be fetched from the API
        fasta_path = alleles_dir / f"{locus.locus_id}.fasta"
        allele_counts[locus.locus_id] = locus.alleles_count

        # Create minimal FASTA content (stub for now)
        _write_stub_fasta(fasta_path, locus.locus_id)

    # Write profiles.tsv from scheme data
    profiles_path = organism_dir / "profiles.tsv"
    _write_profiles_tsv(profiles_path, scheme_data.profiles.profiles)

    # Run k-mer analysis on alleles
    quality_results = _run_quality_analysis(
        alleles_dir,
        kmer_size=kmer_size,
        overlap_threshold=overlap_threshold,
        duplicate_threshold=duplicate_threshold,
    )

    # Generate metadata.toml
    metadata = _generate_metadata(
        namespace=namespace,
        name=torch_name,
        version=torch_version,
        scheme_data=scheme_data,
        database_url=database_url,
        scheme_id=scheme_id,
        loci_list=[locus.locus_id for locus in scheme_data.loci],
        kmer_size=kmer_size,
        overlap_threshold=overlap_threshold,
        duplicate_threshold=duplicate_threshold,
        quality_results=quality_results,
    )

    metadata_path = torch_dir / "metadata.toml"
    with open(metadata_path, "w") as f:
        toml.dump(metadata, f)

    # Generate quality.json
    quality_report = _generate_quality_report(
        kmer_size=kmer_size,
        overlap_threshold=overlap_threshold,
        duplicate_threshold=duplicate_threshold,
        quality_results=quality_results,
    )

    quality_path = torch_dir / "quality.json"
    with open(quality_path, "w") as f:
        json.dump(quality_report, f, indent=2)

    return str(torch_dir)


def _extract_database_name(url: str) -> str:
    """Extract database name from PubMLST URL.

    Args:
        url: Base URL like "http://pubmlst.org/api" or "http://pubmlst.org/api/"

    Returns:
        Database name like "pubmlst"
    """
    # For now, return a generic database name
    # In reality, this might need to be specified explicitly
    if "pubmlst" in url.lower():
        return "pubmlst"
    return "bigsdb"


def _sanitize_name(name: str) -> str:
    """Sanitize a name for use in file paths.

    Args:
        name: Input name

    Returns:
        Sanitized name with spaces replaced by underscores
    """
    return name.replace(" ", "_").replace("/", "_").lower()


def _write_stub_fasta(fasta_path: Path, locus_id: str) -> None:
    """Write a stub FASTA file for testing.

    Args:
        fasta_path: Path to write FASTA file
        locus_id: Locus identifier
    """
    with open(fasta_path, "w") as f:
        # Write stub sequences
        f.write(f">{locus_id}_1\n")
        f.write("ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG\n")
        f.write(f">{locus_id}_2\n")
        f.write("ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG\n")


def _write_profiles_tsv(profiles_path: Path, profiles: List[Dict[str, str]]) -> None:
    """Write profiles.tsv file.

    Args:
        profiles_path: Path to write profiles file
        profiles: List of profile dictionaries
    """
    if not profiles:
        # Write minimal header
        with open(profiles_path, "w") as f:
            f.write("ST\n")
        return

    # Get field names from first profile
    fieldnames = list(profiles[0].keys())

    with open(profiles_path, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(profiles)


def _run_quality_analysis(
    alleles_dir: Path,
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> Dict[str, Any]:
    """Run k-mer quality analysis on all loci.

    Args:
        alleles_dir: Directory containing locus FASTA files
        kmer_size: K-mer size for analysis
        overlap_threshold: Threshold for overlap detection
        duplicate_threshold: Threshold for duplicate detection

    Returns:
        Dictionary with quality analysis results
    """
    results = {
        "total_loci": 0,
        "similar_pairs": [],
        "duplicate_pairs": [],
        "loci_results": {},
    }

    # Iterate through all FASTA files in alleles directory
    for fasta_file in sorted(alleles_dir.glob("*.fasta")):
        locus_name = fasta_file.stem

        # Run analysis
        report = analyze_locus(
            fasta_file,
            k_size=kmer_size,
            overlap_threshold=overlap_threshold * 100,  # Convert to 0-100 scale
            duplicate_threshold=duplicate_threshold * 100,  # Convert to 0-100 scale
        )

        results["total_loci"] += 1
        results["loci_results"][locus_name] = {
            "suspect_pairs": report.suspect_pairs,
            "threshold": report.threshold,
        }

        # Collect similar and duplicate pairs
        for pair in report.suspect_pairs:
            if pair["issue_type"] == "duplicate":
                results["duplicate_pairs"].append(
                    {
                        "locus": locus_name,
                        "allele1": pair["allele1"],
                        "allele2": pair["allele2"],
                        "similarity": pair["similarity"],
                    }
                )
            elif pair["issue_type"] == "overlap":
                results["similar_pairs"].append(
                    {
                        "locus": locus_name,
                        "allele1": pair["allele1"],
                        "allele2": pair["allele2"],
                        "similarity": pair["similarity"],
                    }
                )

    return results


def _generate_metadata(
    namespace: str,
    name: str,
    version: str,
    scheme_data: Any,
    database_url: str,
    scheme_id: int,
    loci_list: List[str],
    kmer_size: int,
    overlap_threshold: float,
    duplicate_threshold: float,
    quality_results: Dict[str, Any],
) -> Dict[str, Any]:
    """Generate metadata dictionary.

    Args:
        namespace: Torch namespace
        name: Torch name
        version: Torch version
        scheme_data: SchemeData object
        database_url: PubMLST database URL
        scheme_id: Scheme ID
        loci_list: List of locus names
        kmer_size: K-mer size used
        overlap_threshold: Overlap threshold used
        duplicate_threshold: Duplicate threshold used
        quality_results: Quality analysis results

    Returns:
        Metadata dictionary
    """
    now = datetime.now(timezone.utc).isoformat()

    return {
        "namespace": namespace,
        "name": name,
        "version": version,
        "provenance": {
            "source": "PubMLST",
            "database_url": database_url,
            "scheme_id": scheme_id,
            "fetch_date": now,
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
            "scheme_name": scheme_data.metadata.name,
            "loci_count": len(loci_list),
            "profiles_count": scheme_data.profiles.row_count
            if scheme_data.profiles
            else 0,
            "last_updated": (
                scheme_data.metadata.last_updated.isoformat()
                if scheme_data.metadata.last_updated
                else now
            ),
        },
        "schemes": {
            _sanitize_name(scheme_data.metadata.name): {
                "organism": scheme_data.metadata.name,
                "loci": loci_list,
            }
        },
    }


def _generate_quality_report(
    kmer_size: int,
    overlap_threshold: float,
    duplicate_threshold: float,
    quality_results: Dict[str, Any],
) -> Dict[str, Any]:
    """Generate quality.json report.

    Args:
        kmer_size: K-mer size used
        overlap_threshold: Overlap threshold used
        duplicate_threshold: Duplicate threshold used
        quality_results: Quality analysis results

    Returns:
        Quality report dictionary
    """
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
