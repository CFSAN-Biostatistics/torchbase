"""Quality report generation from SimilarityReport into formatted outputs.

Supports JSON, human-readable text with ASCII histograms, or both formats.
Generates hierarchical flagging data for metadata.
"""

import json
from typing import Dict, Union, Any, List, Tuple
from dataclasses import dataclass, field

from torchbase.quality.kmer_analysis import SimilarityReport


@dataclass
class QualitySummary:
    """Summary of quality report across all loci."""

    total_loci: int = 0
    total_suspect_allele_pairs: int = 0
    suspect_loci: List[str] = field(default_factory=list)
    suspect_profiles: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert summary to dictionary."""
        return {
            "total_loci": self.total_loci,
            "total_suspect_allele_pairs": self.total_suspect_allele_pairs,
            "suspect_loci": self.suspect_loci,
            "suspect_profiles": self.suspect_profiles,
        }


def _get_suspect_alleles(report: SimilarityReport) -> set:
    """Extract all alleles involved in suspect pairs."""
    suspect_alleles = set()
    for pair in report.suspect_pairs:
        suspect_alleles.add(pair["allele1"])
        suspect_alleles.add(pair["allele2"])
    return suspect_alleles


def _build_histogram(
    similarities: Dict[Tuple[str, str], float], num_bins: int = 10
) -> Tuple[List[int], float, float]:
    """Build histogram data for similarity distribution.

    Args:
        similarities: Dictionary of (allele1, allele2) -> similarity pairs
        num_bins: Number of bins for histogram (default 10)

    Returns:
        Tuple of (bin_counts, min_val, max_val)
    """
    if not similarities:
        return [0] * num_bins, 0.0, 0.0

    values = list(similarities.values())
    if not values:
        return [0] * num_bins, 0.0, 0.0

    min_val = min(values)
    max_val = max(values)

    # Handle case where all values are the same
    if min_val == max_val:
        bins = [0] * num_bins
        bins[0] = len(values)
        return bins, min_val, max_val

    # Create bins
    bins = [0] * num_bins
    bin_width = (max_val - min_val) / num_bins

    for val in values:
        # Calculate which bin this value belongs to
        bin_idx = int((val - min_val) / bin_width)
        # Handle edge case where val == max_val
        if bin_idx >= num_bins:
            bin_idx = num_bins - 1
        bins[bin_idx] += 1

    return bins, min_val, max_val


def _render_histogram_text(
    bins: List[int], min_val: float, max_val: float, width: int = 40
) -> str:
    """Render histogram as ASCII bar chart.

    Args:
        bins: List of bin counts
        min_val: Minimum similarity value
        max_val: Maximum similarity value
        width: Width of the histogram in characters

    Returns:
        ASCII histogram string
    """
    if not bins or all(b == 0 for b in bins):
        return "  (No similarity data)"

    max_count = max(bins) if bins else 1
    if max_count == 0:
        max_count = 1

    lines = []
    num_bins = len(bins)

    # Header
    lines.append(f"  Similarity Distribution ({min_val:.1f} - {max_val:.1f}):")
    lines.append("")

    # Bars
    for i, count in enumerate(bins):
        bin_start = min_val + (i / num_bins) * (max_val - min_val)
        bin_end = min_val + ((i + 1) / num_bins) * (max_val - min_val)
        bar_width = int((count / max_count) * width) if count > 0 else 0
        bar = "|" * bar_width
        lines.append(f"  [{bin_start:5.1f}-{bin_end:5.1f}] {bar} ({count})")

    return "\n".join(lines)


def _format_json_output(
    locus_reports: Dict[str, SimilarityReport],
) -> Dict[str, Any]:
    """Generate JSON-formatted quality report.

    Args:
        locus_reports: Dictionary of locus_name -> SimilarityReport

    Returns:
        Dictionary representing the JSON structure
    """
    # Build the loci section
    loci_data = {}
    for locus_name, report in locus_reports.items():
        loci_data[locus_name] = {
            "similarities": {
                f"{pair[0]}-{pair[1]}": sim
                for pair, sim in report.similarities.items()
            },
            "threshold": report.threshold,
            "statistics": report.statistics,
        }

    # Build the suspect_pairs section
    suspect_pairs_data = {}
    for locus_name, report in locus_reports.items():
        suspect_pairs_data[locus_name] = report.suspect_pairs

    # Build summary
    summary = _build_summary(locus_reports)

    return {
        "loci": loci_data,
        "similarity_stats": {
            locus_name: report.statistics
            for locus_name, report in locus_reports.items()
        },
        "suspect_pairs": suspect_pairs_data,
        "summary": summary.to_dict(),
    }


def _format_text_output(locus_reports: Dict[str, SimilarityReport]) -> str:
    """Generate human-readable text report with ASCII histograms.

    Args:
        locus_reports: Dictionary of locus_name -> SimilarityReport

    Returns:
        Formatted text string
    """
    lines = []
    lines.append("=" * 80)
    lines.append("QUALITY REPORT: Allele Similarity Analysis")
    lines.append("=" * 80)
    lines.append("")

    # Report for each locus
    for locus_name in sorted(locus_reports.keys()):
        report = locus_reports[locus_name]
        lines.append(f"Locus: {locus_name}")
        lines.append("-" * 80)
        lines.append("")

        # Statistics section
        stats = report.statistics
        lines.append("  Statistics:")
        mean = stats.get("mean", 0.0)
        lines.append(f"    Mean similarity:       {mean:.2f}%")
        std_dev = stats.get("std_dev", 0.0)
        lines.append(f"    Std deviation:        {std_dev:.2f}%")
        min_sim = stats.get("min", 0.0)
        lines.append(f"    Min:                  {min_sim:.2f}%")
        max_sim = stats.get("max", 0.0)
        lines.append(f"    Max:                  {max_sim:.2f}%")
        perc_99 = stats.get("percentile_99", 0.0)
        lines.append(f"    99th percentile:      {perc_99:.2f}%")
        lines.append(f"    Threshold:            {report.threshold:.2f}%")
        threshold_type = stats.get("threshold_type", "unknown")
        lines.append(f"    Detection method:     {threshold_type}")
        lines.append("")

        # Histogram
        bins, min_val, max_val = _build_histogram(
            report.similarities, num_bins=10
        )
        histogram_text = _render_histogram_text(
            bins, min_val, max_val, width=40
        )
        lines.append(histogram_text)
        lines.append("")

        # Suspect pairs section
        if report.suspect_pairs:
            lines.append("  Suspect Pairs:")
            for i, pair in enumerate(report.suspect_pairs, 1):
                allele1 = pair["allele1"]
                allele2 = pair["allele2"]
                lines.append(f"    {i}. {allele1} <-> {allele2}")
                sim = pair["similarity"]
                lines.append(f"       Similarity:      {sim:.2f}%")
                cont_1_2 = pair["containment_1_in_2"]
                lines.append(f"       Containment 1→2: {cont_1_2:.2f}%")
                cont_2_1 = pair["containment_2_in_1"]
                lines.append(f"       Containment 2→1: {cont_2_1:.2f}%")
                issue = pair["issue_type"]
                lines.append(f"       Issue type:      {issue}")
            lines.append("")
        else:
            lines.append("  No suspect pairs detected.")
            lines.append("")

    # Summary section
    summary = _build_summary(locus_reports)
    lines.append("=" * 80)
    lines.append("SUMMARY")
    lines.append("=" * 80)
    total = summary.total_loci
    lines.append(f"Total loci analyzed:        {total}")
    total_pairs = summary.total_suspect_allele_pairs
    lines.append(f"Total suspect allele pairs: {total_pairs}")
    suspect_count = len(summary.suspect_loci)
    lines.append(f"Loci with suspects:         {suspect_count}")
    lines.append("")

    if summary.suspect_loci:
        lines.append("Suspect loci:")
        for locus in sorted(summary.suspect_loci):
            lines.append(f"  - {locus}")
        lines.append("")

    if summary.suspect_profiles:
        msg = "Suspect profiles (loci needing review):"
        lines.append(msg)
        for profile in sorted(summary.suspect_profiles):
            lines.append(f"  - {profile}")
        lines.append("")

    lines.append("=" * 80)

    return "\n".join(lines)


def _build_summary(
    locus_reports: Dict[str, SimilarityReport],
) -> QualitySummary:
    """Build summary from locus reports.

    Args:
        locus_reports: Dictionary of locus_name -> SimilarityReport

    Returns:
        QualitySummary object
    """
    summary = QualitySummary(total_loci=len(locus_reports))

    suspect_loci = []
    total_suspect_pairs = 0

    for locus_name, report in locus_reports.items():
        if report.suspect_pairs:
            suspect_loci.append(locus_name)
            total_suspect_pairs += len(report.suspect_pairs)

    summary.suspect_loci = sorted(suspect_loci)
    # Hierarchical flagging: suspect loci propagate to profiles
    summary.suspect_profiles = sorted(suspect_loci)
    summary.total_suspect_allele_pairs = total_suspect_pairs

    return summary


def generate_report(
    locus_reports: Dict[str, SimilarityReport], format: str = "text"
) -> Union[str, Dict[str, str]]:
    """Generate quality report from SimilarityReport data.

    Transforms similarity reports into formatted outputs: JSON, text with
    ASCII histograms, or both. Generates hierarchical flagging data
    for metadata.toml.

    Args:
        locus_reports: Dictionary mapping locus names to SimilarityReport
        format: Output format: 'text' (default), 'json', or 'both'.

    Returns:
        - If format='text': Human-readable string with ASCII histograms
        - If format='json': JSON string with complete report structure
        - If format='both': Dictionary with 'json' and 'text' keys
    """
    if format not in ("text", "json", "both"):
        msg = f"Invalid format: {format}. Must be 'text', 'json', or 'both'."
        raise ValueError(msg)

    json_output = None
    text_output = None

    if format in ("json", "both"):
        json_dict = _format_json_output(locus_reports)
        json_output = json.dumps(json_dict, indent=2)

    if format in ("text", "both"):
        text_output = _format_text_output(locus_reports)

    if format == "text":
        return text_output
    elif format == "json":
        return json_output
    else:  # format == "both"
        return {
            "json": json_output,
            "text": text_output,
        }
