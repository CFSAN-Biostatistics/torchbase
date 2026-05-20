"""K-mer quality analysis module for detecting allele issues in loci."""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
import statistics


@dataclass
class SimilarityReport:
    """Report of k-mer similarity analysis for alleles in a locus."""

    similarities: Dict[Tuple[str, str], float]
    """Pairwise k-mer similarities (0-100 scale)."""

    threshold: float
    """Computed threshold for suspect pair detection."""

    suspect_pairs: List[Dict]
    """List of suspect pairs flagged as overlaps or duplicates."""

    statistics: Dict
    """Statistics about the similarity distribution."""

    def __repr__(self) -> str:
        threshold_type = self.statistics.get('threshold_type',
                                             'unknown')
        return (
            f"SimilarityReport(similarities={len(self.similarities)} "
            f"pairs, "
            f"threshold={self.threshold:.2f}, "
            f"suspect_pairs={len(self.suspect_pairs)}, "
            f"threshold_type={threshold_type})"
        )


def _parse_fasta(fasta_path):
    """Parse a FASTA file and return {sequence_id: sequence}."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]  # Take first word as ID
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)

    return sequences


def _get_kmers(sequence, k):
    """Extract k-mers from a sequence as a set."""
    if len(sequence) < k:
        return set()
    kmers = set()
    for i in range(len(sequence) - k + 1):
        kmers.add(sequence[i:i + k])
    return kmers


def _jaccard_similarity(kmers1, kmers2):
    """Calculate Jaccard similarity between two k-mer sets (0-100 scale)."""
    if not kmers1 and not kmers2:
        return 100.0
    if not kmers1 or not kmers2:
        return 0.0

    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    return (intersection / union * 100) if union > 0 else 0.0


def _containment_similarity(kmers1, kmers2):
    """Calculate one-way containment similarities.

    Returns (containment_1_in_2, containment_2_in_1) as percentages.
    """
    if not kmers1 or not kmers2:
        return 0.0, 0.0

    intersection = len(kmers1 & kmers2)
    containment_1_in_2 = (intersection / len(kmers2) * 100
                          if kmers2 else 0.0)
    containment_2_in_1 = (intersection / len(kmers1) * 100
                          if kmers1 else 0.0)

    return containment_1_in_2, containment_2_in_1


def _detect_gap(sorted_similarities):
    """Detect outliers via gap detection (>2% jumps in sorted similarities).

    Returns (threshold, found_gap) tuple.
    """
    if len(sorted_similarities) < 2:
        return None, False

    max_gap = 0.0
    gap_index = -1

    for i in range(len(sorted_similarities) - 1):
        gap = sorted_similarities[i + 1] - sorted_similarities[i]
        if gap > max_gap:
            max_gap = gap
            gap_index = i

    # Consider it a significant gap if > 2%
    if max_gap > 2.0:
        # Threshold is just above the lower group
        threshold = sorted_similarities[gap_index] + (max_gap / 2)
        return threshold, True

    return None, False


def _get_percentile(values, percentile):
    """Calculate percentile of values."""
    if not values:
        return 0.0
    sorted_vals = sorted(values)
    index = (percentile / 100.0) * (len(sorted_vals) - 1)
    lower_idx = int(index)
    upper_idx = min(lower_idx + 1, len(sorted_vals) - 1)

    if lower_idx == upper_idx:
        return sorted_vals[lower_idx]

    # Linear interpolation
    frac = index - lower_idx
    return sorted_vals[lower_idx] * (1 - frac) + sorted_vals[upper_idx] * frac


def analyze_locus(
        fasta_path,
        k_size=21,
        overlap_threshold=95.0,
        duplicate_threshold=98.0,
):
    """Analyze k-mer similarities within alleles of a locus.

    Auto-tuning detects outliers via gap detection (>2% jumps) or falls back
    to 99th percentile. Flags overlaps (95% one-way containment) and
    duplicates (98% symmetric similarity).

    Args:
        fasta_path: Path to FASTA file with allele sequences.
        k_size: K-mer size (default 21).
        overlap_threshold: Containment threshold for overlap
            detection (default 95).
        duplicate_threshold: Symmetric similarity threshold for
            duplicate detection (default 98).

    Returns:
        SimilarityReport with similarities, threshold, suspect
        pairs, and statistics.
    """
    fasta_path = Path(fasta_path)

    # Parse FASTA
    sequences = _parse_fasta(fasta_path)

    # Handle edge cases
    if len(sequences) < 2:
        return SimilarityReport(
            similarities={},
            threshold=0.0,
            suspect_pairs=[],
            statistics={
                "mean": 0.0,
                "std_dev": 0.0,
                "min": 0.0,
                "max": 0.0,
                "percentile_99": 0.0,
                "threshold_type": "none",
            },
        )

    # Extract k-mers for all sequences
    kmer_sets = {}
    for allele_id, seq in sequences.items():
        kmer_sets[allele_id] = _get_kmers(seq, k_size)

    # Compute pairwise Jaccard similarities
    similarities = {}
    allele_ids = sorted(kmer_sets.keys())

    for i, allele1 in enumerate(allele_ids):
        for allele2 in allele_ids[i + 1:]:
            kmers1 = kmer_sets[allele1]
            kmers2 = kmer_sets[allele2]

            # Use Jaccard similarity
            sim = _jaccard_similarity(kmers1, kmers2)
            similarities[(allele1, allele2)] = sim

    # Calculate statistics
    if not similarities:
        return SimilarityReport(
            similarities={},
            threshold=0.0,
            suspect_pairs=[],
            statistics={
                "mean": 0.0,
                "std_dev": 0.0,
                "min": 0.0,
                "max": 0.0,
                "percentile_99": 0.0,
                "threshold_type": "none",
            },
        )

    similarity_values = list(similarities.values())
    sorted_similarities = sorted(similarity_values)

    mean_sim = statistics.mean(similarity_values)
    if len(similarity_values) > 1:
        std_dev = statistics.stdev(similarity_values)
    else:
        std_dev = 0.0
    min_sim = min(similarity_values)
    max_sim = max(similarity_values)
    percentile_99 = _get_percentile(similarity_values, 99)

    # Auto-tune: try gap detection first
    gap_threshold, found_gap = _detect_gap(sorted_similarities)
    if found_gap:
        threshold_type = "gap_detection"
        threshold = gap_threshold
    else:
        threshold_type = "percentile"
        threshold = percentile_99

    # Identify suspect pairs
    suspect_pairs = []

    for (allele1, allele2), sim in similarities.items():
        kmers1 = kmer_sets[allele1]
        kmers2 = kmer_sets[allele2]

        # Check for overlaps (one-way containment)
        cont_1_in_2, cont_2_in_1 = _containment_similarity(kmers1,
                                                           kmers2)
        is_overlap = (cont_1_in_2 >= overlap_threshold
                      or cont_2_in_1 >= overlap_threshold)

        # Check for duplicates (symmetric similarity)
        is_duplicate = sim >= duplicate_threshold

        if is_overlap or is_duplicate:
            issue_type = "duplicate" if is_duplicate else "overlap"
            suspect_pairs.append(
                {
                    "allele1": allele1,
                    "allele2": allele2,
                    "similarity": sim,
                    "containment_1_in_2": cont_1_in_2,
                    "containment_2_in_1": cont_2_in_1,
                    "issue_type": issue_type,
                }
            )

    # Build statistics dictionary
    stats = {
        "mean": mean_sim,
        "std_dev": std_dev,
        "min": min_sim,
        "max": max_sim,
        "percentile_99": percentile_99,
        "threshold_type": threshold_type,
    }

    return SimilarityReport(
        similarities=similarities,
        threshold=threshold,
        suspect_pairs=suspect_pairs,
        statistics=stats,
    )
