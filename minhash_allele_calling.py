#!/usr/bin/env python3
"""
MinHash allele calling using sourmash.

Takes query sequences and allele database (FASTA format) and produces
allele calls with MinHash-based similarity scores.
"""

import json
import sys
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict


def parse_fasta(fasta_path):
    """Parse FASTA file and return list of (header, sequence) tuples."""
    sequences = []
    with open(fasta_path) as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    return sequences


def extract_locus_and_allele(header):
    """Extract locus and allele from FASTA header.

    Expected format: locus_allele or similar.
    Examples: adk_1, fumC_2, gyrB_1
    """
    parts = header.split('_')
    if len(parts) >= 2:
        # Last part should be allele number
        allele_id = parts[-1]
        locus = '_'.join(parts[:-1])
        return locus, allele_id
    return header, "unknown"


def compute_sequence_similarity(seq1, seq2):
    """Compute sequence similarity as normalized matching k-mers.

    For short sequences, use exact matching ratio.
    For longer sequences, use k-mer overlap (k=31).
    """
    if len(seq1) == 0 or len(seq2) == 0:
        return 0.0

    # For exact match or near-exact sequences
    if seq1 == seq2:
        return 1.0

    # Compute matching positions (simple percentage identity)
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    simple_identity = matches / min_len if min_len > 0 else 0.0

    # Use k-mer based approach for better similarity
    ksize = 31
    if len(seq1) >= ksize and len(seq2) >= ksize:
        # Extract k-mers
        kmers1 = set(seq1[i:i+ksize]
                     for i in range(len(seq1) - ksize + 1))
        kmers2 = set(seq2[i:i+ksize]
                     for i in range(len(seq2) - ksize + 1))

        if len(kmers1) == 0 or len(kmers2) == 0:
            return simple_identity

        # Jaccard similarity
        intersection = len(kmers1 & kmers2)
        union = len(kmers1 | kmers2)
        jaccard_sim = intersection / union if union > 0 else 0.0

        # Return average of simple identity and Jaccard similarity
        return (simple_identity + jaccard_sim) / 2.0
    else:
        # For short sequences, use simple identity
        return simple_identity


def run_sourmash_sketch(input_file, ksize=31, scaled=1000):
    """Create a MinHash sketch for the input file using sourmash.

    Returns the path to the sketch file on success, or None on failure.
    """
    sketch_path = f"{input_file}.k{ksize}.sig"
    try:
        cmd = [
            "sourmash", "sketch", "dna",
            "-p", f"k={ksize},scaled={scaled}",
            "-o", sketch_path,
            input_file
        ]
        result = subprocess.run(cmd, check=False, capture_output=True, timeout=30)
        if result.returncode == 0 and Path(sketch_path).exists():
            return sketch_path
    except Exception:
        pass
    return None


def run_sourmash_compare(query_sketch, allele_sketch):
    """Compare two sketches using sourmash and return similarity score.

    Returns similarity as float in [0, 1], or None on failure.
    """
    try:
        csv_path = "/tmp/comparison.csv"
        cmd = [
            "sourmash", "compare",
            query_sketch,
            allele_sketch,
            "--csv", csv_path
        ]
        result = subprocess.run(cmd, check=False, capture_output=True, timeout=30)
        if result.returncode == 0 and Path(csv_path).exists():
            with open(csv_path) as f:
                lines = f.readlines()
                if len(lines) > 1:
                    # Parse CSV: first row is headers, second row contains values
                    parts = lines[1].split(',')
                    if len(parts) > 1:
                        try:
                            sim = float(parts[1])
                            return max(0.0, min(1.0, sim))  # Clamp to [0, 1]
                        except (ValueError, IndexError):
                            pass
    except Exception:
        pass
    return None


def main():
    if len(sys.argv) < 3:
        usage = ("Usage: minhash_allele_calling.py <query_fasta> "
                 "<allele_fasta> [reads_fasta] [min_coverage]")
        print(usage, file=sys.stderr)
        sys.exit(1)

    query_file = sys.argv[1]
    allele_file = sys.argv[2]
    # reads_file and min_coverage available for future depth filtering
    # reads_file = sys.argv[3] if len(sys.argv) > 3 else None
    # min_coverage = int(sys.argv[4]) if len(sys.argv) > 4 else 3

    # Parse query sequences
    query_sequences = parse_fasta(query_file)
    if not query_sequences:
        print("{}", file=sys.stdout)
        sys.exit(0)

    # Parse allele database and group by locus
    allele_sequences = parse_fasta(allele_file)
    alleles_by_locus = defaultdict(list)
    for header, seq in allele_sequences:
        locus, allele_id = extract_locus_and_allele(header)
        alleles_by_locus[locus].append({
            'allele_id': allele_id,
            'header': header,
            'sequence': seq
        })

    # Results dictionary
    results = {}

    # For each locus, find best match across all query sequences
    for locus, alleles in alleles_by_locus.items():
        best_match = None
        best_similarity = -1.0

        for query_header, query_seq in query_sequences:
            for allele in alleles:
                allele_seq = allele['sequence']

                # Try sourmash first if available
                similarity = None
                try:
                    # Create temporary FASTA files
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as qf:
                        qf.write(f">{query_header}\n{query_seq}\n")
                        query_temp = qf.name

                    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as af:
                        af.write(f">{allele['header']}\n{allele_seq}\n")
                        allele_temp = af.name

                    query_sketch = run_sourmash_sketch(query_temp)
                    allele_sketch = run_sourmash_sketch(allele_temp)

                    if query_sketch and allele_sketch:
                        similarity = run_sourmash_compare(query_sketch, allele_sketch)

                    # Cleanup
                    Path(query_temp).unlink(missing_ok=True)
                    Path(allele_temp).unlink(missing_ok=True)
                    if query_sketch:
                        Path(query_sketch).unlink(missing_ok=True)
                    if allele_sketch:
                        Path(allele_sketch).unlink(missing_ok=True)
                except Exception:
                    similarity = None

                # Fallback to sequence-based similarity
                if similarity is None:
                    similarity = compute_sequence_similarity(query_seq, allele_seq)

                # Track best match
                if similarity > best_similarity:
                    best_similarity = similarity
                    best_match = allele['allele_id']

        # Store result for this locus
        if best_match is not None:
            results[locus] = {
                'allele_id': best_match,
                'similarity': max(0.0, min(1.0, best_similarity)),
                'confidence': best_similarity > 0.9
            }

    # Output JSON
    json.dump(results, sys.stdout, indent=2)


if __name__ == '__main__':
    main()
