version 1.0

task run_alignment {
    input {
        File query_sequences
        File allele_fasta
        String preset = "asm5"
        Float confidence_threshold = 0.95
    }

    command <<<
        set -e
        # Install minimap2
        apt-get update && apt-get install -y minimap2 && rm -rf /var/lib/apt/lists/*

        python3 <<'PYTHON_SCRIPT'
import json
import subprocess
from collections import defaultdict

def parse_fasta(fasta_path):
    """Parse FASTA file into dictionary of header -> sequence."""
    sequences = {}
    with open(fasta_path) as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)
    return sequences

def extract_locus_and_allele(header):
    """Extract locus name and allele ID from FASTA header."""
    parts = header.split('_')
    if len(parts) >= 2:
        allele_id = parts[-1]
        locus = '_'.join(parts[:-1])
        return locus, allele_id
    return header, "unknown"

# Parse inputs
query_seqs = parse_fasta("~{query_sequences}")
allele_seqs = parse_fasta("~{allele_fasta}")

# Group alleles by locus
alleles_by_locus = defaultdict(list)
for header, seq in allele_seqs.items():
    locus, allele_id = extract_locus_and_allele(header)
    alleles_by_locus[locus].append({
        'allele_id': allele_id,
        'header': header,
        'sequence': seq
    })

# Run alignment for each query sequence
alignment_results = {}

for query_name, query_seq in query_seqs.items():
    # Write query to temporary file
    query_file = f"{query_name}_query.fasta"
    with open(query_file, 'w') as f:
        f.write(f">{query_name}\n{query_seq}\n")

    # Run minimap2 with specified preset
    sam_file = f"{query_name}_alignment.sam"
    cmd = [
        'minimap2',
        '-a',
        '-x', '~{preset}',
        '~{allele_fasta}',
        query_file,
        '-o', sam_file
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Warning: minimap2 failed for {query_name}")
        alignment_results[query_name] = {}
        continue

    # Parse SAM file to extract alignment metrics
    alignments_by_ref = {}

    try:
        with open(sam_file) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue

                ref_name = fields[2]
                mapq = int(fields[4])

                if mapq == 0 or ref_name == '*':  # Skip unmapped
                    continue

                # Extract alignment stats
                query_start = int(fields[3])
                match_count = 0
                mismatch_count = 0

                # Look for NM tag (number of mismatches)
                nm = 0
                for tag in fields[11:]:
                    if tag.startswith('NM:i:'):
                        nm = int(tag.split(':')[-1])
                        break

                # Calculate identity
                query_len = len(query_seq)
                identity = max(0.0, (query_len - nm) / query_len)

                if ref_name not in alignments_by_ref:
                    alignments_by_ref[ref_name] = {
                        'identity': identity,
                        'mapq': mapq,
                        'mismatches': nm
                    }
                else:
                    # Keep best match
                    if identity > alignments_by_ref[ref_name]['identity']:
                        alignments_by_ref[ref_name] = {
                            'identity': identity,
                            'mapq': mapq,
                            'mismatches': nm
                        }
    except Exception as e:
        print(f"Warning: Error parsing SAM for {query_name}: {e}")

    alignment_results[query_name] = alignments_by_ref

# Convert alignments to allele calls per locus
allele_calls = {}

for query_name, alignments in alignment_results.items():
    # Group alignments by locus
    by_locus = defaultdict(list)
    for ref_header, metrics in alignments.items():
        locus, allele_id = extract_locus_and_allele(ref_header)
        by_locus[locus].append({
            'allele_id': allele_id,
            'identity': metrics['identity'],
            'mapq': metrics['mapq'],
            'mismatches': metrics['mismatches']
        })

    # Select best match per locus
    for locus, candidates in by_locus.items():
        if candidates:
            best = max(candidates, key=lambda x: x['identity'])
            allele_calls[locus] = {
                'allele_id': best['allele_id'],
                'identity': float(best['identity']),
                'mapq': int(best['mapq']),
                'mismatches': int(best['mismatches']),
                'confidence': best['identity'] >= ~{confidence_threshold}
            }

# Write JSON output
with open('alignment_calls.json', 'w') as f:
    json.dump(allele_calls, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File alignment_calls = "alignment_calls.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 2
        memory: "4 GB"
    }
}
