version 1.0

task align_and_call {
    input {
        File query_sequences
        File allele_fasta
        String input_type = "contigs"
        Float identity_threshold = 0.90
    }

    command <<<
        set -e
        # Install minimap2 if needed
        which minimap2 > /dev/null || (apt-get update && apt-get install -y minimap2 && rm -rf /var/lib/apt/lists/*)

        python3 <<'PYTHON_SCRIPT'
import json
import subprocess
import tempfile
from collections import defaultdict
import os

def parse_fasta(fasta_path):
    """Parse FASTA file into dict of {header: sequence}."""
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

def run_minimap2_alignment(query_sequences, allele_fasta, preset):
    """Run minimap2 to align queries against allele database."""
    # Write queries to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        query_file = f.name
        for header, seq in query_sequences.items():
            f.write(f">{header}\n{seq}\n")

    # Run minimap2
    sam_file = "alignment.sam"
    cmd = [
        'minimap2',
        '-a',
        '-x', preset,
        '--secondary=no',
        allele_fasta,
        query_file,
        '-o', sam_file
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Warning: minimap2 failed: {e}")
        return {}
    finally:
        os.unlink(query_file)

    # Parse SAM file to extract best alignment per query
    results = {}
    query_to_best = {}

    try:
        with open(sam_file) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue

                query_name = fields[0]
                ref_name = fields[2]
                mapq = int(fields[4])

                if mapq == 0 or ref_name == '*':  # Skip unmapped
                    continue

                # Extract identity from NM tag (number of mismatches)
                nm = 0
                for tag in fields[11:]:
                    if tag.startswith('NM:i:'):
                        nm = int(tag.split(':')[-1])
                        break

                # Calculate identity
                query_len = int(fields[8]) if int(fields[8]) > 0 else len(query_sequences.get(query_name, ""))
                identity = max(0.0, (query_len - nm) / query_len) if query_len > 0 else 0.0

                # Track best alignment per query
                if query_name not in query_to_best or identity > query_to_best[query_name]['identity']:
                    locus, allele_id = extract_locus_and_allele(ref_name)
                    query_to_best[query_name] = {
                        'allele_id': allele_id,
                        'identity': identity,
                        'locus': locus
                    }
    except Exception as e:
        print(f"Warning: Error parsing SAM: {e}")

    # Aggregate by locus (take best identity per locus across all queries)
    locus_results = defaultdict(lambda: {'allele_id': None, 'identity': 0.0})
    for query_name, best in query_to_best.items():
        locus = best['locus']
        if best['identity'] > locus_results[locus]['identity']:
            locus_results[locus] = {
                'allele_id': best['allele_id'],
                'identity': best['identity']
            }

    return dict(locus_results)

# Main
query_seqs = parse_fasta("~{query_sequences}")
input_type = "~{input_type}"
preset = "sr" if input_type == "reads" else "asm5"

alignment_results = run_minimap2_alignment(query_seqs, "~{allele_fasta}", preset)

# Format output with confidence
output_results = {}
for locus, result in alignment_results.items():
    identity = result['identity']
    output_results[locus] = {
        'allele_id': result['allele_id'],
        'identity': max(0.0, min(1.0, identity)),
        'confidence': identity >= ~{identity_threshold}
    }

# Write JSON output
with open('alignment_results.json', 'w') as f:
    json.dump(output_results, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File alignment_results = "alignment_results.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 2
        memory: "4 GB"
    }
}
