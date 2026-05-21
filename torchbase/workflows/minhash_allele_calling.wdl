version 1.0

workflow minhash_allele_calling {
    input {
        File query_sequences
        File allele_fasta
        File? reads_file
        Int min_coverage = 3
        Int ksize = 31
        Int sketch_size = 1000
    }

    call minhash_allele_calling_task {
        input:
            query_sequences = query_sequences,
            allele_fasta = allele_fasta,
            reads_file = reads_file,
            min_coverage = min_coverage,
            ksize = ksize,
            sketch_size = sketch_size
    }

    output {
        File results = minhash_allele_calling_task.results
    }
}

task minhash_allele_calling_task {
    input {
        File query_sequences
        File allele_fasta
        File? reads_file
        Int min_coverage = 3
        Int ksize = 31
        Int sketch_size = 1000
    }

    command <<<
        set -e

        # Create output directory
        mkdir -p results

        # Install sourmash if not available
        pip install -q sourmash 2>/dev/null || true

        # Create a Python script to handle the allele calling
        cat > /tmp/minhash_allele_calling.py << 'PYTHON_EOF'
#!/usr/bin/env python3
import json
import sys
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

def parse_fasta(fasta_path):
    sequences = []
    with open(fasta_path) as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    return sequences

def extract_locus_and_allele(header):
    parts = header.split('_')
    if len(parts) >= 2:
        allele_id = parts[-1]
        locus = '_'.join(parts[:-1])
        return locus, allele_id
    return header, "unknown"

def compute_sequence_similarity(seq1, seq2):
    if len(seq1) == 0 or len(seq2) == 0:
        return 0.0
    if seq1 == seq2:
        return 1.0
    min_len = min(len(seq1), len(seq2))
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    simple_identity = matches / min_len if min_len > 0 else 0.0
    ksize = 31
    if len(seq1) >= ksize and len(seq2) >= ksize:
        kmers1 = set(seq1[i:i+ksize] for i in range(len(seq1) - ksize + 1))
        kmers2 = set(seq2[i:i+ksize] for i in range(len(seq2) - ksize + 1))
        if len(kmers1) == 0 or len(kmers2) == 0:
            return simple_identity
        intersection = len(kmers1 & kmers2)
        union = len(kmers1 | kmers2)
        jaccard_sim = intersection / union if union > 0 else 0.0
        return (simple_identity + jaccard_sim) / 2.0
    return simple_identity

def run_sourmash_sketch(input_file, ksize=31, scaled=1000):
    sketch_path = f"{input_file}.k{ksize}.sig"
    try:
        cmd = ["sourmash", "sketch", "dna", "-p", f"k={ksize},scaled={scaled}", "-o", sketch_path, input_file]
        result = subprocess.run(cmd, check=False, capture_output=True, timeout=30)
        if result.returncode == 0 and Path(sketch_path).exists():
            return sketch_path
    except Exception:
        pass
    return None

def run_sourmash_compare(query_sketch, allele_sketch):
    try:
        csv_path = "/tmp/comparison.csv"
        cmd = ["sourmash", "compare", query_sketch, allele_sketch, "--csv", csv_path]
        result = subprocess.run(cmd, check=False, capture_output=True, timeout=30)
        if result.returncode == 0 and Path(csv_path).exists():
            with open(csv_path) as f:
                lines = f.readlines()
                if len(lines) > 1:
                    parts = lines[1].split(',')
                    if len(parts) > 1:
                        try:
                            sim = float(parts[1])
                            return max(0.0, min(1.0, sim))
                        except (ValueError, IndexError):
                            pass
    except Exception:
        pass
    return None

query_file = "~{query_sequences}"
allele_file = "~{allele_fasta}"

query_sequences = parse_fasta(query_file)
if not query_sequences:
    print("{}")
    sys.exit(0)

allele_sequences = parse_fasta(allele_file)
alleles_by_locus = defaultdict(list)
for header, seq in allele_sequences:
    locus, allele_id = extract_locus_and_allele(header)
    alleles_by_locus[locus].append({'allele_id': allele_id, 'header': header, 'sequence': seq})

results = {}

for locus, alleles in alleles_by_locus.items():
    best_match = None
    best_similarity = -1.0

    for query_header, query_seq in query_sequences:
        for allele in alleles:
            allele_seq = allele['sequence']
            similarity = None
            try:
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
                Path(query_temp).unlink(missing_ok=True)
                Path(allele_temp).unlink(missing_ok=True)
                if query_sketch:
                    Path(query_sketch).unlink(missing_ok=True)
                if allele_sketch:
                    Path(allele_sketch).unlink(missing_ok=True)
            except Exception:
                similarity = None

            if similarity is None:
                similarity = compute_sequence_similarity(query_seq, allele_seq)

            if similarity > best_similarity:
                best_similarity = similarity
                best_match = allele['allele_id']

    if best_match is not None:
        results[locus] = {
            'allele_id': best_match,
            'similarity': max(0.0, min(1.0, best_similarity)),
            'confidence': best_similarity > 0.9
        }

json.dump(results, open('results/allele_calls.json', 'w'), indent=2)

PYTHON_EOF

        # Run the Python script
        python /tmp/minhash_allele_calling.py
    >>>

    output {
        File results = "results/allele_calls.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 2
        memory: "4 GB"
    }

    meta {
        author: "MinHash Allele Calling Implementation"
        description: "MinHash-based allele calling using sourmash for sequence comparison"
    }
}
