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

    call sketch_sequences as sketch_queries {
        input:
            sequences = query_sequences,
            ksize = ksize,
            scaled = sketch_size
    }

    call sketch_sequences as sketch_alleles {
        input:
            sequences = allele_fasta,
            ksize = ksize,
            scaled = sketch_size
    }

    call compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = allele_fasta
    }

    call call_alleles {
        input:
            similarity_matrix = compare_sketches.similarity_csv,
            query_sequences = query_sequences,
            allele_fasta = allele_fasta
    }

    output {
        File results = call_alleles.results
        String allele_profile = call_alleles.allele_profile
    }
}

task sketch_sequences {
    input {
        File sequences
        Int ksize = 31
        Int scaled = 1000
    }

    command <<<
        set -e
        sourmash sketch dna \
            -p k=~{ksize},scaled=~{scaled} \
            -o sequences.sig \
            ~{sequences}
    >>>

    output {
        File sketch = "sequences.sig"
    }

    runtime {
        docker: "quay.io/biocontainers/sourmash:4.8.11--hdfd78af_0"
        cpu: 1
        memory: "2 GB"
    }
}

task compare_sketches {
    input {
        File query_sketch
        File allele_sketch
        File allele_fasta
    }

    command <<<
        set -e
        sourmash compare \
            ~{query_sketch} \
            ~{allele_sketch} \
            --csv similarity.csv
    >>>

    output {
        File similarity_csv = "similarity.csv"
    }

    runtime {
        docker: "quay.io/biocontainers/sourmash:4.8.11--hdfd78af_0"
        cpu: 1
        memory: "2 GB"
    }
}

task call_alleles {
    input {
        File similarity_matrix
        File query_sequences
        File allele_fasta
    }

    command <<<
        set -e
        python3 <<CODE
import json
import csv
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

# Parse inputs
query_seqs = parse_fasta("~{query_sequences}")
allele_seqs = parse_fasta("~{allele_fasta}")

# Group alleles by locus
alleles_by_locus = defaultdict(list)
for idx, (header, seq) in enumerate(allele_seqs):
    locus, allele_id = extract_locus_and_allele(header)
    alleles_by_locus[locus].append({
        'allele_id': allele_id,
        'header': header,
        'index': idx
    })

# Read similarity matrix
with open("~{similarity_matrix}") as f:
    reader = csv.reader(f)
    rows = list(reader)

# Extract similarity scores (query vs alleles)
# Matrix format: first row/col are identifiers
# We want row 1 (first query) vs all allele columns
similarities = []
if len(rows) > 1 and len(rows[1]) > 1:
    # Skip first column (label), take rest as similarities
    similarities = [float(x) if x else 0.0 for x in rows[1][1:]]

# Find best match per locus
results = {}
profile_parts = []

for locus, alleles in sorted(alleles_by_locus.items()):
    best_match = None
    best_similarity = -1.0

    for allele in alleles:
        idx = allele['index']
        if idx < len(similarities):
            sim = similarities[idx]
            if sim > best_similarity:
                best_similarity = sim
                best_match = allele['allele_id']

    if best_match is not None:
        results[locus] = {
            'allele_id': best_match,
            'similarity': max(0.0, min(1.0, best_similarity)),
            'confidence': best_similarity > 0.9
        }
        profile_parts.append(f"{locus}_{best_match}")

# Write JSON output
with open('allele_calls.json', 'w') as f:
    json.dump(results, f, indent=2)

# Write profile string output
with open('allele_profile.txt', 'w') as f:
    f.write(','.join(profile_parts))

CODE
    >>>

    output {
        File results = "allele_calls.json"
        String allele_profile = read_string("allele_profile.txt")
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
