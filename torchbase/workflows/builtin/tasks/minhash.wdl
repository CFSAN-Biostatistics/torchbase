version 1.0

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
            ~{sequences} \
            -o sequences.sig \
            --singleton
    >>>

    output {
        File sketch = "sequences.sig"
    }

    runtime {
        cpu: 1
        memory: "2 GB"
    }
}

task compare_sketches {
    input {
        File query_sketch
        File allele_sketch
        File allele_fasta
        Int ksize = 31
        Int scaled = 1000
    }

    command <<<
        set -e
        # Sketch allele FASTA with one signature per sequence
        sourmash sketch dna \
            -p k=~{ksize},scaled=~{scaled} \
            ~{allele_fasta} \
            -o alleles.sig \
            --singleton

        # Compare query against allele sketches, output ANI-based CSV matrix
        sourmash compare \
            ~{query_sketch} alleles.sig \
            --csv similarity.csv \
            --ani \
            -k ~{ksize}
    >>>

    output {
        File similarity_csv = "similarity.csv"
    }

    runtime {
        cpu: 2
        memory: "4 GB"
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

# Handle empty similarity matrix
if len(rows) <= 1:
    # Empty result
    with open('allele_calls.json', 'w') as f:
        json.dump({}, f)
    with open('allele_profile.txt', 'w') as f:
        f.write('')
    exit(0)

# Extract similarity scores (query vs alleles)
# Matrix format with --singleton: all-vs-all NxN matrix
# Row 0: header with N identifiers (query1...queryN, allele1...alleleM)
# Rows 1..N: similarity data (no row labels)
# Query seqs occupy first num_queries rows/cols
# Allele seqs occupy next num_alleles cols

num_queries = len(query_seqs)
num_alleles = len(allele_seqs)

# Validate matrix dimensions
expected_size = num_queries + num_alleles
if len(rows) != expected_size + 1:  # +1 for header
    raise ValueError(f"Matrix size mismatch: expected {expected_size+1} rows, got {len(rows)}")

# Extract max similarity across all queries for each allele
max_similarities = [0.0] * num_alleles

for query_idx in range(num_queries):
    data_row_idx = query_idx + 1  # Skip header row
    if data_row_idx < len(rows):
        row = rows[data_row_idx]
        # Allele columns start at num_queries
        for allele_idx in range(num_alleles):
            col_idx = num_queries + allele_idx
            if col_idx < len(row):
                sim = float(row[col_idx]) if row[col_idx] else 0.0
                max_similarities[allele_idx] = max(max_similarities[allele_idx], sim)

# Find best match per locus
results = {}
profile_parts = []

for locus, alleles in sorted(alleles_by_locus.items()):
    best_match = None
    best_similarity = -1.0

    for allele in alleles:
        idx = allele['index']
        if idx < len(max_similarities):
            sim = max_similarities[idx]
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
        cpu: 2
        memory: "2 GB"
    }
}

task call_alleles_minhash {
    input {
        File similarity_matrix
        File query_sequences
        File allele_fasta
        Float confidence_threshold = 0.85
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

# Handle empty similarity matrix
if len(rows) <= 1:
    # Empty result
    with open('allele_calls.json', 'w') as f:
        json.dump({}, f)
    exit(0)

# Extract similarity scores (query vs alleles)
num_queries = len(query_seqs)
num_alleles = len(allele_seqs)

# Validate matrix dimensions
expected_size = num_queries + num_alleles
if len(rows) != expected_size + 1:
    raise ValueError(f"Matrix size mismatch: expected {expected_size+1} rows, got {len(rows)}")

# Extract max similarity across all queries for each allele
max_similarities = [0.0] * num_alleles

for query_idx in range(num_queries):
    data_row_idx = query_idx + 1
    if data_row_idx < len(rows):
        row = rows[data_row_idx]
        for allele_idx in range(num_alleles):
            col_idx = num_queries + allele_idx
            if col_idx < len(row):
                sim = float(row[col_idx]) if row[col_idx] else 0.0
                max_similarities[allele_idx] = max(max_similarities[allele_idx], sim)

# Find best match per locus
results = {}

for locus, alleles in sorted(alleles_by_locus.items()):
    best_match = None
    best_similarity = -1.0

    for allele in alleles:
        idx = allele['index']
        if idx < len(max_similarities):
            sim = max_similarities[idx]
            if sim > best_similarity:
                best_similarity = sim
                best_match = allele['allele_id']

    if best_match is not None:
        confidence = best_similarity >= ~{confidence_threshold}
        results[locus] = {
            'allele_id': best_match,
            'similarity': max(0.0, min(1.0, best_similarity)),
            'confidence': confidence
        }

# Write JSON output
with open('allele_calls.json', 'w') as f:
    json.dump(results, f, indent=2)

CODE
    >>>

    output {
        File allele_calls = "allele_calls.json"
    }

    runtime {
        cpu: 2
        memory: "2 GB"
    }
}
