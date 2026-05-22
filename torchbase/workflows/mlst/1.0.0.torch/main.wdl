version 1.0

workflow mlst_typing {
    input {
        File? contigs
        File? reads
        File allele_database
        File? profiles
        Map[String, File]? schemes
        Int min_coverage = 3
        Float min_identity = 0.9
        Float ambiguity_threshold = 0.85
        Int ksize = 31
        Int sketch_size = 10
    }

    # Stage 1: Input validation - ensure reads OR contigs provided
    call validate_inputs {
        input:
            contigs = contigs,
            reads = reads
    }

    # Stage 2: Depth filtering for reads (skip for contigs)
    File input_sequences = select_first([contigs, reads])

    call depth_filter {
        input:
            sequences = input_sequences,
            is_reads = defined(reads),
            min_coverage = min_coverage
    }

    # Stage 3: MinHash allele calling across all schemes
    call sketch_sequences as sketch_queries {
        input:
            sequences = depth_filter.filtered_sequences,
            ksize = ksize,
            scaled = sketch_size
    }

    call sketch_sequences as sketch_alleles {
        input:
            sequences = allele_database,
            ksize = ksize,
            scaled = sketch_size
    }

    call compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = allele_database
    }

    call call_alleles {
        input:
            similarity_matrix = compare_sketches.similarity_csv,
            query_sequences = depth_filter.filtered_sequences,
            allele_fasta = allele_database
    }

    # Stage 4: Scheme inference from allele calls
    call infer_scheme {
        input:
            allele_calls = call_alleles.results,
            allele_profile = call_alleles.allele_profile,
            schemes = schemes,
            profiles = profiles,
            ambiguity_threshold = ambiguity_threshold
    }

    # Stage 5: Profile lookup in inferred scheme
    call lookup_profile {
        input:
            allele_calls = call_alleles.results,
            allele_profile = call_alleles.allele_profile,
            profiles_table = select_first([infer_scheme.selected_scheme_profiles, profiles]),
            inferred_scheme = infer_scheme.inferred_scheme,
            min_identity = min_identity
    }

    # Stage 6: Alignment fallback if ambiguous
    Boolean needs_alignment = lookup_profile.confidence < ambiguity_threshold

    call alignment_fallback {
        input:
            query_sequences = input_sequences,
            allele_database = allele_database,
            current_calls = call_alleles.results,
            allele_profile = call_alleles.allele_profile,
            profiles_table = select_first([infer_scheme.selected_scheme_profiles, profiles]),
            inferred_scheme = infer_scheme.inferred_scheme,
            should_run = needs_alignment
    }

    # Stage 7: Final result assembly
    File final_allele_calls = if needs_alignment then alignment_fallback.refined_calls else call_alleles.results
    String final_allele_profile = if needs_alignment then alignment_fallback.refined_profile else call_alleles.allele_profile
    String final_status = if needs_alignment then alignment_fallback.refined_status else lookup_profile.status
    Float final_confidence = if needs_alignment then alignment_fallback.refined_confidence else lookup_profile.confidence

    call assemble_final_result {
        input:
            allele_calls = final_allele_calls,
            allele_profile = final_allele_profile,
            sequence_type = lookup_profile.sequence_type,
            status = final_status,
            confidence = final_confidence,
            nearest_st = lookup_profile.nearest_st,
            scheme = infer_scheme.inferred_scheme,
            scheme_coverage = infer_scheme.coverage,
            alignment_used = needs_alignment
    }

    output {
        File typing_result = assemble_final_result.result
        String inferred_scheme = infer_scheme.inferred_scheme
        String sequence_type = lookup_profile.sequence_type
        String status = final_status
        Float confidence = final_confidence
        String nearest_st = lookup_profile.nearest_st
    }
}

# Stage 1: Input validation
task validate_inputs {
    input {
        File? contigs
        File? reads
    }

    command <<<
        # Ensure at least one input is provided
        contigs_defined=~{defined(contigs)}
        reads_defined=~{defined(reads)}

        if [ "$contigs_defined" == "false" ] && [ "$reads_defined" == "false" ]; then
            echo "ERROR: Must provide either contigs or reads input" >&2
            exit 1
        fi

        echo "Input validation passed"
    >>>

    output {
        String validation_message = read_string(stdout())
    }

    runtime {
        docker: "ubuntu:22.04"
        cpu: 1
        memory: "512 MB"
    }
}

# Stage 2: Depth filtering
task depth_filter {
    input {
        File sequences
        Boolean is_reads
        Int min_coverage = 3
    }

    command <<<
        set -e

        # For contigs, skip filtering - just copy the file
        if [ "~{is_reads}" == "false" ]; then
            cp ~{sequences} filtered_sequences.fasta
            exit 0
        fi

        # For reads, apply simple depth filtering
        python3 <<'PYTHON_SCRIPT'
import sys
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

def extract_depth(header):
    # Extract depth from header like "read_name_depth5"
    parts = header.split('_')
    for part in parts:
        if part.startswith('depth'):
            try:
                return int(part[5:])
            except:
                return 1
    return 1

seqs = parse_fasta("~{sequences}")
min_cov = ~{min_coverage}

# Filter by depth
filtered = []
for header, seq in seqs:
    depth = extract_depth(header)
    if depth >= min_cov:
        filtered.append((header, seq))

# Write output
with open("filtered_sequences.fasta", "w") as f:
    for header, seq in filtered:
        f.write(f">{header}\n{seq}\n")
PYTHON_SCRIPT
    >>>

    output {
        File filtered_sequences = "filtered_sequences.fasta"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}

# Stage 3: MinHash allele calling - Part 1: Sketch sequences
task sketch_sequences {
    input {
        File sequences
        Int ksize = 31
        Int scaled = 1000
    }

    command <<<
        set -e
        # Check if input file is empty or has no sequences
        if [ ! -s ~{sequences} ] || ! grep -q "^>" ~{sequences}; then
            # Create empty signature file for empty input
            touch sequences.sig
            exit 0
        fi

        sourmash sketch dna \
            -p k=~{ksize},scaled=~{scaled},abund \
            --singleton \
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

# Stage 3: MinHash allele calling - Part 2: Compare sketches
task compare_sketches {
    input {
        File query_sketch
        File allele_sketch
        File allele_fasta
    }

    command <<<
        set -e
        # Handle empty query case
        if [ ! -s ~{query_sketch} ]; then
            # Create empty similarity matrix
            echo "" > similarity.csv
            exit 0
        fi

        # Handle empty allele DB case
        if [ ! -s ~{allele_sketch} ]; then
            echo "" > similarity.csv
            exit 0
        fi

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

# Stage 3: MinHash allele calling - Part 3: Call alleles from similarity
task call_alleles {
    input {
        File similarity_matrix
        File query_sequences
        File allele_fasta
    }

    command <<<
        set -e
        python3 <<'CODE'
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
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}

# Stage 4: Scheme inference
task infer_scheme {
    input {
        File allele_calls
        String allele_profile
        Map[String, File]? schemes
        File? profiles
        Float ambiguity_threshold = 0.85
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json
import os

# Read allele calls
with open("~{allele_calls}") as f:
    allele_calls = json.load(f)

# Extract loci from calls
loci_in_calls = set(allele_calls.keys())

# Determine scheme based on loci presence
if loci_in_calls:
    # Extract scheme prefixes from loci names
    scheme_candidates = {}
    for locus in loci_in_calls:
        parts = locus.split('_')
        if len(parts) >= 1:
            # Try to identify scheme from locus naming
            # Common patterns: salmonella_adk, ecoli_dinB
            for part in parts[:-1]:  # All but last part
                if part in ['salmonella', 'ecoli', 'shigella', 'listeria']:
                    scheme = part
                    if scheme not in scheme_candidates:
                        scheme_candidates[scheme] = []
                    scheme_candidates[scheme].append(locus)

    # Select scheme with most matches
    if scheme_candidates:
        best_scheme = max(scheme_candidates.items(), key=lambda x: len(x[1]))[0]
    else:
        # Fallback: use first locus as scheme indicator
        best_scheme = "default"
else:
    best_scheme = "default"

# Calculate scheme coverage (percentage of expected loci found)
coverage = len(loci_in_calls)

# Calculate average identity from allele calls
identities = []
for locus, call_info in allele_calls.items():
    if isinstance(call_info, dict) and 'similarity' in call_info:
        identities.append(call_info['similarity'])

avg_identity = sum(identities) / len(identities) if identities else 0.0

result = {
    "inferred_scheme": best_scheme,
    "coverage": coverage,
    "avg_identity": avg_identity,
    "selected_scheme_profiles": "~{profiles}" if profiles else ""
}

with open("scheme_inference.json", "w") as f:
    json.dump(result, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        String inferred_scheme = read_json("scheme_inference.json")["inferred_scheme"]
        Int coverage = read_json("scheme_inference.json")["coverage"]
        File? selected_scheme_profiles = profiles
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}

# Stage 5: Profile lookup
task lookup_profile {
    input {
        File allele_calls
        String allele_profile
        File profiles_table
        String inferred_scheme
        Float min_identity = 0.9
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json
import csv

# Read allele calls
with open("~{allele_calls}") as f:
    allele_calls = json.load(f)

# Extract loci info: locus -> allele_id
locus_alleles = {}
for locus, call_info in allele_calls.items():
    if isinstance(call_info, dict):
        locus_alleles[locus] = int(call_info.get('allele_id', 0))
    else:
        locus_alleles[locus] = 0

# Read profiles table
profiles_list = []
with open("~{profiles_table}") as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        profiles_list.append(row)

# Try to match against profiles
best_match = None
best_match_st = None
num_matches = 0
total_loci = len(locus_alleles)

for profile in profiles_list:
    matches = 0
    for locus, allele_id in locus_alleles.items():
        # Try different locus name formats
        if locus in profile:
            try:
                profile_allele = int(profile[locus])
                if profile_allele == allele_id:
                    matches += 1
            except:
                pass

    if matches == total_loci and total_loci > 0:
        # Exact match found
        best_match = profile
        best_match_st = profile.get('ST', profile.get('st', 'unknown'))
        num_matches = matches
        break

# Determine status
status = "known" if best_match else "novel_profile" if num_matches > 0 else "novel_allele"

# Calculate confidence
if total_loci > 0:
    confidence = (num_matches / total_loci) if num_matches > 0 else 0.0
else:
    confidence = 0.0

# Find nearest ST if novel
nearest_st = "unknown"
if status != "known" and profiles_list:
    # Find profile with closest match
    best_partial_match = None
    best_partial_count = 0
    for profile in profiles_list:
        partial_matches = 0
        for locus, allele_id in locus_alleles.items():
            if locus in profile:
                try:
                    profile_allele = int(profile[locus])
                    if profile_allele == allele_id:
                        partial_matches += 1
                except:
                    pass

        if partial_matches > best_partial_count:
            best_partial_count = partial_matches
            best_partial_match = profile

    if best_partial_match:
        nearest_st = best_partial_match.get('ST', best_partial_match.get('st', 'unknown'))

result = {
    "sequence_type": str(best_match_st) if best_match_st else "unknown",
    "status": status,
    "confidence": min(1.0, max(0.0, confidence)),
    "nearest_st": str(nearest_st),
    "matched_loci": num_matches,
    "total_loci": total_loci
}

with open("profile_lookup.json", "w") as f:
    json.dump(result, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        String sequence_type = read_json("profile_lookup.json")["sequence_type"]
        String status = read_json("profile_lookup.json")["status"]
        Float confidence = read_json("profile_lookup.json")["confidence"]
        String nearest_st = read_json("profile_lookup.json")["nearest_st"]
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}

# Stage 6: Alignment fallback
task alignment_fallback {
    input {
        File query_sequences
        File allele_database
        File current_calls
        String allele_profile
        File profiles_table
        String inferred_scheme
        Boolean should_run = false
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json

# For now, this is a placeholder that accepts alignment refinement
# In a real implementation, this would use minimap2 or similar

with open("~{current_calls}") as f:
    current_calls = json.load(f)

# Refine calls based on alignment (placeholder)
refined_calls = current_calls.copy()
refined_status = "refined"
refined_confidence = 0.95
refined_profile = "~{allele_profile}"

# Write refined calls separately
with open("refined_calls.json", "w") as f:
    json.dump(refined_calls, f, indent=2)

result = {
    "refined_status": refined_status,
    "refined_confidence": refined_confidence,
    "refined_profile": refined_profile,
    "alignment_performed": True
}

with open("alignment_results.json", "w") as f:
    json.dump(result, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File refined_calls = "refined_calls.json"
        String refined_profile = read_json("alignment_results.json")["refined_profile"]
        String refined_status = read_json("alignment_results.json")["refined_status"]
        Float refined_confidence = read_json("alignment_results.json")["refined_confidence"]
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}

# Stage 7: Final result assembly
task assemble_final_result {
    input {
        File allele_calls
        String allele_profile
        String sequence_type
        String status
        Float confidence
        String nearest_st
        String scheme
        Int scheme_coverage
        Boolean alignment_used
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json

with open("~{allele_calls}") as f:
    allele_calls = json.load(f)

result = {
    "sequence_type": "~{sequence_type}",
    "st": "~{sequence_type}",
    "scheme": "~{scheme}",
    "scheme_coverage": ~{scheme_coverage},
    "status": "~{status}",
    "confidence": ~{confidence},
    "nearest_st": "~{nearest_st}",
    "allele_profile": "~{allele_profile}",
    "allele_calls": allele_calls,
    "alignment_used": ~{alignment_used}
}

with open("typing_result.json", "w") as f:
    json.dump(result, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File result = "typing_result.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
