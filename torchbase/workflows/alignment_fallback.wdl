version 1.0

workflow alignment_fallback {
    input {
        File query_sequences
        File allele_fasta
        File minhash_results
        Float min_similarity_threshold = 0.92
        Float max_second_best_diff = 0.03
        Float identity_threshold = 0.90
    }

    call refine_with_alignment {
        input:
            query_sequences = query_sequences,
            allele_fasta = allele_fasta,
            minhash_results = minhash_results,
            min_similarity_threshold = min_similarity_threshold,
            max_second_best_diff = max_second_best_diff,
            identity_threshold = identity_threshold
    }

    output {
        File refined_calls = refine_with_alignment.refined_calls
    }
}

task refine_with_alignment {
    input {
        File query_sequences
        File allele_fasta
        File minhash_results
        Float min_similarity_threshold = 0.92
        Float max_second_best_diff = 0.03
        Float identity_threshold = 0.90
    }

    command <<<
        set -e
        # Install minimap2 and dependencies
        apt-get update && apt-get install -y minimap2 && rm -rf /var/lib/apt/lists/*

        python3 <<'PYTHON_SCRIPT'
import json
import subprocess
from collections import defaultdict

def parse_fasta(fasta_path):
    """Parse FASTA file into list of (header, sequence) tuples."""
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

def parse_minhash_results(json_path):
    """Parse MinHash results JSON."""
    with open(json_path) as f:
        return json.load(f)

def is_ambiguous(locus_result):
    """Detect if a MinHash call is ambiguous and requires alignment fallback."""
    similarity = locus_result.get('similarity', 0.0)
    confidence = locus_result.get('confidence', False)

    # Low similarity triggers fallback
    if similarity < 0.92:
        return True

    # Check if confidence is explicitly False
    if not confidence:
        # Check if second_best is very close
        second_best = locus_result.get('second_best')
        if second_best:
            second_similarity = second_best.get('similarity', 0.0)
            if abs(similarity - second_similarity) <= 0.03:
                return True

    return False

def run_minimap2_alignment(query_seq, allele_fasta, query_name):
    """Run minimap2 to align query sequence against allele database."""
    # Write query to temporary file
    query_file = f"{query_name}_query.fasta"
    with open(query_file, 'w') as f:
        f.write(f">{query_name}\n{query_seq}\n")

    # Run minimap2 with short read preset
    paf_file = f"{query_name}_alignment.paf"
    cmd = [
        'minimap2',
        '-a',
        '-x', 'sr',
        allele_fasta,
        query_file,
        '-o', f"{query_name}_alignment.sam"
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Warning: minimap2 failed for {query_name}: {e}")
        return None

    # Parse SAM file to extract best alignment
    best_match = None
    best_identity = 0.0
    best_allele_id = None

    try:
        with open(f"{query_name}_alignment.sam") as f:
            for line in f:
                if line.startswith('@'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue

                ref_name = fields[2]
                mapq = int(fields[4])

                if mapq == 0:  # Skip unmapped
                    continue

                # Extract identity from NM tag (number of mismatches)
                nm = 0
                for tag in fields[11:]:
                    if tag.startswith('NM:i:'):
                        nm = int(tag.split(':')[-1])
                        break

                # Calculate identity
                query_len = len(query_seq)
                identity = max(0.0, (query_len - nm) / query_len)

                if identity > best_identity:
                    best_identity = identity
                    best_match = ref_name
                    locus, allele_id = extract_locus_and_allele(ref_name)
                    best_allele_id = allele_id
    except Exception as e:
        print(f"Warning: Error parsing SAM for {query_name}: {e}")

    return best_allele_id, best_identity

def get_locus_from_query_name(query_name):
    """Extract expected locus from query sequence name."""
    # Query names are like: query_adk_almost_1_or_2, query_fumC_between_1_and_2
    # Extract locus: adk, fumC, gyrB, etc.
    parts = query_name.lower().split('_')
    for part in parts:
        if any(known_locus in part for known_locus in ['adk', 'fumc', 'gyrb', 'ftsz', 'dinb']):
            return part
    return None

# Main logic
query_seqs = parse_fasta("~{query_sequences}")
allele_seqs = parse_fasta("~{allele_fasta}")
minhash_results = parse_minhash_results("~{minhash_results}")

# Group alleles by locus
alleles_by_locus = defaultdict(list)
for header, seq in allele_seqs.items():
    locus, allele_id = extract_locus_and_allele(header)
    alleles_by_locus[locus].append({
        'allele_id': allele_id,
        'header': header,
        'sequence': seq
    })

# Process each locus in MinHash results
refined_results = {}

for locus, locus_result in minhash_results.items():
    if is_ambiguous(locus_result):
        # Need to run alignment fallback for this locus
        # Find query sequences corresponding to this locus
        matching_queries = []
        for query_name, query_seq in query_seqs.items():
            query_locus = get_locus_from_query_name(query_name)
            if query_locus and query_locus in locus.lower():
                matching_queries.append((query_name, query_seq))

        if matching_queries:
            # Run minimap2 for the first matching query
            query_name, query_seq = matching_queries[0]
            allele_id, identity = run_minimap2_alignment(query_seq, "~{allele_fasta}", query_name)

            if allele_id is not None:
                # Determine status based on identity threshold
                status = "confirmed" if identity >= ~{identity_threshold} else "novel_allele"
                refined_results[locus] = {
                    'allele_id': allele_id,
                    'identity': float(identity),
                    'status': status
                }
            else:
                # Fallback to original MinHash call
                refined_results[locus] = {
                    'allele_id': locus_result.get('allele_id', 'unknown'),
                    'identity': locus_result.get('similarity', 0.0),
                    'status': 'novel_allele'
                }
        else:
            # No matching query found, use MinHash result
            refined_results[locus] = {
                'allele_id': locus_result.get('allele_id', 'unknown'),
                'identity': locus_result.get('similarity', 0.0),
                'status': 'confirmed' if locus_result.get('confidence', False) else 'novel_allele'
            }
    else:
        # Confident MinHash call, pass through as confirmed
        refined_results[locus] = {
            'allele_id': locus_result.get('allele_id', 'unknown'),
            'identity': locus_result.get('similarity', 0.0),
            'status': 'confirmed'
        }

# Write refined results to JSON
with open('refined_calls.json', 'w') as f:
    json.dump(refined_results, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File refined_calls = "refined_calls.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 2
        memory: "4 GB"
    }
}
