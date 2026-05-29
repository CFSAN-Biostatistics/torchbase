version 1.0

import "tasks/minhash.wdl" as minhash
import "tasks/alignment.wdl" as alignment
import "tasks/profile_lookup.wdl" as profile_lookup
import "tasks/filter_alleles.wdl" as filter

workflow balanced_typing {
    input {
        File query_sequences
        File allele_fasta
        File profiles_table
        String input_type = "contigs"
        Float confidence_threshold = 0.85
        File? quality_json
        Boolean exclude_suspect_alleles = false
        Boolean exclude_suspect_loci = false
        Boolean exclude_suspect_profiles = false
    }

    # Step 0: Filter alleles based on quality.json if provided
    call filter.filter_alleles {
        input:
            allele_fasta = allele_fasta,
            quality_json = quality_json,
            exclude_suspect_alleles = exclude_suspect_alleles,
            exclude_suspect_loci = exclude_suspect_loci,
            exclude_suspect_profiles = exclude_suspect_profiles
    }

    # Use filtered alleles for all downstream tasks
    File working_allele_fasta = filter_alleles.filtered_fasta

    # Step 1: MinHash sketching and comparison
    call minhash.sketch_sequences as sketch_queries {
        input:
            sequences = query_sequences,
            ksize = 31,
            scaled = 1000
    }

    call minhash.sketch_sequences as sketch_alleles {
        input:
            sequences = working_allele_fasta,
            ksize = 31,
            scaled = 1000
    }

    call minhash.compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = working_allele_fasta
    }

    # Step 2: Call alleles using MinHash
    call minhash.call_alleles_minhash {
        input:
            similarity_matrix = compare_sketches.similarity_csv,
            query_sequences = query_sequences,
            allele_fasta = working_allele_fasta,
            confidence_threshold = confidence_threshold
    }

    # Step 3: Determine if alignment fallback is needed
    # We need to check if any locus has confidence < threshold
    call check_confidence_for_alignment {
        input:
            allele_calls = call_alleles_minhash.allele_calls,
            confidence_threshold = confidence_threshold
    }

    # Step 4: Conditional alignment fallback
    # Uses asm5 preset for contigs, sr preset for reads
    call alignment.align_and_call as alignment_fallback {
        input:
            query_sequences = query_sequences,
            allele_fasta = working_allele_fasta,
            input_type = input_type,
            identity_threshold = 0.90
    }

    # Step 5: Choose allele calls based on whether alignment was needed
    call merge_allele_calls {
        input:
            minhash_calls = call_alleles_minhash.allele_calls,
            alignment_calls = alignment_fallback.alignment_results,
            use_alignment = check_confidence_for_alignment.use_alignment,
            input_type = input_type,
            confidence_threshold = confidence_threshold
    }

    # Step 6: Lookup profile and generate final result with method metadata
    call profile_lookup.lookup_profile {
        input:
            allele_calls = merge_allele_calls.final_calls,
            profiles_table = profiles_table,
            strategy = "balanced",
            alignment_used = check_confidence_for_alignment.use_alignment
    }

    # Step 7: Merge exclusion metadata into result
    call add_exclusion_metadata {
        input:
            typing_result = lookup_profile.result,
            exclusions = filter_alleles.exclusions
    }

    output {
        # Output result JSON with standardized format including method metadata and exclusions
        File result = add_exclusion_metadata.final_result
    }
}

task check_confidence_for_alignment {
    input {
        File allele_calls
        Float confidence_threshold = 0.85
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json

def check_alignment_needed(json_path, threshold):
    """Check if any allele call has low confidence."""
    with open(json_path) as f:
        calls = json.load(f)

    for locus, call_data in calls.items():
        if 'confidence' in call_data:
            if not call_data['confidence']:
                return True
        else:
            # No confidence field, check similarity
            if 'similarity' in call_data:
                if call_data['similarity'] < threshold:
                    return True

    return False

needs_alignment = check_alignment_needed("~{allele_calls}", ~{confidence_threshold})

with open('use_alignment.json', 'w') as f:
    json.dump(needs_alignment, f)

PYTHON_SCRIPT
    >>>

    output {
        Boolean use_alignment = read_json("use_alignment.json")
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "1 GB"
    }
}

task merge_allele_calls {
    input {
        File minhash_calls
        File alignment_calls
        Boolean use_alignment
        String input_type = "contigs"
        Float confidence_threshold = 0.85
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json

def merge_calls(minhash_file, alignment_file, use_alignment, threshold):
    """Merge MinHash and alignment calls."""
    with open(minhash_file) as f:
        minhash_calls = json.load(f)

    with open(alignment_file) as f:
        alignment_calls = json.load(f)

    if not use_alignment:
        # Use MinHash calls only
        return minhash_calls

    # Use alignment calls where available, fall back to MinHash
    merged = {}
    for locus in minhash_calls.keys():
        if locus in alignment_calls:
            # Use alignment call if it has good confidence
            alignment_call = alignment_calls[locus]
            if alignment_call.get('identity', 0.0) >= threshold:
                merged[locus] = alignment_call
            else:
                # Fall back to MinHash
                merged[locus] = minhash_calls[locus]
        else:
            # No alignment result, use MinHash
            merged[locus] = minhash_calls[locus]

    return merged

use_align = ~{true} if ~{use_alignment} else ~{false}
final_calls = merge_calls("~{minhash_calls}", "~{alignment_calls}", use_align, ~{confidence_threshold})

with open('final_calls.json', 'w') as f:
    json.dump(final_calls, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File final_calls = "final_calls.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "1 GB"
    }
}

task add_exclusion_metadata {
    input {
        File typing_result
        File exclusions
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json

# Load typing result
with open("~{typing_result}") as f:
    result = json.load(f)

# Load exclusions
with open("~{exclusions}") as f:
    exclusions = json.load(f)

# Add exclusion metadata to result
if 'notes' not in result:
    result['notes'] = {}

result['notes']['exclusions'] = {
    'excluded_alleles': exclusions['excluded_alleles'],
    'excluded_loci': exclusions['excluded_loci'],
    'num_excluded_alleles': exclusions['num_excluded_alleles'],
    'num_excluded_loci': exclusions['num_excluded_loci']
}

# Write final result
with open('final_result.json', 'w') as f:
    json.dump(result, f, indent=2)

PYTHON_SCRIPT
    >>>

    output {
        File final_result = "final_result.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "1 GB"
    }
}
