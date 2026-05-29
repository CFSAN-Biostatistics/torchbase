version 1.0

import "tasks/minhash.wdl" as minhash
import "tasks/alignment.wdl" as alignment
import "tasks/profile_lookup.wdl" as profile_lookup
import "tasks/filter_alleles.wdl" as filter

workflow sensitive_typing {
    input {
        File query_sequences
        File allele_database
        File profiles
        String preset = "asm5"
        Float confidence_threshold = 0.95
        File? quality_json
        Boolean exclude_suspect_alleles = false
        Boolean exclude_suspect_loci = false
        Boolean exclude_suspect_profiles = false
    }

    # Step 0: Filter alleles if quality.json provided
    call filter.filter_alleles {
        input:
            allele_fasta = allele_database,
            quality_json = quality_json,
            exclude_suspect_alleles = exclude_suspect_alleles,
            exclude_suspect_loci = exclude_suspect_loci,
            exclude_suspect_profiles = exclude_suspect_profiles
    }

    File working_allele_fasta = filter_alleles.filtered_fasta

    # Step 1: Sketch query sequences with MinHash (for guidance only)
    call minhash.sketch_sequences as sketch_queries {
        input:
            sequences = query_sequences,
            ksize = 31,
            scaled = 1000
    }

    # Step 2: Sketch allele database with MinHash (for guidance only)
    call minhash.sketch_sequences as sketch_alleles {
        input:
            sequences = working_allele_fasta,
            ksize = 31,
            scaled = 1000
    }

    # Step 3: Compare sketches (guidance only)
    call minhash.compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = working_allele_fasta
    }

    # Step 4: ALWAYS run full alignment with strict parameters using minimap2
    # In sensitive mode, alignment is not conditional - it always runs
    # Uses minimap2 with asm5 or asm5+eqx preset for high accuracy
    call alignment.align_and_call as alignment_call {
        input:
            query_sequences = query_sequences,
            allele_fasta = working_allele_fasta,
            input_type = "contigs",
            identity_threshold = confidence_threshold
    }

    # Step 5: Lookup profile from alignment-based allele calls
    call profile_lookup.lookup_profile as profile_call {
        input:
            allele_calls = alignment_call.alignment_results,
            profiles_table = profiles,
            strategy = "sensitive",
            alignment_used = true
    }

    # Step 6: Add exclusion metadata
    call add_exclusion_metadata {
        input:
            typing_result = profile_call.result,
            exclusions = filter_alleles.exclusions
    }

    output {
        File typing_result = add_exclusion_metadata.final_result
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

with open("~{typing_result}") as f:
    result = json.load(f)

with open("~{exclusions}") as f:
    exclusions = json.load(f)

if 'notes' not in result:
    result['notes'] = {}

result['notes']['exclusions'] = {
    'excluded_alleles': exclusions['excluded_alleles'],
    'excluded_loci': exclusions['excluded_loci'],
    'num_excluded_alleles': exclusions['num_excluded_alleles'],
    'num_excluded_loci': exclusions['num_excluded_loci']
}

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
}
