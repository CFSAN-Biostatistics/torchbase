version 1.0

import "tasks/minhash.wdl" as minhash_tasks
import "tasks/profile_lookup.wdl" as profile_tasks
import "tasks/filter_alleles.wdl" as filter

workflow fast_typing {
    input {
        File query_sequences
        File allele_database
        File profiles_table
        Int ksize = 31
        Int sketch_size = 1000
        File? quality_json
        Boolean exclude_suspect_alleles = false
        Boolean exclude_suspect_loci = false
        Boolean exclude_suspect_profiles = false
    }

    # Filter alleles if quality.json provided
    call filter.filter_alleles {
        input:
            allele_fasta = allele_database,
            quality_json = quality_json,
            exclude_suspect_alleles = exclude_suspect_alleles,
            exclude_suspect_loci = exclude_suspect_loci,
            exclude_suspect_profiles = exclude_suspect_profiles
    }

    File working_allele_fasta = filter_alleles.filtered_fasta

    call minhash_tasks.sketch_sequences as sketch_queries {
        input:
            sequences = query_sequences,
            ksize = ksize,
            scaled = sketch_size
    }

    call minhash_tasks.sketch_sequences as sketch_alleles {
        input:
            sequences = working_allele_fasta,
            ksize = ksize,
            scaled = sketch_size
    }

    call minhash_tasks.compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = working_allele_fasta
    }

    call minhash_tasks.call_alleles {
        input:
            similarity_matrix = compare_sketches.similarity_csv,
            query_sequences = query_sequences,
            allele_fasta = working_allele_fasta
    }

    call profile_tasks.lookup_profile {
        input:
            allele_calls = call_alleles.results,
            profiles_table = profiles_table,
            strategy = "fast",
            alignment_used = false
    }

    call add_exclusion_metadata {
        input:
            typing_result = lookup_profile.result,
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
