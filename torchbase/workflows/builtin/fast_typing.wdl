version 1.0

import "tasks/minhash.wdl" as minhash_tasks
import "tasks/profile_lookup.wdl" as profile_tasks

workflow fast_typing {
    input {
        File query_sequences
        File allele_database
        File profiles_table
        Int ksize = 31
        Int sketch_size = 1000
    }

    call minhash_tasks.sketch_sequences as sketch_queries {
        input:
            sequences = query_sequences,
            ksize = ksize,
            scaled = sketch_size
    }

    call minhash_tasks.sketch_sequences as sketch_alleles {
        input:
            sequences = allele_database,
            ksize = ksize,
            scaled = sketch_size
    }

    call minhash_tasks.compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = allele_database
    }

    call minhash_tasks.call_alleles {
        input:
            similarity_matrix = compare_sketches.similarity_csv,
            query_sequences = query_sequences,
            allele_fasta = allele_database
    }

    call profile_tasks.lookup_profile {
        input:
            allele_calls = call_alleles.results,
            profiles_table = profiles_table,
            strategy = "fast",
            alignment_used = false
    }

    output {
        File typing_result = lookup_profile.result
    }
}
