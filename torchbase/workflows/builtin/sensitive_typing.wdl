version 1.0

import "tasks/minhash.wdl" as minhash
import "tasks/alignment.wdl" as alignment
import "tasks/profile_lookup.wdl" as profile_lookup

workflow sensitive_typing {
    input {
        File query_sequences
        File allele_database
        File profiles
        String preset = "asm5"
        Float confidence_threshold = 0.95
    }

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
            sequences = allele_database,
            ksize = 31,
            scaled = 1000
    }

    # Step 3: Compare sketches (guidance only)
    call minhash.compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = allele_database
    }

    # Step 4: ALWAYS run full alignment with strict parameters using minimap2
    # In sensitive mode, alignment is not conditional - it always runs
    # Uses minimap2 with asm5 or asm5+eqx preset for high accuracy
    call alignment.align_and_call as alignment_call {
        input:
            query_sequences = query_sequences,
            allele_fasta = allele_database,
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

    output {
        File typing_result = profile_call.result
    }
}
