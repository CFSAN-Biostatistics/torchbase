version 1.0

# MinHash screen: sketch an assembly or read set and a torch's allele database,
# then compare them. This is the `fast` typing strategy's whole compute step and
# the first step of `balanced`.
#
# The workflow emits sourmash's similarity matrix and stops there. Turning
# similarities into allele calls is a threshold decision, not compute, and
# happens at the package layer (torchbase/allele_calls.py) — which is also
# where the strategy decides whether the result is confident enough to skip
# alignment.

import "tasks/minhash.wdl" as minhash

workflow minhash_screen {
    input {
        File query_sequences
        File allele_fasta
        Int ksize = 31
        Int scaled = 1000
    }

    call minhash.sketch_sequences as sketch_queries {
        input:
            sequences = query_sequences,
            ksize = ksize,
            scaled = scaled
    }

    call minhash.compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_fasta = allele_fasta,
            ksize = ksize,
            scaled = scaled
    }

    output {
        File similarity_matrix = compare_sketches.similarity_csv
    }
}
