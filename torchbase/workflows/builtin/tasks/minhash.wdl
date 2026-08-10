version 1.0

# MinHash sketching and comparison with sourmash.
#
# Compute only: sketches go in, a similarity matrix comes out. The threshold
# that turns a similarity into an allele call lives at the package layer
# (torchbase/allele_calls.py), so `fast` and `balanced` cannot disagree about
# what "confident" means.

task sketch_sequences {
    input {
        File sequences
        Int ksize = 31
        Int scaled = 1000
        String image = "quay.io/biocontainers/sourmash:4.8.11--hdfd78af_0"
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
        docker: image
        cpu: 1
        memory: "2 GB"
    }
}

task compare_sketches {
    input {
        File query_sketch
        File allele_fasta
        Int ksize = 31
        Int scaled = 1000
        String image = "quay.io/biocontainers/sourmash:4.8.11--hdfd78af_0"
    }

    command <<<
        set -e
        # One signature per allele sequence.
        sourmash sketch dna \
            -p k=~{ksize},scaled=~{scaled} \
            ~{allele_fasta} \
            -o alleles.sig \
            --singleton

        # ANI-based similarity of every query against every allele. Row and
        # column order is queries-then-alleles, which is the layout
        # torchbase.allele_calls expects.
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
        docker: image
        cpu: 2
        memory: "4 GB"
    }
}
