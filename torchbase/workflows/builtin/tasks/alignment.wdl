version 1.0

# Nucleotide alignment of query sequences against an allele database.
#
# Compute only: minimap2 runs, SAM comes out. Identity computation, best-allele
# selection and thresholding live at the package layer in
# torchbase/allele_calls.py, so the same numbers are used whether the calls come
# from an alignment or a MinHash screen.
task align_sequences {
    input {
        File query
        File reference
        # Passed verbatim to minimap2 -x. asm5/asm20 for assemblies, sr for
        # short reads; an invalid preset makes minimap2 exit non-zero rather
        # than silently degrading.
        String preset = "asm5"
        Int num_threads = 2
        String image = "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
    }

    command <<<
        set -e
        minimap2 \
            -a \
            -x ~{preset} \
            --secondary=no \
            -t ~{num_threads} \
            ~{reference} \
            ~{query} \
            -o alignments.sam
    >>>

    output {
        File alignments = "alignments.sam"
    }

    runtime {
        docker: image
        cpu: num_threads
        memory: "4 GB"
    }
}
