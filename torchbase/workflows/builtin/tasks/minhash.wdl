version 1.0

# MinHash sketching and comparison tasks for allele matching
# These tasks provide sourmash-based sketch generation and comparison
# for fast allele calling and similarity estimation

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

task compare_sketches {
    input {
        File query_sketch
        File allele_sketch
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
