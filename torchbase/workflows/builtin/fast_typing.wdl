version 1.0

workflow fast_typing {
    input {
        File query_sequences
        File allele_fasta
        Int ksize = 31
        Int sketch_size = 1000
    }

    call sketch_sequences as sketch_queries {
        input:
            sequences = query_sequences,
            ksize = ksize,
            scaled = sketch_size
    }

    call sketch_sequences as sketch_alleles {
        input:
            sequences = allele_fasta,
            ksize = ksize,
            scaled = sketch_size
    }

    call compare_sketches {
        input:
            query_sketch = sketch_queries.sketch,
            allele_sketch = sketch_alleles.sketch,
            allele_fasta = allele_fasta
    }

    call call_alleles {
        input:
            similarity_matrix = compare_sketches.similarity_csv,
            query_sequences = query_sequences,
            allele_fasta = allele_fasta
    }

    output {
        File results = call_alleles.results
        String allele_profile = call_alleles.allele_profile
    }
}

task sketch_sequences {
    input {
        File sequences
        Int ksize = 31
        Int scaled = 1000
    }

    command <<<
        set -e
        if [ ! -s ~{sequences} ] || ! grep -q "^>" ~{sequences}; then
            touch sequences.sig
            exit 0
        fi
        sourmash sketch dna -p k=~{ksize},scaled=~{scaled},abund --singleton -o sequences.sig ~{sequences}
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
        File allele_fasta
    }

    command <<<
        set -e
        if [ ! -s ~{query_sketch} ]; then
            echo "" > similarity.csv
            exit 0
        fi
        if [ ! -s ~{allele_sketch} ]; then
            echo "" > similarity.csv
            exit 0
        fi
        sourmash compare ~{query_sketch} ~{allele_sketch} --csv similarity.csv
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

task call_alleles {
    input {
        File similarity_matrix
        File query_sequences
        File allele_fasta
    }

    command <<<
        echo "Calling alleles"
    >>>

    output {
        File results = "allele_calls.json"
        String allele_profile = ""
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
