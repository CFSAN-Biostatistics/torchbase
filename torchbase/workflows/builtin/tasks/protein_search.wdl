version 1.0

# Translated protein search: tblastn of an operon's subunit reference set
# against nucleotide contigs. Mirrors StxTyper's own invocation
# (stxtyper.cpp:887 plus Hsp::blastp_fast in seq.hpp) so that hit sets, and
# therefore calls, match: composition-based statistics off, no SEG filtering,
# fixed -dbsize so scores do not drift with assembly size, word size 5,
# gapextend 2.
#
# The task emits BLAST's tabular output unchanged. Parsing it is the package
# layer's job (torchbase.operon.normalize_hsps), which keeps this container
# holding nothing but BLAST+ and keeps the column list and the code that reads
# it in one place.
task search_subunits {
    input {
        File contigs
        File subunit_reference
        Int gapextend = 2
        # String, not Float: WDL renders Float with %f, which silently turns
        # 1e-10 into "0.000000" and BLAST rejects it.
        String evalue = "1e-10"
        Int word_size = 5
        Int dbsize = 10000
        Int max_target_seqs = 10000
        Int db_gencode = 11
        Int num_threads = 4
        String image = "ncbi/blast:2.16.0"
    }

    command <<<
        set -e
        makeblastdb -in ~{contigs} -dbtype nucl -out contigs_db

        tblastn \
            -query ~{subunit_reference} \
            -db contigs_db \
            -comp_based_stats 0 \
            -seg no \
            -max_target_seqs ~{max_target_seqs} \
            -dbsize ~{dbsize} \
            -evalue ~{evalue} \
            -word_size ~{word_size} \
            -gapextend ~{gapextend} \
            -db_gencode ~{db_gencode} \
            -num_threads ~{num_threads} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen sstrand qseq sseq" \
            -out hits.tsv

        # tblastn writes nothing when there are no hits; downstream code reads
        # an empty table as "no operons", so make sure the file exists.
        touch hits.tsv
    >>>

    output {
        File hits = "hits.tsv"
    }

    runtime {
        docker: image
        cpu: num_threads
        memory: "4 GB"
    }
}
