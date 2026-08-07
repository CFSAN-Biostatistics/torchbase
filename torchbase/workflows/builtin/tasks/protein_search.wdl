version 1.0

# Protein-space search: tblastn of the operon's subunit reference set against
# translated nucleotide contigs (docs/operon-strategy-plan.md §4, mirrors
# stxtyper.cpp:887 including -gapextend 2). Normalizes BLAST tabular output to
# the internal HSP JSON list at this task boundary — nothing downstream ever
# sees BLAST format, so swapping the search engine later (e.g. for a Phraya
# protein-mode aligner once it gains affine gaps) only touches this task.
#
# Normalization itself lives in torchbase/operon.py, localized into the
# container as `operon_module`, so the algorithm that runs here is the one the
# unit tests cover.
task search_subunits {
    input {
        File contigs
        File subunit_reference
        File operon_config_json
        File operon_module
        Int gapextend = 2
        # String, not Float: WDL renders Float with %f, which silently turns
        # 1e-10 into "0.000000" and BLAST rejects it.
        String evalue = "1e-10"
        Int word_size = 5
        Int dbsize = 10000
        Int max_target_seqs = 10000
        Int db_gencode = 11
        Int num_threads = 4
    }

    command <<<
        set -e
        makeblastdb -in ~{contigs} -dbtype nucl -out contigs_db

        # Search parameters mirror StxTyper's own tblastn invocation
        # (stxtyper.cpp:887 + Hsp::blastp_fast in seq.hpp): composition-based
        # statistics off, no SEG filtering, fixed -dbsize so bitscores and
        # e-values do not drift with assembly size, word size 5, gapextend 2.
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

        python3 <<'PYTHON_SCRIPT'
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname("~{operon_module}"))
from operon import normalize_hsps, parse_fasta

with open("~{operon_config_json}") as f:
    cfg = json.load(f)
header_format = cfg.get("reference", {}).get(
    "header_format", "accession|subunit_role|reference_subtype|class"
)

with open("~{subunit_reference}") as f:
    references = parse_fasta(f.read())

with open("hits.tsv", newline="") as f:
    rows = [row for row in csv.reader(f, delimiter="\t") if row]

hsps = normalize_hsps(rows, references, header_format)

with open("hsps.json", "w") as f:
    json.dump([h.to_dict() for h in hsps], f, indent=2)
PYTHON_SCRIPT
    >>>

    output {
        File hsps = "hsps.json"
    }

    runtime {
        docker: "ncbi/blast:2.16.0"
        cpu: num_threads
        memory: "4 GB"
    }
}
