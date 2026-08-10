version 1.0

# Compute step for the "operon" typing model (docs/operon-strategy-plan.md §4):
# translated protein search of an operon's subunit reference set against an
# assembly. That is the only expensive part of operon typing, and the only part
# that needs a tool container.
#
# Everything downstream — frameshift stitching, locus reduction, synteny
# pairing, thresholds, residue tables, the status ladder — is interpretation,
# not compute, and runs at the package layer in torchbase/operon.py. The
# workflow therefore emits raw BLAST tabular output; `torchbase run` turns it
# into operon calls. Keeping the search behind one task boundary is also what
# makes swapping the search engine later (e.g. a protein-mode Phraya) a
# single-task change.

import "tasks/protein_search.wdl" as search

workflow operon_typing {
    input {
        File contigs
        File subunit_reference
    }

    call search.search_subunits {
        input:
            contigs = contigs,
            subunit_reference = subunit_reference
    }

    output {
        File hits = search_subunits.hits
    }
}
