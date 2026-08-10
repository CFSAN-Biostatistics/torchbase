version 1.0

# Alignment of query sequences against a torch's allele database with minimap2.
# This is the expensive step of the `sensitive` strategy, and the fallback the
# `balanced` strategy dispatches when a MinHash screen is not confident.
#
# The workflow emits SAM and stops there. Deriving per-allele identity from the
# SAM `NM` tag, choosing the best allele per locus and applying an identity
# threshold are decisions, not compute, and happen at the package layer
# (torchbase/allele_calls.py).

import "tasks/alignment.wdl" as alignment

workflow allele_alignment {
    input {
        File query_sequences
        File allele_fasta
        # minimap2 preset, passed through verbatim: asm5/asm20 for assemblies,
        # sr for short reads.
        String preset = "asm5"
    }

    call alignment.align_sequences {
        input:
            query = query_sequences,
            reference = allele_fasta,
            preset = preset
    }

    output {
        File alignments = align_sequences.alignments
    }
}
