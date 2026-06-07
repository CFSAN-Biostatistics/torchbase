version 1.0

import "../../builtin/balanced_typing.wdl" as balanced

workflow mlst {
    input {
        File query_sequences
        File allele_fasta
        File profiles_table
        String input_type = "contigs"
        Float confidence_threshold = 0.85
        Int min_depth = 3
        File? quality_json
        Boolean exclude_suspect_alleles = false
        Boolean exclude_suspect_loci = false
        Boolean exclude_suspect_profiles = false
    }

    call balanced.balanced_typing {
        input:
            query_sequences = query_sequences,
            allele_fasta = allele_fasta,
            profiles_table = profiles_table,
            input_type = input_type,
            confidence_threshold = confidence_threshold,
            min_depth = min_depth,
            quality_json = quality_json,
            exclude_suspect_alleles = exclude_suspect_alleles,
            exclude_suspect_loci = exclude_suspect_loci,
            exclude_suspect_profiles = exclude_suspect_profiles
    }

    output {
        File result = balanced_typing.result
    }
}
