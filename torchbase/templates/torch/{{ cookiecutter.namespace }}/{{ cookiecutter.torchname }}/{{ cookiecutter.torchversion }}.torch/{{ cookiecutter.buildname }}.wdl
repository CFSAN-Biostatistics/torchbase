version 1.0


import "https://github.com/CFSAN-Biostatistics/wdl-commons/blob/main/identify.wdl" as identify
import "https://github.com/CFSAN-Biostatistics/wdl-commons/blob/main/compression.wdl" as compression
import "https://github.com/CFSAN-Biostatistics/wdl-commons/blob/main/paired_reads.wdl" as paired_reads


workflow {{ cookiecutter.torchname }} {

    # Take reads (paired/interleaved) or contigs

    input {
        File? interlaced_or_single_reads
        Boolean? are_interlaced = false
        File? contigs
        File? longreads
        String? user_provided_identifier

        File profiles
        Array[File] references
    }

    # 

    if (defined(interlaced_or_single_reads)) {
        call identify_reads { input:reads=reads }
        String reads_name = identify_reads.name
        Map[String, String] qualities =  fastp.output
    }

    if (defined(contigs)){
        call identify_contigs { input:contigs=contigs }
        String contigs_name = identify_contigs.name
    }

    if (defined(longreads)){

    }

    File query = select_first([contigs, reads])
    String name = select_first([user_provided_identifier, contigs_name, reads_name])
    # The rest of the pipeline

}


# Tasks