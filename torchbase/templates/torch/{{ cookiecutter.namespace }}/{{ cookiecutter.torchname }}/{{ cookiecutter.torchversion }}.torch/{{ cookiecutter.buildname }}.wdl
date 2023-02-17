version 1.0


import "https://github.com/CFSAN-Biostatistics/wdl-commons/blob/main/identify.wdl" as identify
import "https://github.com/CFSAN-Biostatistics/wdl-commons/blob/main/compression.wdl" as compression
import "https://github.com/CFSAN-Biostatistics/wdl-commons/blob/main/paired_reads.wdl" as paired_reads


workflow {{ cookiecutter.torchname }} {

    # Take reads (paired/interleaved) or contigs

    input {
        Pair[File, File]? paired_reads
        File? interlaced_or_single_reads
        File? contigs
        String? user_provided_identifier

        File profiles
        Array[File] references
    }

    # Identify and interleave the reads if necessary

    if (defined(paired_reads) || defined(interlaced_or_single_reads)) {
        if (defined(paired_reads)){
            call interlace { input:forward=paired_reads.left, reverse=paired_reads.right }
            File reads = interlace.reads
        }
        if (!defined(paired_reads)){
            File reads = interlaced_or_single_reads
        }
        call identify_reads { input:reads=reads }
        String reads_name = identify_reads.name
        Map[String, String] qualities =  fastp.output
    }

    if (defined(contigs)){
        call identify_contigs { input:contigs=contigs }
        String contigs_name = identify_contigs.name
    }

    File query = select_first([contigs, reads])
    String name = select_first([user_provided_identifier, contigs_name, reads_name])


    # The rest of the pipeline

}