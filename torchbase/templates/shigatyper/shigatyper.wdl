version 1.0

# This project constitutes a work of the United States Government and is not
# subject to domestic copyright protection under 17 USC § 105. No Rights Are 
# Reserved.

# This program is distributed in the hope that it will be useful. Responsibility
# for the use of the system and interpretation of documentation and results lies
# solely with the user. In no event shall FDA be liable for direct, indirect,
# special, incidental, or consequential damages resulting from the use, misuse,
# or inability to use the system and accompanying documentation. Third parties'
# use of or acknowledgment of the system does not in any way represent that
# FDA endorses such third parties or expresses any opinion with respect to their
# statements. 

# This program is free software: you can redistribute it and/or modify it.

import "https://github.com/biowdl/tasks/blob/develop/minimap2.wdl" as minimap2

workflow shigatyper {

    input {
        Pair[File, File]? paired_reads
        File? interlaced_or_single_reads
        File? contigs
        String? user_provided_identifier
        Array[File] references
    }

    if (defined(paired_reads) || defined(interlaced_or_single_reads)) {
        if (defined(paired_reads)){
            call interlace { input:forward=paired_reads.left, reverse=paired_reads.right }
            File reads = interlace.reads
        }
        if (!defined(paired_reads)){
            File reads = interlaced_or_single_reads
        }
        call identify_reads { input:reads=reads }
        call fastp { input:reads=reads }
        String reads_name = identify_reads.name
        Map[String, String] qualities =  fastp.output
    }

    if (defined(contigs)){
        call identify_contigs { input:contigs=contigs }
        String contigs_name = identify_contigs.name
    }

    File query = select_first([contigs, reads])
    String name = select_first([user_provided_identifier, contigs_name, reads_name])

    scatter (reference in references) {
        call uncompress { input:input=reference }
        call minimap2.Indexing { 
            input:
                outputPrefix=name, 
                referenceFile=uncompress.output 
        }
        call minimap2.Mapping { 
            input:
                referenceFile=Indexing.indexFile,
                outputPrefix=name,
                queryFile=query
        }
    }

    call callFromAlignments {input:alignments=Mapping.alignmentFile}

}

task uncompress {

}

task fastp {

}

task callFromAlignments {

}

