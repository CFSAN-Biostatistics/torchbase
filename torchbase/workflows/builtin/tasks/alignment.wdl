version 1.0

# Minimap2 alignment task for allele sequence alignment
# Supports multiple preset configurations for different accuracy/speed tradeoffs
#
# Available presets:
#   - asm20: Fast assembly-to-assembly, suitable for quick screening (fast strategy)
#   - asm5: Balanced assembly-to-assembly with moderate accuracy (balanced strategy)
#   - asm5+eqx: Sensitive assembly-to-assembly with equivalent-to operations (sensitive strategy)
#
# The preset parameter can be adjusted per invocation to support different typing strategies

task align_sequences {
    input {
        File query
        File reference
        String preset = "asm5"
    }

    command <<<
        set -e

        # Create alignment output directory and files
        minimap2 \
            -a \
            -x ~{preset} \
            ~{reference} \
            ~{query} \
            > alignment.sam

        # Convert SAM to sorted BAM for further processing if needed
        samtools view -bh alignment.sam | samtools sort -o alignment.sorted.bam -
        samtools index alignment.sorted.bam
    >>>

    output {
        File alignment = "alignment.sam"
        File bam = "alignment.sorted.bam"
        File bam_index = "alignment.sorted.bam.bai"
    }

    runtime {
        docker: "quay.io/biocontainers/minimap2:2.28--he4a0461_0"
        cpu: 4
        memory: "8 GB"
    }
}
