#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input
params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container "community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867"

    publishDir 'results', mode: 'copy'

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    """
    samtools index '$input_bam'
    """

}

workflow {

    // Create input channel

    // Create index file for input BAM file
    
}
