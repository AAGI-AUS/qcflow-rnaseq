#!/usr/bin/env nextflow
nextflow.enable.dsl=2

genes           = params.genes

process run_rnaseqc {

    tag { "rnaseqc: ${sample_id}" }
    publishDir "${outdir}/align-qc/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.metrics.tsv"), emit: rnaseqc_matrix

    """
    rnaseqc ${genes} ${bam} ${sample_id} --sample=${sample_id}
    """
}
