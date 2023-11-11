#!/usr/bin/env nextflow

nextflow.enable.dsl=2

genes            = params.genes
outdir           = params.output_dir

process convert_bed {
    tag {"convert gtf to bed"}
    publishDir "${outdir}/align-qc/", mode: 'copy', overwrite: true

    input:
    path genes
    
    output:
    path("${genes.simpleName}.bed12")

    script:
    """
    $PWD/bin/gtf2bed.py $genes | awk -F'\t' 'NF==12 {print}' > ${genes.simpleName}.bed12
    """
}

process run_rnaseqc {

    tag {"rnaseqc: ${sample_id}"}
    publishDir "${outdir}/align-qc/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)
    path(genes)

    output:
    tuple val(sample_id), path("infer_experiment/${sample_id}_infer-experiment.out")
    tuple val(sample_id), path("bam_stat/${sample_id}_bam-stats.out")
    tuple val(sample_id), path("junctions_annotation/${sample_id}_junction-annotation.r")

    script:
    """
    mkdir infer_experiment
    mkdir bam_stat
    mkdir junction_annotation

    infer_experiment.py -i $bam -r $genes > infer_experiment/${sample_id}_infer-experiment.out
    bam_stat.py -i $bam > bam_stat/${sample_id}_bam-stats.out
    junction_annotation.py -i $bam -o junction_annotation/${sample_id}_junction-annotation -r $genes
    """
}
