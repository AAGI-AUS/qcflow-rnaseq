#!/usr/bin/env nextflow

nextflow.enable.dsl=2

genes            = params.genes
outdir           = params.output_dir

process convert_bed {
    tag {"convert gtf to bed"}
    publishDir "${outdir}/align_rseqc/", mode: 'copy', overwrite: true

    input:
    path genes
    
    output:
    path("${genes.simpleName}.bed12"), emit: bedfile

    script:
    """
    $PWD/bin/gtf2bed.py $genes | awk -F'\t' 'NF==12 {print}' > ${genes.simpleName}.bed12
    """
}


process run_infer_experiment {
    tag {"rnaseqc - infer experiment: ${sample_id}"}
    publishDir "${outdir}/align_rseqc/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)
    path(genes)

    output:
    tuple val(sample_id), path("infer_experiment/${sample_id}_infer-experiment.out")

    script:
    """
    mkdir infer_experiment

    infer_experiment.py -i $bam -r $genes > infer_experiment/${sample_id}_infer-experiment.out
    """
}

process run_bam_stats {
    
    tag {"rnaseqc - bam stats: ${sample_id}"}
    publishDir "${outdir}/align_rseqc/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)
    path(genes)

    output:
    tuple val(sample_id), path("bam_stat/${sample_id}_bam-stats.out")

    script:
    """
    mkdir bam_stat

    bam_stat.py -i $bam > bam_stat/${sample_id}_bam-stats.out
    """
}


process run_junction_annotation {

    tag {"rnaseqc - junction annotation: ${sample_id}"}
    publishDir "${outdir}/align_rseqc/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)
    path(genes)

    output:
    tuple val(sample_id), path("junction_annotation/${sample_id}/")

    script:
    """
    mkdir -p junction_annotation/${sample_id}

    junction_annotation.py -i $bam -o junction_annotation/${sample_id}/${sample_id} -r $genes
    """
}
