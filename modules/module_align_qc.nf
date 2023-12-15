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

process run_infer_experiment_rseqc {
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

process run_check_strandedness {

    tag { "check_strandedness: ${sample_id}"  }    
    publishDir "${outdir}/check_strandedness/", mode: 'copy', overwrite: true

    input:
    val(genes)
    val(cdna)
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_check-strandedness.out")

    script:
    """
    check_strandedness -g $genes -fa $cdna -r1 ${reads[0]} -r2 ${reads[1]} > ${sample_id}_check-strandedness.out
    """

}

process run_bam_stats {
    
    tag {"rnaseqc - bam stats: ${sample_id}"}
    publishDir "${outdir}/align_rseqc/", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)
    path(genes)

    output:
    tuple val(sample_id), path("bam_stat/${sample_id}_bam-stats.out"), emit: bam_stats

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

process combine_bam_stats {

    tag {"combine bam stats"}
    publishDir "${outdir}/align_rseqc/", mode: 'copy', overwrite: true

    input:
    path(bam_stats)

    output:
    path("StatsBam.tsv"), emit: bamstats

    script:
    """
    $PWD/bin/Combine-bam-stats.py --input $bam_stats --output StatsBam.tsv
    """

}
