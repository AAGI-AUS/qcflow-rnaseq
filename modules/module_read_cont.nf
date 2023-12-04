#!/usr/bin/env nextflow
nextflow.enable.dsl=2

outdir           = params.output_dir

process combine_cont_bbt {

    tag {"combine cont bbt"}
    publishDir "${outdir}/bbt-contamination/", mode: 'copy', overwrite: true

    input:
    path(rates)

    output:
    path("ContaminationRate.tsv"), emit: rates_cont

    script:
    """
    $PWD/bin/CombineCont.py --input $rates --output ContaminationRate.tsv
    """

}


process run_cont {

    label 'bbt'
    tag { "bbt-categorizer: ${sample}" }
    
    publishDir "${outdir}/bbt-contamination", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)
    val(bbt_filters)
    
    output:
    tuple val(sample_id), path("${sample_id}_cont_BBT_summary.tsv"), emit: cont_out
    
    script:
    """
    biobloomcategorizer -p ${sample_id}_cont_BBT \
	-t ${task.cpus} \
        -e \
        -i \
        -f "${bbt_filters.join(" ")}" ${reads[0]} ${reads[1]}
    """
}
