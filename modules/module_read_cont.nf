#!/usr/bin/env nextflow
nextflow.enable.dsl=2

outdir           = params.output_dir

process run_cont {

    label 'bbt-cont'
    tag { "bbt-categorizer: ${sample}" }
    
    input:
    tuple val(sample_id), path(reads)
    val(bbt_filters)
    
    output:
    path("{sample_id}_count_BBT.txt"), emit: cont_counts
    
    script:
    """
    biobloomcategorizer -p cont_BBT \
	-t ${tasks.cpu} \
        -e \
        -i \
        -f ${bbt_filters} ${reads[0]} ${reads[1]}
    """
}
