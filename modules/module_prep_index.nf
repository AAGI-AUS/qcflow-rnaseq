#!/usr/bin/env nextflow

nextflow.enable.dsl=2

genome           = params.genome
genes            = params.genes
index_dir        = file(genome, type: 'file', checkIfExists: true).getParent()

//STAR
sjOverhang       = params.sjOverhang

//HISAT
//exons            = params.exons_tsv
//splicesites      = params.splicesites_tsv

cpus             = params.max_cpus
max_memory       = params.max_memory

process run_star_index {

    label 'star'
    tag { "star: index" }
    memory = max_memory
    publishDir "${index_dir}", mode: 'copy', overwrite: true
    
    output:
    tuple val("starIndex"), path("*"), emit: star_index
    
    script:
    """
    STAR --runMode genomeGenerate \
         --genomeDir star_index \
	 --genomeFastaFiles ${genome} \
         --runThreadN ${cpus} \
         --sjdbGTFfile ${genes} \
         --sjdbOverhang ${sjOverhang}
    """
}

process run_hisat_index {
    label 'hisat2'
    tag { "hisat2: index" }
    publishDir "${index_dir}", mode: 'copy', overwrite: true  

    output:
    path("hisat_index"), emit: hisat_index
    
    script:
    """
    mkdir hisat_index

    hisat2_extract_exons.py ${genes} > exons.tsv
    hisat2_extract_splice_sites.py ${genes} > splicesites.tsv
    hisat2-build -p ${cpus} --ss splicesites.tsv --exon exons.tsv ${genome} hisat_index/hisat_index
    """
}
