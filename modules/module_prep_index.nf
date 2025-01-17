#!/usr/bin/env nextflow

nextflow.enable.dsl=2

genome           = params.genome
genes            = params.genes
index_dir        = file(genome, type: 'file', checkIfExists: true).getParent()
output_dir       = params.output_dir
input_snps       = params.snps

//STAR
sjOverhang       = params.sjOverhang

//HISAT

process run_star_index {

    label 'star'
    tag { "star: index" }

    publishDir "${index_dir}/${output_dir}", mode: 'copy', overwrite: true

    input:
    path(genome)
    path(genes)
    val(gen_value)
	    
    output:
    tuple val("starIndex"), path("*"), emit: star_index
    
    script:
    """
    mkdir star_index

    STAR --runMode genomeGenerate \
         --genomeDir star_index \
	 --genomeFastaFiles ${genome} \
         --runThreadN ${task.cpus} \
         --sjdbGTFfile ${genes} \
         --sjdbOverhang ${sjOverhang} \
	 --genomeSAindexNbases ${gen_value}
    """
}

process run_star_index_snps {
    label "star"

    tag { "star_snps: index" }
    publishDir "${index_dir}/${output_dir}", mode: 'copy', overwrite: true

    input:
    path(genome)
    path(genes)
    val(gen_value)

    output:
    tuple val("starIndex"), path("*"), emit: star_index

    script:
    """
    STAR --runMode genomeGenerate \
         --genomeDir star_index_snps \
         --genomeFastaFiles ${genome} \
         --genomeTransformVCF ${input_snps} \
         --genomeTransformType Haploid \
         --runThreadN ${task.cpus} \
         --sjdbGTFfile ${genes} \
         --sjdbOverhang ${sjOverhang} \
	 --genomeSAindexNbases ${gen_value}
    """
}


process run_hisat_index {
    label 'hisat'
    tag { "hisat: index" }
    publishDir "${index_dir}/${output_dir}", mode: 'copy', overwrite: true  

    output:
    path("hisat_index"), emit: hisat_index

    input:
    path(genome)
    val(genes)
    
    script:
    """
    mkdir hisat_index
    
    hisat2_extract_exons.py $genes > exons.tsv
    hisat2_extract_splice_sites.py $genes > splicesites.tsv
    hisat2-build -p ${task.cpus} --ss splicesites.tsv --exon exons.tsv $genome hisat_index/hisat_index
    """
}

process run_hisat_index_high_mem {
    label 'hisat_highmem'
    tag { "hisat: index" }

    publishDir "${index_dir}/${output_dir}", mode: 'copy', overwrite: true

    output:
    path("hisat_index"), emit: hisat_index

    input:
    path(genome)
    val(genes)

    script:
    """
    mkdir hisat_index

    hisat2_extract_exons.py $genes > exons.tsv
    hisat2_extract_splice_sites.py $genes > splicesites.tsv
    
    hisat2-build -p ${task.cpus} --large-index --noauto --bmaxdivn 8 --dcv 4096 \
	--ss splicesites.tsv --exon exons.tsv \
	$genome hisat_index/hisat_index
    """
}
