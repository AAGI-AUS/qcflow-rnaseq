#!/usr/bin/env nextflow

nextflow.enable.dsl=2

genome           = params.genome
genes            = params.genes
index_dir        = params.index_dir
output_dir       = params.output_dir

sjOverhang       = params.sjOverhang

process run_star_align {
     label 'star'     
     tag "Star align reads for ${sample_id}"

     publishDir "${output_dir}/alignment_star", mode: 'copy'

     input:
        tuple val(sample_id), path(reads1), path(reads2)

     output:
	tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_Aligned.out.bam"), emit: star_alignments
    	tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_Log.final.out"), emit: star_reports

script:
        """
        STAR --runThreadN ${task.cpus} \
        --runMode alignReads \
        --readFilesCommand zcat \
        --sjdbScore 2 \
        --sjdbOverhang ${sjOverhang} \
        --limitSjdbInsertNsj 1000000 \
        --outFilterMultimapNmax 10 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMunmapped Within \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --outTmpDir star_aligned/${sample_id}/_STARtmp \
        --outFileNamePrefix star_aligned/${sample_id}/${sample_id}_ \
        --genomeDir ${index_dir} \
        --sjdbGTFfile ${genes} \
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")}
        """
}

process run_multiqc_star {
    tag { 'multiqc: star' }
    publishDir "${output_dir}/alignment_star", mode: 'copy', overwrite: false

    input:
    path(dir) 
    
    output:
    path("qc_star"), emit: multiqc_report
    
    """
    multiqc . --force --outdir qc_star
    """
}
