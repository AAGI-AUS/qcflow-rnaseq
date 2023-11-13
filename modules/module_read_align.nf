#!/usr/bin/env nextflow

nextflow.enable.dsl=2

genome           = params.genome
genes            = params.genes
index_dir        = params.index_dir
output_dir       = params.output_dir
workflow         = params.workflow


sjOverhang       = params.sjOverhang

process run_star_align_plants {

     // Function uses specific parameters for large and gappy plant genomes (>3Gbp)
     label 'star'     
     tag "Star align reads for ${sample_id}"

     publishDir "${output_dir}/alignements", mode: 'copy'

     input:
        tuple val(sample_id), path(reads1), path(reads2)

     output:
	tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam"), emit: star_alignments
    	tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_Log.final.out"), emit: star_reports
	tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_ReadsPerGene.out.tab"), emit: star_counts

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

process run_hisat {

     label 'hisat_highmem'
     tag "Star align reads for ${sample_id}"
     publishDir "${output_dir}/alignements", mode: 'copy'

     input:
     tuple val(sample_id), path(reads1), path(reads2)
     
     output:
     tuple val(meta), path("hisat_aligned/${sample_id}/${sample_id}_Aligned.sortedByCoord.bam"), emit: bam
     tuple val(meta), path("hisat_aligned/${sample_id}/${sample_id}.hisat2.summary.log"), emit: summary
     tuple val(meta), path("hisat_aligned/${sample_id}/*splicesite.txt"), emit: splicesites

     script:
     """
     hisat2 \
     	-x ${index_dir} \
     	-1 ${reads1.join(",")} \
     	-2 ${reads2.join(",")} \\
     	--summary-file ${sample_id}.hisat2.summary.log \\
     	--threads ${task.cpus} \
     	| samtools sort --threads ${task.cpus} -o ${sample_id}/${sample_id}_Aligned.sortedByCoord.bam -
     """
}



process run_multiqc {
    tag { 'multiqc run' }
    publishDir "${output_dir}/", mode: 'copy', overwrite: true

    input:
    path(dir)
    
    """
    multiqc . --force --outdir qc_align
    """
}
