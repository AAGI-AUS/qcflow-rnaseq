#!/usr/bin/env nextflow

nextflow.enable.dsl=2

index_dir        = params.index_dir
output_dir       = params.output_dir
hisat_prefix     = params.hisat_prefix
algner           = params.aligner

strandedness     = params.strandedness 
sjOverhang       = params.sjOverhang

process run_star_align {
     
     label 'star'
     tag "Star align reads for ${sample_id}"

     publishDir "${output_dir}/alignements", mode: 'copy'

     output:
        tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam"), emit: alignements
        tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_Log.final.out"), emit: reports
	tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_ReadsPerGene.out.tab"), emit: counts
	tuple val(sample_id), path("star_aligned/${sample_id}/${sample_id}_SJ.out.tab"), emit: splicejunctions

     input:
        tuple val(sample_id), path(reads1), path(reads2)
	val(index_dir)
        val(genes)
     
     script:
	"""
	STAR --runThreadN ${task.cpus} \
	--runMode alignReads \
        --readFilesCommand zcat \
        --sjdbOverhang ${sjOverhang} \
        --outSAMunmapped Within \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outTmpDir star_aligned/${sample_id}/_STARtmp \
        --outFileNamePrefix star_aligned/${sample_id}/${sample_id}_ \
        --genomeDir $index_dir \
        --sjdbGTFfile $genes \
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")}
     """
}


process run_star_align_plants {

     // Function uses specific parameters for large and gappy plant genomes (>3Gbp)
     label 'star'     
     tag "Star align reads for ${sample_id}"
     

     publishDir "${output_dir}/alignements", mode: 'copy'

     input:
        tuple val(sample_id), path(reads1), path(reads2)
	val(index_dir)
	val(genes)

     output:
	tuple val(sample_id), path("star_aligned_plants/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam"), emit: alignements
	tuple val(sample_id), path("star_aligned_plants/${sample_id}/${sample_id}_Log.final.out"), emit: reports
	tuple val(sample_id), path("star_aligned_plants/${sample_id}/${sample_id}_ReadsPerGene.out.tab"), emit: counts
	tuple val(sample_id), path("star_aligned_plants/${sample_id}/${sample_id}_SJ.out.tab"), emit: splicejunctions

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
        --outTmpDir star_aligned_plants/${sample_id}/_STARtmp \
        --outFileNamePrefix star_aligned_plants/${sample_id}/${sample_id}_ \
        --genomeDir $index_dir \
        --sjdbGTFfile $genes \
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")}
        """
}


process run_hisat_align {

     label 'hisat'
     tag "Hisat align reads for ${sample_id}"
     
     publishDir "${output_dir}/alignements", mode: 'copy'

     input:
     tuple val(sample_id), path(reads1), path(reads2)
     val(index_dir)
     val(genes)

     output:
     tuple val(sample_id), path("hisat_aligned/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam"), emit: alignements
     tuple val(sample_id), path("hisat_aligned/${sample_id}/${sample_id}_hisat.summary.log"), emit: reports
     tuple val(sample_id), path("hisat_aligned/${sample_id}/${sample_id}_splicejunctions.txt"), emit: splicejunctions

     script:
        """
        mkdir -p hisat_aligned/${sample_id}
        hisat2_extract_splice_sites.py $genes > splicesites.tsv

        hisat2 \
        -x $index_dir/${hisat_prefix} \
        -1 ${reads1.join(",")} \
        -2 ${reads2.join(",")} \
        --known-splicesite-infile splicesites.tsv \
        --summary-file hisat_aligned/${sample_id}/${sample_id}_hisat.summary.log \
        --rna-strandness ${strandedness} \
        --dta \
        --threads ${task.cpus} \
        --novel-splicesite-outfile hisat_aligned/${sample_id}/${sample_id}_splicejunctions.txt \
        | samtools view -bS -F 4 -F 256 - | samtools sort - -o hisat_aligned/${sample_id}/${sample_id}_Aligned.sortedByCoord.out.bam
        """
}

process run_multiqc_align {
    tag { 'multiqc run' }
    publishDir "${output_dir}", mode: 'copy', overwrite: true

    input:
    path(dir)
    val(aligner)

    output:
    path("multi_qc-${aligner}"), emit: multiQC
    
    """
    multiqc . --force --outdir multi_qc-${aligner}
    """
}
