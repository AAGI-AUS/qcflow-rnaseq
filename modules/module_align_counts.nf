#!/usr/bin/env nextflow
nextflow.enable.dsl=2

genome             = params.genome
genes              = params.genes
strandedness       = params.strandedness
outdir             = params.output_dir

process combine_counts_star {
    
    tag {"combine counts"}
    publishDir "${outdir}/counts/", mode: 'copy', overwrite: true

    input:
    path(counts)

    output:
    path("RawCountsStar.tsv"), emit: rawcounts

    script:
    """
    CombineCounts-star.py --input $counts --strandedness $strandedness --output RawCountsStar.tsv
    """
}

process combine_counts_featurecounts {

    tag {"combine counts"}
    publishDir "${outdir}/counts/", mode: 'copy', overwrite: true

    input:
    path(counts)

    output:
    path("RawCountsFeatureCounts.tsv"), emit: counts

    script:
    """
    CombineCounts-featurecounts.py --input $counts --output RawCountsFeatureCounts.tsv
    """

}

process run_htseq_counts {

    label 'htseqcount'
    tag { "htseqcount: ${sample_id}" }
    publishDir "${outdir}/counts/htseqCounts", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam)
    path(genes)
    
    output:
    tuple val(sample_id), path("${sample_id}_counts.txt"), emit: htseqcounts_raw
    
    script:
    def stranded = (strandedness == "RF") ? "reverse" : ((strandedness == "FR") ? "yes" : "no")
    """
    htseq-count --format=bam \
        --order=name \
        --idattr=gene_id \
        --stranded=${stranded} \
        --mode=union \
        --type=transcript \
        $bam $genes > ${sample_id}_counts.txt
    """
}

process run_feature_counts {
    
    label 'featurecounts'
    tag { "featureCounts: ${sample_id}" }
    publishDir "${outdir}/counts/featureCounts/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam)
    path(genome)
    path(genes)

    output:
    tuple val(sample_id), path("${sample_id}_counts"), emit: counts
    tuple val(sample_id), path("${sample_id}_counts.summary"), emit: summary
    tuple val(sample_id), path("${sample_id}_counts.jcounts"), emit: jcounts

    script:
    def stranded = (strandedness == "RF") ? 2 : ((strandedness == "FR") ? 1 : 0)
    """
    featureCounts -p -B	-C -P -J -s ${stranded} \
        -G ${genome} \
	-J \
        -t exon \
        -g gene_id \
        -a ${genes} \
        -T ${task.cpus} \
        -o ${sample_id}_counts \
        $bam
    """
}
