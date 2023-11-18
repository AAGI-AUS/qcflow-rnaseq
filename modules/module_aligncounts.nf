#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
    $PWD/bin/CombineCounts-star.py --input $counts --strandedness $strandedness > RawCountsStar.tsv
    """

}
