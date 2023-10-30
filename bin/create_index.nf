/*
 * Index genomes
 */

params.input_dir     = "${baseDir}/data/*"
params.sjOverhang  = 149 //read length -1

/*
* Create the genome index file for STAR
*/

if( params.help ) {

log.info """
Indexing - supports STAR aligner
=============================================
Usage:
    nextflow run -resume -profile pawsey_nimbus,c16r64 create_index.nf \
	-with-singularity docker://0000000000777/reads_mapping_rnaseq:v0.0.2 \
	--input_dir "dir_name/*" \
	--sjOverhang 149 
Input:
    * --input_dir: Path to directory with fasta and gtf files for indexing [${params.input_dir}]
    * --sjOverhang: Length of sj overhang - "sjdbOverhang" parameter (see STAR manual for more details). Default [${params.sjOverhang}]
"""
    exit 0
}


process prepare_star_genome_index {
    
    cpus 16

    input:
    tuple val(genome_id), path(genome_files)

    output:
    path genome_id

    script:
    """
    mkdir ${genome_id}

    STAR --runMode genomeGenerate \\
         --genomeDir ${genome_id} \\
         --genomeFastaFiles ${genome_files[0]} \\
         --runThreadN ${task.cpus} \\
         --sjdbGTFfile ${genome_files[1]} \\
         --sjdbOverhang ${params.sjOverhang}
    """
}

workflow {

    Channel
        .fromPath(params.input_dir)
        .map { file -> tuple(file.baseName, file) }
        .groupTuple()
	.set {file_pairs}

    prepare_star_genome_index(file_pairs)
}
