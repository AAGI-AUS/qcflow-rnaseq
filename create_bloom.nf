#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
line="=".multiply(100)
ver="qcflow-rnaseq v1.0.0"

// General params
help                  = params.help          
outdir                = params.output_dir 
fasta_dir             = params.input_fasta

//-------------------------------------------------------------------------------------------------------
// Workflow

samples = Channel.fromPath("${fasta_dir}")
	.map { tuple( it.getBaseName(), it ) }

samples.view()

process create_bf {
   
   label "bbt"
   
   tag {"bbt-make: ${sample_id}"}
   publishDir "${outdir}/biobloom-filters/", mode: 'copy', overwrite: true

   input:
   tuple val(sample_id), path(fasta_file)

   output:
   path("${sample_id}.bf")
   path("${sample_id}.txt")

   script:
   """
   biobloommaker -d -p ${sample_id} -f 0.001 ${fasta_file}
   """
}

workflow {

    create_bf(samples)

}
