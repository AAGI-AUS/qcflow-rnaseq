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

/*

fq1=$1
fq2=$2

echo $fq2

cont_p=../../../analysis/preliminary_reads/contamination/cont_filters/contamination_plants

biobloomcategorizer -p cont_BBT -t 24 \
	-e \
	-i \
	-f "$cont_p/aphids/aphids.bf $cont_p/mites/mites.bf $cont_p/fungi/fungi.bf $cont_p/archaea/archaea.bf $cont_p/bacteria/bacteria.bf $cont_p/protozoa/protozoa.bf $cont_p/thrips/thrips.bf $cont_p/univec/univec.bf $cont_p/viral/viral.bf" \
	$fq1 $fq2
*/
