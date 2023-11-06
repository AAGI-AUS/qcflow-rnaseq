#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
line="=".multiply(100)
ver="qcflow-rnaseq v0.2"

params.help                  = null
params.genome                = null
params.genes                 = null

params.max_cpus              = 1
params.max_memory            = 1

config_profile_description   = null 
config_profile_contact       = null
config_profile_url           = null

params.workflow              = "all"

//STAR
params.sjOverhang            = 149
//params.index_dir           = "star_index"

//HISAT2
//params.exons_tsv             = null
//params.splicesites_tsv       = null 

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

if (params.help) {
   log.info paramsHelp("nextflow run main.nf --outdir dir_name --inputdir_fastq dir_name/*_R{1,2}.fastq.gz")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

workflow_input         = params.workflow
sjOverhang             = params.sjOverhang
index_dir              = params.index_dir
splicesites            = params.splicesites_tsv
exons                  = params.exons_tsv

println(workflow_input)

switch (workflow_input) {
    case ["genome-index"]:
        include { run_star_index; run_hisat_index } from './modules/module_prep_index.nf'
        genome = params.genome
        genes = params.genes
	//exons = params.exons_tsv
	//splicesites = params.splicesites_tsv
        //index_dir = params.index_dir
        break
}


workflow GENOME_INDEXING {
    main:
    //run_star_index()
    run_hisat_index()
}

workflow {
    switch (workflow_input) {
        case ["genome-index"]:
            GENOME_INDEXING()
            break
     }
}
