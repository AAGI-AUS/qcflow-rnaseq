#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
line="=".multiply(100)
ver="qcflow-rnaseq v0.0.5"

// General params
params.help                  = null
params.output_dir            = "results"
params.genome                = null
params.genes                 = null
params.data_dir              = "data"

// system specific parameters
params.max_cpus              = 1
params.max_memory            = 1

config_profile_description   = null 
config_profile_contact       = null
config_profile_url           = null

params.workflow              = "all"
params.aligner 	             = "star"

// STAR
params.sjOverhang            = 149

// Fastp params
params.adapters              = "data/truseq_adapters.fasta"
params.qual_phred            = 20
params.min_read_length       = 50

//--------------------------------------------------------------------------------------------------------
// Replace values

workflow_input         = params.workflow
aligner                = params.aligner

sjOverhang             = params.sjOverhang

genome                 = params.genome
genes                  = params.genes
output_dir             = params.output_dir
data_dir               = params.data_dir

adapters               = params.adapters
qual_phred             = params.qual_phred
min_len                = params.min_read_length

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

if (params.help) {
   log.info paramsHelp("nextflow run main.nf --outdir dir_name --inputdir_fastq dir_name/*_R{1,2}.fastq.gz")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

switch (workflow_input) {
    case ["genome-index"]:
        include { run_star_index; run_hisat_index; run_hisat_index_high_mem } from './modules/module_prep_index.nf'
        aligner = params.aligner
	genome = params.genome
        genes = params.genes
	output_dir = params.output_dir
        break;
    case ["read-trimming"]:
	include { process_reads_fastp; run_fastqc; run_multiqc } from './modules/module_read_trimming.nf'
	adapters = params.adapters
	qual_phred = params.qual_phred
	min_len = params.min_read_length
	break;
     case ["reads-qc"]:
	include { run_fastqc; run_multiqc } from './modules/module_read_qc.nf'
	data_dir = params.data_dir
        output_dir = params.output_dir
        samples = Channel.fromFilePairs("${data_dir}", type: 'file')
                    .ifEmpty { exit 1, data_dir }
        break;
}


workflow GENOME_INDEXING_STAR {
    main:
    run_star_index()
}

workflow GENOME_INDEXING_HISATHIGHMEM {
    main:
    run_hisat_index_high_mem()
}

workflow READ_QC {
    take:
    samples
    
    main:
    run_fastqc(samples)
    run_fastqc.out.qc_html
        .map { it -> it[1] }
        .collect()
        .set { fastqc_out }
    run_multiqc(fastqc_out)
}


workflow {
	switchVariable = 0
	
	if (workflow_input == "genome-index" && aligner == "star") {
    		switchVariable = 1;
	} else if (workflow_input == "genome-index" && aligner == "hisat") {
   		switchVariable = 2;
	} else if (workflow_input == "genome-index" && aligner == "hisat-highmem") {
    		switchVariable = 3;
	} else if (workflow_input == "reads-qc") {
		switchVariable = 4
	}

	switch (switchVariable) {
    	case 1:
        	GENOME_INDEXING_STAR();
        	break;
    	case 2:
        	GENOME_INDEXING_HISAT();
        	break;
    	case 3:
        	GENOME_INDEXING_HISATHIGHMEM();
        	break;
	case 4:
		READ_QC(samples);
		break;
    	default:
        	println("Please provide the correct input options")
		break;
}
}
