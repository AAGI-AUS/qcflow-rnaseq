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

// STAR index
params.sjOverhang            = 149
params.snps                  = null

// Fastp params
params.adapters              = "data/truseq_adapters.fasta"
params.qual_phred            = 20
params.min_read_length       = 50

//--------------------------------------------------------------------------------------------------------
// Replace values

workflow_input         = params.workflow
aligner                = params.aligner

sjOverhang             = params.sjOverhang
snp                    = params.snps

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

//-----------------------------------------------------------------------------------------------------------------------
// Functions

def getLibraryId( prefix ){
  // Return the ID number, you can change for other file formats, here it just takes the first part before "_"
  prefix.split("_")[0]
}

def is_null(var) {
  if (var == null) {
    return true
  } else {
    return false
  }
}


switch (workflow_input) {
    case ["genome-index"]:
        include { run_star_index; run_star_index_snps; run_hisat_index; run_hisat_index_high_mem } from './modules/module_prep_index.nf'
        aligner = params.aligner
	genome = params.genome
        genes = params.genes
	output_dir = params.output_dir
	snp = params.snps
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
     case ["align"]:
	include { run_star_align; run_multiqc } from './modules/module_read_align.nf'
	include { convert_bed; run_rnaseqc } from './modules/module_align_qc.nf'
	genes = params.genes
	data_dir = params.data_dir
	index_dir = params.index_dir
	output_dir = params.output_dir
	sjOverhang = params.sjOverhang
	samples = Channel.fromFilePairs("${data_dir}", flat: true)
                .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
                .groupTuple(sort: true)
		.ifEmpty { exit 1, data_dir }
	break;
}


workflow GENOME_INDEXING_STAR {
    main:
    run_star_index()
}

workflow GENOME_INDEXING_STAR_SNPS {
    main:
    run_star_index_snps()
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

workflow ALIGN_STAR {
     take:
     samples

     main:
     run_star_align(samples)
     run_star_align.out.star_reports
        .map { it -> it[1]}
        .flatten()
        .collect()
        .set { reports }
     convert_bed(genes)
     run_rnaseqc(run_star_align.out.star_alignments, genes)
     run_multiqc(reports)
}

workflow {
	switchVariable = 0
	
	if (workflow_input == "genome-index" && aligner == "star") {
    		switchVariable = 1;
	} else if (workflow_input == "genome-index" && aligner == "star-snps" && !is_null(snps)){
                switchVariable = 2;
	} else if (workflow_input == "genome-index" && aligner == "hisat") {
   		switchVariable = 3;
	} else if (workflow_input == "genome-index" && aligner == "hisat-highmem") {
    		switchVariable = 4;
	} else if (workflow_input == "reads-qc") {
		switchVariable = 5;
	} else if (workflow_input == "align" && (aligner == "star" || aligner == "star-snps")) {
		switchVariable = 6;
	}

	switch (switchVariable) {
    	case 1:
        	GENOME_INDEXING_STAR();
        	break;
    	case 2:
		GENOME_INDEXING_STAR_SNPS();
		break;
	case 3:
        	GENOME_INDEXING_HISAT();
        	break;
    	case 4:
        	GENOME_INDEXING_HISATHIGHMEM();
        	break;
	case 5:
		READ_QC(samples);
		break;
	case 6:
		ALIGN_STAR(samples)
		break;
    	default:
        	println("Please provide the correct input options")
		break;
	}
}
