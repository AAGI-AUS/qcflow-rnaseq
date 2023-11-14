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
params.fastq_dir             = "data"

// system specific parameters
params.max_cpus              = 1
params.max_memory            = 1

config_profile_description   = null 
config_profile_contact       = null
config_profile_url           = null

params.workflow              = "all"
params.aligner	             = "star"

// STAR index
params.sjOverhang            = 149
params.snps                  = null

// STAR align 
params.library_name          = 0
params.index_dir             = null

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
library_name           = params.library_name
index_dir              = params.index_dir

genome                 = params.genome
genes                  = params.genes
output_dir             = params.output_dir
fastq_dir              = params.fastq_dir

adapters               = params.adapters
qual_phred             = params.qual_phred
min_len                = params.min_read_length

//--------------------------------------------------------------------------------------------------------
// Validation

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

def extractCharacters(String inputString, library_name) {

    def intervals = library_name.split(',').collect { it.toInteger() }

    def adjustedStart
    def adjustedEnd

    if (intervals.size() == 1) {
        // Single value provided, create 'end'
        adjustedStart = 0
        adjustedEnd = inputString.length()
    } else if (intervals.size() == 2) {
        // Two values provided, use them as 'start' and 'end'
        adjustedStart = intervals[0] - 1
        adjustedEnd = intervals[1]
    } else {
        // Invalid input, handle appropriately (e.g., raise an error)
        throw new IllegalArgumentException("Invalid input provided: $library_name")
    }
    return inputString.substring(adjustedStart, adjustedEnd)
}

def is_null(var) {
  if (var == null) {
    return true
  } else {
    return false
  }
}

def selectTool(inputParameter) {
   selectedTool = ""
   if (inputParameter in ["star", "star-plants", "star-snps"]) {
        selectedTool = "star"
   } else if (inputParameter in ["hisat", "hisat_highmem"]) {
        selectedTool = "hisat"
   } else {
        selectedTool = "reads"
   }
   return selectedTool
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
    case ["trim"]:
	include { run_fastp; run_fastqc; run_multiqc } from './modules/module_read_trimming.nf'
	//adapters = params.adapters
	//qual_phred = params.qual_phred
	//min_len = params.min_read_length
	adapters = Channel.FromPath(params.adapters)
	samples = Channel.fromFilePairs("${fastq_dir}", type: 'file')
                    .ifEmpty { exit 1, fastq_dir }
	break;
     case ["reads-qc"]:
	include { run_fastqc; run_multiqc_reads } from './modules/module_read_qc.nf'
	fastq_dir = params.fastq_dir
        //output_dir = params.output_dir
        samples = Channel.fromFilePairs("${fastq_dir}", type: 'file')
                    .ifEmpty { exit 1, fastq_dir }
        break;
     case ["align"]:
	include { run_star_align_plants; run_hisat_align; run_multiqc } from './modules/module_read_align.nf'
	//include { convert_bed; run_bam_stats; run_infer_experiment; run_junction_annotation } from './modules/module_align_qc.nf'
	//genes = params.genes
	fastq_dir = params.fastq_dir
	//index_dir = params.index_dir
	//output_dir = params.output_dir
	//sjOverhang = params.sjOverhang
        library_name = params.library_name
	genes = Channel.fromPath(params.genes)
	samples = Channel.fromFilePairs("${fastq_dir}", flat: true)
                .map { prefix, file1, file2 -> tuple(extractCharacters(prefix,"${library_name}"), file1, file2) }
                .groupTuple(sort: true)
                .ifEmpty { exit 1, fastq_dir }
	break;
}


workflow GENOME_INDEX {
    main:

    if (aligner == "hisat-highmem") {
	run_hisat_index_high_mem()
     } else if (aligner == "hisat") {
	run_hisat_index_()
     } else if (aligner == "star")  {
	run_star_index()
     } else if (aligner == "star-snps") {
	run_star_index_snps()
     }
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
    run_multiqc_reads(fastqc_out)
}

workflow TRIM_READS {
    take:
    samples
    adapters

    main:
    fastp_out = run_fastp(samples, adapters)
    fastqc_trimmed_out = run_fastqc(fastp_out.trimmed_reads)
    run_multiqc(fastp_out.json.mix(fastqc_trimmed_out).collect())
}

workflow ALIGN_READS {
     take:
     samples
     aligner
     genes

     main:
     Channel.fromPath( params.index_dir ).set{ index }

     if (aligner == "star-plants") {
	output_align = run_star_align_plants(samples, index, genes)

     } else if (aligner == "hisat") {
	output_align = run_hisat_align(samples, index, genes)
     }
     
     output_align.reports
        .map { it -> it[1]}
        .flatten()
        .collect()
        .set { reports }
     convert_bed(genes)
     // QC stages
     run_bam_stats(output_align.alignements, convert_bed.out)
     run_junction_annotation(output_align.alignements, convert_bed.out)
     run_infer_experiment(output_align.alignements, convert_bed.out)
     run_multiqc(reports, selectTool(params.aligner))
}


workflow {
	switchVariable = 0
	
	if (workflow_input == "genome-index") {
		switchVariable = 1;
	} else if (workflow_input == "reads-qc") {
		switchVariable = 2;
	} else if (workflow_input == "trim") {
		switchVariable = 3;
	} else if (workflow_input == "align") {
		switchVariable = 4;
	}
	switch (switchVariable) {
	case 1:
		GENOME_INDEX();
		break;
	case 2:
		READ_QC(samples);
		break;
	case 3:
		TRIM_READS(samples, adapters)
		break;
	case 4:
		ALIGN_READS(samples, aligner, genes)
		break;
	default:
		println("Please provide the correct input options")
		break;
	}
}
