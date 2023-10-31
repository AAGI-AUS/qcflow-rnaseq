// qcflow-rnaseq pipeline - an updated and faster pipeline based on qcflow

/*
NXF ver 23.03 needed because DSL2 language properties
*/

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.outdir = "results"
params.inputdir_fastq = "input/*_R{1,2}.fastq.gz"

params.multi_qc_conf = "conf/multiqc/multiqc_config.yml"
params.min_read_length = 50
params.qual_phred = 20
params.adapters = "data/truseq_adapters.fasta"

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run main.nf --outdir dir_name --inputdir_fastq dir_name/*_R{1,2}.fastq.gz")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// END OF VALIDATION

// fastp trimmed files are published, json are only sent in the channel and used by multiqc
process fastp {

    tag "filter $sample_id"
    //echo true
    publishDir params.outdir, mode: 'copy', pattern: 'trimmed_reads/*fastq.gz', overwrite: false  // publish trimmed fastq files
    
    cpus 16

    input:
	tuple val(sample_id), path(reads)
	path "truseq_adapters.fasta"
    
    output:
	tuple val(sample_id), path("trimmed_reads/${sample_id}_filt_R*.fastq.gz"), emit: trimmed_reads
	path("trimmed_reads/${sample_id}.fastp.json"), emit: json

    script:
    """
    mkdir trimmed_reads
    fastp -i ${reads[0]} -I ${reads[1]} \\
      -o 'trimmed_reads/${sample_id}_filt_R1.fastq.gz' -O 'trimmed_reads/${sample_id}_filt_R2.fastq.gz' \\
      -q $params.qual_phred \\
      -l $params.min_read_length \\
      --adapter_fasta "truseq_adapters.fasta" \\
      -w ${task.cpus} \\
      -j 'trimmed_reads/${sample_id}.fastp.json'
      2>&1 | tee > fastp.log
    """
}

process fastqc {
    
    tag "Fastqc on $sample_id"

    cpus 16

    input:
	tuple val(sample_id), path(reads)

    output:
	path "fastqc/fastqc_${sample_id}_logs"
    
    script:
    """
    mkdir -p 'fastqc/fastqc_${sample_id}_logs'
    fastqc -t ${task.cpus} \\
	-o 'fastqc/fastqc_${sample_id}_logs' \\
	-f fastq -q ${reads[0]} ${reads[1]}
    """
}

process multiqc {
    
    publishDir params.outdir, mode: 'copy'

    input:
        path x
        //path "multiqc.config"

    output:
        file("multiqc/multiqc_report.html")

    script:
    """
    multiqc ${x} --filename "multiqc/multiqc_report.html" 
    """
} 


workflow {
	Channel.fromFilePairs(params.inputdir_fastq, checkIfExists: true, size: 2).set {fastqPairs_ch}
	Channel.fromPath(params.adapters, checkIfExists: true, type: "file").set {adapters}
	fastp_out = fastp(fastqPairs_ch, adapters)
	//fastqc_raw_out = fastqc(fastqPairs_ch)
	fastqc_trimmed_out = fastqc(fastp_out.trimmed_reads)
	multiqc(fastp_out.json.mix(fastqc_trimmed_out).collect())
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc/multiqc_report.html\n" : "Oops .. something went wrong" )
}
