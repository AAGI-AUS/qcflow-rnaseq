// qcflow-rnaseq pipeline - an updated and faster pipeline based on qcflow

/*
NXF ver 23.03 needed because DSL2 language properties
*/

def helpMessage() {
    log.info"""
    # qcflow-rnaseq

    A [nextflow](www.nextflow.io) pipeline for running multiple illumina
    short-read filtering, trimming, and QC steps. QC is performed at each stage.

    ## Usage

    ```bash
    nextflow run -resume darcyabjones/qcflow \
      --fastq "fastq/*_R{1,2}.fastq.gz" \
      --adapters "data/truseq_adapters.fasta" \
      --references "my_genome*.fasta"
    ```

    Note that the quotes around globbing patterns are important for some
    shells that automatically expand globs (e.g. zsh).

    ## Mandatory Arguments

    ```
    param                   | description
    ---------------------------------------------------------------------------
    `--fastq <fastq pairs>` | The fastq pairs to process. This must be
                            | provided as a glob pattern to capture the pairs.
                            | E.G. `sample{1,2}.fastq` will capture
                            | `sample1.fastq sample2.fastq`.
    ```

    ## Options

    ```
    param                      | default          | description
    ---------------------------------------------------------------------------
    `--outdir <path>`          | ./               | The base directory to
                               |                  | publish the results in.

    `--references <*fasta>`    | none             | A glob pattern of reference
                               |                  | genomes to use for mapping
                               |                  | and for masking contaminant
                               |                  | filtering.

    `--adapters <*fasta>`      | data/            | Adapter sequences to trim
                               | truseq_fwd.fasta | from the 5' end.

    `--scontaminants <*fasta>` | data/            | Synthetic contaminants to
                               | synth_cont.fasta | filter from the reads.
                               |                  | This is for things like
                               |                  | PHiX or primer dimer, or
                               |                  | common lab vectors that are
                               |                  | likely contaminants of
                               |                  | sequencing rather than of
                               |                  | the samples themselves.

    `--contaminants <*fasta>   | none             | Filter reads that match
                               |                  | these sequences out.
                               |                  | This is for filtering out
                               |                  | sample contaminants, e.g.
                               |                  | bacteria or endophytes.
                               |                  | If a reference genome is
                               |                  | provided, regions in this
                               |                  | database matching the
                               |                  | reference will be masked.
                               |                  | Generally I would run kraken
                               |                  | to check for contamination
                               |                  | before using this.

    `--filter_phred <int>`     | 5                | Filter out reads that have
                               |                  | lower average phred scores
                               |                  | than this. Keep this low.
                               |                  | bbduk seems to be a bit
                               |                  | "filter happy" for this.

    `--trim_phred <int>`       | 2                | Trim bases with phred
                               |                  | qualities lower than this
                               |                  | off the end of reads.
                               |                  | Generally quality trimming
                               |                  | is not useful unless you
                               |                  | have very poor data.

    `--min_read_length <int>`  | 50               | Filter out reads with
                               |                  | lengths less than this
                               |                  | after trimming.

    `--map`                    | false            | Align the raw reads to the
                               |                  | reference genomes and get
                               |                  | qc stats. Useful for
                               |                  | estimating insert/fragment
                               |                  | size, or error rates.

    `--biobloomtoolsDb <dir>`  | none             | Search for matches to
                               |                  | potential contaminants in
                               |                  | this Biobloomtools mRNA database.
                               |                  | This should be seen as a
                               |                  | first pass, check. If you
                               |                  | see something more
                               |                  | substantial you might want
                               |                  | to do a more accurate
                               |                  | alignment.

    ## Running Biobloomtools

    Biobloomtools is useful as a first pass to detect potential contaminants in your
    sequencing.The reads can be then removed with fastp

    There is a slurm script to download and prepare a Biobloom database in the
    `batch_scripts` directory in the repo.

    ## Outputs

    TBD

    ## Requirements

    * `BBMap` <https://sourceforge.net/projects/bbmap/>.
      Developed with v38.39.
    * `fastqc` <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>.
      Developed with v0.11.8.
    * `fastp` <https://github.com/OpenGene/fastp>.
      Developed with vxx.
      Required by default.
    * `multiqc` <https://multiqc.info/>
      Developed with v1.xx.
    """.stripIndent()
}

params.help = "NA"
if (params.help){
    helpMessage()
    exit 0
}

// Default parameter values
params.outdir = "."
params.fastq = false
params.references = false
params.adapters = "data/truseq_adapters.fasta"
params.scontaminants = "data/synth_cont.fasta"
params.qscore_cutoff = 20
//params.filter_phred = 5
//params.trim_phred = 2
params.min_read_length = 50
params.contaminants = false
params.map = false
//params.krakendb = false

// INPUT VALIDATION

if ( params.fastq ) {
    fastqPairs_ch = Channel.fromFilePairs(
        params.fastq,
        checkIfExists: true,
        size: 2,
        type: "file"
    )
} else {
    log.info "Hey I need some fastq files to look at please."
    exit 1
}

if ( params.references ) {
    references = Channel.fromPath(
        params.references,
        checkIfExists: true,
        type: "file"
    )
}

// END OF INPUT VALIDATION

//fastqPairs.into {
//    fastqPairs4QC;
//    fastqPairs4AdapterTrimming;
//    fastqPairs4Alignment;
//}

// fastp trimmed files are published, json are only sent in the channel and used only by multiqc
process fastp {

    tag "filter $sample_id"
    //echo true
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy', pattern: 'trimmed_reads/*' // publish only trimmed fastq files

    input:
	tuple val(sample_id), path(reads)
    
    output:
	tuple val(sample_id), path("${sample_id}_filt_R*.fastq.gz"), emit: trimmed_reads
	path("${sample_id}.fastp.json"), emit: json

    script:
    """
    mkdir trimmed_reads

    fastp -i ${reads[0]} -I ${reads[1]} \\
      -o 'trimmed_reads/${sample_id}_filt_R1.fastq.gz' -O 'trimmed_reads/${sample_id}_filt_R2.fastq.gz' \\
      -q $params.qscore_cutoff \\
      -l $params.min_read_length \\
      -w ${task.cpus} \\
      -j 'trimmed_reads/${sample_id}.fastp.json'    
    """
}

workflow {
    fastqPairs_ch.view()
    fastp(fastqPairs_ch)
}
