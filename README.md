# qcflow-rnaseq

QC and trimming pipeline for RNA-seq reads.      
The pipeline is based on the [qcflow](https://github.com/ccdmb/qcflow/tree/master) but includes newer software and is also specifically focused on RNAseq.      

The pipeline runs the following:    
* FastQC for reads quality control
* MultiQC for a prettier visualization of reads and alignments
* Fastp for general reads statistics and trimming
* Star or Hisat for alignment
* RSeQC for alignment quality control
* how_are_we_stranded_here for library strandedness 
* Biobloomtools for contamination screening

The pipeline is meant to be executed on the Pawsey-Setonix supercomputer and Nimbus VM. Please get in touch if you want more configurations.       

## Prepare bloom filters

Contamination screening relies on separate fasta files as input. The screening is executed using cds/cdna sequences from protein coding genes. You need to create a separte cds/fasta file for the sequences you want to check. Please check ```docs/cds_fast_download.md``` to find out how efficiently ypu can download fasta files by taxon. See ```runme_create-bloom.sh``` for example run.    

## Workflow options

A script with an example run can be found in ```runme.sh``` for local machines (i.e. Nimbus); the ```batch_scripts/run_main.sbatch``` contains a batch script for the Pawsey HPC.      
The main script involves the following workflow options:    

| Option    | Description | Software  |
| ----------- | ----------- |-----------|
| genome-index | Generate genomic index | star,hisat |
| read-trimming | Trim input reads and output QC | fastp,fastqc,multiqc |
| reads-qc | Perform QC on input reads | fastqc,multiqc |
| reads-qc-cont | Perform QC and contamination screening on reads | fastqc,multiqc,biobloomtools |
| align | Align RNAseq short reads and output counts | star or hisat,multiqc,rnaseqc,featureCount|
| infer-strandedness | Infer RNAseq library preparation strandedness | how_are_we_stranded_here | 

| Aligner     | Description | Software  |
| ----------- | ----------- |-----------|
| star | 2 pass star aligner default | star |
| star-plants | 2 pass star aligner tuned for plant genomes | star |
| star-snps | [Star consensus](https://github.com/alexdobin/STAR/blob/master/CHANGES.md#star-277a-----20201228) aligner | star |
| hisat | Hisat aligner default | hisat |
| hisat-highmem | Hisat high-memory aligner | hisat |


## reads-qc-cont

Example command run
```
nextflow run -resume -profile local -r main ccdmb/qcflow-rnaseq \
  --workflow reads-qc-cont \
  --genome "$PWD/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "$PWD/genes/Morex.gtf" \
  --output_dir results \
  --fastq_dir "$PWD/reads/*_{R1,R2}.test.fastq.gz" \
  --bbt_filters "$PWD/biobloom-filters/*"
```


## Multi library alignment

The aligner allows running multiple library alignments. Is samples are provided as ```NAMELIB1_file.fastq.gz```, ```NAMELIB2_file.fastq.gz```.... it's possible to align the reads by grouping them by sample name, provided by the user through the ```---library_name``` option.      
Example: If you want to run the samples as ```LIB1_LIB001.fastq.gz``` and ```LIB1_LIB002.fastq.gz``` together in the alignemnt, you can use ```---library_name 1,4```. The pipeline will use the characters from 1 to 4 (inclusive) to determine the samples that will be grouped together. If the user leaves the default option, each pair of reads will be processed separately, even if they belong to the same sample.      
