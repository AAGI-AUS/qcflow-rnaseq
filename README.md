# qcflow-rnaseq

QC and trimming pipeline for RNA-seq reads.      
The pipeline is based on the [qcflow](https://github.com/ccdmb/qcflow/tree/master) but including newer software and also specifically focused to RNAseq.      

The pipeline runs the follwing:    
* FastQC for reads quality control
* MultiQC for a prettier visualizatio of reads and alignements
* Fastp for general reads statistics and trimming
* Star or Hisat for alignment
* RSeQC for alignment quality control
* Biobloomtools for contamination screening (TBD)

The pipeline is meant to be executed on the Pawsey-Setonix supercomputer and Nimbus VM. Please get in touch if you want more configurations.       

## Workflow options

| Option    | Description | Software  |
| ----------- | ----------- |-----------|
| genome-index | Generate genomic index | star,hisat |
| read-trimming | Trim input reads and output QC | fastp,fastqc,multiqc |
| reads-qc | Perform QC on input reads | fastqc, multiqc |
| align | Align RNAseq short reads and output counts | [star|hisat],multiqc,rnaseqc,featureCounts|

| Aligner     | Description | Software  |
| ----------- | ----------- |-----------|
| star | 2 pass star aligner (default) | star |
| star-plants | 2 pass star aligner tuned for plant genomes | star |
| star-snps | [Star consensus](https://github.com/alexdobin/STAR/blob/master/CHANGES.md#star-277a-----20201228) aligner | star |
| hisat | Hisat aligner (default) | fastqc, multiqc |
| hisat-highmem | Hisat high-memory aligner | hisat |

## Multi library alignment

The aligner allows to run multiple library alignements. Is samples are provided as ```NAMELIB1_file.fastq.gz```, ```NAMELIB2_file.fastq.gz```.... it's possible to align the reads grouping them by sample name, provided by the user through the ```---library_name``` option.      
Example: If you want to run the samples as ```LIB1_LIB001.fastq.gz``` and ```LIB1_LIB002.fastq.gz``` together in the alignemnt, you can use ```---library_name 1,4```. The pipeline will use the characters from 1 to 4 (inclusive) to determine the samples that will be grouped together. If the user leaves the default option, each pair of reads will ne processed separately, even if they belong to the same sample.      

A script with an example run can be found in ```runme.sh``` for local machines (i.e. Nimbus); the ```batch_scripts/run_main.sbatch``` contains a batch script for the Pawsey HPC.    
