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
* Biobloomtools for contamination screening (TBD)

The pipeline is meant to be executed on the Pawsey-Setonix supercomputer and Nimbus VM. Please get in touch if you want more configurations.       

## Workflow options

| Option    | Description | Software  |
| ----------- | ----------- |-----------|
| genome-index | Generate genomic index | star,hisat |
| read-trimming | Trim input reads and output QC | fastp,fastqc,multiqc |
| reads-qc | Perform QC on input reads | fastqc,multiqc |
| align | Align RNAseq short reads and output counts | star or hisat,multiqc,rnaseqc,(HTseq)|
| infer-strandedness | Infer RNAseq library preparation strandedness | how_are_we_stranded_here | 

| Aligner     | Description | Software  |
| ----------- | ----------- |-----------|
| star | 2 pass star aligner default | star |
| star-plants | 2 pass star aligner tuned for plant genomes | star |
| star-snps | [Star consensus](https://github.com/alexdobin/STAR/blob/master/CHANGES.md#star-277a-----20201228) aligner | star |
| hisat | Hisat aligner default | fastqc,multiqc |
| hisat-highmem | Hisat high-memory aligner | hisat |


## Multi library alignment

The aligner allows running multiple library alignments. Is samples are provided as ```NAMELIB1_file.fastq.gz```, ```NAMELIB2_file.fastq.gz```.... it's possible to align the reads by grouping them by sample name, provided by the user through the ```---library_name``` option.      
Example: If you want to run the samples as ```LIB1_LIB001.fastq.gz``` and ```LIB1_LIB002.fastq.gz``` together in the alignemnt, you can use ```---library_name 1,4```. The pipeline will use the characters from 1 to 4 (inclusive) to determine the samples that will be grouped together. If the user leaves the default option, each pair of reads will be processed separately, even if they belong to the same sample.      

A script with an example run can be found in ```runme.sh``` for local machines (i.e. Nimbus); the ```batch_scripts/run_main.sbatch``` contains a batch script for the Pawsey HPC.    
