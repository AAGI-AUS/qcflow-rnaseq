# qcflow-rnaseq

QC and trimming pipeline for RNA-seq reads.      
The pipeline is based on the [qcflow](https://github.com/ccdmb/qcflow/tree/master) but including newer software and also specifically focused to RNAseq.      

The pipeline runs the follwing:    
* FastQC for quality control
* MultiQC for a prettier visualizatio of all the samples
* Fastp for general reads statistics and trimming
* Biobloomtools for contamination screening

The pipeline is meant to be executed on the Pawsey-Setonix supercomputer but other settings are comings soon. 

How to run the pipeline on VM or locally     
```
# in this case we use the nimbus profile - contact me if you want other running options
nextflow run -resume \
	-with-singularity "docker://0000000000777/qcflow-rnaseq:v0.0.5" \
	-profile pawsey_nimbus \
	 main.nf --outdir test \
	--inputdir_fastq "input/*R{1,2}.fastq.gz"
```
