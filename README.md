# qcflow-rnaseq

QC and trimming pipeline for RNA-seq reads.      
The pipeline is based on the [qcflow](https://github.com/ccdmb/qcflow/tree/master) but including newer software and also specifically focused to RNAseq.      

The pipeline runs the follwing:       
* FastQC for quality control
* MultiQC for a prettier visualizatio of all the samples
* Fastp for general reads statistics and trimming
* Biobloomtools for contamination screening

The pipeline is meant to be executed on the Pawsey-Setonix supercomputer but other settings are comings soon. 


