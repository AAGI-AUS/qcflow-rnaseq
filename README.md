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

In order to see the full input options, please use the following
```
nextflow run -r main ccdmb/qcflow-rnaseq --help
...
Typical local pipeline command:

  nextflow run ./main.nf ...

Input/output options
  --workflow        [string]  Preferred pipeline to run (accepted: genome-index, reads-qc, reads-qc-cont, trim, align, infer-strandedness) [default:align]
  --aligner         [string]  Preferred aligner (accepted: star, star-snps, star-plants, hisat, hisat_highmem) [default: star]
  --output_dir      [string]  Output directory [default: results]
  --genome          [string]  Reference genome
  --genes           [string]  Reference annotation
  --fastq_dir       [string]  Input directory with fastq files

Fastp parameters
  --adapters        [string]  File containing adaters to trim [default: ./data/truseq_adapters.fasta]
  --min_read_length [number]  Minimum read length (-l in fastp) [default: 50.0]
  --qual_phred      [number]  Phred score used for trimming (-q in fastp) [default: 30.0]

STAR index
  --sjOverhang      [integer] sjDbOverhang value to use

Other parameters
  --snps            [string]  Vcf file (uncompressed) with snps to include to indexing. Must be the same regerence as genome
  --library_name    [string]  String intervals to extract library name (ex. LIB1_name - 1,4 extracts LIB1)
  --index_dir       [string]  Path to index (star or hisat)
  --hisat_prefix    [string]  Prefix of hisat index [default: hisat_index]
  --cdna            [string]  Fasta file with cDNA sequences [default: None]
  --strandedness    [string]  RNAseq library strandedness (accepted: RF, FR, unstranded) [default: RF]
  --bbt_filters     [string]  null [default: None]
```

## genome-index

This command can be run with any of the aligner types. ```star``` and ```star-plants``` will generate a standard STAR index, ```star-snps``` generates consensus index. All the hisat aligner options generate the same type of index. Please consider your genome size and chose accordingly between ```hisat``` and ```hisat-highmem```. Large genomes such as human or larger require the adequate memory so please use ```hisat-highmem```.      
Example run command
```
nextflow run -resume -profile singularity,pawsey_setonix -r main ccdmb/qcflow-rnaseq \
  --workflow genome-index \
  --aligner [type aligner] \
  --genome "$PWD/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "$PWD/genome/Morex.gtf" \
  --output_dir results
```
Output: returns the index directory in the same directory of the input genome.       

## reads-qc

Example run command
```
nextflow run -resume -profile singularity,pawsey_setonix -r main ccdmb/qcflow-rnaseq \
  --workflow reads-qc \
  --output_dir results \
  --fastq_dir "$PWD/test/*R{1,2}.fastq.gz"
```

Output
```
results
   |-reads-qc
   |---multiqc_data
```

## reads-qc-cont

Example run command
```
nextflow run -resume -profile local -r main ccdmb/qcflow-rnaseq \
  --workflow reads-qc-cont \
  --output_dir results \
  --fastq_dir "$PWD/test/*_R{1,2}.fastq.gz" \
  --bbt_filters "$PWD/biobloom-filters/*"
```

Output
```
results
   |-bbt-contamination
   |-reads-qc
   |---multiqc_data
```

## trim

Example run command
```
nextflow run -resume -profile singularity,pawsey_setonix -r main ccdmb/qcflow-rnaseq \
  --workflow trim \
  --output_dir results \
  --fastq_dir "$PWD/test/*R{1,2}.fastq.gz"
```

Output
```
results
   |-multi_qc-reads
   |---multiqc_data
   |-trimmed_reads
```

## infer-strandedness

Example run command
```
nextflow run -resume -profile singularity,pawsey_setonix -r main ccdmb/qcflow-rnaseq \
  --workflow infer-strandedness \
  --genome "$PWD/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "$PWD/genome/Morex.gtf" \
  --output_dir results \
  --fastq_dir "$PWD/test/*R{1,2}.fastq.gz" \
  --cdna "$PWD/genome/Morex.cds.fasta"
```

Output
```
 results
   |-check_strandedness
```

## align

Example run command
```
nextflow run -resume -profile singularity,pawsey_setonix -r main ccdmb/qcflow-rnaseq \
  --workflow align \
  --aligner star \
  --genome "$PWD/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "$PWD/genome/Morex.gtf" \
  --output_dir results \
  --library_name 1,6 \
  --index_dir "$PWD/index/Morex/" \
  --fastq_dir "$PWD/test/*R{1,2}.fastq.gz" \
  --strandedness RF
```

In case you are using hisat for alignment, please add the parameter ```--hisat_prefix [prefix_used]```.    

Output
```
results
   |-align_rseqc
   |---bam_stat
   |---junction_annotation
   |-----sample
   |-alignements
   |---[aligner]_aligned
   |-----sample
   |-counts
   |-multi_qc-[aligner]
   |---multiqc_data
``` 

## Multi library alignment

The aligner allows running multiple library alignments. Is samples are provided as ```NAMELIB1_file.fastq.gz```, ```NAMELIB2_file.fastq.gz```.... it's possible to align the reads by grouping them by sample name, provided by the user through the ```---library_name``` option.      
Example: If you want to run the samples as ```LIB1_LIB001.fastq.gz``` and ```LIB1_LIB002.fastq.gz``` together in the alignemnt, you can use ```---library_name 1,4```. The pipeline will use the characters from 1 to 4 (inclusive) to determine the samples that will be grouped together. If the user leaves the default option, each pair of reads will be processed separately, even if they belong to the same sample.      

## Memory usage
Hisat and STAR has been tested with several types of genomes, including large genomes (genome size > human genome size), especially looking at memory usage.      
In the current setup of the pipleine, the requested memory usage is calculated from the ```genome_size``` parameter, provided by the user.          

### Examples
* Lupin (*Lupinus albus sp.*) ~50Mb, works well with both Hisat and STAR aligners.     
* Barley (*Hordeum vulgare*) ~5 GB, the genome size memory usage works well with STAR aligner (mapping and index) and default parameters. It requires a high memory usage for Hisat, at the index generation stage.     
* Wheat (*Triticum aestivum*) 17GB, STAR mapping and index have been optimized with high memory options (400 GB) and has been used on High Performance Computing clusters. Hisat is not recommended because it requires large computing memory, which may exceed 1TB.        

| Genome | Size                | Star index         | Star map           | Hisat index        | Hisat map          |
|--------|---------------------|--------------------|--------------------|--------------------|--------------------|
| Lupin  | _Lupinus albus sp._ | :white_check_mark: | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| Barley | _Hordeum vulgare_   | :white_check_mark: | :white_check_mark: | :warning:          | :white_check_mark: |
| Wheat  | _Triticum aestivum_ | :warning:          | :warning:          | :x:                | :x:                |
