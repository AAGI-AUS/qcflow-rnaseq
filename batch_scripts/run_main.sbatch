#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2GB
#SBATCH --partition=work
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load nextflow/22.10.0
module load singularity/3.11.4-slurm

nextflow run -resume -profile singularity,pawsey_setonix -r main ccdmb/qcflow-rnaseq \
  --workflow align \
  --aligner star \
  --sjOverhang 149 \
  --genome "$PWD/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "$PWD/Morex/annotation/Morex.gtf" \
  --output_dir results \
  --library_name 1,6 \
  --index_dir "$PWD/indices/star/Mor/star_index/Morex/" \
  --hisat_prefix hisat_index \
  --fastq_dir "$PWD/rawreads/*R{1,2}.fastq.gz" \
  --strandedness RF
