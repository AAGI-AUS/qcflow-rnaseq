#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=work
#SBATCH --time=01:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load nextflow/23.04.3
module load singularity/3.11.4-slurm

nextflow run -resume -profile singularity,pawsey_setonix ./main.nf \
  --workflow genome-index \
  --aligner hisat-highmem \
  --genome "$PWD/genomes/Akashinriki_pseudomolecule_v1.fasta" \
  --genes "$PWD/genomes/Akashinriki.gtf" \
  --output_dir results 
