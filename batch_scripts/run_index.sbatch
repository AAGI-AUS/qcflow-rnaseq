#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=work
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load nextflow/22.10.0
module load singularity/3.11.4-slurm

nextflow run -resume -profile singularity,pawsey_setonix -r main ccdmb/qcflow-rnaseq \
  --workflow genome-index \
  --aligner hisat-highmem \
  --genome "$PWD/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "$PWD/genome/Morex.gtf" \
  --output_dir results
