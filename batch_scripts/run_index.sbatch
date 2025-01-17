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
module load singularity/4.1.0-slurm

nextflow run -resume -profile singularity,pawsey_setonix -r main AAGI-AUS/qcflow-rnaseq \
  --workflow genome-index \
  --aligner star \
  --genome "$PWD/references/GCA_009771035.1_CNRS_Lalb_1.0_genomic.fna" \
  --genes "$PWD/references/GCA_009771035.1_CNRS_Lalb_1.0_genomic.gtf" \
  --genome_size 420000000 \
  --output_dir results
