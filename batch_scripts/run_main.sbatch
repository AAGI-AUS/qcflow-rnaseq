#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --partition=debug
#SBATCH --time=01:00:00
#SBATCH --account=y95
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load nextflow/23.04.3
module load singularity/3.11.4-slurm

nextflow run -resume -profile pawsey_setonix ./main.nf \
  -with-singularity "docker://0000000000777/qcflow-rnaseq:v0.0.1" \
  --outdir results \
  --fastq "reads/*_R{1,2}.fastq.gz"
