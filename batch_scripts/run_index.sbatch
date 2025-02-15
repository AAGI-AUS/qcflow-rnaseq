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

# use container structure
nextflow run -resume -profile singularity,pawsey_setonix ./main.nf \
  --workflow genome-index \
  --aligner star-plants \
  --genome "/path/to/genome/gemome.fa" \
  --genes "/path/to/genes/genes.gtf" \
  --genome_size 14000000000 \
  --output_dir results

# genome size must be included here
