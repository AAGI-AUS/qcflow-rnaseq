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
module load singularity/4.1.0-slurm

nextflow run -resume -profile singularity,pawsey_setonix ./main.nf \
  --workflow align \
  --aligner star-plants \
  --sjOverhang 149 \
  --genome "/scratch/y95/kgagalova/wheatLeon/genome/Triticum_aestivum.IWGSC.dna.toplevel.fa" \
  --genes "/scratch/y95/kgagalova/wheatLeon/genome/Triticum_aestivum.IWGSC.60.gtf" \
  --output_dir results \
  --index_dir "/scratch/y95/kgagalova/wheatLeon/genome/results/star_index/" \
  --fastq_dir '/scratch/y95/kgagalova/wheatLeon/reads/*R{1,2}.trimmed.fastq.gz' \
  --strandedness RF \
  --genome_size 14000000000
