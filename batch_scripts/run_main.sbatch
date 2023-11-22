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

nextflow run -resume -profile singularity,pawsey_setonix ./main.nf \
  --workflow align \
  --aligner star \
  --sjOverhang 149 \
  --genome "/scratch/y95/kgagalova/barley/barley_rnaseq/barley_reference/Morex/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "/scratch/y95/kgagalova/barley/barley_rnaseq/barley_reference/Morex/annotation/Morex.gtf" \
  --output_dir results \
  --library_name 1,6 \
  --index_dir "/scratch/y95/kgagalova/barley/barley_rnaseq/alignments/indices/star/Mor/star_index/Morex/" \
  --hisat_prefix hisat_index \
  --fastq_dir "/scratch/y95/kgagalova/barley/barley_rnaseq/alignments/alignments/qcflow-rnaseq/test/*R{1,2}.fastq.gz" \
  --strandedness RF
