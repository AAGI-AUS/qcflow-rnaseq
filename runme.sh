#!/bin/env bash

if [ "$1" == "-h" ]; then
  echo "Usage: bash `basename $0` -- bash script to run locally the nextflow pipeline"
  exit 0
fi

# alignment example
nextflow run -resume -profile nextflow run -resume -profile singularity,pawsey_nimbus,c16r64 ./main.nf \
  --workflow align \
  --aligner star-plants \
  --sjOverhang 149 \
  --genome "/scratch/y95/kgagalova/barley/barley_rnaseq/barley_reference/Morex/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "/scratch/y95/kgagalova/barley/barley_rnaseq/barley_reference/Morex/annotation/Morex.gtf" \
  --output_dir results \
  --library_name 1,6 \
  --index_dir "/scratch/y95/kgagalova/barley/barley_rnaseq/alignments/indices/star/Mor/star_index/Morex" \
  --fastq_dir "/scratch/y95/kgagalova/barley/barley_rnaseq/trimmed_reads/reads_Mor/*filt_R{1,2}.fastq.gz"
