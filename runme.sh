#!/bin/env bash

if [ "$1" == "-h" ]; then
  echo "Usage: bash `basename $0` -- bash script to run locally the nextflow pipeline"
  exit 0
fi

# alignment example
nextflow run -resume -profile singularity,pawsey_nimbus,c16r64 ./main.nf \
  --workflow infer-strandedness \
  --aligner hisat \
  --sjOverhang 149 \
  --genome "/home/ubuntu/Documents/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "/home/ubuntu/Documents/genome/Morex.gtf" \
  --output_dir results \
  --library_name 1,6 \
  --index_dir "/home/ubuntu/Documents/barley/rnaseq/alignments/hisat_indices/Morex/" \
  --hisat_prefix "hisat_index" \
  --fastq_dir "/home/ubuntu/Documents/barley/rnaseq/raw_reads/test_runs/HPR453_227325LT3_GAACTGAGCG-TCGTGGAGCG_L001_{R1,R2}.fastq.gz" \
  --cdna /home/ubuntu/Documents/genome/Morex.cds.fasta
