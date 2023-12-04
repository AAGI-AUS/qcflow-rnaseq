#!/bin/env bash

if [ "$1" == "-h" ]; then
  echo "Usage: bash `basename $0` -- bash script to run locally the nextflow pipeline"
  exit 0
fi

# alignment example
nextflow run -resume -with-singularity /home/294622J/singularity_images/docker/qcflow-rnaseq/qcflow-rnaseq_v1.0.0.sif -profile local ./main.nf \
  --workflow reads-qc-cont \
  --aligner hisat \
  --sjOverhang 149 \
  --genome "../Morex_pseudomolecules_v2.fasta" \
  --genes "../Morex.gtf" \
  --output_dir results \
  --library_name 1,6 \
  --hisat_prefix "hisat_index" \
  --fastq_dir "$PWD/../reads/*_{R1,R2}.test.fastq.gz" \
  --bbt_filters "$PWD/results/bbt-filters/*"
