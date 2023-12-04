#!/bin/env bash

if [ "$1" == "-h" ]; then
  echo "Usage: bash `basename $0` -- bash script to run locally the nextflow pipeline"
  exit 0
fi

# alignment example
nextflow run -resume -profile singularity,pawsey_nimbus,c16r64 -r main ccdmb/qcflow-rnaseq \
  --workflow reads-qc-cont \
  --aligner hisat \
  --sjOverhang 149 \
  --genome "$PWD/genome/Morex_pseudomolecules_v2.fasta" \
  --genes "$PWD/genome/Morex.gtf" \
  --output_dir results \
  --library_name 1,6 \
  --index_dir "$PWD/hisat_indices/Morex/" \
  --hisat_prefix "hisat_index" \
  --fastq_dir "$PWD/rawreads/*_{R1,R2}.fastq.gz" \
  --cdna "$PWD/Morex.cds.fasta" \
  --strandedness RF
