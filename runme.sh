#!/bin/env bash

if [ "$1" == "-h" ]; then
  echo "Usage: bash `basename $0` -- bash script to run locally the nextflow pipeline"
  exit 0
fi

# input files do not work with symlinks
nextflow run -resume -profile singularity,pawsey_nimbus,c16r64 ./main.nf \
  --workflow genome-index \
  --genome "$PWD/genome/Akashinriki_pseudomolecule_v1.fasta" \
  --genes "$PWD/genome/Akashinriki.gtf" \
  --output_dir results
