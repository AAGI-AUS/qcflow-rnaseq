#!/bin/env bash

if [ "$1" == "-h" ]; then
  echo "Usage: bash `basename $0` -- bash script to create biobloom tools used for contamination screening"
  exit 0
fi

nextflow run -resume -profile singularity ./create_bloom.nf \
	--input_fasta "$PWD/bloom-filters/*fasta" \
	--output_dir results
