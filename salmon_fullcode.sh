#!/bin/bash

# Runs salmon on files in a data directory, matching forward and reverse reads within the
# directory. Assumes forward and reverse reads end in _1.fq.gz and _2.fq.gz. Takes three 
# arguments: path to the data directory, and desired output directory path and name 
# (should end in .salmon) and the salmon index location. 

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <data_directory> <salmon_output_dir> <salmon_index>"
    exit 1
fi

data_directory="$1"
salmon_out="$2"
salmon_index="$3"

mkdir -p "$salmon_out"

eval "$(conda shell.bash hook)"
conda activate salmon

find "$data_directory" -name '*_clean_1.fq.gz' -print0 |
	while IFS= read -r -d '' line; do
		echo "running salmon analysis on $line"
		sample_name=$(basename "$line" _clean_1.fq.gz)
		in1="$line"
		in2="${line/_clean_1.fq.gz/_clean_2.fq.gz}"
		sample_out="$salmon_out/$sample_name"
		mkdir -p "$sample_out"
		salmon quant \
			-i "$salmon_index" \
			-l A \
			-1 "$in1" \
			-2 "$in2" \
			-o "$sample_out" \
			--seqBias \
			--useVBOpt \
			--validateMappings
	done


