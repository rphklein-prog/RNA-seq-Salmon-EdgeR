#!/bin/bash

# takes two arguments: path to the data directory, and desired output directory path and name. 
# searches through given data_directory and identifies paired end reads, runs bbduk
# trimming with given adapt_seqs. Assumes forward and reverse reads end in _1.fq.gz and 
# _2.fq.gz. 

set -euo pipefail

data_directory=$1
trimmed_out=$2
adapt_seqs="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT,GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG"

	
mkdir -p "$trimmed_out"


find "$data_directory" -name '*_1.fq.gz' -print0 |
	while IFS= read -r -d '' line; do
		echo "running bbduk analysis on $line"
		sample_name=$(basename "$line" _1.fq.gz)
		
		in1="$line"
		in2="${line/_1.fq.gz/_2.fq.gz}"
		out1="${trimmed_out}/${sample_name}_clean_1.fq.gz"
		out2="${trimmed_out}/${sample_name}_clean_2.fq.gz"
		
		if [[ ! -f "$in2" ]]; then
        	echo "$in2 not found! Skipping"
        	continue
    	fi
		
		bbduk.sh \
			in1="$in1" \
			in2="$in2" \
			out1="$out1" \
			out2="$out2" \
			literal="$adapt_seqs" \
			ktrim=r \
			k=23 \
			mink=11 \
			hdist=1 \
			tpe \
			tbo
	done
