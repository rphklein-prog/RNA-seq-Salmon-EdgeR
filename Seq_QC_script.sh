#!/bin/bash

#Sequencing QC script

# Takes two arguments: first the data directory where the fastq files are stored, second
# the directory to output fastqc results into. Assumes files are in .fq.gz format, and 
# F and R sequences are located together in a subfolder within the data directory.

set -euo pipefail

data_directory=$1
fastqc_out=$2

mkdir -p "$fastqc_out"

find "$data_directory" -name '*.fq.gz' -print0 |
	while IFS= read -r -d '' line; do
		echo "running fastqc analysis on $line"
		fastqc -o "$fastqc_out" -f fastq "$line"
	done

