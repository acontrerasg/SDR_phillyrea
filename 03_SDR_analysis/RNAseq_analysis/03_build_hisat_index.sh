#!/usr/bin/env bash

# Module load
module load hisat2/2.2.1
module load samtools/1.16.1

# Threads
threads=8

# Move to directory and create output folder
dir="PATH/TO/FOLDER/"
mkdir -p ${dir}
cd ${dir}
# Create output directory
mkdir -p "./index_hisat/"

# Path to the reference genome
ref="PATH/TO/FILE"
ref_name=$(basename $ref .fasta)

# Copy reference genome to index folder:
cp ${ref} ./index_hisat/

# Index the reference genome
hisat2-build -p ${threads} -f ./index_hisat/${ref_name}.fasta ./index_hisat/${ref_name}
