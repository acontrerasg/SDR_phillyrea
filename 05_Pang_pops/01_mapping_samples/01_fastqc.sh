#!/usr/bin/env bash

# Load module
module load fastqc/0.11.9
module load multiqc/1.21

# Threads and paths
threads=12
folder="PATH/TO/FOLDER/"
outdir="PATH/TO/FOLDER/"
# Make sure outdir is created
mkdir -p ${outdir}
# Read quality control
fastqc -t ${threads} ${folder}*.fastq.gz -o ${outdir}

#$ Multiqc
cd ${outdir}
multiqc ${outdir}
