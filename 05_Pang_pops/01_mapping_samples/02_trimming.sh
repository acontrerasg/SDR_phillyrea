#!/usr/bin/env bash

# Load module
module load fastqc/0.11.9
module load trim_galore/0.6.6
module load multiqc/1.21

# Threads and paths
threads=12
inputdir="PATH/TO/FOLDER/"
outdir="PATH/TO/FOLDER/"
sample_list="PATH/TO/FILE"

# Move and create paths
cd ${inputdir}
mkdir -p ${outdir}

# Iterate over each FASTQ file in the input directory
for sample in $(cat ${sample_list} ) ; do
     # Trim adapters and perform quality control using trim_galore
     trim_galore --illumina -j ${threads} --output_dir "$outdir" \
     --paired ${sample}_1.fastq.gz ${sample}_2.fastq.gz \
     -fastqc_args "-t ${threads} --outdir ${outdir}"
done 
# Do a multiqc 
cd ${outdir}-illumina || exit
multiqc .
