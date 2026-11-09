#!/usr/bin/env bash

# Module load 
module load hisat2/2.2.1
module load samtools/1.16.1
threads=24
# Move to directory and create output folder
dir="PATH/TO/FOLDER/"
mkdir -p ${dir} 
cd ${dir} 

# Create output directory
mkdir -p "./RNA_map_Haplo/"

# Reference genome path
ref="PATH/TO/FILE"
sample_names="PATH/TO/FOLDER/"
# Loop through samples
for sample in $(cat sample_names.txt); do
    # Assuming paired-end reads with *R1 and *R2 suffixes
    input1="${dir}trimmed_reads/${sample}_val_1.fq.gz"
    input2="${dir}trimmed_reads/${sample}_val_2.fq.gz"
    
    # Output SAM file
    output="${dir}RNA_map_Haplo/${sample}_aligned.sam"
    
    # HISAT2 mapping command
    hisat2 \
        -x "$ref" \
        -1 "$input1" \
        -2 "$input2" \
        -S "$output" \
        -p ${threads}
    
    # Samtools compress and sort
    samtools view -@ ${threads} -b "$output" | samtools sort -@ ${threads} > "${dir}RNA_map_Haplo/${sample}_sorted.bam"
    
    # Index the BAM file
    samtools index -@ ${threads} "${dir}RNA_map_Haplo/${sample}_sorted.bam" 
    
    # Remove the sam
    rm "${output}"
done
