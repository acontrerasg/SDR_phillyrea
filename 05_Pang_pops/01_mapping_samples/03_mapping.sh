#!/usr/bin/env bash

# Module load
module load bwa/0.7.17
module load samtools/1.16.1
threads=36

# Path and variables
SampleList="PATH/TO/FILE"
Maindir="PATH/TO/FOLDER/"
Refdir="PATH/TO/FOLDER/"
ReadPath="PATH/TO/FOLDER/"
MappingPath="PATH/TO/FOLDER"
# Reference
reference="PATH/TO/FILE"
# Got the reference to the Refdir
cp ${reference} ${Refdir}

# Move to directory
cd "${Maindir}"

# Function to check if bwa index files exist
index_exists() {
    if [ -e "$1.sa" ] && [ -e "$1.amb" ] && [ -e "$1.ann" ] && [ -e "$1.bwt" ] && [ -e "$1.pac" ]; then
        return 0  # Index files exist
    else
        return 1  # Index files do not exist
    fi
}

# Function to check if BAM file exists (assuming this function is defined elsewhere)
bam_exists() {
    [ -f "$1" ]
}

# Check if indexes for genomes exist
for Reference in $reference ; do
    if [ -z "$Reference" ]; then
        echo "Warning: Reference variable is empty"
        continue
    fi
    
    if index_exists "$Reference"; then
        echo "Index for $Reference already exists, skipping indexing."
    else
        echo "Index for $Reference does not exist, creating index..."
        if ! bwa index "$Reference"; then
            echo "Error: Failed to create index for $Reference"
            exit 1
        fi
    fi
done

# Ensure MappingPath exists and ends with /
MappingPath="${MappingPath%/}/"
mkdir -p ${MappingPath}

# Step 1: Main loop
echo "Starting alignment process..."
while IFS= read -r sample; do
    [ -z "$sample" ] && continue  # Skip empty lines
    
    echo "Processing sample: $sample"
    
    # Check if input files exist
    read1="${ReadPath}${sample}_1_val_1.fq.gz"
    read2="${ReadPath}${sample}_2_val_2.fq.gz"
    
    if [[ ! -f "$read1" ]]; then
        echo "Warning: Read 1 file not found: $read1"
        continue
    fi
    if [[ ! -f "$read2" ]]; then
        echo "Warning: Read 2 file not found: $read2"
        continue
    fi
    
    # Hap1 alignment
    output_bam="${MappingPath}${sample}_maleh1.bam"
    if bam_exists "$output_bam"; then
        echo "BAM file already exists: $output_bam, skipping hap1..."
    else
        echo "Aligning to hap1: $sample"
        if ! bwa mem -t "$threads" "$reference" "$read1" "$read2" | \
             samtools view -@ "$threads" -b -o "$output_bam" -; then
            echo "Error: Failed to align $sample to $reference"
            [ -f "$output_bam" ] && rm "$output_bam"  # Remove incomplete file
            continue
        fi
    fi
done < "$SampleList"

echo "Alignment completed."
