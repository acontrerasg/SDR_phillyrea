#!/usr/bin/env bash

# This script maps the  trimmed reads to both the male (040) and hermaphrodite (01N25) P.angustifolia genomes. 

# Module load
module load bwa/0.7.17
module load samtools/
module load bedtools/

# Path and variables
Maindir="PATH/TO/FOLDER/"
GmapM=${MainDir}"gmap_M.csv" # List of names of the samples
read_path=${MainDir}"LinkageMap/" # Path were the trimmed reads are stored
cd ${read_path}

# Function to check if bwa index files exist
index_exists() {
    if [ -e "$1.sa" ] && [ -e "$1.amb" ] && [ -e "$1.ann" ] && [ -e "$1.bwt" ] && [ -e "$1.pac" ]; then
        return 0  # Index files exist
    else
        return 1  # Index files do not exist
    fi
}

# Function to check if BAM files exist
bam_exists() {
    if [ -e "$1" ] && [ -e "$2" ]; then
        return 0  # BAM files exist
    else
        return 1  # BAM files do not exist
    fi
}

# Reference locations
genome_paths=(
    "PATH/TO/FILE"
    "PATH/TO/FILE"
    "PATH/TO/FILE"
    "PATH/TO/FILE"
)

# Create output directory if it doesn't exist
mkdir -p bam_output

## MAIN 

# Loop through the array and process each genome
for path in "${genome_paths[@]}"; do
    # Extract just the filename from the path
    filename=$(basename "$path" .fasta)
    
    echo "Processing $filename..."
    
    # Create the symlink, with error handling
    if [ -e "$filename" ]; then
        echo "Warning: $filename already exists in current directory, skipping symlink creation..."
    else
        if ln -s "$path" "$filename"; then
            echo "Created symlink for $filename"
        else
            echo "Error creating symlink for $filename. Skipping this genome."
            continue
        fi
    fi

    # Check if genome file exists and is readable
    if [ ! -r "$filename" ]; then
        echo "Error: Cannot read $filename. Skipping this genome."
        continue
    fi  # Fixed syntax: changed } to fi

    # Check if indexes for genomes exists
    if index_exists "$filename"; then
        echo "Index for $filename already exists, skipping indexing."
    else
        echo "Index for $filename does not exist, creating index..."
        if ! bwa index "$filename"; then
            echo "Error: BWA indexing failed for $filename. Skipping this genome."
            continue
        fi
    fi

    # Define output files
    temp_bam="bam_output/${filename}_markers.bam"
    output_bam="bam_output/${filename}_markers.sorted.bam"
    
    echo "Running BWA alignment for $filename..."
    
    # Run BWA alignment and process with samtools in steps with error checking
    bwa mem -t 48 "$filename" mapped_markers.fasta |
    samtools view -@ 48 -b -o "$temp_bam" - ;
    samtools sort -@ 48 -o "$output_bam" "$temp_bam";
    samtools index "$output_bam" ; 
    # Clean up temporary BAM file
    rm -f "$temp_bam"
    
    echo "Successfully processed $filename"
    echo "----------------------------------------"
done

# Get bed positions
for bam in *sorted.bam ; do
    name=$(basename ${bam} .bam) ;
    bedtools  bamtobed -i ${bam} > ${name}.contigpos.bed ;
done 

echo "All processing complete!"
