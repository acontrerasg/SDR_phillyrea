#!/usr/bin/env bash

# Load needed modules
module load bwa/0.7.17
module load samtools/1.16.1
module load parallel

# Specify the threads:
threads_per_job=4
jobs=4

# Variables and path
dir="PATH/TO/FOLDER/"
mkdir -p ${dir}
cd ${dir}
# Force temporary directory to current working directory
export TMPDIR=$(pwd)


# Step 1: Create an array of BAM files in the current directory
bam_files=(${dir}*_maleh1.bam)
# Step 2: Define the processing function
process_bam() {
    bam_file="$1"
    counter="$2"
    threads="$3"
    # Force temporary directory to current working directory
    export TMPDIR=$(pwd)

    SM=$(basename "${bam_file}" | cut -f1 -d '_' )
    ID="ID_${counter}"
    final_bam="${bam_file%.bam}_filt.bam"
    temp_prefix="${bam_file%.bam}_temp_collate"
    stats_file="${bam_file%.bam}_filt.stats"

    echo "Processing ${bam_file} → ${final_bam}"

    samtools collate -@ "${threads}" -O "${bam_file}" "${temp_prefix}" | \
    samtools fixmate -@ "${threads}" - - -m -u | \
    samtools sort -@ "${threads}" - -u | \
    samtools markdup -@ "${threads}" -u -s - - | \
    samtools addreplacerg -@ "${threads}" -r "@RG\tID:${ID}\tSM:${SM}" - -o "${final_bam}"

    samtools index "${final_bam}"
    samtools flagstats -@ "${threads}" "${final_bam}" > "${stats_file}"

    echo "Done: ${bam_file}"
}
export -f process_bam
# Step 3: Run the function in parallel
echo "Running in parallel using ${threads_per_job} threads per job..."
parallel --jobs "${jobs}" process_bam {} {#} "${threads_per_job}" ::: "${bam_files[@]}"
