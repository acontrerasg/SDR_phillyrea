#!/usr/bin/env bash


# Module load
module load parallel
module load meryl

# Count Kmers fun
Kmer_count() {
    local sample="$1"
    local ReadPath="$2"
    local kmer="$3"
    local threads="$4"

    local memory=60
    
    # Skip empty lines
    [ -z "$sample" ] && return 0
    
    echo "Processing sample: $sample"
    
    # Check if input files exist
    local read1="${ReadPath}${sample}_1_val_1.fq.gz"
    local read2="${ReadPath}${sample}_2_val_2.fq.gz"
    
    if [[ ! -f "$read1" ]]; then
        echo "Warning: Read 1 file not found: $read1" >&2
        return 1
    fi
    
    if [[ ! -f "$read2" ]]; then
        echo "Warning: Read 2 file not found: $read2" >&2
        return 1
    fi
    
    # Run meryl
    meryl count k=${kmer} memory=${memory} threads=${threads} \
        "$read1" \
        "$read2" \
        output "${sample}_k${kmer}.meryl"
    
    if [[ $? -eq 0 ]]; then
        echo "Successfully processed: $sample"
        return 0
    else
        echo "Error processing: $sample" >&2
        return 1
    fi
}

# Export the function so GNU parallel can use it
export -f Kmer_count

# Path and variables
SampleList="PATH/TO/FILE"
Kdir="PATH/TO/FOLDER/"
ReadPath="PATH/TO/FOLDER/"
kmer=21

# Total threads=32
threads_per_job=8  # Threads for each meryl job
parallel_jobs=2    # Number of parallel jobs to run. Each parallel job takes 16 gb memory

# Set  dir
mkdir -p ${Kdir} && cd ${Kdir}

# Run with GNU parallel
cat "$SampleList" | parallel -j ${parallel_jobs} \
     Kmer_count {} "$ReadPath" "$kmer" "$threads_per_job"
# Check exit status
if [[ $? -eq 0 ]]; then
    echo "All samples processed successfully!"
else
    echo "Some samples failed. Check the output above." >&2
fi
