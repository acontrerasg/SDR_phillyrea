#!/usr/bin/env bash

set -u
set -o pipefail

# Module load
module load parallel
module load meryl

# Path and variables
Sample_metadata="PATH/TO/FILE"
Kmer_dir="PATH/TO/FOLDER/"
ReadPath="PATH/TO/FOLDER/"
kmer=21

# Total threads=32
threads_per_job=8  # Threads for each meryl job
parallel_jobs=2    # Number of parallel jobs to run. Each parallel job takes 16 gb memory

# Set  dir
mkdir -p ${Kmer_dir} && cd ${Kmer_dir}

# Get bathch samples
# Herma
cat ${Sample_metadata} | 
grep Fabregues | 
awk '{if($2 == "H") print $1}' > Hermaphrodite_samples.txt
# Male
cat ${Sample_metadata} |
grep Fabregues |
awk '{if($2 == "M") print $1}' > Male_samples.txt


HermaList="Hermaphrodite_samples.txt"     # File with Hermaphrodites sample names (one per line)
MaleList="Male_samples.txt"     # File with Males sample names (one per line)
kmer=21
threads=16

# Output directories
mkdir -p Herma_kmers Male_kmers

echo "=== Processing Hermaphrodites ==="

# Read Hermaphrodites samples into array
mapfile -t Herma_samples < "$HermaList"

if [ ${#Herma_samples[@]} -eq 0 ]; then
    echo "Error: No samples in Hermaphrodites" >&2
    exit 1
fi

echo "Hermaphrodites contains ${#Herma_samples[@]} samples"

# Build meryl union-sum command for Hermaphrodites
Herma_dbs=()
for sample in "${Herma_samples[@]}"; do
    db="${Kmer_dir}/${sample}_k${kmer}.meryl"
    if [[ -d "$db" ]]; then
        Herma_dbs+=("$db")
        echo "  - Added: $sample"
    else
        echo "  - Warning: Database not found: $db" >&2
    fi
done

if [ ${#Herma_dbs[@]} -eq 0 ]; then
    echo "Error: No valid databases found for Hermaphrodites" >&2
    exit 1
fi

# Union all Hermaphrodites k-mers (combine all k-mers from all samples)
echo "Creating union of Hermaphrodites k-mers..."
meryl union-sum threads=${threads} output Herma_kmers/union.meryl "${Herma_dbs[@]}"

# Intersect all Hermaphrodites databases to find common k-mers
echo "Finding common k-mers in Hermaphrodites..."
meryl intersect-sum threads=${threads} output Herma_kmers/common.meryl "${Herma_dbs[@]}"

# Get statistics
echo "Hermaphrodites statistics:"
echo "  Union k-mers:"
meryl statistics Herma_kmers/union.meryl | grep "Found"
echo "  Common k-mers (present in ALL samples):"
meryl statistics Herma_kmers/common.meryl | grep "Found"

echo ""
echo "=== Processing Males ==="

# Read Males samples into array
mapfile -t Male_samples < "$MaleList"

if [ ${#Male_samples[@]} -eq 0 ]; then
    echo "Error: No samples in Males" >&2
    exit 1
fi

echo "Males contains ${#Male_samples[@]} samples"

# Build meryl command for Males
Male_dbs=()
for sample in "${Male_samples[@]}"; do
    db="${Kmer_dir}/${sample}_k${kmer}.meryl"
    if [[ -d "$db" ]]; then
        Male_dbs+=("$db")
        echo "  - Added: $sample"
    else
        echo "  - Warning: Database not found: $db" >&2
    fi
done

if [ ${#Male_dbs[@]} -eq 0 ]; then
    echo "Error: No valid databases found for Males" >&2
    exit 1
fi

# Union all Males k-mers
echo "Creating union of Males k-mers..."
meryl union-sum threads=${threads} output Male_kmers/union.meryl "${Male_dbs[@]}"

# Intersect all Males databases to find common k-mers
echo "Finding common k-mers in Males..."
meryl intersect-sum threads=${threads} output Male_kmers/common.meryl "${Male_dbs[@]}"

# Get statistics
echo "Males statistics:"
echo "  Union k-mers:"
meryl statistics Male_kmers/union.meryl | grep "distinct"
echo "  Common k-mers (present in ALL samples):"
meryl statistics Male_kmers/common.meryl | grep "distinct"

# Find k-mers unique to Hermaphrodites (present in Hermaphrodites but not in Males)
echo "Finding k-mers unique to Hermaphrodites..."
meryl difference threads=${threads} \
    Herma_kmers/common.meryl \
    Male_kmers/union.meryl \
    output Herma_kmers/unique.meryl

# Find k-mers unique to Males (present in Males but not in Hermaphrodites)
echo "Finding k-mers unique to Males..."
meryl difference threads=${threads} \
    Male_kmers/common.meryl \
    Herma_kmers/union.meryl \
    output Male_kmers/unique.meryl

echo ""
echo "=== Final Statistics ==="
echo "Hermaphrodites unique k-mers:"
meryl statistics Herma_kmers/unique.meryl | grep "distinct"
echo "Males unique k-mers:"
meryl statistics Male_kmers/unique.meryl | grep "distinct"
echo "K-mers common to both batches:"
meryl statistics common_to_both_batches.meryl | grep "distinct"

echo "Exporting Unique k-mers to text..."
meryl print Herma_kmers/unique.meryl > Herma_unique_kmers.txt
meryl print Male_kmers/unique.meryl > Male_unique_kmers.txt

echo "Done! Results saved in:"
echo "  - Herma_kmers/common.meryl (k-mers present in ALL Hermaphrodites samples)"
echo "  - Male_kmers/common.meryl (k-mers present in ALL Males samples)"
echo "  - Herma_kmers/unique.meryl (k-mers unique to Hermaphrodites)"
echo "  - Male_kmers/unique.meryl (k-mers unique to Males)"