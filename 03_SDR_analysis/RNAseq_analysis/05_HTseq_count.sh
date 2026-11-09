#!/usr/bin/env bash

# Load module
module load htseq/2.0.5
module load samtools/1.16.1
module load parallel/

# Create feautre files:
hap1_anno="PATH/TO/FILE"
hap2_anno="PATH/TO/FILE"

feature_folder="PATH/TO/FOLDER/"
mkdir -p ${feature_folder}
cut -f1 -d ';' ${hap1_anno}  >  ${feature_folder}Pan_040_hap1.noinfo.genes.gff3
cut -f1 -d ';' ${hap2_anno}  >  ${feature_folder}Pan_040_hap2.noinfo.genes.gff3

# Feature files
hap1_file="PATH/TO/FILE"
hap2_file="PATH/TO/FILE"

# Move to directory
dir="PATH/TO/FOLDER/"
cd ${dir}

# Function to run HTSeq-count
run_htseq_count() {
    local bam_file="$1"
    local feature_file="$2"
    local output_dir="$3"
    
    # Extract sample name
    local sample_name=$(basename "$bam_file" "_sorted.bam")
    local output_file="${output_dir}/${sample_name}_counts.txt"
    
    htseq-count \
        -f bam \
        -r pos \
        -s no \
        -t gene \
        -i ID \
        --nonunique none \
        "$bam_file" \
        "$feature_file" \
        > "$output_file"
}
# Export the function so it can be used by parallel
export -f run_htseq_count
# HAP1
# Create output directory
outdir="./Counts_h1/"
mkdir -p "${outdir}"
parallel -j 12 run_htseq_count {} "$hap1_file" "$outdir" ::: *_sorted.bam
# Make the count matrix
count_matrix="PATH/TO/FILE"
python3 "${count_matrix}" -i "${outdir}" -o "${outdir}/Count_matrix_h1.tsv"
# HAP2
# Create output directory
outdir="./Counts_h2/"
mkdir -p "${outdir}"
parallel -j 12 run_htseq_count {} "$hap2_file" "$outdir" ::: *_sorted.bam
# Make the count matrix for HAP2 (this was missing in the original)
python3 "${count_matrix}" -i "${outdir}" -o "${outdir}/Count_matrix_h2.tsv"
