#!/usr/bin/env bash

# Module load
module load parallel
module load meryl
module load bwa
module load samtools
module load bedtools

# Create main path and set global varaibles
# Path and variables
Kdir="PATH/TO/FOLDER/"
# Set dir
mkdir -p ${Kdir} && cd ${Kdir}
# Threads
threads=16
# Kmer lists
Herma_unique_kmers="PATH/TO/FILE"
Male_unique_kmers="PATH/TO/FILE"

#  Mapping kmers
run_kmer_mapping() {
    local reference="$1"
    local kmer_txt="$2"
    local threads="${3:-1}"
    local outdir="$4"
    local prefix="$5"
    local window="${6:-100}"

    if [[ -z "$reference" || -z "$kmer_txt" || -z "$outdir" || -z "$prefix" ]]; then
        echo "Usage: run_kmer_mapping <reference.fa> <kmer_file.txt> [threads] <outdir> <output_prefix> [window_size]"
        return 1
    fi

    mkdir -p "$outdir"

    local fasta="${outdir}/${prefix}.fa"
    local output_prefix="${outdir}/${prefix}"

    # Check inputs
    [[ -f "$reference" ]] || { echo "Reference not found: $reference"; return 1; }
    [[ -f "$kmer_txt" ]] || { echo "Kmer file not found: $kmer_txt"; return 1; }

    echo "Step 1: Convert kmers to FASTA"
    awk '{print ">ID_"NR"_kmerfreq_"$2"\n"$1}' "$kmer_txt" > "$fasta"

    echo "Step 2: Mapping with bwa"
    bwa mem -t "$threads" -k 21 -T 21 -a -c 1000 "$reference" "$fasta" | \
        samtools view -@ "$threads" -b -o "${output_prefix}.bam" -

    echo "Step 3: Sorting and indexing"
    samtools sort -@ "$threads" "${output_prefix}.bam" -o "${output_prefix}.sorted.bam"
    samtools index "${output_prefix}.sorted.bam"

    echo "Step 4: Convert BAM to BED"
    bedtools bamtobed -i "${output_prefix}.sorted.bam" > "${output_prefix}.sorted.bed"

    echo "Step 5: Create genome windows (${window} bp)"
    awk '{print $1"\t0\t"$2}' "${reference}.fai" | \
        bedtools makewindows -b - -w "$window" > "${outdir}/windows_${window}bp.bed"

    echo "Step 6: Intersect with windows"
    bedtools intersect \
        -a "${outdir}/windows_${window}bp.bed" \
        -b "${output_prefix}.sorted.bed" \
        -c > "${output_prefix}.sorted_intersection_counts.bed"

    echo "Done! Output prefix: ${output_prefix}"
}

# Mapping unique kmers to hap1
hap1="PATH/TO/FILE"

run_kmer_mapping ${hap1} ${Male_unique_kmers} ${threads} ./kmer_mapping Male_unique_hap1 100
run_kmer_mapping ${hap1} ${Herma_unique_kmers} ${threads} ./kmer_mapping Herma_unique_hap1 100

# Mapping unique kmers to hap2
hap2="PATH/TO/FILE"

run_kmer_mapping ${hap2} ${Male_unique_kmers} ${threads} ./kmer_mapping Male_unique_hap2 100
run_kmer_mapping ${hap2} ${Herma_unique_kmers} ${threads} ./kmer_mapping Herma_unique_hap2 100
