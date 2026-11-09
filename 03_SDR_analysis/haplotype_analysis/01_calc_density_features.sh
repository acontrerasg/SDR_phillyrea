#!/usr/bin/env bash

# Modules
module load bedtools
# Functions

# Function to calculate genomic feature coverage in sliding windows
# Usage: calculate_feature_coverage <fai> <genes_gff> <satellites_gff> <TEs_gff> <outdir> <name> [winsize]
calculate_feature_coverage() {
    local fai="$1"
    local genes="$2"
    local sat="$3"
    local TEs="$4"
    local outdir="$5"
    local name="$6"
    local winsize="${7:-500000}"  # Default window size: 500kb
    
    # Validate inputs
    if [[ -z "$fai" || -z "$genes" || -z "$sat" || -z "$TEs" || -z "$outdir" || -z "$name" ]]; then
        echo "Error: Missing required arguments"
        echo "Usage: calculate_feature_coverage <fai> <genes_gff> <satellites_gff> <TEs_gff> <outdir> <name> [winsize]"
        return 1
    fi
    
    # Ensure output directory exists and ends with /
    mkdir -p "${outdir}"
    [[ "${outdir}" != */ ]] && outdir="${outdir}/"
    
    echo "===== Genomic Feature Coverage Analysis ====="
    echo "Window size: ${winsize} bp"
    echo "Output directory: ${outdir}"
    echo "Output prefix: ${name}"
    echo ""
    
    # Step 1: Create genomic windows
    echo "[Step 1/3] Creating ${winsize}bp genomic windows..."
    bedtools makewindows -g "${fai}" -w "${winsize}" -s 10000 | \
    sort -k1,1V -k2,2n -k3,3n > "${outdir}${name}.${winsize}.bed"
    
    local window_count=$(wc -l < "${outdir}${name}.${winsize}.bed")
    echo "Created ${window_count} windows"
    echo ""
    
    # Step 2: Concatenate and sort all GFF files
    echo "[Step 2/3] Concatenating and sorting feature annotations..."
    cat "${genes}" "${sat}" "${TEs}" | \
    grep -v "^#" | \
    sort -k1,1V -k4,4n -k5,5n > "${outdir}${name}.features.gff3"
    
    local feature_count=$(wc -l < "${outdir}${name}.features.gff3")
    echo "Total features: ${feature_count}"
    echo ""
    
    # Step 3: Calculate coverage for each feature type
    echo "[Step 3/3] Calculating coverage for each feature type..."
    
    # Get unique feature types
    local feature_types=$(cut -f3 "${outdir}${name}.features.gff3" | sort | uniq)
    local type_count=$(echo "${feature_types}" | wc -l)
    echo "Found ${type_count} feature types"
    echo ""
    
    # Process each feature type
    local counter=0
    for item in ${feature_types}; do
        counter=$((counter + 1))
        local item_name=$(echo "${item}" | sed 's|/|_|g')  
        # Extract features of this type and calculate coverage
        awk -v item="${item}" '$3 == item {print $0}' "${outdir}${name}.features.gff3" | \
        bedtools coverage -a "${outdir}${name}.${winsize}.bed" -b - | \
        awk -v item="${item}" '{print $0"\t"item}'
    done > "${outdir}${name}.feat_cov.${winsize}.bed"
    
    echo ""
    echo "===== Analysis Complete ====="
    echo "Output files:"
    echo "  - ${outdir}${name}.${winsize}.bed (genomic windows)"
    echo "  - ${outdir}${name}.features.gff3 (combined features)"
    echo "  - ${outdir}${name}.feat_cov.${winsize}.bed (feature coverage)"
    echo ""
    
    
    return 0
}

# Example usage (commented out):
# calculate_feature_coverage \
#     genome.fa.fai \
#     genes.gff3 \
#     satellites.gff3 \
#     TEs.gff3 \
#     ./output/ \
#     my_genome \
#     5000

# directory
dir="PATH/TO/FOLDER/"
cd ${dir}
# Create output directory
outdir="./calc_density/"
mkdir -p ${outdir}

# Calculate  density across the 4 haplotypes:
anno_dir="PATH/TO/FOLDER/"
win_size=50000
# 500 Mbp
for haplotype in Pan_01N25.hap1-rt Pan_01N25.hap2-rt Pan_040.hap1-rt Pan_040.hap2-rt ; do
  genome=${anno_dir}${haplotype}.fasta
  fai=${genome}.fai
  genes=${anno_dir}${haplotype}.complete.FunAnno.gff3
  TE=${anno_dir}${haplotype}_Final_TEanno_Chrom.gff3
  Sat=${anno_dir}${haplotype}_Satellite_Chrom.gff3
  
  echo "Processing haplotype: ${haplotype}"
  calculate_feature_coverage \
      ${fai} \
      ${genes} \
      ${Sat} \
      ${TE} \
      ${outdir} \
      ${haplotype} \
      ${win_size}
done

# Now annotate the  intervals as part of the genomic backhround or belonging to the SDR. 
SDR_dir="PATH/TO/FOLDER/"

for haplotype in Pan_01N25.hap1-rt Pan_01N25.hap2-rt Pan_040.hap1-rt Pan_040.hap2-rt ; do
  SDR_file=${SDR_dir}${haplotype}_annotated_SDR_synteny.gff3
  
  echo "Adding SDR information for: ${haplotype}"
  
  # Add SDR information to the WG density calculation 
  bedtools intersect -wao  -f 1.0 \
    -a ${outdir}${haplotype}.feat_cov.${win_size}.bed \
    -b ${SDR_file} | \
  cut -f1-8,17 > ${outdir}${haplotype}.feature_coverage.${win_size}.SDR.bed
done






