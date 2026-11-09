#!/usr/bin/env bash

# After running 01 and 01b retreive SDR coordinates. 

# Main dir
hap_minimap_dir="PATH/TO/FOLDER/"
mkdir -p ${hap_minimap_dir}
cd ${hap_minimap_dir}
# Module load
module load minimap2/2.17
module load seqkit
module load bedtools
# Genomes
Pan_01N25_hap1="PATH/TO/FILE"
Pan_01N25_hap2="PATH/TO/FILE"
Pan_040_hap1="PATH/TO/FILE"
Pan_040_hap2="PATH/TO/FILE"
#####

# Function
adjust_bed_boundaries() {
    local bed_file=$1
    local fai_file=$2
    local chrom_name=$3
    
    # Get chromosome size from fai file
    local chrom_size=$(awk -v chrom="$chrom_name" '$1 == chrom {print $2}' "$fai_file")
    
    if [ -z "$chrom_size" ]; then
        echo "Error: Chromosome $chrom_name not found in $fai_file" >&2
        return 1
    fi
    
    # Use sed to modify the bed file
   sed "1s/\t[0-9]\+\t/\t0\t/; \$s/\t[0-9]\+\$/\t$chrom_size/" "$bed_file"
}
# Usage:
# adjust_bed_boundaries input.bed reference.fai chromosome_name > output.bed

# Straigthen syntenic files
cat Pang_01N25_h1_syntenic.coord.bed  |
bedtools  sort |
bedtools  merge  -d 1000000 |
adjust_bed_boundaries - Pan_01N25.h1_chr12.fa.fai Pan_01N25.h1-12 > Pan_01N25.h1-12_syntenic_region.bed

# Straigthen syntenic files
cat Pang_01N25_h2_syntenic.coord.bed  |
bedtools  sort |
bedtools  merge  -d 1000000 |
adjust_bed_boundaries - Pan_01N25.h2_chr12.fa.fai Pan_01N25.h2-12 > Pan_01N25.h2-12_syntenic_region.bed

# Straigthen syntenic files
cat Pang_040_h1_syntenic.coord.bed  |
bedtools  sort |
bedtools  merge  -d 1000000 |
adjust_bed_boundaries - Pan_040.h1_chr12.fa.fai Pan_040.h1-12 > Pan_040.h1-12_syntenic_region.bed

# Straigthen syntenic files
cat Pang_040_h2_syntenic.coord.bed  |
bedtools  sort |
bedtools  merge  -d 1000000 |
adjust_bed_boundaries - Pan_040.h2_chr12.fa.fai Pan_040.h2-12 > Pan_040.h2-12_syntenic_region.bed

# Retrieve the SDR region from the synteic portion:
bedtools complement -i Pan_01N25.h1-12_syntenic_region.bed -g Pan_01N25.h1_chr12.fa.fai > Pan_01N25.h1-12_SDR_m_allele.bed
bedtools complement -i Pan_01N25.h2-12_syntenic_region.bed -g Pan_01N25.h2_chr12.fa.fai > Pan_01N25.h2-12_SDR_m_allele.bed
bedtools complement -i Pan_040.h1-12_syntenic_region.bed -g Pan_040.h1_chr12.fa.fai  > Pan_040.h1-12_SDR_M_allele.bed
bedtools complement -i Pan_040.h2-12_syntenic_region.bed -g Pan_040.h2_chr12.fa.fai > Pan_040.h2-12_SDR_m_allele.bed
