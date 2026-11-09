#!/usr/bin/env bash

# Main dir
hap_minimap_dir="PATH/TO/FOLDER/"
mkdir -p ${hap_minimap_dir}
cd ${hap_minimap_dir}

# Module load
module load minimap2/2.17
module load seqkit
module load samtools
# Genomes
Pan_01N25_hap1="PATH/TO/FILE"
Pan_01N25_hap2="PATH/TO/FILE"
Pan_040_hap1="PATH/TO/FILE"
Pan_040_hap2="PATH/TO/FILE"

# subset hermaphrodite haps
seqkit grep -p "Pan_01N25.h1-12" $Pan_01N25_hap1 > Pan_01N25.h1_chr12.fa ;
seqkit grep -p "Pan_01N25.h2-12" $Pan_01N25_hap2 > Pan_01N25.h2_chr12.fa ;
# subset Male haps
seqkit grep -p "Pan_040.h1-12" $Pan_040_hap1 > Pan_040.h1_chr12.fa ; 
seqkit grep -p "Pan_040.h2-12" $Pan_040_hap2 > Pan_040.h2_chr12.fa ; 

# Create faidx index 
for i in *fa ; do
	samtools faidx $i ;
done

# MAIN 
minimap2  -x asm5 Pan_01N25.h1_chr12.fa Pan_01N25.h2_chr12.fa > Pan_01N25_asm5.paf
minimap2  -x asm5 Pan_040.h1_chr12.fa Pan_040.h2_chr12.fa > Pan_040_asm5.paf

minimap2  -x asm5 Pan_040.h1_chr12.fa Pan_01N25.h1_chr12.fa > Pan_040_h1_01N25_h1_asm5.paf
minimap2  -x asm5 Pan_040.h1_chr12.fa Pan_01N25.h2_chr12.fa > Pan_040_h1_01N25_h2_asm5.paf

minimap2  -x asm5 Pan_040.h2_chr12.fa Pan_01N25.h1_chr12.fa > Pan_040_h2_01N25_h1_asm5.paf
minimap2  -x asm5 Pan_040.h2_chr12.fa Pan_01N25.h2_chr12.fa > Pan_040_h2_01N25_h2_asm5.paf

# Merge all
cat *_asm5.paf > Pang_haplotypes_asm5.paf

# Clean other results
rm Pan_0*paf




