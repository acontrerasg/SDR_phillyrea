#!/usr/bin/env bash

# Load modules
module load bedtools

# This script uses the SDR annotation to retrieve genes and TEs belonging to each SDR block.
# directory
dir="PATH/TO/FOLDER/"
cd ${dir}

# Annotation directory
anno_dir="PATH/TO/FOLDER/"
# SDR directory 
SDR_dir="PATH/TO/FOLDER/"

for haplotype in Pan_01N25.hap1-rt Pan_01N25.hap2-rt Pan_040.hap1-rt Pan_040.hap2-rt ; do
    # Features
    genes=${anno_dir}${haplotype}.complete.FunAnno.gff3 ;
    TE=${anno_dir}${haplotype}_Final_TEanno_Chrom.gff3 ;
    # Extract genes
    bedtools intersect  -f 0.5 ${genes} -b ${SDR_dir}${haplotype}_annotated_SDR_synteny.gff3 -wb |
    cut -f1-9,18 | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9";"$10}'  >  ${SDR_dir}${haplotype}_SDR_genes.gff3 ;

    bedtools intersect  -f 0.5 -a ${TE} -b ${SDR_dir}${haplotype}_annotated_SDR_synteny.gff3 -wb |
    cut -f1-9,18 | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9";"$10}'  >   ${SDR_dir}${haplotype}_SDR_TEs.gff3 ;
done 
