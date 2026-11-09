#!/usr/bin/env bash

# Load modules
module load bedtools

# GET SDR files"
Pan_040_h1_SDR="PATH/TO/FILE"
Pan_040_h2_SDR="PATH/TO/FILE"

# Get fai files:
Pan_040_h1_fai="PATH/TO/FILE"
Pan_040_h2_fai="PATH/TO/FILE"

# Get gene annotation files:
Pan_040_h1_gene="PATH/TO/FILE"
Pan_040_h2_gene="PATH/TO/FILE"
# Get TE annotation files:
Pan_040_h1_TE="PATH/TO/FILE"
Pan_040_h2_TE="PATH/TO/FILE"

# Get path:
path="PATH/TO/FOLDER/"
mkdir -p ${path}
cd ${path}

# SET win size
win_size=100000

# HAP1
# Create windows
cat $Pan_040_h1_fai  | awk '{print $1"\t"$2}' > Pan_040_h1_chromsize.txt ;
bedtools makewindows -g Pan_040_h1_chromsize.txt -w ${win_size} -i winnum |
sort -k1,1V -k2,2n -k3,3n  > Pan_040_h1_${win_size}win.bed ;
# Associate with SDR
bedtools intersect -f 0.30 -wao -a Pan_040_h1_${win_size}win.bed  -b ${Pan_040_h1_SDR} |
awk 'BEGIN{OFS=FS="\t"} {if($13 == ".") print $1,$2,$3,"ID="$4"_""background" ;
    else print $1,$2,$3,"ID="$4"_"$13}' >  Pan_040_h1_${win_size}winID.bed
# calculate density per window for gene
cat $Pan_040_h1_gene | awk 'BEGIN{OFS=FS="\t"}{if($3 == "gene") print $0}' |
cut -f1 -d ';' | sort -k1,1V -k4,4n -k5,5n  |
bedtools coverage -a Pan_040_h1_${win_size}winID.bed -b - |
awk 'BEGIN{OFS=FS="\t"} {print $1,$2,$3,$4,$8,"Gene"}' > Pan_040_h1_${win_size}win_gene_density.bed
# For TEs
grep -v "^#"  $Pan_040_h1_TE |  cut -f3 | sort | uniq > TE_list.txt

for TE in $(cat TE_list.txt) ; do 
    awk -v item=$TE '{if($3 == item) print $0 }' $Pan_040_h1_TE |
    cut -f1 -d ';' | sort -k1,1V -k4,4n -k5,5n  |
    bedtools coverage -a Pan_040_h1_${win_size}winID.bed -b - |
    awk -v item=$TE 'BEGIN{OFS=FS="\t"} {print $1,$2,$3,$4,$8,item}'
done > Pan_040_h1_${win_size}win_TE_density.bed

# Merge both:
cat Pan_040_h1_${win_size}win_gene_density.bed Pan_040_h1_${win_size}win_TE_density.bed | 
sort -k1,1V -k2,2n -k3,3n > Pan_040_h1_${win_size}_SDR_density_results.tsv

# HAP2
# Create windows
cat $Pan_040_h2_fai  | awk '{print $1"\t"$2}' > Pan_040_h2_chromsize.txt ;
bedtools makewindows -g Pan_040_h2_chromsize.txt -w ${win_size} -i winnum |
sort -k1,1V -k2,2n -k3,3n  > Pan_040_h2_${win_size}win.bed ;
# Associate with SDR
bedtools intersect -f 0.30 -wao -a Pan_040_h2_${win_size}win.bed  -b ${Pan_040_h2_SDR} |
awk 'BEGIN{OFS=FS="\t"} {if($13 == ".") print $1,$2,$3,"ID="$4"_""background" ;
    else print $1,$2,$3,"ID="$4"_"$13}' >  Pan_040_h2_${win_size}winID.bed


# calculate density per window for gene
cat $Pan_040_h2_gene | awk 'BEGIN{OFS=FS="\t"}{if($3 == "gene") print $0}' |
cut -f1 -d ';' | sort -k1,1V -k4,4n -k5,5n  |
bedtools coverage -a Pan_040_h2_${win_size}winID.bed -b - |
awk 'BEGIN{OFS=FS="\t"} {print $1,$2,$3,$4,$8,"Gene"}' > Pan_040_h2_${win_size}win_gene_density.bed

# For TEs
grep -v "^#"  $Pan_040_h2_TE |  cut -f3 | sort | uniq > TE_list.txt

for TE in $(cat TE_list.txt) ; do 
    awk -v item=$TE '{if($3 == item) print $0 }' $Pan_040_h2_TE |
    cut -f1 -d ';' | sort -k1,1V -k4,4n -k5,5n  |
    bedtools coverage -a Pan_040_h2_${win_size}winID.bed -b - |
    awk -v item=$TE 'BEGIN{OFS=FS="\t"} {print $1,$2,$3,$4,$8,item}'
done > Pan_040_h2_${win_size}win_TE_density.bed

# Merge both:
cat Pan_040_h2_${win_size}win_gene_density.bed Pan_040_h2_${win_size}win_TE_density.bed | 
sort -k1,1V -k2,2n -k3,3n > Pan_040_h2_${win_size}_SDR_density_results.tsv
