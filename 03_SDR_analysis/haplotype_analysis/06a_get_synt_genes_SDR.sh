#!/usr/bin/env bash
# Env and module
module load diamond
module load seqkit
module load bedops

# Dir
main_dir="PATH/TO/FOLDER/"
mkdir -p ${main_dir}
cd ${main_dir}
threads=4
###
# Files
gene_hap1_SDR_gff="PATH/TO/FILE"
gene_hap2_SDR_gff="PATH/TO/FILE"

protein_hap1_SDR_fa="PATH/TO/FILE"
protein_hap2_SDR_fa="PATH/TO/FILE"

# Get bed file with genes at the SDR
cat $gene_hap1_SDR_gff | awk '{if($3=="gene") print $0}' |
grep -v "ID=CA" | sed 's/;.*ID=/;ID=/' |
gff2bed  | cut -f1-3,10 | sed 's/;/\t/' |
sed 's/ID=//' > gene_h1.bed

cat $gene_hap2_SDR_gff | awk '{if($3=="gene") print $0}' |
grep -v "ID=CA" | sed 's/;.*ID=/;ID=/' |
gff2bed  | cut -f1-3,10 | sed 's/;/\t/' |
sed 's/ID=//' > gene_h2.bed

# Get SDR proteins
cut -f4 gene_h1.bed | seqkit grep -nrf - ${protein_hap1_SDR_fa} > protein_h1.fa
cut -f4 gene_h2.bed | seqkit grep -nrf - ${protein_hap2_SDR_fa} > protein_h2.fa


# Create diamond db diamond all vs all
for fasta in protein_h1.fa protein_h2.fa ; do 
    name=$(basename $fasta .fa) ;
    diamond makedb \
        --in $fasta \
        --db ${name}| \
        --threads "$threads" ;
done 

# Running DIAMOND blastp: 
diamond blastp \
    --query protein_h1.fa \
    --db    protein_h2 \
    --out   protein_h1vsh2.tsv \
    --outfmt  6 \
    --header \
    --ultra-sensitive \
    --max-target-seqs 1 \
    --threads $threads

diamond blastp \
    --query protein_h2.fa \
    --db    protein_h1 \
    --out   protein_h2vsh1.tsv \
    --outfmt  6 \
    --header \
    --ultra-sensitive \
    --max-target-seqs 1 \
    --threads $threads
##########################

# Format outputs

tail -n +3  protein_h1vsh2.tsv  | sed 's/# Fields: //;s/, /\t/g'  > tmp
mv tmp protein_h1vsh2.tsv

tail -n +3  protein_h2vsh1.tsv  | sed 's/# Fields: //;s/, /\t/g'  > tmp
mv tmp protein_h2vsh1.tsv
