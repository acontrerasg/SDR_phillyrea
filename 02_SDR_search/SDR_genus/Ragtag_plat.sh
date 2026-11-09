#!/usr/bin/env bash
# Load needed modules
module load ragtag/
# Parallelization
threads=23
# Set directory
MainDir="PATH/TO/FOLDER/"
mkdir -p ${MainDir} && cd ${MainDir}

# assemblies
pang_mal_hap1="PATH/TO/FILE"
pang_mal_hap2="PATH/TO/FILE"

# Pla haplotypes
plat_mal_hap1="PATH/TO/FILE"
plat_mal_hap2="PATH/TO/FILE"
#########

# scaffold a query assembly
ragtag.py scaffold -r -t $threads $pang_mal_hap1 $plat_mal_hap2_flipped -o ./plat_hap2_ragtag/
ragtag.py scaffold -r -t $threads $pang_mal_hap2 $plat_mal_hap1_flipped -o ./plat_hap1_ragtag

# -r infer gap sizes 

# Rename ragtag  contigs and remove unplaced bits:
# Hap1
sed 's/Pan_040.h2/Pla.h1/' plat_hap1_ragtag/ragtag.scaffold.fasta > Plat.hap1.ragtag.all.fasta
seqkit grep -vnrp "Pla_h1-16a|h1tg" Plat.hap1.ragtag.all.fasta > Plat.hap1.ragtag.fasta

# Hap2
sed 's/Pan_040.h1/Pla.h2/' plat_hap2_ragtag/ragtag.scaffold.fasta > Plat.hap2.ragtag.all.fasta
seqkit grep -vnrp "Pla_h2-3_2|h2tg" Plat.hap2.ragtag.all.fasta > Plat.hap2.ragtag.fasta


