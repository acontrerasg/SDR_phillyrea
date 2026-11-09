#!/usr/bin/env bash

# This  script retrieves the CDS of the public anotations of the 
# Load module 
module load seqkit/ 
module load mafft/ 
module load revtrans/
module load trimal/1.5.1/
module load iqtree/3.0.1

# Set path and dir
dir="PATH/TO/FOLDER/"
mkdir -p ${dir} ; 
cd ${dir} ; 
#### 
orthogroups="PATH/TO/FILE"
ortho_proteomes="PATH/TO/FOLDER/"
#######################################
# UIF-like orthogroup:
ortho_name=$(grep "augustus_masked-chr_M_12-processed-gene-95.2" ${orthogroups} | cut -f1)
grep "${ortho_name}" ${orthogroups} | tr -s '\t' | cut -f2- |  sed 's/\t/, /g;s/, /\n/g' >  ${ortho_name}.txt
# Collect with seqkit
seqkit  grep -nrf ${ortho_name}.txt  *fasta  > ${ortho_name}.CDS.fa
#############################################3
# Collect proteomes from Orthofinder folder
seqkit  grep -nrf ${ortho_name}.txt  ${ortho_proteomes}*fasta  > ${ortho_name}.aa.fa
# Merge all proteins and make an aliment with mafft:
mafft --maxiterate 1000 --localpair ${ortho_name}.aa.fa  > ${ortho_name}.aln.fa
# Get the correctly CDS
revtrans.py  ${ortho_name}.CDS.fa ${ortho_name}.aln.fa -match name >  ${ortho_name}.CDS.aln.fa
# Cleanup aligment with trimal
trimal -in ${ortho_name}.CDS.aln.fa -out ${ortho_name}.CDS.aln.trimmed.fa  -automated1 ;

# IQTREE After manual curation 
iqtree3 -s ${ortho_name}.CDS.aln.trimmed.fa -m TEST -B 1000 --mem 30GB -T AUTO -st CODON
