#!/usr/bin/env bash

# threads used
threads=24
# module load
module load seqkit
module load diamond
# Get paths to complete haplotype annotations
source ./00_paths.conf
# DIR
dir=${Main_anno_dir}${homolog_dir}
mkdir -p ${dir}
cd ${dir}


# Format Ara11 genes so  we can parse them better:
cat $Ara11 | sed 's/ | /;/g;s/; Chr.* /;/'  | tr -s ';'  > Araport11.genes.fasta

# RUN DIAMOND:
for proteome in  Pan_01N25.hap1_genes.aa.fasta  Pan_01N25.hap2_genes.aa.fasta  Pan_040.hap1_genes.aa.fasta  Pan_040.hap2_genes.aa.fasta ; do
  name=$(basename $proteome .fasta) ;
  bash ${run_diamond_script} -q $proteome -s Araport11.genes.fasta -t $threads -m 1 -e 0.001 -o ${name}_Araport11.dmd.out ;
done

for proteome in *genes.aa.fasta ; do
   name=$(basename $proteome _genes.aa.fasta) ; 
   tail -n +2 ${name}_genes.aa_Araport11.dmd.out |
   awk 'BEGIN{FS=OFS="\t"}{print $1,"Ara11Hom="$2";Ara11Ident="$3";Ara11Mlen="$4}' |
   sed 's/-mRNA-1//;s/Symbols:/Ara11Desc=/' >  ${name}_Ara11.dmd.out2parse.tsv ;
done