#!/usr/bin/env bash

# threads used
threads=24
# Get paths to complete haplotype annotations
source ./00_paths.conf
# DIR
dir=${Main_anno_dir}${homolog_dir}
mkdir -p ${dir}
cd ${dir}
# module load 
module load seqkit
module load diamond

# Because this new proteome has proteins already lifted from OLEA genome, we are going to remove:
for proteome in *genes.aa.fasta ; do
  name=$(basename $proteome _genes.aa.fasta) ;
  seqkit grep -vnrp "^MP" ${proteome}  > ${name}_clean_maker.aa.fasta ;
done

# RUN DIAMOND:
for proteome in  *_clean_maker.aa.fasta ; do
    name=$(basename $proteome _clean_maker.aa.fasta) ; 
    bash ${run_diamond_script} -q $proteome -s ${OE9D} -t $threads -m 1 -e 0.001 -o ${name}_OLEA9.dmd.out ;
done

# Format the file to be able to be parsed into a final GFF:
# Column1: Pang gene name, Column2: Dimaondhitatrributes:sseqid=<val>;pident=<val>;matchlen=<val>
for proteome in *genes.aa.fasta ; do
   name=$(basename $proteome _genes.aa.fasta) ; 
   tail -n +2 ${name}_OLEA9.dmd.out |
   awk 'BEGIN{FS=OFS="\t"}{print $1,"OE9Hom="$2";OE9Ident="$3";OE9Mlen="$4}' |
   sed 's/-mRNA-1//' >  ${name}_OLEA9.dmd.out2parse.tsv ;
done

# Add OLEA9 description column on those homologs which have it.
cat $OE9_anno  | grep Description | cut -f9 | sed 's/^.*product=/OE9Hom=/' > OLEA9_descriptions.tsv

for table in  *_OLEA9.dmd.out2parse.tsv ; do 
  name=$(basename $table _OLEA9.dmd.out2parse.tsv)
  python3 $run_add_description $table OLEA9_descriptions.tsv ${name}_OLEA9.dmd.out2parse.descrip.tsv ; 
  sed -i 's/Description=/OE9Desc=/;s/$/;/' ${name}_OLEA9.dmd.out2parse.descrip.tsv ; 
done 

# Final files:
# Pan_01N25.hap1_OLEA9.dmd.out2parse.descrip.tsv
# Pan_01N25.hap2_OLEA9.dmd.out2parse.descrip.tsv
# Pan_040.hap1_OLEA9.dmd.out2parse.descrip.tsv
# Pan_040.hap2_OLEA9.dmd.out2parse.descrip.tsv
