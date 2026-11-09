#!/usr/bin/env bash

# Module load
module load miniprot
module load bedtools
# Source the paths used in this script
source ./00_paths.conf
# Threads
threads=24
# Create used dir
out_dir=${Main_anno_dir}${miniprot_subdir}
input_dir=${Main_anno_dir}${gene_anno_subdir}
mkdir -p ${out_dir}
cd ${out_dir}
# Get anntations
cp ${input_dir}*.liftoff.gff ${out_dir}

# Mask annotated loci
bedtools maskfasta -fi $Pan_040_hap1   -bed Pan_040.hap1-rt.liftoff.gff   -fo Pan_040.hap1-rt.masked.fasta
bedtools maskfasta -fi $Pan_040_hap2   -bed Pan_040.hap2-rt.liftoff.gff   -fo Pan_040.hap2-rt.masked.fasta
bedtools maskfasta -fi $Pan_01N25_hap1 -bed Pan_01N25.hap1-rt.liftoff.gff -fo Pan_01N25.hap1-rt.masked.fasta
bedtools maskfasta -fi $Pan_01N25_hap2 -bed Pan_01N25.hap2-rt.liftoff.gff -fo Pan_01N25.hap2-rt.masked.fasta

# Run miniprot
for genome in *.masked.fasta ; do
  miniprot_output=$(basename $genome .fasta) ;
  miniprot -t $threads --gff $genome $OE9D > ${miniprot_output}.aln
  # Split results in GFF and PAF
  cat ${miniprot_output}.aln | grep "^##PAF" | cut -f2- > ${miniprot_output}.paf ;
  cat ${miniprot_output}.aln | grep -v "^##PAF" > ${miniprot_output}.gff3 ;
done
