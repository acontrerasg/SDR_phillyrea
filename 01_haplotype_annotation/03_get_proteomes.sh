#!/usr/bin/env bash


# Get paths to complete haplotype annotations
source ./00_paths.conf

dir=${Main_anno_dir}${homolog_dir}
mkdir -p ${dir}
cd ${dir}

# Extract CDS and Genes from files:
for assembly in Pan_040 Pan_01N25 ; do
  for hap in hap1 hap2 ; do
    genome_var="${assembly}_${hap}"
    genome=${!genome_var}
    anno_var="${assembly}_${hap}_complete"
    annotation=${!anno_var}
    name="${assembly}.${hap}"
    # MAIN
    bash $translate_gff -g $annotation \
                    -f ${genome} \
                    -c ${name}_genes.cds.fasta \
                    -p ${name}_genes.aa.fasta
  done
done
