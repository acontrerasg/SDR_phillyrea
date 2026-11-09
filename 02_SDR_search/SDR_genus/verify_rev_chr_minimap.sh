#!/usr/bin/env bash
# Load module and set threads
module load minimap2

threads=8
dir="PATH/TO/FOLDER/"
cd ${dir}

# Path to haplotype resolved latifolia genomes
Pla_hap1="PATH/TO/FILE"
Pla_hap2="PATH/TO/FILE"

# Path to good old Philly
Pang_hap1="PATH/TO/FILE"
Pang_hap2="PATH/TO/FILE"


for asm in asm5 ; do
    minimap2 -t ${threads} -x ${asm} ${Pang_hap1} ${Pla_hap1} > lati.flip.h1_pangh1_${asm}.paf ;
    minimap2 -t ${threads} -x ${asm} ${Pang_hap1} ${Pla_hap2} > lati.flip.h2_pangh1_${asm}.paf ;
    minimap2 -t ${threads} -x ${asm} ${Pang_hap2} ${Pla_hap1} > lati.flip.h1_pangh2_${asm}.paf ;
    minimap2 -t ${threads} -x ${asm} ${Pang_hap2} ${Pla_hap2} > lati.flip.h2_pangh2_${asm}.paf ;
done
