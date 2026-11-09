#!/usr/bin/env bash

threads=12
module load agat 
module load seqkit
module load minimap2

# PAth and threads
dir="/PATH/TO/FOLDER/"
mkdir -p ${dir}
cd ${dir}


# Get genomes here:
# Pang
ln -s PATH/TO/FILE ./
ln -s PATH/TO/FILE ./
# Plat
ln -s PATH/TO/FILE ./
ln -s PATH/TO/FILE ./
# Ofa
ln -s PATH/TO/FILE ./

# MINIMAP2
for genome in *fasta ; do 
    name=$(basename ${genome} .fasta)
    minimap2 -x asm10 OFR.fa ${genome} -t $threads > ${name}.vsOFR.asm10.paf ;
done


