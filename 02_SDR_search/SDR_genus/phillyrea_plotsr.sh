#!/usr/bin/env bash
## Load environment where plotsr is installed
module load plotsr/0.5.2 syri/1.6 minimap2/2.29 seqkit/2.9.0 samtools/
# Set directory
MainDir="PATH/TO/FOLDER/"
mkdir -p ${MainDir} && cd ${MainDir}


# Create aux files for plotsr
# Get the length of the  chrom files
for i in *fasta ; do 
       samtools faidx $i ;
done

for i in *fai ; do echo $i ; 
       name=$(basename $i .fasta.fai) ; 
       cut -f1-2 $i > ${name}.chrlen ;
done 
# Get chromsome order:
cut -f1 pang_mal_hap1.chrlen | sort -V > chrord.txt
# Chromosome file
cat <<EOF > genomes_pan.txt
#file  name   tags
pang_her_hap1.chrlen AngHerm.h1    ft:cl;lw:1.5
pang_her_hap2.chrlen AngHerm.h2    ft:cl;lw:1.5
pang_mal_hap2.chrlen AngMale.h2    ft:cl;lw:1.5
pang_mal_hap1.chrlen AngMale.h1    ft:cl;lw:1.5
Plat.hap1.mod.chrlen LatMale.h1    ft:cl;lw:1.5
Plat.hap2.mod.chrlen LatMale.h2    ft:cl;lw:1.5
EOF

##############################################################################

# MAIN
plotsr --sr A_syri_syri.out \
       --sr B_syri_syri.out \
       --sr C_syri_syri.out \
       --sr D_syri_syri.out \
       --sr E_syri_syri.out \
       --genomes genomes.txt \
       -R -S 0.75  -W 7 -H 20 -f 8 \
       -o output_plot_local.pdf \
       --chrord chrord.txt


# Plot chrom 12 only 
plotsr --sr A_syri_syri.out \
       --sr B_syri_syri.out \
       --sr C_syri_syri.out \
       --sr D_syri_syri.out \
       --sr E_syri_syri.out \
       --genomes genomes.txt \
       -R -S 0.75  -W 6 -H 4 -f 8 \
       -o output_plot_chr12_local.pdf \
       --chr chr12
