#!/usr/bin/env bash

# In this script we are goign to prepare the protein sets that we are goign to use in our synteny analysis 
# between haplotype 1 and Haplotype 2 of P. angustifolia.

# Needs bedop and seqkit
module load seqkit bedops

# Prepare directory:
MainDir="PATH/TO/FOLDER/"
mkdir -p ${MainDir}"genespace_run/"
cd ${MainDir}"genespace_run/"

# Get anntations
cat PATH/TO/FILE  |
grep "^Pan_040.h1"   > Pan_040_Hap1.gff
cat PATH/TO/FILE  |
grep "^Pan_040.h2"   > Pan_040_Hap2.gff
cat PATH/TO/FILE  |
grep "^Pan_01N25.h1" > Pan_01N25_Hap1.gff
cat PATH/TO/FILE  |
grep "^Pan_01N25.h2" > Pan_01N25_Hap2.gff

# Get assemblies
cat PATH/TO/FILE |
seqkit grep -nrp "Pan_040.h1" > Pan_040_Hap1.fasta
cat PATH/TO/FILE |
seqkit grep -nrp "Pan_040.h2" > Pan_040_Hap2.fasta
cat PATH/TO/FILE |
seqkit grep -nrp "Pan_01N25.h1" > Pan_01N25_Hap1.fasta
cat PATH/TO/FILE |
seqkit grep -nrp "Pan_01N25.h2" > Pan_01N25_Hap2.fasta

# Extract CDS and Genes from files:
# Main script 
translate_gff="PATH/TO/FILE"

# Extract CDS and Genes from files:
for assembly in Pan_040 Pan_01N25 ; do
  for hap in Hap1 Hap2 ; do
    genome="${assembly}_${hap}.fasta"
    annotation="${assembly}_${hap}.gff"
    name="${assembly}.${hap}"
    # MAIN
    bash $translate_gff -g $annotation \
                    -f ${genome} \
                    -c ${name}_genes.cds.fasta \
                    -p ${name}_genes.aa.fasta
  done
done

# Change the GFF with the annotaiton to a bed file with only the ID name on it so it
# matches the gene anotation.
for gff in *gff; do 
    echo $gff ;
    name=$(basename $gff .gff) ;
    gff2bed < $gff |
    awk 'BEGIN{OFS="\t"}{if($8 == "mRNA") print $1,$2,$3,$10}' |
    sort -V -k1 | 
    sed 's/ID=//;s/;.*$//' > ${name}.bed ;
done


# Clean bed file:
for type in 040 01N25 ; do
    echo ${type} ;
    for hap in 1 2 ; do
        echo $hap ;

        grep ">" Pan_${type}.Hap${hap}_genes.aa.fasta | sort | sed 's/>//' > namesfa ;
        grep -wf namesfa Pan_${type}_Hap${hap}.bed > allbed ;
        mv allbed Pan_${type}_Hap${hap}.bed ;
    done ;
done


for type in 040 01N25 ; do
    echo ${type} ;
    for hap in 1 2 ; do
        echo $hap
        grep ">" Pan_${type}.Hap${hap}_genes.aa.fasta | wc -l 
        wc -l Pan_${type}_Hap${hap}.bed ;
    done ;
done

# Compare names in gene file and in bed file so they are the same
for type in 040 01N25 ; do
        echo ${type} ;
        for hap in 1 2 ; do
                grep ">" Pan_${type}.Hap${hap}_genes.aa.fasta  | sort | sed 's/>//' > namesfa
                cut -f4 Pan_${type}_Hap${hap}.bed | sort > namesbed
                echo $hap ;
                diff namesbed namesfa ;
                rm namesbed namesfa ;
        done
done

# Make a new directory so names are clearer:
mkdir -p ./{peptide,bed}

# Bed files
cp Pan_040_Hap1.bed ./bed/Male_hap1.bed
cp Pan_040_Hap2.bed ./bed/Male_hap2.bed
cp Pan_01N25_Hap1.bed ./bed/Herm_hap1.bed
cp Pan_01N25_Hap2.bed ./bed/Herm_hap2.bed

# Modify the chromsome names so they are more readable in the final files:
sed -i 's/Pan.*.h[1,2]-//' ./bed/*.bed 

# Fasta proteins
# Change the fasta end to fa (very important for geenspace, it does not recognize .fasta endings).
cp Pan_040.Hap1_genes.aa.fasta ./peptide/Male_hap1.fa
cp Pan_040.Hap2_genes.aa.fasta ./peptide/Male_hap2.fa
cp Pan_01N25.Hap1_genes.aa.fasta ./peptide/Herm_hap1.fa
cp Pan_01N25.Hap2_genes.aa.fasta ./peptide/Herm_hap2.fa

# That is it
echo "DONE!"

#### RUN GENESPACE
module load orthofinder
module load mcscanx
module load R/4.3.1 
module load diamond

Rscript ./genespace.R