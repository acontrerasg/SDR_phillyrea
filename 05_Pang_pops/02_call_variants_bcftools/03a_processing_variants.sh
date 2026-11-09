#!/usr/bin/env bash

# module load
module load bcftools/
module load plink2/2.0a1
module load samtools
# Set env
threads=4
memory=120000

# Path
dir="PATH/TO/FOLDER/"
cd ${dir}
#####

# Set variables
vcf="PATH/TO/FILE"
name=$(basename $vcf .vcf.gz)

# Reheader to remove suffix
bcftools reheader \
    --samples <(bcftools query -l ${vcf} | sed 's/_maleh1.*$//') \
    ${vcf} \
    -o ./${name}.reheader.vcf.gz

# Index with tabix
tabix -p vcf ./${name}.reheader.vcf.gz
# Add header line for DP format field
bcftools annotate \
    --header-lines <(echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">') \
    --threads ${threads} \
    -Oz -o ${name}.reheader.h.vcf.gz \
    ${name}.reheader.vcf.gz

# Index with tabix
tabix -p vcf ${name}.reheader.h.vcf.gz

# Add variant IDs using bcftools (format: CHROM:POS)
bcftools annotate \
  --set-id '%CHROM:%POS' \
  --threads ${threads} \
  -Oz \
  -o ./${name}.annotated.vcf.gz \
   ${name}.reheader.h.vcf.gz

# Index the output file
tabix -p vcf ${name}.annotated.vcf.gz

# Apply a  quality filter and set final name:
bcftools filter -i 'QUAL>=20' ${name}.annotated.vcf.gz |
bcftools +setGT -O z -o Final.hap1.popvars.vcf.gz -- -t q -n . -i 'FMT/DP<5'


# Extract chromosome 12
bcftools view \
    -r Pan_040.h1-12 \
    --threads ${threads} \
    -Oz -o Final.hap1.popvars.chr12.vcf.gz \
   Final.hap1.popvars.vcf.gz

# Index with tabix
tabix -p vcf Final.hap1.popvars.chr12.vcf.gz


cp Final.hap1.popvars.* ../