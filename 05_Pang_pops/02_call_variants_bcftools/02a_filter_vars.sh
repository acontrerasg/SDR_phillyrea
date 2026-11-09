#!/usr/bin/env bash

# Step 0: Load needed modules and set working path
module load samtools/1.16.1
module load bcftools/1.12
module load vcftools/0.1.16
module load parallel
module load R/4.3.1
threads=12
# Variables
# Path and variables
vardir="PATH/TO/FOLDER/"
mkdir -p ${vardir}
cd ${vardir}

# Reference
reference="PATH/TO/FILE"
# Visu script
visu_script="PATH/TO/FILE"
# Atomize, normalize and split multiallelic sites

vcf=${vardir}full_results.raw.vcf.gz
# Get base name without extension
name=$(basename ${vcf} .vcf.gz)

echo "===== VCF Filtering Pipeline ====="
echo "Input VCF: ${vcf}"
echo "Reference: ${reference}"
echo "Threads: ${threads}"
echo ""

# Step 1: Stats on raw VCF
echo "[Step 1/5] Calculating raw VCF statistics..."
bcftools stats --threads ${threads} ${vcf} > ${name}.raw.stats
echo "Raw stats saved to: ${name}.raw.stats"
echo ""

# Step 2: Normalize VCF (left-align and split multiallelic sites)
echo "[Step 2/5] Atomize and Normalizing VCF (left-align and split multiallelic)..."
bcftools norm \
    	-m -any \
    	-f ${reference} \
    	-c x \
    	--threads ${threads} \
    	${vcf} \
    	-Oz \
    	-o ${name}.normalized.vcf.gz

bcftools index --threads ${threads} ${name}.normalized.vcf.gz
echo "Normalized VCF saved to: ${name}.normalized.vcf.gz"
echo ""

# Step 3: Stats on normalized VCF
echo "[Step 3/5] Calculating normalized VCF statistics..."
bcftools stats --threads ${threads} ${name}.normalized.vcf.gz > ${name}.normalized.stats
echo "Normalized stats saved to: ${name}.normalized.stats"
echo ""

# Step 4: Filter by MAF and missingness
echo "[Step 4/5] Filtering by MAF (>=0.05) and missingness (<=0.2)..."
vcftools \
    --gzvcf "${name}".normalized.vcf.gz \
    --maf 0.05 \
    --max-missing 0.2 \
    --recode \
    --recode-INFO-all \
    --stdout | bgzip -c > ${name}.temp.variantonly.vcf.gz

# create a filtered VCF containing only invariant sites
vcftools --gzvcf ${name}.normalized.vcf.gz \
        --max-maf 0 \
        --recode \
        --stdout | bgzip -c > ${name}.temp.invariantsites.vcf.gz

# index both vcfs using tabix
tabix ${name}.temp.variantonly.vcf.gz
tabix ${name}.temp.invariantsites.vcf.gz

# combine the two VCFs using bcftools concat
bcftools concat \
    --allow-overlaps \
    ${name}.temp.variantonly.vcf.gz ${name}.temp.invariantsites.vcf.gz \
    -O z -o ${name}.filtered.vcf.gz

tabix "${name}".filtered.vcf.gz
echo "Filtered VCF saved to: ${name}.filtered.vcf.gz"
echo ""

# Step 5: Stats on final filtered VCF
echo "[Step 5/5] Calculating final VCF statistics..."
bcftools stats --threads ${threads} ${name}.filtered.vcf.gz > ${name}.filtered.stats
echo "Final stats saved to: ${name}.filtered.stats"
echo ""

# Summary
echo "===== Pipeline Complete ====="
echo "Output files:"
echo "  - ${name}.normalized.vcf.gz (normalized)"
echo "  - ${name}.filtered.vcf.gz (final filtered)"
echo ""
echo "Statistics files:"
echo "  - ${name}.raw.stats"
echo "  - ${name}.normalized.stats"
echo "  - ${name}.filtered.stats"
echo ""

# Print variant counts summary
echo "===== Variant Count Summary ====="
echo -n "Raw variants: "
bcftools view -H ${vcf} | wc -l
echo -n "Normalized variants: "
bcftools view -H ${name}.normalized.vcf.gz | wc -l
echo -n "Filtered variants: "
bcftools view -H ${name}.filtered.vcf.gz | wc -l
echo ""
