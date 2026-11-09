#!/usr/bin/env bash

# Step 0: Load needed modules and set working path
module load samtools/1.16.1
module load bcftools/1.12
module load vcftools/0.1.16
module load parallel

threads=46
threads_chrom=2
# Variables
# Path and variables
MappingPath="PATH/TO/FOLDER/"
outdir="PATH/TO/FOLDER/"
mkdir -p ${outdir}
cd ${outdir} || exit
# Reference
reference="PATH/TO/FILE"


call_variants_by_chromosome() {
    local chrom=$1
    local path=$2
    local outdir=$3
    local threads_per_chrom=${4:-1}
    
    echo "Processing chromosome ${chrom}..."
    
    # Create output directory if it doesn't exist
    mkdir -p "${outdir}"
    
    # Variant calling for specific chromosome using chromosome-specific reference
    bcftools mpileup \
        --threads ${threads_per_chrom} \
        --redo-BAQ \
        --min-BQ 20 \
        -d 200 \
        -r ${chrom} \
        -f chromosomes/${chrom}.fasta \
        ${path}/*_maleh1_filt.bam | \
        bcftools call \
        --threads ${threads_per_chrom} \
        -Oz \
        -c \
        -o ${outdir}${chrom}_results.raw.vcf.gz
    
    # Index the output VCF
    bcftools index ${outdir}${chrom}_results.raw.vcf.gz
    
    echo "Done with chromosome ${chrom}!"
}
# Export to be used by parallel
export -f call_variants_by_chromosome



# MAIN 

# Pre-extract all chromosomes from reference
echo "Extracting chromosomes from reference..."
mkdir -p ${outdir}chromosomes/
cut -f1 ${reference}.fai > chromosomes.txt
cut -f1 ${reference}.fai | while read chr; do
    samtools faidx ${reference} ${chr} > ${outdir}chromosomes/${chr}.fasta
    samtools faidx chromosomes/${chr}.fasta
done

echo "Starting parallel variant calling..."

# Run parallel variant calling - one job per chromosome
parallel -j ${threads} call_variants_by_chromosome {} ${MappingPath} ${outdir} ${threads_chrom} :::: chromosomes.txt

echo "All chromosomes processed!"

# Merge all chromosome VCFs into one file
echo "Merging chromosome VCFs..."
bcftools concat \
    -Oz \
    -o ${outdir}full_results.raw.vcf.gz \
    ${outdir}*_results.raw.vcf.gz

bcftools index ${outdir}full_results.raw.vcf.gz

echo "Merging complete! Final file: ${outdir}full_results.raw.vcf.gz"
