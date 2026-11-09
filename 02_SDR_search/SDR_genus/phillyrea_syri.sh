#!/bin/bash

# Load needed modules
module load minimap2
module load syri/1.6 
module load samtools
# Parallelization
threads=23
# Set directory
MainDir="PATH/TO/FOLDER/"
mkdir -p ${MainDir} && cd ${MainDir}

# Genome paths
# Pang haplotypes
pang_her_hap1="PATH/TO/FILE"
pang_her_hap2="PATH/TO/FILE"
pang_mal_hap1="PATH/TO/FILE"
pang_mal_hap2="PATH/TO/FILE"
# Pla haplotypes
plat_mal_hap1="PATH/TO/FILE"
plat_mal_hap2="PATH/TO/FILE"
#########

# Function to run minimap2 and Syri

run_minimap2_syri() {
    local threads="$1"
    local reference="$2"
    local query="$3"
    local out="$4"
    
    # Validate required arguments
    if [[ $# -ne 4 ]]; then
        echo "Error: Missing required arguments" >&2
        echo "Usage: run_minimap2_syri <threads> <hap1_file> <hap2_file> <output_prefix>" >&2
        return 1
    fi
    
    # Check if input files exist
    if [[ ! -f "$reference" ]]; then
        echo "Error: Input file '$reference' not found" >&2
        return 1
    fi
    
    if [[ ! -f "$query" ]]; then
        echo "Error: Input file '$query' not found" >&2
        return 1
    fi
    
    # Check if BAM file already exists
    if [[ -f "${out}.bam" ]]; then
        echo "BAM file '${out}.bam' already exists. Skipping minimap2 alignment step..."
    else
        # Using minimap2 for generating alignment in sequential order
        echo "Running minimap2 alignment..."
        minimap2 -t "${threads}" -ax asm5 --eqx "${reference}" "${query}" | \
            samtools view -b -o "${out}.bam"
        # Index
        samtools index ${out}.bam ;
        # Check if BAM file was created successfully
        if [[ ! -f "${out}.bam" ]]; then
            echo "Error: Failed to create BAM file" >&2
            return 1
        fi
    fi
    
    echo "Running SyRI analysis..."
    syri -c "${out}.bam" -r "${reference}" -q "${query}" -k -F B --nc "${threads}" --prefix "${out}"_syri_
    
    echo "Analysis complete!"
}
# Example usage:
# run_minimap2_syri 8 reference_hap1.fa query_hap2.fa output_prefix

# Normalizing naming schemes for syri
sed 's/>Pan_01N25.h1-/>chr/' $pang_her_hap1 > pang_her_hap1.fasta
sed 's/>Pan_01N25.h2-/>chr/' $pang_her_hap2 > pang_her_hap2.fasta
sed 's/>Pan_040.h1-/>chr/' $pang_mal_hap1 > pang_mal_hap1.fasta
sed 's/>Pan_040.h2-/>chr/' $pang_mal_hap2 > pang_mal_hap2.fasta
pang_her_hap1="./pang_her_hap1.fasta"
pang_her_hap2="./pang_her_hap2.fasta"
pang_mal_hap1="./pang_mal_hap1.fasta"
pang_mal_hap2="./pang_mal_hap2.fasta"
# Pla haplotypes
sed 's/_RagTag//;s/>Pla.h1-/>chr/' Plat.hap1.ragtag.fasta > Plat.hap1.mod.fasta
sed 's/_RagTag//;s/>Pla.h2-/>chr/' Plat.hap2.ragtag.fasta > Plat.hap2.mod.fasta
plat_mal_hap1="./Plat.hap1.mod.fasta"
plat_mal_hap2="./Plat.hap2.mod.fasta"
#########
# perform whole genome aligments.
# Comparitions are goign to be made in the follwoing order:
# her hap1 vs her hap2 ; her hap2 vs mal hap2 ; mal hap2 vs male hap1 ; 
# male hap1 vs pla hap1 ; pla hap1 vs pla hap2 ;

# Using minimap2 for generating alignment in sequential order

run_minimap2_syri ${threads} ${pang_her_hap1}  ${pang_her_hap2}  A
run_minimap2_syri ${threads} ${pang_her_hap2}  ${pang_mal_hap2}  B
run_minimap2_syri ${threads} ${pang_mal_hap2}  ${pang_mal_hap1}  C
run_minimap2_syri ${threads} ${pang_mal_hap1}  ${plat_mal_hap1}  D
run_minimap2_syri ${threads} ${plat_mal_hap1}  ${plat_mal_hap2}  E

# Alternative comparitions
run_minimap2_syri ${threads} ${pang_mal_hap2}  ${plat_mal_hap2}  D_1
run_minimap2_syri ${threads} ${pang_mal_hap1}  ${plat_mal_hap2}  D_2
run_minimap2_syri ${threads} ${plat_mal_hap2}  ${plat_mal_hap1}  E_2
