#!/usr/bin/env bash

# modules
module load seqkit

# Directory
dir="PATH/TO/FOLDER/"
cd ${dir}
# Input and output files
reference="PATH/TO/FILE"
# output:
output="Plat.hap1.flipped.fasta"
# List of chromosomes to flip (reverse complement)
chromosomes_to_flip=(
    "Pla_h1-1"
    "Pla_h1-2"
    "Pla_h1-3"
    "Pla_h1-4"
    "Pla_h1-5"
    "Pla_h1-6"
    "Pla_h1-7"
    "Pla_h1-8"
    "Pla_h1-9"
    "Pla_h1-10"
    "Pla_h1-11b"
    "Pla_h1-12"
    "Pla_h1-13"
    "Pla_h1-15"
    "Pla_h1-16b"
    "Pla_h1-19"
    "Pla_h1-20"
    "Pla_h1-21a"
    "Pla_h1-22"
)

# Create a pattern file for seqkit
pattern_file=$(mktemp)
printf "%s\n" "${chromosomes_to_flip[@]}" > "$pattern_file"

# Extract and reverse complement the specified chromosomes
seqkit grep -f "$pattern_file" "$reference" | seqkit seq -rp > flipped_temp.fasta

# Extract chromosomes NOT in the list (keep as is)
seqkit grep -v -f "$pattern_file" "$reference" > unflipped_temp.fasta

# Combine both files
cat flipped_temp.fasta unflipped_temp.fasta > "$output"

# Clean up temporary files
rm "$pattern_file" flipped_temp.fasta unflipped_temp.fasta

echo "Flipped chromosomes saved to: $output"
# remove unscaffolded contings
#seqkit grep -vnrp "h1tg"  $output" > temph1
#mv temph1 $output
