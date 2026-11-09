#!/usr/bin/env bash

# modules
module load seqkit

# Directory
dir="PATH/TO/FOLDER/"
cd ${dir}
# Input and output files
reference="PATH/TO/FILE"
# output:
output="Plat.hap2.flipped.fasta"
# List of chromosomes to flip (reverse complement)
chromosomes_to_flip=(
    "Pla_h2-2"
    "Pla_h2-3_1"
    "Pla_h2-4"
    "Pla_h2-9"
    "Pla_h2-10_1"
    "Pla_h2-13"
    "Pla_h2-15"
    "Pla_h2-22"
    "Pla_h2-23_2"
    "Pla_h2-6_2"
    "Pla_h2-8_1"
    "Pla_h2-8_2"
    "Pla_h2-22"

)

# Create a pattern file for seqkit
pattern_file=$(mktemp)
printf "%s\n" "${chromosomes_to_flip[@]}" > "$pattern_file"

# Extract and reverse complement the specified chromosomes
seqkit grep -f "$pattern_file" "$reference" | seqkit seq -rp > flipped_temp_h2.fasta

# Extract chromosomes NOT in the list (keep as is)
seqkit grep -v -f "$pattern_file" "$reference" > unflipped_temp_h2.fasta

# Combine both files
cat flipped_temp_h2.fasta unflipped_temp_h2.fasta > "$output"

# Clean up temporary files
rm "$pattern_file" flipped_temp_h2.fasta unflipped_temp_h2.fasta

echo "Flipped chromosomes saved to: $output"

# Change the name of a chromsome to the correct match:
sed -i 's/>Pla_h2-6_2/>Pla_h2-3_4/' $output

# Remove unscafooled contigs
#seqkit grep -vnrp "h2tg" $output > temph2
#mv temph2 $output
