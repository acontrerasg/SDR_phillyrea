#!/usr/bin/env bash
# Env and module
module load bedops
module load bedtools
module load minimap2
module load seqkit
threads=16
#  Extract fasta
gff_to_fasta() {
    local gff="$1"
    local reference="$2"
    local output="${3:-output.fasta}"
    
    # Check if required arguments are provided
    if [[ -z "$gff" ]] || [[ -z "$reference" ]]; then
        echo "Usage: gff_to_fasta <gff_file> <reference_fasta> [output_fasta]"
        return 1
    fi
    
    # Check if input files exist
    if [[ ! -f "$gff" ]]; then
        echo "Error: GFF file '$gff' not found"
        return 1
    fi
    
    if [[ ! -f "$reference" ]]; then
        echo "Error: Reference file '$reference' not found"
        return 1
    fi
    
    # Create temporary BED file
    local temp_bed=$(mktemp)
    
    # Convert GFF to BED using bedops
    echo "Converting GFF to BED format..."
    gff2bed < "$gff"  | cut -f1-3,10 > "$temp_bed"
    
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to convert GFF to BED"
        rm -f "$temp_bed"
        return 1
    fi
    
    # Extract sequences using bedtools getfasta
    echo "Extracting sequences from reference..."
    bedtools getfasta -name -fi "$reference" -bed "$temp_bed" -fo "$output"
    
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to extract sequences"
        #rm -f "$temp_bed"
        return 1
    fi
    
    # Clean up temporary file
    #rm -f "$temp_bed"
    
    echo "Done! Sequences written to $output"
    return 0
}

# Fasta mapping
align_sequence_by_id() {
    local seq_id="$1"
    local fasta_file1="$2"
    local fasta_file2="$3"
    local out_name="$4"
    local threads="${5:-8}"
    
    # Check if required arguments are provided
    if [[ -z "$seq_id" ]] || [[ -z "$fasta_file1" ]] || [[ -z "$fasta_file2" ]] || [[ -z "$out_name" ]]; then
        echo "Usage: align_sequence_by_id <sequence_id> <fasta_file1> <fasta_file2> <out_name> [threads]"
        return 1
    fi
    
    # Check if input files exist
    if [[ ! -f "$fasta_file1" ]]; then
        echo "Error: FASTA file '$fasta_file1' not found"
        return 1
    fi
    
    if [[ ! -f "$fasta_file2" ]]; then
        echo "Error: FASTA file '$fasta_file2' not found"
        return 1
    fi
    
    echo "Processing sequence: $seq_id (using $threads threads)"
    
    # Extract sequences with matching ID from both files
    seqkit grep -nrp "$seq_id" "$fasta_file1" > temp1.fa
    seqkit grep -nrp "$seq_id" "$fasta_file2" > temp2.fa
    
    # Align sequences with minimap2
    minimap2 -t "$threads" -c --secondary=no -x asm20 temp1.fa temp2.fa > "${seq_id}_${out_name}_align.paf"
    
    # Clean up temporary files
    rm temp1.fa temp2.fa
    
    echo "Done! Alignment written to ${seq_id}_${out_name}_align.paf"
    return 0
}

# Set  vars
dir="PATH/TO/FOLDER/"
coordinates_SDR_dir="PATH/TO/FOLDER/"
reference_dir="PATH/TO/FOLDER/"


# MAIN
mkdir -p ${dir} && cd ${dir}

# Get fasta files
for gff in ${coordinates_SDR_dir}*annotated_SDR_synteny.gff3 ; do
    name=$(basename $gff _annotated_SDR_synteny.gff3) ;
    gff_to_fasta ${gff} ${reference_dir}${name}.fasta ${name}_annotated_SDR_synteny.fasta ;
done 

# Extract syntenic IDs from the 040 hap1"
cat Pan_040.hap1-rt_annotated_SDR_synteny.fasta |
grep ">" | sed 's/>//'  | grep -v "ID=M.*" >  syntenic_hap1.txt

for ID in $(cat syntenic_hap1.txt) ; do
    align_sequence_by_id ${ID} \
                         Pan_040.hap1-rt_annotated_SDR_synteny.fasta \
                         Pan_040.hap2-rt_annotated_SDR_synteny.fasta \
                         "040_hap_comp" \
                         ${threads} ;
done

# Remove the ID= from the files:
rename.ul 'ID=' '' *paf