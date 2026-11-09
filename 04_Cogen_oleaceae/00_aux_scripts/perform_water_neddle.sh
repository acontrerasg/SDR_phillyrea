#!/bin/bash
# Sylvain Legrand 01 mars 2019
# Translated from perl to bash by Adrian Contreras 2025.
# It needs the embo package installed to work. 
# Script to work with orthologous pairs (tabulated file)
# 1-search for cds sequences
# 2-perform alignment with needle/faster
# 3-transform results to provide to KAKA_calculator

# Function to display usage information
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -1, --cds1 FILE     First CDS tabulated file (required)"
    echo "  -2, --cds2 FILE     Second CDS tabulated file (required)"
    echo "  -p, --pairs FILE    Gene pairs file (required)"
    echo "  -o, --output FILE   Output file prefix (default: OUT)"
    echo "  -a, --aligner TOOL  Alignment tool to use: water or needle (default: water)"
    echo "  -g, --gap FLOAT     Gap open penalty (default: 10.0)"
    echo "  -e, --extend FLOAT  Gap extension penalty (default: 0.5)"
    echo "  -h, --help          Display this help message and exit"
    exit 1
}

# Default values
CDS_tab_1=""
CDS_tab_2=""
listpairs=""
output_prefix="OUT"
aligner="water"
gap_open="10.0"
gap_extend="0.5"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -1|--cds1)
            CDS_tab_1="$2"
            shift 2
            ;;
        -2|--cds2)
            CDS_tab_2="$2"
            shift 2
            ;;
        -p|--pairs)
            listpairs="$2"
            shift 2
            ;;
        -o|--output)
            output_prefix="$2"
            shift 2
            ;;
        -a|--aligner)
            aligner="$2"
            shift 2
            ;;
        -g|--gap)
            gap_open="$2"
            shift 2
            ;;
        -e|--extend)
            gap_extend="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [ -z "$CDS_tab_1" ] || [ -z "$CDS_tab_2" ] || [ -z "$listpairs" ]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if files exist
if [ ! -f "$CDS_tab_1" ]; then
    echo "Error: Cannot open input file $CDS_tab_1"
    exit 1
fi

if [ ! -f "$CDS_tab_2" ]; then
    echo "Error: Cannot open input file $CDS_tab_2"
    exit 1
fi

if [ ! -f "$listpairs" ]; then
    echo "Error: Cannot open input file $listpairs"
    exit 1
fi

# Check aligner validity
if [ "$aligner" != "water" ] && [ "$aligner" != "needle" ]; then
    echo "Error: Aligner must be either 'water' or 'needle'"
    exit 1
fi

# Part1: Search for sequences
function process_sequences {
    # Load sequence files into associative arrays for both genomes
    declare -A listeGenome1
    declare -A listeGenome2

        echo "Loading sequences from $CDS_tab_2..."
    while IFS=$'\t' read -r id seq; do
        listeGenome2["$id"]="$seq"
    done < "$CDS_tab_2"

    echo "Loading sequences from $CDS_tab_1..."
    while IFS=$'\t' read -r id seq; do
        listeGenome1["$id"]="$seq"
    done < "$CDS_tab_1"

    echo "Processing gene pairs from $listpairs..."
    outfile=0
    # For each pair, create a temporary file for genome 1 and genome 2
    while IFS=$'\t' read -r id1 id2; do
    # Trim any whitespace from the IDs
    id1=$(echo "$id1" | tr -d '[:space:]')
    id2=$(echo "$id2" | tr -d '[:space:]')
    
    # Debug output to help with troubleshooting
    echo "Processing pair: '$id1' and '$id2'"
    
    # Create temporary files for genome 1
    if [[ ${listeGenome1[$id1]+_} ]]; then
        echo ">$id1" > Temp_genome1.fasta
        echo "${listeGenome1[$id1]}" >> Temp_genome1.fasta
        echo "Found sequence for $id1 in $CDS_tab_1"
    else
        echo ">Discard_Genome1" > Temp_genome1.fasta
        echo "TTTTTTT" >> Temp_genome1.fasta
        echo "Warning: Sequence for $id1 not found in $CDS_tab_1"
    fi
    
    # Create temporary files for genome 2
    if [[ ${listeGenome2[$id2]+_} ]]; then
        echo ">$id2" > Temp_genome2.fasta
        echo "${listeGenome2[$id2]}" >> Temp_genome2.fasta
        echo "Found sequence for $id2 in $CDS_tab_2"
    else
        echo ">Discard_Genome2" > Temp_genome2.fasta
        echo "AAAAAAAA" >> Temp_genome2.fasta
        echo "Warning: Sequence for $id2 not found in $CDS_tab_2"
    fi
    
    ((outfile++))
    outfilename="${output_prefix}_${outfile}"
    echo "Aligning pair $outfile: $id1 vs $id2..."
    # Launch aligner (water or needle) and directly output in fasta format
    $aligner -asequence Temp_genome1.fasta -bsequence Temp_genome2.fasta -gapopen $gap_open -gapextend $gap_extend -outfile "${outfilename}.fasta" -aformat3 fasta
    # Clean up temporary files
    rm -f Temp_genome1.fasta Temp_genome2.fasta
done < "$listpairs"
}

# Function to parse alignment output
function parse_alignments {
    output_file="${output_prefix}_general.txt"
    # Make sure this file doesnt exits previously
    if [ -f "${output_file}" ]; then
    rm "${output_file}"
    fi
    echo "Creating summary file $output_file..."
    touch "$output_file"
    # Process all output files
    for file in ${output_prefix}_[0-9]*.fasta; do
        cat ${file} | grep ">" | sed 's/>//' | tr '\n' '_'  | sed 's/_$/\n/' >> "$output_file" ;
        seqkit seq -w 0 ${file} | grep -v ">" >> "$output_file" ;
        echo -e "\n" >> "$output_file" ; 
    done 
    # remove double spaces:
    sed -i '$!N; /^\(.*\)\n\1$/!P; D' "$output_file"
    # Clean up
    echo "Cleaning up temporary files..."
    rm  ${output_prefix}_[0-9]*.txt.fasta
    
    echo "Process completed successfully!"
}

# Main execution
echo "Starting orthologous sequence processing..."
process_sequences
parse_alignments