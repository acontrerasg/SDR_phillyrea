#!/bin/bash
#SBATCH --partition=COMPUTE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --job-name=OsmaMini
#SBATCH -o PATH/TO/FILE
#SBATCH -e PATH/TO/FILE

# Threads
threads=12
module load agat seqkit
# PAth and threads
dir="PATH/TO/FOLDER/"
genome_dir="PATH/TO/FOLDER/"
output_dir="PATH/TO/FOLDER/"
# Script
miniprot_genomes_run="PATH/TO/FILE"

# Main  folder
cd ${dir}
proteins="PATH/TO/FILE"

# Find MYB genes in osmanthus genomes:
bash $miniprot_genomes_run -q ${proteins} -d ${genome_dir} -o ${output_dir}  -t $threads

# Extract  found proteins:
cd miniprot_output
module load agat/

#genomes
OFL_genome="PATH/TO/FILE"
OFR_genome="PATH/TO/FILE"


echo "Extracting protein sequences (merged)..."

agat_sp_extract_sequences.pl \
    -g  OFL_vs_query.paf.gff \
    -f  $OFL_genome \
    -t cds \
    -p \
    --merge \
    -o OFL_proteins.fasta

agat_sp_extract_sequences.pl \
    -g  OFR_vs_query.paf.gff \
    -f  $OFR_genome \
    -t cds \
    -p \
    --merge \
    -o OFR_proteins.fasta


# Combine output wiith query and align
sed -i 's/ .*//;s/>/>OFL_/' OFL_proteins.fasta 
sed -i 's/ .*//;s/>/>OFR_/' OFR_proteins.fasta 

# Merge 
cat $proteins OFL_proteins.fasta  OFR_proteins.fasta > combined_proteins.fa 

align_and_tree() {
    # load modules
    module load iqtree/3.0.1 mafft/7.475 
    # Check if input file is provided
    if [ -z "$1" ]; then
        echo "Error: No input file provided"
        echo "Usage: align_and_tree <input.fa> [threads]"
        return 1
    fi
    # Set input file and threads
    local input_file="$1"
    local threads="${2:-4}"  # Default to 4 threads if not specified
    # Check if input file exists
    if [ ! -f "$input_file" ]; then
        echo "Error: Input file '$input_file' not found"
        return 1
    fi
    # Get basename without extension
    local basename="${input_file%.*}"
    
    # Define output files
    local aln_file="${basename}.aln"
    
    echo "Starting alignment with MAFFT..."
    echo "Input: $input_file"
    echo "Threads: $threads"
    echo "Output alignment: $aln_file"
    
    # Run MAFFT
    mafft --thread "$threads" --auto "$input_file" > "$aln_file"
    
    if [ $? -ne 0 ]; then
        echo "Error: MAFFT alignment failed"
        return 1
    fi
    
    echo "Alignment completed successfully"
    echo ""
    echo "Starting phylogenetic tree construction with IQ-TREE..."
    
    # Run IQ-TREE
    iqtree3 -s "$aln_file" -m MFP -bb 1000 -nt "$threads"
    
    if [ $? -ne 0 ]; then
        echo "Error: IQ-TREE analysis failed"
        return 1
    fi
    
    echo "Analysis completed successfully"
    echo "Output files: ${basename}.*"
}

export -f align_and_tree

# MAIN
align_and_tree combined_proteins.fa  $threads
