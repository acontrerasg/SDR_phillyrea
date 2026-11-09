#!/usr/bin/env bash

# Function to run InterProScan analysis
run_interproscan() {
    local protein_fasta=$1
    local cpu_cores=$2
    local iprscan_path=$3
    # Extract basename without extension for naming output directory and log files
    local basename_no_ext=$(basename "$protein_fasta" | sed 's/\.[^.]*$//')
    # Create output directory based on input filename
    local output_dir="${basename_no_ext}_interproscan_results"
    mkdir -p "$output_dir"
    # Run InterProScan
    ${iprscan_path} -i ${protein_fasta} -d ${output_dir} -f TSV,GFF3 -goterms -cpu ${cpu_cores} -appl CDD,PANTHER,Pfam,SUPERFAMILY,TIGRFAM
}

export -f run_interproscan

# Main
threads=8
interproscan_path="PATH/TO/FILE"

for file in *fasta ; do
        name=$(basename $file .fasta) ;
        sed 's/*/X/g' $file > ${name}.X.fasta ;
        run_interproscan ${name}.X.fasta  $threads $interproscan_path ;
done

# Once done  sumarize the reuslts using bedtools groupby handy function:
for result in *_interproscan_results/*_genes.aa.X.fasta.tsv ; do
    name=$(basename "$result" _genes.aa.X.fasta.tsv) ;
    cut -f1,4,5,6,12,13,14 "$result" |
    sort -Vk1 | uniq  |
    sed 's/,//g' > temp_sort ; 
    bedtools  groupby -i temp_sort -g 1 -c 2,3,4,5,6,7 -o distinct  > temp_distc ;
    awk 'BEGIN{FS=OFS="\t"}{print $1"\tIP_DB="$2";IP_SigAcc="$3";IP_SigDesc="$4";IP_Acc="$5";IP_Anno="$6";IP_GOterm="$7}' temp_distc |
    sed 's/;IP_GOterm=-$//;s/-mRNA-1//' > "${name}"_interpro2parse.tsv ;
    rm temp_sort temp_distc ;
done
