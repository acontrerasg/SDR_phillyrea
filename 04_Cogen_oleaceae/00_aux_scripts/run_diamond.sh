#!/usr/bin/env bash

# Get module
module load diamond/2.1.9

# Function to display usage
usage() {
  echo "Usage: $0 [-q <query_file>] [-s <subject_file>] [-o <output_file>] [-t <threads>] [-m <max_target_seqs>] [-e <evalue_threshold>]"
  exit 1
}

# Parse arguments using getopt
while getopts ":q:s:o:t:m:e:" opt; do
  case $opt in
    q) query_file="$OPTARG"
    ;;
    s) subject_file="$OPTARG"
    ;;
    o) output_file="$OPTARG"
    ;;
    t) threads="$OPTARG"
    ;;
    m) max_target_seqs="$OPTARG"
    ;;
    e) evalue_threshold="$OPTARG"
    ;;
    \?) usage
    ;;
  esac
done

# Check if required arguments are provided
if [ -z "$query_file" ] || [ -z "$subject_file" ]; then
  usage
fi

# Set default output file if not provided
if [ -z "$output_file" ]; then
  output_file="diamond_output.m8"
fi

# Set default number of threads if not provided
if [ -z "$threads" ]; then
  threads=4
fi

# Set default max target sequences if not provided
if [ -z "$max_target_seqs" ]; then
  max_target_seqs=1
fi

# Set default e-value threshold if not provided
if [ -z "$evalue_threshold" ]; then
  evalue_threshold=0.001
fi

# Get the current date
current_date=$(date +"%Y-%m-%d")

# Set the log file name
log_file="${output_file%.*}_$current_date.log"

# Run DIAMOND alignment and capture the log
diamond blastp \
  --query "$query_file" \
  --db "$subject_file" \
  --out "$output_file" \
  --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --max-target-seqs "$max_target_seqs" \
  --evalue "$evalue_threshold" \
  --threads "$threads" \
  2> "$log_file"

# Add header to the output file
sed -i '1i qseqid\tqlen\tsseqid\tslen\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' "$output_file"

echo "DIAMOND alignment complete."
echo "Output saved to: $output_file"
echo "Log saved to: $log_file"
echo "DIAMOND alignment complete. Output saved to: $output_file"