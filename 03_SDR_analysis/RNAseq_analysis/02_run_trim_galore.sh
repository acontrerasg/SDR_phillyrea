#!/usr/bin/env bash

# Module load 
module load trim_galore/
module load fastqc/0.12.1
module load multiqc/
# Threads to use
threads=8

# FUnction for checking compression of files
process_file () {
  infile="$1"

  if [[ "$infile" == *.bz2 ]]; then
    base=$(basename "$infile" .bz2)
    out="$new_compress_path/${base}.gz"

    if [[ ! -f "$out" ]]; then
      echo "Converting: $infile -> $out"
      bzcat "$infile" | bgzip -c > "$out"
    else
      echo "Using cached: $out"
    fi

    echo "$out"
  else
    echo "$infile"
  fi
}

# Move to path and create output dir
dir="PATH/TO/FOLDER/"
cd ${dir}

# Make output directory for trim_galore
output_dir=${dir}/trimmed_reads/
new_compress_path=${dir}/gz_compressed_reads/
mkdir -p "$new_compress_path"


# MAIN LOOP: Process each sample, find corresponding R1 and R2, convert if needed, and run trim_galore
while read -r sample; do
  echo "Processing sample: $sample"

  read1=$(grep -E "/${sample}.*(R1\.fastq\.gz|D_1\.fastq\.bz2)" read_path_list.txt)
  read2=$(grep -E "/${sample}.*(R2\.fastq\.gz|D_2\.fastq\.bz2)" read_path_list.txt)

  if [[ -z "$read1" || -z "$read2" ]]; then
    echo "WARNING: missing reads for $sample"
    continue
  fi

  read1_new=$(process_file "$read1")
  read2_new=$(process_file "$read2")

  trim_galore --paired \
    --length 50 \
    --output_dir "$output_dir" \
    --basename "$sample" \
    "$read1_new" "$read2_new"

done < sample_names.txt

## Run FASTQC one by one
for fq in ls ${output_dir}/*.gz ; do
        fastqc -o ${output_dir} -t ${threads} $fq
done
# MULTIQC
multiqc ${output_dir} -o ${output_dir} \
    --filename multiqc_report \
    --title "P.ANG SET 2 and SET 3 QC After Trimming" \
    -v
