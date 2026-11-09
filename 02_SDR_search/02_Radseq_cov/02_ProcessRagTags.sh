#!/usr/bin/env bash

# Load module and variables 
# module load stacks/2.60
module load parallel/v20200322
module load fastqc/0.12.1
module load multiqc/1.17

# Directories
rawdir="PATH/TO/FOLDER/"
outdir="PATH/TO/FOLDER/"

# Move to outdir
cd ${outdir}

# Number of threads
threads=12

# Function for parallelizing
func_process_radtag() {
    module load stacks/2.60 ;

    local input_r1=$1
    local input_r2=$(echo ${input_r1} | sed 's/_1.fastq.gz/_2.fastq.gz/')
    local output_dir=$2
    local log_file=$(basename ${input_r1} | sed 's/_1.fastq.gz/_process_radtags.log/')
    echo $input_r1 ;
    echo $input_r2 ;
    echo $log_file ;
    
    process_radtags -1 $input_r1 \
                    -2 $input_r2 \
                    -o ${output_dir} \
                    --renz-1 'pstI' \
                    --renz-2 'mseI' \
                    -E phred33 \
                    -c \
                    -q 2>&1 | tee -a ${output_dir}"/"${log_file}
}

# Export function
export -f func_process_radtag

# Run process_radseq function using GNU Parallel
ls ${rawdir}/*_1.fastq.gz |
parallel -j ${threads} ' func_process_radtag {} ./ '


# We ignore the rem files as they are:
# The process_radtags program wants to keep the reads in phase, 
# so that the first read in the sample_XXX.1.fq file is the mate of 
# the first read in the sample_XXX.2.fq file. Likewise for the second
# pair of reads being the second record in each of the two files and so on. 
# When one read in a pair is discarded due to low quality or a missing
# restriction enzyme cut site, the remaining read can't simply be output to 
# the sample_XXX.1.fq or sample_XXX.2.fq files as it would cause the remaining
# reads to fall out of phase. Instead, this read is considered a remainder
# read and is output into the sample_XXX.rem.1.fq file if the paired-end was
# discarded, or the sample_XXX.rem.2.fq file if the single-end was discarded.

# Also our rem files had mostly 0 and tope 400 reads only. 

# Merge output files of both pair end reads into the  single, merged file
sample_list="PATH/TO/FILE" 

# Loop through each entry in the entries file
while IFS= read -r entry; do
    zcat "${entry}_1.1.fq.gz" "${entry}_2.2.fq.gz" > ${entry}_merged.fq
    # Compress the merged file
    gzip  ${entry}_merged.fq
    fastqc -o ./ -t ${threads} ${entry}_merged.fq.gz ;

done < ${sample_list}

# Done
echo "Processing radtags complete."
