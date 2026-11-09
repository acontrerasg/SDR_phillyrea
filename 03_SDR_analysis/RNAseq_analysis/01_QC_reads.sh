#!/usr/bin/env bash
# Module load 
module load fastqc/0.12.1
module load multiqc/

# Threads to use
threads=8
# Move to path and create output dir
dir="PATH/TO/FOLDER/"
mkdir -p ${dir} 
cd ${dir}  
outdir=${dir}/reads_qc/
mkdir -p ${outdir}

## Run FASTQC one by one
for fq in $(cat read_path_list.txt) ; do
        fastqc -o ${outdir} -t ${threads} $fq
done
# MULTIQC
multiqc ${outdir} -o ${outdir} \
    --filename multiqc_report \
    --title "P.ANG SET 2  -SET 3 QC" \
    -v
