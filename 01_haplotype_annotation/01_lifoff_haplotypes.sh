#!/usr/bin/env bash


# Modules
module load liftoff/1.6.3
module load seqkit
# Threads
threads=24
# Liftoff function wrapper
liftoff_wrapper() {
    local target="$1"
    local reference="$2"
    local annotation="$3"
    local threads="${4:-1}"
    local outname=$(basename "${target%.*}")

    liftoff "${target}" "${reference}" \
        -g "${annotation}" \
        -o "${outname}.liftoff.gff" \
        -u "${outname}.unmapped.features.txt" \
        -p "$threads"
}
export -f liftoff_wrapper

# Source the paths used in this script 
source ./00_paths.conf
# Create used dir
dir=${Main_anno_dir}${gene_anno_subdir}
mkdir -p ${dir}
cd ${dir}

# LOOP MALE
grep "^chr_M" ${Pan_040_anno} > Pan_040_anno_chrom.gff

for genome in "$Pan_040_hap1" "$Pan_040_hap2"; do
    liftoff_wrapper "$genome" "$Pan_040_genome" Pan_040_anno_chrom.gff  "$threads"
done
# LOOP HERMA
for genome in "$Pan_01N25_hap1" "$Pan_01N25_hap2"; do
    liftoff_wrapper "$genome" "$Pan_01N25_genome" "$Pan_01N25_anno" "$threads"
done