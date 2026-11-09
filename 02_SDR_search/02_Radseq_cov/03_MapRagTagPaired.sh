#!/usr/bin/env bash


# Module load
module load bwa/0.7.17
module load samtools/1.16.1
module load bedtools
module load parallel

# Threads used
threads=48

# Path and variables
sample_list="PATH/TO/FILE"
ragtagdir="PATH/TO/FOLDER"
MainDir="PATH/TO/FOLDER/"
outdir="${MainDir}/Male_haplo_mapping/"
mkdir -p ${outdir}
cd ${outdir}

# Reference genome paths
Male_hap1="PATH/TO/FILE"
Male_hap2="PATH/TO/FILE"

# Make symlinks to genomes
ln -s ${Male_hap1} .
ln -s ${Male_hap2} .

# Get base names of reference files
hap1_base=$(basename ${Male_hap1})
hap2_base=$(basename ${Male_hap2})

# Function to check if bwa index files exist
index_exists() {
    if [ -e "$1.sa" ] && [ -e "$1.amb" ] && [ -e "$1.ann" ] && [ -e "$1.bwt" ] && [ -e "$1.pac" ]; then
        return 0
    else
        return 1
    fi
}

# Index both reference genomes if needed
for genome in "$hap1_base" "$hap2_base"; do
    if index_exists "$genome"; then
        echo "Index for $genome already exists, skipping indexing."
    else
        echo "Index for $genome does not exist, creating index..."
        bwa index "$genome"
    fi
done

# Counter for read group IDs
counter=0

# Process each sample against both references
for sample in $(cat ${sample_list} | tail -n +2); do
    for genome in "$hap1_base" "$hap2_base"; do
        hap_suffix=$(basename "$genome" -renamed.fasta)
        
        # Step 1: Map with bwa mem
        output_bam1="${sample}_${hap_suffix}_radtag_paired.bam"
        bwa mem -t ${threads} "$genome" ${ragtagdir}/${sample}_1.1.fq.gz ${ragtagdir}/${sample}_2.2.fq.gz |
        samtools view -@ ${threads} -b -o "$output_bam1" -

        # Step 2: Process the generated BAM file
        ((counter++))
        sorted_bam="${output_bam1%.bam}_sort.bam"
        name=$(basename "${output_bam1}" .bam)

        # Samtools pipeline
        samtools sort -@ "${threads}" "${output_bam1}" > "${sorted_bam}"
        samtools addreplacerg -r "ID:${counter}" -r "SM:${name}" -o "${name}_tagged.bam" "${sorted_bam}"
        samtools index "${name}_tagged.bam"
        samtools stats -@ "${threads}" "${name}_tagged.bam" > "${name}_tagged.bamstats"

        # Clean up temporary files
        rm "${output_bam1}" "${sorted_bam}"
    done
done

# Create indexes and windows for both references
for genome in "$hap1_base" "$hap2_base"; do
    # Create fasta index if needed
    if [ ! -f "${genome}.fai" ]; then
        echo "Creating index for ${genome}..."
        samtools faidx "${genome}"
    fi

    # Create 10kb windows
    win_file="${genome%.fasta}.win10kb.bed"
    bedtools makewindows -g "${genome}.fai" -w 10000 > "$win_file"
done

# Function to calculate coverage
calc_cov() {
    local bam="$1"
    local winfile="$2"
    local name=$(basename "${bam}" | sed 's/_Pan_040.hap[1,2]-rt.fasta_radtag_paired_tagged.bam//')
    bedtools coverage -nonamecheck -a "${winfile}" -b "${bam}" |
    awk -v name="$name" 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,$5,$7,name}' > "${name}.bamcov.bed"
    gzip "${name}.bamcov.bed"
}
export -f calc_cov

# Process coverage for each haplotype
for genome in "$hap1_base" "$hap2_base"; do
    hap_suffix=$(basename "$genome" .fasta)
    win_file="${genome%.fasta}.win10kb.bed"
    
    if ls *${hap_suffix}*_tagged.bam 1> /dev/null 2>&1; then
        # Calculate coverage in parallel
        parallel -j ${threads} calc_cov {} "$win_file" ::: *${hap_suffix}*_tagged.bam

        # Merge coverage files
        zcat *${hap_suffix}*.bamcov.bed.gz > ALL.${hap_suffix}.bamcov.bed ;
        rm *${hap_suffix}*.bamcov.bed.gz ;       
        # Add header
        sed -i '1 i\Chrom\tStart\tEnd\tReadNum\tBasesNum\tFracCov\tSample' ALL.${hap_suffix}.bamcov.bed      
        # Compress final file
        gzip ALL.${hap_suffix}.bamcov.bed
    else
        echo "No BAM files found for ${hap_suffix}"
    fi
done

echo "Processing complete for both haplotypes"
