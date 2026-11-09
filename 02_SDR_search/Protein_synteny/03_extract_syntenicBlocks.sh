#!/usr/bin/env bash
# Here I will describe how I extracted the syntenic blocks from the 
# output of Genespace from the reference Genomes Pan_040 and Pan_0125

# Function ot modify  gff start and end coordinates.
# Start for the first entry, end for the last entry:
function modify_ends_gff() {
    local file="$1"
    local start_val="$2"
    local end_val="$3"
    # main 
    awk -F '\t' -v start="$start_val" -v end="$end_val" '
    BEGIN {OFS="\t"}
    {
        lines[NR] = $0;
    }
    END {
        for (i=1; i<=NR; i++) {
            $0 = lines[i];
            if (i == 1) $4 = start;
            if (i == NR) $5 = end;
            print;
        }
    }' "$file"
}

# Path of the genespace results: 
path="PATH/TO/FOLDER/"
cd ${path}
# Create output folder
out_dir=${path}"synteny_chrom12/"
mkdir -p ${out_dir}
# Get the syntenicBlck coordinats file and transform it to a tsv file for ease of handling
coordinate_file=${path}"syntenicBlock_coordinates.csv"
cat ${coordinate_file} | tr ',' '\t' > ${out_dir}"coordinate_genespace_file.tsv"
# move to the output folder
cd ${out_dir}
# Modify the format of the synteny block name and remofe the self blocks:
cat coordinate_genespace_file.tsv | grep -v "selfBlk" | sed 's/: /_/'  > temp
mv temp coordinate_genespace_file.tsv
# Get chromsome 12 only and Male haplotype synteny:
awk 'NR == 1 || $3 == "12" || $4 == "12"' coordinate_genespace_file.tsv | 
awk '{if($1 ~ "Male_hap" && $2 ~ "Male_hap") print $0}'  > coordinate_genespace_file_chrom12.tsv 

 # Get chromsome 12 only and Male haplotype synteny:
awk 'NR == 1 || $3 == "12" || $4 == "12"' coordinate_genespace_file.tsv | 
awk '{if($1 ~ "Herm_hap" && $2 ~ "Male_hap1") print $0}'  > coordinate_genespace_file_herm_chrom12.tsv 

# For each haplotype get the follwoing columns:
# chr1,startBp1,endBp1,orientation,blkID ;
# Modify back the chromosome name and add transform to a GFF file 

# Male hap 1 
awk '$1 == "Male_hap1"' coordinate_genespace_file_chrom12.tsv |
awk 'BEGIN{OFS=FS="\t"}{print "Pan_040.h1-"$3,"GENESPACE","SynBlk",$6+1,$7,".",$20,".","ID="$5}' |
bedtools  sort |
bedtools  merge -s -i - -c 2,3,6,7,8,9 -o distinct |
awk 'BEGIN{OFS=FS="\t"}{print $1,$4,$5,$2,$3,$6+1,$7,$8,$9}' |
sed 's/,ID=/,/g' > Pan_040.hap1-rt_SynBlk_chr12.gff3
# Now extend the  first and last record to match with the start and the end of the chromosome
# using the fai file
fai="PATH/TO/FILE"
start=1
end=$(grep "Pan_040.h1-12" $fai  | cut -f2)
# Apply function to change end and start 
modify_ends_gff Pan_040.hap1-rt_SynBlk_chr12.gff3  $start $end > mod.temp
mv mod.temp Pan_040.hap1-rt_SynBlk_chr12.gff3
## Now extract Nonsyntenic regions for M and m: 
bedtools complement -i Pan_040.hap1-rt_SynBlk_chr12.gff3 -g $fai |
grep "Pan_040.h1-12" |
awk 'BEGIN{OFS=FS="\t"}{print $1,"GENESPACE","NonSynBlk",$2,$3,".",".",".","ID=NonSynBlk_M_"NR}'  > Pan_040.hap1-rt_NonSynBlk.gff3
############################################

# Male hap2 
awk '$1 == "Male_hap2" && $2 == "Male_hap1"' coordinate_genespace_file_chrom12.tsv |
awk 'BEGIN{OFS=FS="\t"}{print "Pan_040.h2-"$3,"GENESPACE","SynBlk",$6+1,$7,".",$20,".","ID="$5}' |
bedtools sort |
bedtools merge -s -i - -c 2,3,6,7,8,9 -o distinct |
awk 'BEGIN{OFS=FS="\t"}{print $1,$4,$5,$2,$3,$6+1,$7,$8,$9}' |
sed 's/,ID=/,/g' > Pan_040.hap2-rt_SynBlk_chr12.gff3 
# Now extend the  first and last record to match with the start and the end of the chromosome
# using the fai file
fai="PATH/TO/FILE"
start=1
end=$(grep "Pan_040.h2-12" $fai  | cut -f2)
# Apply function to change end and start 
modify_ends_gff Pan_040.hap2-rt_SynBlk_chr12.gff3 $start $end > mod.temp
mv mod.temp Pan_040.hap2-rt_SynBlk_chr12.gff3
## Now extract Nonsyntenic regions for M and m: 
bedtools complement -i Pan_040.hap2-rt_SynBlk_chr12.gff3 -g $fai |
grep "Pan_040.h2-12" |
awk 'BEGIN{OFS=FS="\t"}{print $1,"GENESPACE","NonSynBlk",$2,$3,".",".",".","ID=NonSynBlk_m_"NR}' \
> Pan_040.hap2-rt_NonSynBlk.gff3
############################################

# Hermaphrodite hap1
awk '$1 == "Herm_hap1" && $2 == "Male_hap1"' coordinate_genespace_file_herm_chrom12.tsv |
awk 'BEGIN{OFS=FS="\t"}{print "Pan_01N25.h1-"$3,"GENESPACE","SynBlk",$6+1,$7,".",$20,".","ID="$5}' |
bedtools sort |
bedtools merge -s -i - -c 2,3,6,7,8,9 -o distinct |
awk 'BEGIN{OFS=FS="\t"}{print $1,$4,$5,$2,$3,$6+1,$7,$8,$9}' |
sed 's/,ID=/,/g' > Pan_01N25.hap1-rt_SynBlk_chr12.gff3 
# Now extend the  first and last record to match with the start and the end of the chromosome
# using the fai file
fai="PATH/TO/FILE"
start=1
end=$(grep "Pan_01N25.h1-12" $fai  | cut -f2)
# Apply function to change end and start 
modify_ends_gff Pan_01N25.hap1-rt_SynBlk_chr12.gff3 $start $end > mod.temp
mv mod.temp Pan_01N25.hap1-rt_SynBlk_chr12.gff3 
## Now extract Nonsyntenic regions for M and m: 
bedtools complement -i Pan_01N25.hap1-rt_SynBlk_chr12.gff3 -g $fai |
grep "Pan_01N25.h1-12" |
awk 'BEGIN{OFS=FS="\t"}{print $1,"GENESPACE","NonSynBlk",$2,$3,".",".",".","ID=NonSynBlk_m_"NR}' \
> Pan_01N25.hap1-rt_NonSynBlk.gff3
############################################

# Hermaphrodite hap2
awk '$1 == "Herm_hap2" && $2 == "Male_hap1"' coordinate_genespace_file_herm_chrom12.tsv |
awk 'BEGIN{OFS=FS="\t"}{print "Pan_01N25.h2-"$3,"GENESPACE","SynBlk",$6+1,$7,".",$20,".","ID="$5}' |
bedtools sort |
bedtools merge -s -i - -c 2,3,6,7,8,9 -o distinct |
awk 'BEGIN{OFS=FS="\t"}{print $1,$4,$5,$2,$3,$6+1,$7,$8,$9}' |
sed 's/,ID=/,/g' > Pan_01N25.hap2-rt_SynBlk_chr12.gff3 
# Now extend the  first and last record to match with the start and the end of the chromosome
# using the fai file
fai="PATH/TO/FILE"
start=1
end=$(grep "Pan_01N25.h2-12" $fai  | cut -f2)
# Apply function to change end and start 
modify_ends_gff Pan_01N25.hap2-rt_SynBlk_chr12.gff3 $start $end > mod.temp
mv mod.temp Pan_01N25.hap2-rt_SynBlk_chr12.gff3 
## Now extract Nonsyntenic regions for M and m: 
bedtools complement -i Pan_01N25.hap2-rt_SynBlk_chr12.gff3 -g $fai |
grep "Pan_01N25.h2-12" |
awk 'BEGIN{OFS=FS="\t"}{print $1,"GENESPACE","NonSynBlk",$2,$3,".",".",".","ID=NonSynBlk_m_"NR}' \
> Pan_01N25.hap2-rt_NonSynBlk.gff3

####################################################
# Combine syntenic and non syntenic blocks
for haplotype in Pan_01N25.hap1 Pan_01N25.hap2 Pan_040.hap1 Pan_040.hap2 ; do
    cat ${haplotype}-rt_SynBlk_chr12.gff3  ${haplotype}-rt_NonSynBlk.gff3 |
    bedtools sort  > ${haplotype}-rt_combined_synteny_chr12.gff3
done 

# Make an associative array to add the new names to the GFF3 files
declare -A map=(
  ["ID=Male_hap1_vs_Male_hap2_43"]="ID=CA5p"
  ["ID=NonSynBlk_M_1"]="ID=M1"
  ["ID=Male_hap1_vs_Male_hap2_47"]="ID=TIB1"
  ["ID=NonSynBlk_M_2"]="ID=M2"
  ["ID=Male_hap1_vs_Male_hap2_46"]="ID=B2"
  ["ID=NonSynBlk_M_3"]="ID=M3"
  ["ID=Male_hap1_vs_Male_hap2_45"]="ID=TIB3"
  ["ID=NonSynBlk_M_4"]="ID=M4"
  ["ID=Male_hap1_vs_Male_hap2_44"]="ID=TB4"
  ["ID=NonSynBlk_M_5"]="ID=M5"
  ["ID=Male_hap1_vs_Male_hap2_48"]="ID=I5"
  ["ID=NonSynBlk_M_6"]="ID=M6"
  ["ID=Male_hap1_vs_Male_hap2_49"]="ID=CA3p"

  ["ID=Male_hap1_vs_Male_hap2_43"]="ID=CA5p"
  ["ID=NonSynBlk_m_1"]="ID=m1"
  ["ID=Male_hap1_vs_Male_hap2_44"]="ID=TB4"
  ["ID=NonSynBlk_m_2"]="ID=m2"
  ["ID=Male_hap1_vs_Male_hap2_45"]="ID=TIB3"
  ["ID=NonSynBlk_m_3"]="ID=m3"
  ["ID=Male_hap1_vs_Male_hap2_46"]="ID=B2"
  ["ID=NonSynBlk_m_4"]="ID=m4"
  ["ID=Male_hap1_vs_Male_hap2_47"]="ID=TIB1"
  ["ID=NonSynBlk_m_5"]="ID=m5"
  ["ID=Male_hap1_vs_Male_hap2_48"]="ID=I5"
  ["ID=NonSynBlk_m_6"]="ID=m6"
  ["ID=Male_hap1_vs_Male_hap2_49"]="ID=CA3p"

  ["ID=Male_hap1_vs_Herm_hap1_42"]="ID=CA5p"
  ["ID=NonSynBlk_m_1"]="ID=m1"
  ["ID=Male_hap1_vs_Herm_hap1_43"]="ID=TB4"
  ["ID=NonSynBlk_m_2"]="ID=m2"
  ["ID=Male_hap1_vs_Herm_hap1_44"]="ID=TIB3"
  ["ID=NonSynBlk_m_3"]="ID=m3"
  ["ID=Male_hap1_vs_Herm_hap1_45"]="ID=B2"
  ["ID=NonSynBlk_m_4"]="ID=m4"
  ["ID=Male_hap1_vs_Herm_hap1_46"]="ID=TIB1"
  ["ID=NonSynBlk_m_5"]="ID=m5"
  ["ID=Male_hap1_vs_Herm_hap1_47"]="ID=I5"
  ["ID=NonSynBlk_m_6"]="ID=m6"
  ["ID=Male_hap1_vs_Herm_hap1_48"]="ID=CA3p"

  ["ID=Male_hap1_vs_Herm_hap2_35"]="ID=CA5p"
  ["ID=NonSynBlk_m_1"]="ID=m1"
  ["ID=Male_hap1_vs_Herm_hap2_36"]="ID=TB4"
  ["ID=NonSynBlk_m_2"]="ID=m2"
  ["ID=Male_hap1_vs_Herm_hap2_37"]="ID=TIB3"
  ["ID=NonSynBlk_m_3"]="ID=m3"
  ["ID=Male_hap1_vs_Herm_hap2_38"]="ID=B2"
  ["ID=NonSynBlk_m_4"]="ID=m4"
  ["ID=Male_hap1_vs_Herm_hap2_39"]="ID=TIB1"
  ["ID=NonSynBlk_m_5"]="ID=m5"
  ["ID=Male_hap1_vs_Herm_hap2_40"]="ID=I5"
  ["ID=NonSynBlk_m_6"]="ID=m6"
  ["ID=Male_hap1_vs_Herm_hap2_41"]="ID=CA3p"
)

# apply the changes to every syntenic file
for haplotype in Pan_01N25.hap1 Pan_01N25.hap2 Pan_040.hap1 Pan_040.hap2 ; do
    while IFS=$'\t' read -r c1 c2 c3 c4 c5 c6 c7 c8 c9 rest; do
        # If column 9 matches one of the keys, replace it
        if [[ -n "${map[$c9]}" ]]; then
            c9="${map[$c9]}"
        fi
        # Rebuild the line and print
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" \
               "$c1" "$c2" "$c3" "$c4" "$c5" "$c6" "$c7" "$c8" "$c9"
        # Print the rest of the columns if any
        if [[ -n "$rest" ]]; then
            printf "\t%s" "$rest"
        fi
        printf "\n"
    done < ${haplotype}-rt_combined_synteny_chr12.gff3 > ${haplotype}-rt_annotated_SDR_synteny.gff3 ;
done 
