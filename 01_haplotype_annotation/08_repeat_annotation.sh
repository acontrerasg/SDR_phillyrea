#!/usr/bin/env bash

# This script sets the Satellite Repeat Annotation subfolder and main files used for the repeat annotation of P. angustifolia.
# Source the paths used in this script
source ./00_paths.conf


mkdir -p ${Repeat_anno_dir}
cd ${Repeat_anno_dir}

# Load needed modules
module load trf/4.09.1 
module load trash/1.2

# Step 1: Set variables
trf2bed="PATH/TO/FILE"
# EDTA settings for TRF:
align_score=1000 # -e para, dft:1000
max_seed=2000  # maximum period size to report
threads=24

# Step 2: Run TRASH with default settings
output_dir_trash=${Repeat_anno_dir}"/01_TRASH_anno/"
mkdir -p ${output_dir_trash}
cd ${output_dir_trash}

# Primary Assemblies
for genome in ${Repeat_anno_dir}*fasta ; do
	name=$(basename ${genome} .fasta) ;
 	mkdir -p ${name} ;
  	cd ${name}
	TRASH_run.sh ${genome} --par ${threads} ;
 	cd - ;
done
# Fix funny oversigth with the produced GFF file by trash. If you cat the result gff file,
# it will be all in the same line with the "^M" carrdige return delimiter between lines.
# To fix it, run sed 's/^M/\n/g' over the gff file. (to add the ^M itr needs to be added amnually 
# as ctr+v and ctr+M. Ref: https://stackoverflow.com/questions/5843495/what-does-the-m-character-mean-in-vim/

for i in  */*gff ; do
	name=$(basename $i) ;
 	sed 's/^M/\n/g' $i > $name ;
done 

# Step 3: Run TRF
output_dir_trf=${Repeat_anno_dir}"/02_TRF_anno/"
mkdir -p ${output_dir_trf}
cd ${output_dir_trf}

# Primary assemblies
for genome in ${Repeat_anno_dir}*fasta ; do
	# TRF cannot deal with a compress genome	
	genome_name=$(basename ${genome} .fasta) ;
	# MAIN TRF RUN 
	trf ${genome} 2 7 7 80 10 $align_score $max_seed -h -l 6 -m -ngs -d > ${genome_name}.2.7.7.80.10.1000.2000.dat
	python3 ${trf2bed} ${genome_name}.2.7.7.80.10.1000.2000.dat |
	sort -k1,1V -k4,4n -k5,5n > ${genome_name}.2.7.7.80.10.1000.2000.gff
done


for i in *gff ; do 
	name=$(basename $i .gff)
 	bedtools merge -i $i -d 50 |
	awk 'BEGIN{OFS="\t"} {print $1,"TRF","Satellite",$2+1,$3,".",".",".","ID=TRF_Sat_"NR}' > ${name}_Satellite.gff ; 
done

# Merge both anotations in the Main directory
cd ${Repeat_anno_dir}

for genome in *fasta ; do 
	name=$(basename $genome .fasta) ;
 	cat 01_TRASH_anno/TRASH_${name}.fasta.gff 02_TRF_anno/${name}.2.7.7.80.10.1000.2000_Satellite.gff |
  	sort -k1,1V -k4,4n -k5,5n > ${name}_Satellite.gff ;
done 