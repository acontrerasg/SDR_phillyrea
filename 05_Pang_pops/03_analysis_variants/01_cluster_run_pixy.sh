#!/bin/bash
#SBATCH --partition=COMPUTE
#SBATCH -c 8            # number of cores
#SBATCH --mem=40G   # requested memory
#SBATCH --job-name=RunPixy_VCF
#SBATCH -o PATH/TO/FILE               # File to which standard output will be written
#SBATCH -e PATH/TO/FILE               # File to which standard error will be written

# Path
module load pixy/2.0.0.beta14
dir="PATH/TO/FOLDER/"
mkdir -p ${dir} && cd ${dir}

# Variables
run_pixy="PATH/TO/FILE"
popfile="PATH/TO/FILE"
win_size=100000
threads=8

# Create the pop file and select only Fabregues pop:
popfile="PATH/TO/FILE"
cat ${popfile} | grep "Fabregues" | tail -n +2 | awk '{print $1"\t"$2}' > pop_file_pixy.txt

popfile="pop_file_pixy.txt"

# Run pixy all
vcf_dir="PATH/TO/FOLDER/"

for vcf in ${vcf_dir}Pan_040.h2-12.snps.vcf.gz ; do 
    echo ${vcf} ;
    name=$(basename ${vcf} .vcf.gz) ;
	for win_size in  10000 50000 100000 ; do
		echo $win_size ;
		pixy --stats pi dxy fst watterson_theta tajima_d \
     		--vcf $vcf \
     		--populations $popfile \
     		--window_size $win_size \
     		--n_cores $threads \
     		--output_prefix pixy_${name}_${win_size} \
     		--output_folder ./ ;
	done ;
done

for vcf in ${vcf_dir}Pan_040.h*.vcf.gz ; do 
    echo ${vcf} ;
    name=$(basename ${vcf} .vcf.gz) ;
	for win_size in  10000 50000 100000 ; do
		echo $win_size ;
		pixy --stats pi dxy fst watterson_theta tajima_d \
     		--vcf $vcf \
     		--populations $popfile \
     		--window_size $win_size \
     		--n_cores $threads \
     		--output_prefix pixy_${name}_${win_size} \
     		--output_folder ./ ;
	done ;
done
