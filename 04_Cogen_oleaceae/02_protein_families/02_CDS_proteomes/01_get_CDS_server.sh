#!/usr/bin/env bash

# This  script retrieves the CDS of the public anotations of the 
# Load module 
module load seqkit/ 
# Set path and dir
dir="PATH/TO/FOLDER/"
mkdir -p ${dir} ; 
cd ${dir} ; 

# GET CDS from the server and format them to be same as orthofinder run:
# Abeliophyllum_distichum
Abe_di="PATH/TO/FILE"
# Forsythia_ovata
For_ov="PATH/TO/FILE"
# Forsythia_suspensa
For_sus="PATH/TO/FILE"
# Fraxinus pennsylvatica
Fra_pen="Fpennsylvanica_v1.4_genes_CDS.fa.gz"
# Syringa Oblata
Syr_ob="PATH/TO/FILE"
# Fraxinus excelsior
Fra_ex="PATH/TO/FILE"
# Osmanthus fragrans L
Osm_fr="PATH/TO/FILE"
# Oleae europeaea
Ole_eu="PATH/TO/FILE"
# Phillyrea angustifolia
Phi_an_dir="PATH/TO/FOLDER/"
# Phillyrea latifolia
Phi_la_dir="PATH/TO/FOLDER/"

# Public CDS formatting 
# Abeliophyllum_distichum
cat ${Abe_di} |
sed 's/>/>Adi_/;s/ .*$//' > Abeliophyllum_distichum.fasta
# Forsythia_ovata
cat ${For_ov} |
sed 's/>/>Fov_/;s/ .*$//' > Forsythia_ovata.fasta
# Forsythia_suspensa
zcat ${For_sus} |
sed 's/>.*Protein=/>Fsu_/;s/\t.*$//' > Forsythia_suspensa.fasta
# Fraxinus pennsylvanica
zcat ${Fra_pen} |
sed 's/>/>Fpe_/;s/ .*$//' > Fraxinus_pennsylvanica.fasta
# Syringa_oblata
cat ${Syr_ob} |
seqkit grep -nrp ".1 " |
sed 's/>/>Sob_/;s/ .*$//' > Syringa_oblata.fasta
# Fraxinus excelsior
cat ${Fra_ex} |
seqkit grep -nrp ".1 " |
sed 's/>/>Fex_/;s/ .*$//' > Fraxinus_excelsior.fasta
# Osmanthus fragrans L 
cat ${Osm_fr} |
sed 's/>/>Ofr_/;s/\t.*$//' > Osmanthus_fragrans.fasta
# Olea_europaea_unver
cat ${Ole_eu}  |
seqkit grep -nrp ".1 " |
sed 's/>/>Oeu_unv_/;s/ .*$//' > Olea_europaea.fasta
# Phillyrea angustifolia
cat ${Phi_an_dir}Pan_040.hap1-rt.genes.cds.fasta |
sed 's/>/>Pan_/;s/ .*$//' > Phillyrea_angustifolia.fasta
# Phillyrea latifolia
cat ${Phi_la_dir}Pla.HiFiasm.hap1_rv_n.liftoff.cds.fasta  | 
sed 's/>H._/>Pla_/;s/ .*$//'  > temp1
cat ${Phi_la_dir}Pla.HiFiasm.hap2_rv_n.liftoff.cds.fasta | 
sed 's/>H._/>Pla_/;s/ .*$//'  > temp2
cat temp1 temp2 | seqkit rmdup -n | seqkit rmdup -s  >  Phillyrea_latifolia.fasta
rm temp1 temp2
