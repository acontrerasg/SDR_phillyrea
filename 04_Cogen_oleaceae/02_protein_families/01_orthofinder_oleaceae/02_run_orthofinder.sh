#!/usr/bin/env bash

# This scirpt  run orthofinder agaisnt a set of oleaceae proteomes using Dupin tree
# Load module
module load orthofinder/3.1.3
module load seqkit/

threads=12
# Env:
public_proteomes="PATH/TO/FOLDER/"
main_dir="PATH/TO/FOLDER/"
mkdir -p ${main_dir} && cd ${main_dir}
# Manually edit all  proteoms so they start with the initials of the species;
#  ;s/ .*$//;s/\*//'
#  ;s/ .*$//'
# Phillyrea
# Phillyrea angustifolia
#Phi_an_dir="PATH/TO/FOLDER/"
# Phillyrea latifolia
#Phi_la_dir="PATH/TO/FOLDER/"
# Phillyrea angustifolia
#cat ${Phi_an_dir}Pan_040.hap1-rt.genes.aa.fasta |
#sed 's/>/>Pan_/;s/ .*$//' > Phillyrea_angustifolia.fasta
# Phillyrea latifolia
#cat ${Phi_la_dir}Pla.HiFiasm.hap1_rv_n.liftoff.genes.fasta |
#sed 's/>H._/>Pla_/;s/ .*$//'  > temp1
#cat ${Phi_la_dir}Pla.HiFiasm.hap2_rv_n.liftoff.genes.fasta |
#sed 's/>H._/>Pla_/;s/ .*$//'  > temp2
#cat temp1 temp2 | seqkit rmdup -n | seqkit rmdup -s  >  Phillyrea_latifolia.fasta
#rm temp1 temp2
##########################################
# Public proteomes
# Abeliophyllum_distichum
cat ${public_proteomes}Abeliophyllum_distichum.Genemodel.prot.fasta |
sed 's/>/>Adi_/;s/ .*$//' > Abeliophyllum_distichum.fasta
# Forsythia_ovata
cat ${public_proteomes}Forsythia_ovata.Genemodel.prot.fasta  |
sed 's/>/>Fov_/;s/ .*$//' > Forsythia_ovata.fasta
# Fraxinus pennsylvanica
cat ${public_proteomes}Fpennsylvanica_v1.4_genes_AA.fa |
sed 's/>/>Fpe_/;s/ .*$//' > Fraxinus_pennsylvanica.fasta
# Syringa_oblata
cat ${public_proteomes}oblata_prot.fasta |
seqkit grep -nrp ".1 " |
sed 's/>/>Sob_/;s/ .*$//' > Syringa_oblata.fasta
# Fraxinus excelsior
cat ${public_proteomes}Fraxinus_excelsior_38873_TGAC_v2.longestCDStranscript.gff3.pep.fa |
seqkit grep -nrp ".1 " |
sed 's/>/>Fex_/;s/ .*$//' > Fraxinus_excelsior.fasta
# Osmanthus fragrans R OFR (Rixianggui) as it is the male
cat ${public_proteomes}OFR.pep |
sed 's/>/>Ofr_/;s/\t.*$//' > Osmanthus_fragrans.fastaForsythia_suspensa
# Olea_europaea_unver
cat ${public_proteomes}Oeuropaea_451_v1.0.protein.fa  |
seqkit grep -nrp ".1 " |
sed 's/>/>Oeu_unv_/;s/ .*$//' > Olea_europaea.fasta
# Forsythia_suspensa
cat ${public_proteomes}GWHERDT00000000.Protein.faa |
sed 's/>/>Fsu_/;s/ .*$//' > Forsythia_suspensa.fasta
#############################################3
# Not used olive genomes:
# Arabidopsis thaliana
#cat ${public_proteomes}Araport11_pep_20250411.fa |
#seqkit grep -nrp ".1 " |
#sed 's/>/>Ath_/;s/ .*$//' > Arabidopsis_thaliana.fasta
# Olea europaea leccino
#cat ${public_proteomes}GWHEUUU00000000.1.Protein.faa  |
#sed 's/>/>Oeu_lec_/;s/ .*$//' > Oeu_lec.fasta
# Olea_europaea_picual
#cat ${public_proteomes}Olea_europaea_proteins_v061.fasta  |
#sed 's/>/>Oeu_pic_/;s/ .*$//'\ > Oeu_pic.fasta
# Olea europaea arbequina
#cat ${public_proteomes}GWHAOPM00000000.Protein.faa |
#sed 's/>/>Oeu_arb_/;s/ .*$//' > Oeu_arb.fasta
# Olea_farga
#cat ${public_proteomes}OE9D.longest_peptide.fa |
#sed 's/>/>Ofa_/;s/ .*$//' > Ofa.fasta
# Osmanthus fragans L
#cat ${public_proteomes}OFL.pep |
#sed 's/>/>Ofr_/;s/ .*$//'\ > Ofr.fasta

##############################################################

# RUn Orthofinder
species_tree="PATH/TO/FILE"
orthofinder -X -f ${main_dir} -s $species_tree -n "w_dupin_tree" -S diamond_nosensitive
