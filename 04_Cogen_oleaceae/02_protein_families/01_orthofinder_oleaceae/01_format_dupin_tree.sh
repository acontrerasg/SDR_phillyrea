#!/usr/bin/env bash

# This script gets the tree from Dupin et al 2024 for orthofinder run. 

# We downloaded the supplementary material of Dupin et al 2024.
# From the file "v4b_fixedtree.xml" we extracted the newik tree 
# and name it species_tree.txt
# Env:
public_proteomes="PATH/TO/FOLDER/"
main_dir="PATH/TO/FOLDER/"
mkdir -p ${main_dir} && cd ${main_dir}

# modules
module load newick-utils/
# We need to modify the following entries so  they match between the newkick tree and the fasta files:
# Species names are between "C" and ":"
# Extract species names:
grep -o '[^(),:]\+:'  species_tree.txt | sed 's/:$//' > tree.names.txt
# store new names (same as proteome fasta names) in a file 
echo "Abeliophyllum_distichum
Forsythia_ovata
Forsythia_suspensa
Fraxinus_pennsylvanica
Syringa_oblata
Fraxinus_excelsior
Osmanthus_fragrans
Olea_europaea
Phillyrea_angustifolia
Phillyrea_latifolia" > proteome.names.txt
# Modify names in the tree:
while read sp; do
    sed -i -E "s/${sp}_[^,:)]*/${sp}/g" species_tree.txt
done < proteome.names.txt

# Verify extracting the names of the new tree
grep -o '[^(),:]\+:'  species_tree.txt |
sed 's/:$//' |
grep -f proteome.names.txt 

# Make sure species names are not repeated:
# so from :
# Fraxinus_excelsior
# Fraxinus_excelsior
# Fraxinus_excelsior
# to : 
# Fraxinus_excelsior
# Fraxinus_excelsior_2
# Fraxinus_excelsior_3
awk '
BEGIN { FS=":"; OFS=":" }
{
    for (i=1; i<NF; i++) {
        if (match($i, /[^(),;]+$/)) {
            name = substr($i, RSTART, RLENGTH)
            count[name]++
            if (count[name] > 1) {
                $i = substr($i, 1, RSTART-1) name "_" count[name]
            }
        }
    }
    print
}
' species_tree.txt > species_tree_renamed.nwk

# Keep only the branches belonigng to the species used:
nw_prune species_tree_renamed.nwk -vf proteome.names.txt > pruned_tree.nwk