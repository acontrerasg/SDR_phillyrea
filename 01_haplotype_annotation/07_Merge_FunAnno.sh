#!/usr/bin/env bash

# Get paths
source ./00_paths.conf
# DIR
dir=${Main_anno_dir}${merge_anno_subdir}
mkdir -p ${dir}
cd ${dir}

# Araport11 Results
Pan_01N25_h1_Ara11="PATH/TO/FILE"
Pan_01N25_h2_Ara11="PATH/TO/FILE"
Pan_040_h1_Ara11="PATH/TO/FILE"
Pan_040_h2_Ara11="PATH/TO/FILE"

# OLEA9
Pan_01N25_h1_OLEA9="PATH/TO/FILE"
Pan_01N25_h2_OLEA9="PATH/TO/FILE"
Pan_040_h1_OLEA9="PATH/TO/FILE"
Pan_040_h2_OLEA9="PATH/TO/FILE"

# Interpro Results
Pan_01N25_h1_IP="PATH/TO/FILE"
Pan_01N25_h2_IP="PATH/TO/FILE"
Pan_040_h1_IP="PATH/TO/FILE"
Pan_040_h2_IP="PATH/TO/FILE"

# Pan 01N25
python3 ${run_merge_sourc} Pan_01N25_h1_FunAnno.tsv ${Pan_01N25_h1_Ara11} ${Pan_01N25_h1_OLEA9} ${Pan_01N25_h1_IP}
python3 ${run_merge_sourc} Pan_01N25_h2_FunAnno.tsv ${Pan_01N25_h2_Ara11} ${Pan_01N25_h2_OLEA9} ${Pan_01N25_h2_IP}
# Pan 040
python3 ${run_merge_sourc} Pan_040_h1_FunAnno.tsv ${Pan_040_h1_Ara11} ${Pan_040_h1_OLEA9} ${Pan_040_h1_IP}
python3 ${run_merge_sourc} Pan_040_h2_FunAnno.tsv ${Pan_040_h2_Ara11} ${Pan_040_h2_OLEA9} ${Pan_040_h2_IP}

# Alter the Inital ID so it matches with the Annotation
sed -i 's/^/ID=/' *_FunAnno.tsv 

# Merge Funcitonal Annotation with the GFF 

Pan_01N25_hap1_complete="PATH/TO/FILE"
Pan_01N25_hap2_complete="PATH/TO/FILE"
Pan_040_hap1_complete="PATH/TO/FILE"
Pan_040_hap2_complete="PATH/TO/FILE"

# Pan_01N25
cat ${Pan_01N25_hap1_complete} > local_copy.gff3 
python3 ${run_make_gff} local_copy.gff3 Pan_01N25_h1_FunAnno.tsv Pan_01N25.hap1-rt.complete.FunAnno.gff3 ;
rm local_copy.gff3 ;
cat ${Pan_01N25_hap2_complete} > local_copy.gff3 
python3 ${run_make_gff} local_copy.gff3 Pan_01N25_h2_FunAnno.tsv Pan_01N25.hap2-rt.complete.FunAnno.gff3
rm local_copy.gff3 ;
# Pan_040
cat ${Pan_040_hap1_complete} > local_copy.gff3 
python3 ${run_make_gff} local_copy.gff3 Pan_040_h1_FunAnno.tsv Pan_040.hap1-rt.complete.FunAnno.gff3
rm local_copy.gff3 ;
cat ${Pan_040_hap2_complete} > local_copy.gff3 
python3 ${run_make_gff} local_copy.gff3 Pan_040_h2_FunAnno.tsv Pan_040.hap2-rt.complete.FunAnno.gff3
rm local_copy.gff3 ;