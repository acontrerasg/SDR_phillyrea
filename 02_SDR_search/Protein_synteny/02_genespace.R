#!/usr/bin/env Rscript
# I ran GENESPACE locally
# Load required library
library(GENESPACE)

wd="PATH/TO/FOLDER/"
path2mcscanx="PATH/TO/FOLDER"
path2diamond="PATH/TO/FOLDER/"
path2orthofinder="PATH/TO/FOLDER/"
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  path2diamond = path2diamond,
  path2orthofinder = path2orthofinder)
# Run the synteny analysis
gs <- run_genespace(gpar)