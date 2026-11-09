#!/usr/bin/env Rscript
library(readr)
library(pafr)
library(dplyr)
# read the linkage map
remote="PATH/TO/FOLDER"
path_linkage="all.markers.contigpos.bed"
linkage <- read_table(paste0(remote,path_linkage), col_names = c("chr", "start","end", "marker"))
path_linkage="markers_sex_link.bed"
linkage_sex <- read_table(paste0(remote,path_linkage), col_names = c("chr", "start","end", "marker"))
# read the paf file (standard 12-column minimap2 output)
path_paf="path/to/file.paf"
alignment <- pafr::read_paf(file_name =paste0(remote,path_paf))

chr_to_plot <- "chr1"
linkage_chr <- linkage %>% filter(chr == "Pan_040.h1-12") %>%
linkage_sex_chr <- linkage_sex  %>% filter(chr == "Pan_040.h1-12")
linkage_sex_chr2 <- linkage_sex  %>% filter(chr == "Pan_040.h2-12")
paf_chr <- alignment %>%  filter(qname =="Pan_040.h2-12") %>%  
  filter(tname == "Pan_040.h1-12") %>%  filter(mapq == 60 & tp == "P") 

start=0
endt=paf_chr %>%  as_tibble() %>%  .$tlen %>%  max
endq=paf_chr %>%  as_tibble() %>%  .$qlen %>%  max

# Plot the y-axis (ylim) coordinates.
plot(NULL, xlim = c(start, endt), ylim = c(0,4), axes=F, ylab = '', xlab = '')
# Draw two simple rectangles.
rect(0, 3.3, end, 3.5, col =  "#5ab4ac", border = NA)     # blue
rect(0, 2.3, endq, 2.5, col = "#d8b365", border = NA)    # vermillion
# Add info
text(start, 3.7, substitute(paste(italic('Chr12-h1 "M" '))), pos = 4)
text(start, 2, substitute(paste(italic('Chr12-h2 "m" '))), pos = 4)
# Plot the polygons representing the alignments.
for(e in 1:nrow(paf_chr)){
  polygon(x = c(paf_chr$tstart[e], paf_chr$tend[e], paf_chr$qend[e], paf_chr$qstart[e]), 
          y = c(3.3,3.3,2.5,2.5),
          col = adjustcolor('grey40', alpha.f = 0.2), border = FALSE)
}
# Add the markers
points(linkage_sex_chr$start, rep(3.4, nrow(linkage_sex_chr)),
       pch = 6, cex =1, col = "red")

points(linkage_sex_chr2$start, rep(2.4, nrow(linkage_sex_chr2)),
       pch = 2, cex =1, col = "blue")
