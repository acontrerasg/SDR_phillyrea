# Load required libraries
library(tidyverse)
library(pafr)

# Run this in  the shell 
# Define file names
fasta_file1 <- "Pan_040.h1-rt_ALL_SDR.fasta"
fasta_file2 <- "Pan_040.h2-rt_ALL_SDR.fasta"

# Construct the shell command
cmd <- sprintf(
  'for i in $(grep ">" %s | sed "s/>//"); do
    seqkit grep -nrp $i %s > temp1.fa
    seqkit grep -nrp $i %s > temp2.fa
    minimap2 -c --secondary=no -x asm20 temp1.fa temp2.fa > ${i}_align.paf
    rm temp1.fa temp2.fa
done', fasta_file1, fasta_file1, fasta_file2)

# Run the shell command from R
setwd("PATH/TO/FOLDER")
system(cmd, intern = TRUE)

# Calculate  average compress gap divergence
# de column contains this metric, which is explained here:
# https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
process_paf <- function(name, input_paf) {
  read_paf(input_paf) %>%
    as_tibble() %>%
    summarise(
      mean_de = mean(de, na.rm = TRUE),
      median_de = median(de, na.rm = TRUE),
      sd_de = sd(de, na.rm = TRUE),
      count_ali = n(),
      count_na_de = sum(is.na(de))  # Count number of NAs in 'de'
    ) %>%
    mutate(name = name) %>%
    relocate(name)  # Move 'name' column to the first position
}

# Read sequences
file_path="PATH/TO/FOLDER/"
results_divergence <- bind_rows(
  process_paf("CA3p", paste0(file_path,"ID=SynBlk_1_align.paf")),
  process_paf("TIB1", paste0(file_path,"ID=SynBlk_2_align.paf")),
  process_paf("B2", paste0(file_path,"ID=SynBlk_3_align.paf")),
  process_paf("TIB3", paste0(file_path,"ID=SynBlk_4_align.paf")),
  process_paf("TB4", paste0(file_path,"ID=SynBlk_5_align.paf")),
  process_paf("I5", paste0(file_path,"ID=SynBlk_6_align.paf")),
  process_paf("CA5p", paste0(file_path,"ID=SynBlk_7_align.paf"))
)

# Set factor levels to control order
results_divergence$name <- factor(results_divergence$name, levels = c("CA3p", "TIB1", "B2", "TB3", "IB4", "CA5p"))
#http://127.0.0.1:11751/graphics/plot_zoom_png?width=922&height=411

# Plot using ggplot2
results_divergence  %>%
  ggplot( aes(x = name, y = median_de)) +
  geom_point(size = 2.5) +
  geom_segment(aes(x = name, xend = name, 
                   y = median_de - sd_de, 
                   yend = median_de + sd_de), 
               color = "black", linewidth = 1) +  # Straight line for SD
  ggthemes::theme_few(base_size = 16) +
  labs(
    x = "",
    y = "gap compress divergence")