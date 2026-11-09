#!/usr/bin/env Rscript
# Load required libraries
library(tidyverse)
library(pafr)

# Calculate  average compress gap divergence
# de column contains this metric, which is explained here:
# https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
process_paf <- function(name, input_paf) {
  vroom(input_paf, col_names = F, delim = "\t") %>%
    filter(X12 == 60) %>% 
    select(name = X1, de = X21) %>%
    mutate(de = as.numeric(str_remove(de, "de:f:"))) %>% 
    summarise(
      mean_de = mean(de, na.rm = TRUE),
      median_de = median(de, na.rm = TRUE),
      sd_de = sd(de, na.rm = TRUE),
      count_ali = n(),
      count_na_de = sum(is.na(de))  # Count number of NAs in 'de',
    ) %>%
    mutate(name = name) %>%
    relocate(name)  # Move 'name' column to the first position
}

# Read sequences
remote_path="PATH/TO/FOLDER"
divergence_path="PATH/TO/FOLDER/" # nolint

results_divergence <- bind_rows(
  process_paf("B2", paste0(remote_path,divergence_path,"B2_040_hap_comp_align.paf")),
  process_paf("CA3p", paste0(remote_path,divergence_path,"CA3p_040_hap_comp_align.paf")),
  process_paf("CA5p", paste0(remote_path,divergence_path,"CA5p_040_hap_comp_align.paf")),
  process_paf("I5", paste0(remote_path,divergence_path,"I5_040_hap_comp_align.paf")),
  process_paf("TB4", paste0(remote_path,divergence_path,"TB4_040_hap_comp_align.paf")),
  process_paf("TIB1", paste0(remote_path,divergence_path,"TIB1_040_hap_comp_align.paf")),
  process_paf("TIB3", paste0(remote_path,divergence_path,"TIB3_040_hap_comp_align.paf"))
)

# Set factor levels to control order
results_divergence$name <- factor(results_divergence$name, levels = c("CA3p", "TIB1", "B2", "TIB3", "TB4","I5", "CA5p"))

# Plot using ggplot2

# Plot using ggplot2
results_divergence %>%
  ggplot(aes(x = name, y = median_de)) +
  geom_segment(aes(x = name, xend = name, 
                   y = median_de - sd_de, 
                   yend = median_de + sd_de), 
               color = "grey40", linewidth = 0.5) +
  geom_point(size = 3, color = "#2E86AB", alpha = 1) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text( margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "040 'M' vs 'm' divergence \n between syntenic blocks",
    x = "",
    y = "Haplotype Divergence"
  )