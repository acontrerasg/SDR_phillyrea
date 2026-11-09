#!/usr/bin/env Rscript
library(tidyverse)
library(ggplot2)
library(patchwork)

# Paths
remote="PATH/TO/FOLDER"
path="PATH/TO/FOLDER"

# Find all files ending in .500000.SDR.bed
files <- list.files(path = paste0(remote,path), 
                    pattern = "\\.50000\\.SDR\\.bed$", 
                    full.names = TRUE)
# Read all files into a list
data_list <- lapply(files, function(file) {
  cat("Reading:", basename(file), "\n")
  vroom(file, col_names = c("chr", "start", "end", "num","bases","winsize", "freq","type","group")) %>% 
    mutate(group = if_else(group == ".", "background", group)) %>% 
    mutate(genome = str_remove(chr, pattern = "-.*"))
})



final_df <- bind_rows(data_list)

final_df <- final_df %>%  mutate(group = if_else(group == "ID=CA3p", "background", group)) %>%
  mutate(group = if_else(group == "ID=CA5p", "background", group)) %>%
  mutate(type = if_else(type == "LTR/Gypsy", "LTR/TY3", type)) %>%
  mutate(type = if_else(type == "LTR/Copia", "LTR/TY1", type)) %>%
  mutate(type = if_else(type == "gene", "GENE", type)) %>%
  mutate(group = sub("^ID=", "", group)) %>% 
  mutate(group = factor(
    group,
    levels = c(
      "background", 
      "M1", "m1", "TIB1",
      "M2", "m2", "B2",
      "M3", "m3", "TIB3",
      "M4", "m4", "TB4",
      "M5", "m5", "I5",
      "M6", "m6"
    )
  ))

color_scheme <- c("CA3p" = "#666666", 
                  "TIB1" = "#b7b3ca",
                  "M1" = "#5ab4ac",
                  "M2" = "#5ab4ac",
                  "M3" = "#5ab4ac",
                  "M4" = "#5ab4ac",
                  "M5" = "#5ab4ac",
                  "M6" = "#5ab4ac",
                  "m1" = "#d8b365",
                  "m2" = "#d8b365",
                  "m3" = "#d8b365",
                  "m4" = "#d8b365",
                  "m5" = "#d8b365",
                  "m6" = "#d8b365",
                  "B2" = "#729fcf", 
                  "TIB3" = "#b7b3ca", 
                  "TB4" = "#729fcf", 
                  "I5" ="#9999cc",
                  "CA5p" = "#666666",
                  "background" = "#666666")


# Calculate stats 
summary <- final_df %>%
  group_by(group, type, genome) %>%
  summarise(
    mean_num   = mean(num, na.rm = TRUE),
    sd_num     = sd(num, na.rm = TRUE),
    median_num = median(num, na.rm = TRUE),
    
    mean_bases   = mean(bases, na.rm = TRUE),
    sd_bases     = sd(bases, na.rm = TRUE),
    median_bases = median(bases, na.rm = TRUE),
    sum_bases    = sum(bases, na.rm = TRUE),
    
    mean_freq   = mean(freq, na.rm = TRUE),
    sd_freq     = sd(freq, na.rm = TRUE),
    median_freq = median(freq, na.rm = TRUE),
    .groups = "drop"
  )
write.table(summary,
            file = "PATH/TO/FILE",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)
####################################################


# PLot most relevant bits:
exclude_types <- c("Satellite", "satellite_DNA", "SINE", "five_prime_UTR",
                   "three_prime_UTR", "DNA", "exon", "mRNA", 
                   "LTR/Caulimovirus", "LTR", "CDS")

include_types <- c("DNA/DTM", "GENE", "LTR/TY1", "LTR/TY3")

## HAP1 ##
final_df %>% 
  filter(genome == "Pan_040.h1") %>%
  filter(!type %in% exclude_types) %>% 
  ggplot(aes(x=group,y=freq, color=group)) + 
  geom_boxplot(outliers = F) + 
  scale_color_manual(values = color_scheme) +
  facet_wrap(~type)  + 
  coord_flip() +
  xlab("") +
  ylab("Density in 0.5 Mbp") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "none")


a <- final_df %>%
  filter(genome == "Pan_040.h1") %>%
  filter(type %in% include_types) %>% 
  ggplot(aes(x=group,y=freq, color=group)) + 
  geom_boxplot(outliers = F) + 
  scale_color_manual(values = color_scheme) +
  facet_wrap(~type)  + 
  coord_flip() +
  xlab("") +
  ylab("Density in 0.5 Mbp") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "none")


## HAP2 ##

final_df %>% 
  filter(genome == "Pan_040.h2") %>%
  filter(!type %in% exclude_types) %>% 
  ggplot(aes(x=group,y=freq, color=group)) + 
  geom_boxplot(outliers = F) + 
  scale_color_manual(values = color_scheme) +
  facet_wrap(~type)  + 
  coord_flip() +
  xlab("") +
  ylab("Density in 0.5 Mbp") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "none")

b <- final_df %>%
  filter(genome == "Pan_040.h2") %>%
  filter(type %in% include_types) %>% 
  ggplot(aes(x=group,y=freq, color=group)) + 
  geom_boxplot(outliers = F) + 
  scale_color_manual(values = color_scheme) +
  facet_wrap(~type)  + 
  coord_flip() +
  xlab("") +
  ylab("Density in 0.5 Mbp") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "none")


a | b
