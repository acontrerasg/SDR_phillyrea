#!/usr/bin/env Rscrip
library(patchwork)
library(tidyverse)
library(vroom)
library(gtools)
library(ggpubr)
# Load up  Radseq coverage files:
# Path
remote="PATH/TO/FOLDER"
local="/PATH/TO/FOLDER/"
path="PATH/TO/FOLDER/"

# hap 1 coverage
hap1_radseq_coverage <- vroom(paste0(remote,local,"ALL.Pan_040.hap1-rt.bamcov.bed.gz"), col_names = T)

# Cleanup Sample colunn
hap1_radseq_coverage <- hap1_radseq_coverage %>%
  mutate(Sample = str_remove(Sample, "_Pan_040.hap1-rt.fasta_radtag_paired"))

# Load up phenotypic data:

phenotype <- vroom(paste0(remote,local,"Metadata/samples_pheno_parents_Sex.tsv"),
                   col_names = c("Sample", "Phenotype")) %>%   mutate(Phenotype = ifelse(Phenotype == "ParentF", "Hermaphrodite",
                                                                                         ifelse(Phenotype == "G1", "Hermaphrodite",
                                                                                                ifelse(Phenotype == "ParentM", "Male", "Male"))))

hap1_radseq_coverage <- dplyr::left_join(hap1_radseq_coverage,phenotype,by="Sample")
hap1_radseq_coverage <- hap1_radseq_coverage %>%  filter(!is.na(Phenotype))

nrow(hap1_radseq_coverage)
# Remove rows with any NA values
hap1_radseq_coverage <- hap1_radseq_coverage[complete.cases(hap1_radseq_coverage),]
nrow(hap1_radseq_coverage)


hap1_radseq_coverage_means <- hap1_radseq_coverage %>%  
  group_by(Chrom, Start, Phenotype) %>%
  summarize(mean_cov = mean(as.numeric(FracCov))) %>% 
  filter(!is.na(mean_cov))
rm(hap1_radseq_coverage)


sorted_chroms <- mixedsort(unique(hap1_radseq_coverage_means$Chrom), 
                           decreasing = TRUE)

all_chrom_h1 <- hap1_radseq_coverage_means %>%   
  mutate(Chrom = factor(Chrom, levels = sorted_chroms)) %>%  
  ggplot(aes(x=as.numeric(Start)/10^6, y=mean_cov, color=Phenotype, group=Phenotype)) + 
  geom_smooth() + 
  scale_color_manual(values = color_scheme <- c(Hermaphrodite="#d8b365", Male="#5ab4ac")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  scale_y_continuous(transform ="log10") + 
  facet_wrap(~Chrom) + 
  xlab("Position in Mb") +
  ylab("Log10 of mean coverage") + 
  ggthemes::theme_few(base_size = 24) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save plot 
ggsave(paste0(path,"Cosegregation_h1_allchrom.pdf"), 
       all_chrom_h1, 
       width = 24,      
       height = 16,     
       units = "in",
       limitsize = FALSE)
############# ZOOM 12 #####################
chrom12_zoom_hap1 <- hap1_radseq_coverage_means %>%   
  mutate(Chrom = factor(Chrom, levels = sorted_chroms)) %>% 
  filter(Chrom == "Pan_040.h1-12") %>% 
  ggplot(aes(x=as.numeric(Start)/10^6, y=mean_cov, color=Phenotype, group=Phenotype)) + 
  geom_smooth(linewidth = 2) +  
  scale_color_manual(values = color_scheme <- c(Hermaphrodite="#d8b365", Male="#5ab4ac")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  scale_y_continuous(transform ="log10") + 
  facet_wrap(~Chrom) + 
  xlab("Position in Mb") +
  ylab("Log10 of mean coverage \n per phenotype") + 
  ggthemes::theme_few(base_size = 20) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  )

# Save plot
ggsave(paste0(path,"Cosegregation_h1_12.pdf"), 
       chrom12_zoom_hap1, 
       width = 20,      
       height = 12,     
       units = "in",
       limitsize = FALSE)

####################### HAP2

# hap 2 coverage
hap2_radseq_coverage <- vroom(paste0(remote,local,"ALL.Pan_040.hap2-rt.bamcov.bed.gz"), col_names = T)


# Cleanup Sample colunn
hap2_radseq_coverage <- hap2_radseq_coverage %>%
  mutate(Sample = str_remove(Sample, "_Pan_040.hap2-rt.fasta_radtag_paired"))

hap2_radseq_coverage <- dplyr::left_join(hap2_radseq_coverage,phenotype,by="Sample")
hap2_radseq_coverage <- hap2_radseq_coverage %>%  filter(!is.na(Phenotype))

hap2_radseq_coverage_means <- hap2_radseq_coverage %>%  
  drop_na() %>% 
  group_by(Chrom, Start, Phenotype) %>%
  summarize(mean_cov = mean(as.numeric(FracCov)))
rm(hap2_radseq_coverage)

sorted_chroms2 <- mixedsort(unique(hap2_radseq_coverage_means$Chrom), 
                            decreasing = TRUE)

############# ZOOOM
all_chrom_h2 <- hap2_radseq_coverage_means %>%   
  mutate(Chrom = factor(Chrom, levels = sorted_chroms2)) %>%  
  ggplot(aes(x=Start/10^6, y=mean_cov, color=Phenotype, group=Phenotype)) + 
  geom_smooth() + 
  scale_color_manual(values = color_scheme <- c(Hermaphrodite="#d8b365", Male="#5ab4ac")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  scale_y_continuous(transform ="log10") + 
  facet_wrap(~Chrom) + 
  xlab("Position in Mb") +
  ylab("Log10 of mean coverage per phenotype") + 
  ggthemes::theme_few(base_size = 20) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Save plot 
ggsave(paste0(path,"Cosegregation_h2_allchrom.pdf"), 
       all_chrom_h2, 
       width = 24,      
       height = 16,     
       units = "in",
       limitsize = FALSE)

chrom12_zoom_hap2 <- hap2_radseq_coverage_means %>% 
  mutate(Chrom = factor(Chrom, levels = sorted_chroms2)) %>% 
  filter(Chrom == "Pan_040.h2-12") %>% 
  ggplot(aes(x=Start/10^6, y=mean_cov, color=Phenotype, group=Phenotype)) + 
  geom_smooth(linewidth = 2) + 
  scale_color_manual(values = color_scheme <- c(Hermaphrodite="#d8b365", Male="#5ab4ac")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  scale_y_continuous(transform ="log10") + 
  facet_wrap(~Chrom) + 
  xlab("Position in Mb") +
  ylab("Log10 of mean coverage \n per phenotype") +  
  ggthemes::theme_few(base_size = 20) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )
# Save plot
ggsave(paste0(path,"Cosegregation_h2_12.pdf"), 
       chrom12_zoom_hap2, 
       width = 24,      
       height = 16,     
       units = "in",
       limitsize = FALSE)

# Retain bins belonging to the SDR:
# SDR coordinates
# Pan_040.h1-12   7566688 14277721
# Pan_040.h2-12   7686336 12870949

hap1_radseq_coverage_means  <- hap1_radseq_coverage_means %>%  
  mutate(Start = as.numeric(Start)) %>% 
  mutate(region = ifelse(Chrom == "Pan_040.h1-12" & Start >= 7566688 & Start <= 14277721, 
                         "SDR", 
                         "Background"))

hap1_radseq_coverage_means <- hap1_radseq_coverage_means %>%
  filter(!is.na(Phenotype), !is.na(region)) %>%
  mutate(
    Phenotype = factor(Phenotype),
    region = factor(region)
  )

# T test for Background
background_test <- hap1_radseq_coverage_means %>% filter(region == "Background") %>% 
  t.test(mean_cov ~ Phenotype, data = .)
# T test for SDR
sdr_test <- hap1_radseq_coverage_means %>% filter(region == "SDR") %>% 
  t.test(mean_cov ~ Phenotype, data = .)

# View results
background_test
sdr_test

# Plot test
a <- hap1_radseq_coverage_means %>%   filter(mean_cov < 0.03) %>% 
  ggplot(aes(x = Phenotype, y = mean_cov, fill = Phenotype)) +
  geom_boxplot(width = 0.5) + 
  ggtitle("040A11-H1") +
  ggthemes::theme_few(base_size = 20) + 
  theme(strip.background = element_blank(), 
        legend.position = "none",
        plot.title =  element_text(size = 16)) +  
  scale_fill_manual(values = c(Hermaphrodite = "#d8b365", Male = "#5ab4ac")) + 
  labs(y = "Mean Coverage",
       x = "") +
  facet_wrap(~region) +
  stat_compare_means(data = hap1_radseq_coverage_means,  comparisons = list(c("Hermaphrodite", "Male")),
                     method = "t.test",
                     label = "p.signif", label.y = 0.02)  # Changed from "p.format"



# H2
# Retain bins belonging to the SDR:
# SDR coordinates
# Pan_040.h2-12   7686336 12870949

hap2_radseq_coverage_means  <- hap2_radseq_coverage_means %>%  
  mutate(Start = as.numeric(Start)) %>% 
  mutate(region = ifelse(Chrom == "Pan_040.h2-12" & Start >= 7686336 & Start <= 12870949, 
                         "SDR", 
                         "Background"))

hap2_radseq_coverage_means <- hap2_radseq_coverage_means %>%
  filter(!is.na(Phenotype), !is.na(region)) %>%
  mutate(
    Phenotype = factor(Phenotype),
    region = factor(region)
  )

# T test for Background
background_test <- hap2_radseq_coverage_means %>% filter(region == "Background") %>% 
  t.test(mean_cov ~ Phenotype, data = .)
# T test for SDR
sdr_test <- hap2_radseq_coverage_means %>% filter(region == "SDR") %>% 
  t.test(mean_cov ~ Phenotype, data = .)

# View results
background_test
sdr_test

# Plot test
b <- hap2_radseq_coverage_means %>%   filter(mean_cov < 0.03) %>% 
  ggplot(aes(x = Phenotype, y = mean_cov, fill = Phenotype)) +
  geom_boxplot(width = 0.5) + 
  ggthemes::theme_few(base_size = 20) + 
  ggtitle("040A11-H2") +
  theme(strip.background = element_blank(), 
        legend.position = "none", 
        plot.title =  element_text(size = 16)) +  
  scale_fill_manual(values = c(Hermaphrodite = "#d8b365", Male = "#5ab4ac")) + 
  labs(y = "Mean Coverage",
       x = "") +
  facet_wrap(~region) +
  stat_compare_means(data = hap2_radseq_coverage_means,  comparisons = list(c("Hermaphrodite", "Male")),
                     method = "t.test",
                     label = "p.signif", label.y = 0.02)  # Changed from "p.format"


a | b
