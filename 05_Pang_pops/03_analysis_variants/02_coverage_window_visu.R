library(vroom)
library(tidyverse)
library(patchwork)
library(gtools)

# Remote path
remote="PATH/TO/FOLDER"
local="PATH/TO/FOLDER/"

# Sample metadata
metadata <-  vroom(paste0(remote,local,"metadata_samples.tsv"))
# Get big file
WINCOV <- vroom(paste0(remote,local,"ALL_BOTH_10000.regions.bed.gz"), 
                col_names = c("sample", "chrom", "start", "end","coverage"))

WINCOV <-left_join(WINCOV, metadata, by=c("sample"="SAMPLE_ID"))

# Remove regions with aberrant coverage  across all samples:
aberrant_windows <- WINCOV  %>%  filter(coverage > 50) %>%  select(chrom,start,end)
WINCOV_filtered <- WINCOV  %>%
  anti_join(aberrant_windows, by = c("chrom", "start", "end"))

# Remove Camargue individuals
WINCOV_filtered <-  WINCOV_filtered %>%  filter(Location != "Camargue")

# Normalize with median coverage
median_cov_sample <- WINCOV_filtered %>%  
  group_by(sample) %>% 
  summarize(median_cov = median(coverage)) 

WINCOV_filtered <- left_join(WINCOV_filtered,median_cov_sample, by="sample")
WINCOV_filtered <-  WINCOV_filtered %>% 
  mutate(norm_coverage=coverage/median_cov)

write_delim(x = WINCOV_filtered,
            file = paste0(remote,local,"ALL_windowed_10000_coverage_normalized.tsv"), 
            delim = "\t", 
            quote = "none", 
            col_names = TRUE)
#####################################

# Get averages across groups:
WINCOV_filtered_averaged <- WINCOV_filtered %>%
  group_by(chrom, start, end, Location, Sex_phenotype) %>%
  summarize(group_mean_cov = mean(norm_coverage))

# Sort chromosomes
sorted_chroms <- mixedsort(unique(WINCOV_medians_sex_phenotype$chrom), 
                           decreasing = TRUE)

# SDR Rearrangend loci cordinates: 
# Pan_040.h1-12	7470730	14307655
# Pan_040.h2-12	7588215	12900934
WINCOV_filtered_averaged  %>%   
  filter(grepl("Pan_040.h1-12",chrom)) %>% 
  mutate(chrom = factor(chrom, levels = sorted_chroms)) %>%  
  ggplot(aes(x=as.numeric(start)/10^6, y=group_mean_cov, color=Sex_phenotype, group=Sex_phenotype)) + 
  geom_point() + 
  geom_vline(xintercept = 7470730/10^6, linetype="dotted", color="black", size=1.5) + 
  geom_vline(xintercept = 14307655/10^6, linetype="dotted", color="black", size=1.5) + 
  scale_color_manual(values = color_scheme <- c(H="#d8b365", M="#5ab4ac")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  facet_wrap(~chrom+Location) + 
  xlab("Position in Mb") +
  ylab("Group mean coverage") + 
  theme_bw(base_size = 24) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
# SDR Rearrangend loci cordinates: 
# Pan_040.h1-12	7470730	14307655
# Pan_040.h2-12	7588215	12900934
WINCOV_filtered_averaged  %>%   
  filter(grepl("Pan_040.h2-12",chrom)) %>% 
  mutate(chrom = factor(chrom, levels = sorted_chroms)) %>%  
  ggplot(aes(x=as.numeric(start)/10^6, y=group_mean_cov, color=Sex_phenotype, group=Sex_phenotype)) + 
  geom_point() + 
  scale_color_manual(values = color_scheme <- c(H="#d8b365", M="#5ab4ac")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  facet_wrap(~chrom+Location) + 
  geom_vline(xintercept = 7588215/10^6, linetype="dotted", color="black", size=1.5) + 
  geom_vline(xintercept = 12900934/10^6, linetype="dotted", color="black", size=1.5) + 
  xlab("Position in Mb") +
  ylab("Group mean coverage") + 
  theme_bw(base_size = 24) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

WINCOV_filtered_averaged_fabr <- WINCOV_filtered_averaged %>%  filter(Location == "Fabregues")


#############################################
WINCOV_filtered_averaged_fabr <- WINCOV_filtered_averaged_fabr  %>%
  group_by(chrom,start,end) %>% 
  # Pivot wider to get H and M as separate columns
  pivot_wider(
    names_from = Sex_phenotype,
    values_from = group_mean_cov,
    names_prefix = "Sex_"
  ) %>%
  # Calculate the difference (H - M)
  mutate(
    H_minus_M = Sex_H - Sex_M
  )

a <- WINCOV_filtered_averaged_fabr %>%  filter(grepl("Pan_040.h1-12",chrom)) %>% 
  mutate(chrom = factor(chrom, levels = sorted_chroms)) %>%  
  ggplot(aes(x=as.numeric(start)/10^6, y=H_minus_M)) + 
  geom_point(color="grey") +
  geom_vline(xintercept = 7470730/10^6, linetype="dotted", color="black", size=1.5) + 
  geom_vline(xintercept = 14307655/10^6, linetype="dotted", color="black", size=1.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  facet_wrap(~chrom+Location) + 
  xlab("Position in Mb") +
  ylab("Herm-Male") + 
  theme_bw(base_size = 24) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

b <- WINCOV_filtered_averaged_fabr %>%  filter(grepl("Pan_040.h2-12",chrom)) %>% 
  mutate(chrom = factor(chrom, levels = sorted_chroms)) %>%  
  ggplot(aes(x=as.numeric(start)/10^6, y=H_minus_M)) + 
  geom_point(color="grey") +
  geom_vline(xintercept = 7588215/10^6, linetype="dotted", color="black", size=1.5) + 
  geom_vline(xintercept = 12900934/10^6, linetype="dotted", color="black", size=1.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  facet_wrap(~chrom+Location) + 
  xlab("Position in Mb") +
  ylab("Herm-Male") +
  theme_bw(base_size = 24) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )


#############################################
WINCOV_filtered_averaged_cadix <- WINCOV_filtered_averaged %>%  filter(Location == "Cadix")

WINCOV_filtered_averaged_cadix  <- WINCOV_filtered_averaged_cadix  %>%
  group_by(chrom,start,end) %>% 
  # Pivot wider to get H and M as separate columns
  pivot_wider(
    names_from = Sex_phenotype,
    values_from = group_mean_cov,
    names_prefix = "Sex_"
  ) %>%
  # Calculate the difference (H - M)
  mutate(
    H_minus_M = Sex_H - Sex_M
  )

c <- WINCOV_filtered_averaged_cadix %>%  filter(grepl("Pan_040.h1-12",chrom)) %>% 
  mutate(chrom = factor(chrom, levels = sorted_chroms)) %>%  
  ggplot(aes(x=as.numeric(start)/10^6, y=H_minus_M)) + 
  geom_point(color="grey") +
  geom_vline(xintercept = 7470730/10^6, linetype="dotted", color="black", size=1.5) + 
  geom_vline(xintercept = 14307655/10^6, linetype="dotted", color="black", size=1.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  facet_wrap(~chrom+Location) + 
  xlab("Position in Mb") +
  ylab("Herm-Male") + 
  theme_bw(base_size = 24) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

d <- WINCOV_filtered_averaged_cadix %>%  filter(grepl("Pan_040.h2-12",chrom)) %>% 
  mutate(chrom = factor(chrom, levels = sorted_chroms)) %>%  
  ggplot(aes(x=as.numeric(start)/10^6, y=H_minus_M)) + 
  geom_point(color="grey") +
  geom_vline(xintercept = 7588215/10^6, linetype="dotted", color="black", size=1.5) + 
  geom_vline(xintercept = 12900934/10^6, linetype="dotted", color="black", size=1.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n=12) ) +
  facet_wrap(~chrom+Location) + 
  xlab("Position in Mb") +
  ylab("Herm-Male") +
  theme_bw(base_size = 24) +
  theme(
    strip.background = element_rect(fill = "white", color="white"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

c | d 



##############################################################

WINCOV_filtered_averaged_cadix %>% 
  filter(Sex_H == 0) %>%  
  filter(Sex_M > 0.2) %>% 
  ggplot(aes(x=start/10^6,y=H_minus_M)) + 
  geom_point() + 
  facet_wrap(~chrom)

write_delim(x = WINCOV_filtered_averaged_cadix,
            file = paste0(remote,local,"WINCOV_filtered_averaged_cadix.tsv"), 
            delim = "\t", 
            quote = "none", 
            col_names = TRUE)



# Cadix FAB differences

# Get averages across groups:
WINCOV_filtered_popdiff<- WINCOV_filtered %>%
  group_by(chrom, start, end, Location) %>%
  summarize(group_mean_cov = mean(norm_coverage))
