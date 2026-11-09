#!/usr/bin/env Rscript
library(tidyverse)
library(vroom)
library(patchwork)

# remote path
remote="PATH/TO/FOLDER"
# path
path="PATH/TO/FOLDER"
# Genomic block size:
genomic_size <- data.frame(
  Block = c("CA5p", "M1", "TIB1", "M2", "B2", "M3", "TIB3", "M4", "TB4", "M5", "I5", "M6", "CA3p"),
  Size_Mbp = c(7.47073, 0.038302, 0.906352, 1.43419, 0.84567, 0.143925, 0.069353, 0.081672, 0.881481, 
               0.123004, 1.02496, 1.28801, 23.5961)
)
genomic_size %>%
  filter(!Block %in% c("M6","M1", "M2", "M3", "M4", "M5", "M6", "Background","CA3p","CA5p")) %>% 
  mutate(Block = factor(Block, levels = c("TIB1", "B2", "TIB3", "TB4", "I5"))) %>% 
  ggplot(aes(x = Block, y = Size_Mbp)) +
  geom_col() +
  theme_bw() +
  labs(x = "", y = "Size in Mbp") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

# KAKS
kaks_GY <- vroom(file= paste0(remote,path,"KaKs_haplotype/KaKs_final.ALL.tsv"))

# Overview:
# Calculate median before filtering
median_ks <- median(kaks_GY$Ks, na.rm = TRUE)
median_MYA <- median(kaks_GY$Mya, na.rm = TRUE)
# Now create the plot with the pre-calculated median
kaks_GY %>%
  filter(Ks < 0.2) %>%
  ggplot(aes(x=Ks)) +
  geom_histogram(bins = 50, color="grey20", fill="dodgerblue") +
  geom_vline(xintercept = median_ks, 
             color = "grey20", 
             linetype = "dashed", 
             linewidth = 1) +
  annotate("text", 
           x = median_ks, 
           y = Inf, 
           label = paste("Median =", round(median_ks, 3)),
           vjust = 10, 
           hjust = -0.2,
           color = "grey20", 
           size = 4) +
  theme_bw(base_size = 12)
#######

# Prepare Kaks fle for merging
kaks_GY <- kaks_GY %>%
  separate(Sequence, into = c("hap1", "hap2"), 
           sep = "&",
           remove = FALSE) %>% 
  mutate(hap1 == str_remove(string = hap1, pattern = "_.$"))

# Remove the _1 pattern
kaks_GY <- kaks_GY %>% 
  mutate(hap1 = str_remove(string = hap1, pattern = "_\\d+$"),
         hap2 = str_remove(string = hap2, pattern = "_\\d+$"))



# Join it to the  SDR annotation:
remote="PATH/TO/FOLDER"
path="PATH/TO/FOLDER"
SDR_anno_hap1 <- vroom(file= paste0(remote,path,"Pan_040.hap1-rt_SDR_genes.gff3"), col_names = F)
# Calculate median before filtering
SDR_anno_hap1 <- SDR_anno_hap1 %>%
  filter(X3== 'gene') %>%
  mutate(hap1 = str_remove(string = X9, pattern = ";.*"))  %>% 
  mutate(hap1 =str_remove(string = hap1, pattern = "ID=")) %>%
  mutate(Block= str_remove(string = X9, pattern = ".*;ID="))

# Join it to the Kaks dataset
kaks_GY_annotated <- left_join(kaks_GY, SDR_anno_hap1, by = "hap1")
kaks_GY_annotated <- kaks_GY_annotated  %>% 
  mutate(Block = if_else(is.na(Block), "Background", Block))

# Add per block syntenic information
remote="PATH/TO/FOLDER"
path="PATH/TO/FOLDER/"
# Load Syntenic lassification file from disk 
Syntenic_classification_table <- vroom(paste0(remote,path,"syntenic_classification_SDR.tsv")) %>% 
  select(hap1,hap2, Synteny, Identity)
# Join it to the Kaks dataset
kaks_GY_annotated <- left_join(kaks_GY_annotated, Syntenic_classification_table , by = c("hap1","hap2"))
# Filter out NA values
kaks_GY_annotated <- kaks_GY_annotated %>% filter(!is.na(Ks))


# Histogram with counts
p1 <- plot_data %>%
  ggplot(aes(x = Block, fill = Block)) +
  geom_bar(width = 0.6, color="grey20") +
  ylim(c(0,1200)) + 
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("#808080", "#9999cc", "#2eb4e6", "#9999cc","#2eb4e6", "#9999cc", "#808080")) +
  labs(y = "Allele\nNumber") +
  ylim(c(0,1500)) + 
  theme_bw(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# Boxplot
p2 <- plot_data %>%
  filter(Ks < 0.05) %>% 
  ggplot(aes(x = Block, y = Ks, fill = Block)) +
  geom_boxplot(width = 0.2, outliers = F) + 
  scale_fill_manual(values = c("#808080", "#9999cc", "#2eb4e6", "#9999cc","#2eb4e6", "#9999cc", "#808080")) +
  labs(y = "Ds") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  geom_hline(yintercept = median_ks, 
             color = "grey20", 
             linetype = "dashed", 
             linewidth = 1) +
  annotate("text", 
           x = Inf, 
           y = median_ks, 
           label = paste("Genome median =", round(median_ks, 4)),
           vjust = 2, 
           hjust = 1.7,
           color = "grey20", 
           size = 4)

# Combine plots
p1 / p2 + plot_layout(heights = c(1, 3))


kaks_GY_annotated  <- kaks_GY_annotated %>%
  mutate(type = if_else(Block %in% c("Background", "CA3p", "CA5p"), 
                        "genomic", 
                        "SDR"))

# Calcualte overall Ds

kaks_GY_annotated %>% filter( Ks < 0.2 ) %>% 
  ggplot(aes(Ks, fct_rev(type), fill = type)) +
  scale_fill_manual(values = c("grey60", "darkred")) + 
  ggridges::geom_density_ridges(alpha = 0.7, color="grey20") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") + 
  labs(y="", x="Ds")

medians <- kaks_GY_annotated %>%
  filter(Ks < 0.2) %>%
  group_by(type) %>%
  summarise(median_Ks = median(Ks, na.rm = TRUE))

kaks_GY_annotated %>%
  filter(Ks < 0.2) %>%
  ggplot(aes(x = Ks, fill = type)) +
  geom_histogram(alpha = 0.7, color = "grey20", bins = 30) +
  geom_vline(
    data = medians,
    aes(xintercept = median_Ks),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_text(
    data = medians,
    aes(
      x = median_Ks,
      y = Inf,
      label = sprintf("%.3f", median_Ks)
    ),
    vjust = 1.3,
    hjust=-0.2,
    angle = 0,
    size = 4,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("grey60", "darkred")) +
  facet_wrap(~ fct_rev(type), ncol = 1, scales = "free_y") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(colour = "black")
  ) +
  labs(
    x = "Ds",
    y = ""
  )

kaks_GY_annotated %>% filter( Ks < 0.2 ) %>% 
  ggplot(aes(x=Ks)) + geom_histogram(color="black", fill="dodgerblue", bins = 40) + 
  geom_vline(xintercept = median_ks, 
             color = "grey20", 
             linetype = "dashed", 
             linewidth = 1) +
  annotate("text", 
           x = Inf, 
           y = median_ks, 
           label = paste("Genome median =", round(median_ks, 4)),
           vjust = -15, 
           hjust = 2.3,
           color = "grey20", 
           size = 4) + 
  theme_bw(base_size = 16) +
  theme(legend.position = "none") + 
  labs(y="", x="Genome-wide Ds")


# Estimate diverngce time between alleles at the SDR:
syntenic_block <- c("TIB1", "B2", "TIB3", "TB4", "I5")
SDR_syntenic_genes <- kaks_GY_annotated %>%  filter(type == "SDR") %>%  filter(Block %in% syntenic_block)
mu_rate=7.77e-09
SDR_syntenic_genes %>%  mutate(time=Ks/(2*mu_rate)/10^6) %>%  
  group_by(Block) %>%  
  summarise(block_age=median(time), gene_pairs=n())
# Given t=d/2r
# D is ks and r is mutation rate
to_disk <- SDR_syntenic_genes %>%  
  mutate(time=Ks/(2*mu_rate)/10^6) %>%  
  select(Sequence, hap1, hap2, Method, Ka, Ks, `Ka/Ks` , Length, `S-Sites`, `N-Sites`, Substitutions, type, Block, X9) 

write.table(to_disk , file = "PATH/TO/FILE", 
            sep = "\t", quote = F, row.names = F,col.names = T)

## Test whether the Ks atthe SDR is significative higher: 
# Calculate statistics for each type

syntenic_block_background <- c("TIB1", "B2", "TIB3", "TB4", "I5", "Background")
stats_by_type <- kaks_GY_annotated %>%
  filter(Block %in% syntenic_block_background) %>% 
  group_by(type) %>%
  summarise(
    n = n(),
    median = median(Ks, na.rm = TRUE),
    mean = mean(Ks, na.rm = TRUE),
    sd = sd(Ks, na.rm = TRUE),
    min = min(Ks, na.rm = TRUE),
    max = max(Ks, na.rm = TRUE)
  )

print(stats_by_type)

# Extractng data
SDR_ks <- SDR_syntenic_genes %>%  pull(Ks)  # 55 Ks values
genomic_ks <- kaks_GY_annotated %>%  filter(type != "SDR") %>% filter(!is.na(Ks)) %>% pull(Ks)   # 33,139 Ks values
# Set parameters
n_small <- length(SDR_ks)
n_permutations <- 10000
set.seed(93565)

# Calculate observed statistics for SDR group
observed_median <- median(SDR_ks)
# Perform permutation test with replacement
permuted_medians <- numeric(n_permutations)
sampled_values_results <- vector("list", n_permutations) 
shapiro_results <- vector("list", n_permutations)
wilcox_results <- vector("list", n_permutations)

## Loop
for(i in 1:n_permutations) {
  # Sample with replacement from genomic group
  sampled_values <- sample(genomic_ks, size = n_small, replace = TRUE)
  # Store sampled values
  sampled_values_results[[i]] <- sampled_values
  # Calculate statistics
  permuted_medians[i] <- median(sampled_values)
  # Perform Shapiro-Wilk test for normality
  shapiro_results[[i]] <- shapiro.test(sampled_values)
  # Perform Mann-Whitney U test comparing sampled values to SDR_ks
  wilcox_results[[i]] <- wilcox.test(SDR_ks, sampled_values)
}

# Extract Wilcox p-values and plot them
wilcox_results_pvalues <- sapply(wilcox_results, function(x) x$p.value)

df_pvals <- data.frame(pval_wilcox = wilcox_results_pvalues)

df_pvals  %>%  ggplot(aes(x=pval_wilcox)) + 
  geom_histogram() + 
  geom_vline(xintercept = 0.05, 
             color = "red", 
             linetype = "dashed", 
             linewidth = 1) +
  theme_bw(base_size = 12)

# Plot distribution of permuted medians with observed median
tibble(permuted_median = permuted_medians) %>%
  ggplot(aes(x = permuted_median)) +
  geom_histogram(bins = 50, fill = "dodgerblue", color = "grey20") +
  geom_vline(xintercept = observed_median, color = "grey10", linetype = "dashed", linewidth = 1) +
  annotate("text", x = observed_median, y = Inf, 
           label = "Observed SDR\n median: 0.0283", 
           vjust = 2, hjust = 1.1, color = "grey10", size = 4) +
  labs(title = "Distribution of boostrapped genomic Ks medians",
       x = "Median Ks",
       y = "Count") +
  theme_bw(base_size = 12)