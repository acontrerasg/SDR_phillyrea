#!/usr/bin/env Rscript
#Input
path_synt="PATH/TO/FOLDER/"
chr_12_M <- vroom(paste0(path_synt,"Pan_040.h1-rt_ALL_SDR.gff3"),
                  comment = "#", col_names = FALSE)
chr_12_m <- vroom(paste0(path_synt,"Pan_040.h2-rt_ALL_SDR.gff3"),
                  comment = "#", col_names = FALSE)

# Add a y-axis level to separate chromosomes
chr_12_M <- chr_12_M %>% mutate(ymin = 0.5, ymax = 0.7)  # Chromosome 1
chr_12_m <- chr_12_m %>% mutate(ymin = 0.1, ymax = 0.3)

# Combine both chromosomes
chromosomes <- bind_rows(chr_12_M, chr_12_m)
# Colors:
colors_synblk <- c(
  "ID=SynBlk_1" ="#8B8680", # Gray
  "ID=SynBlk_2" = "#B22222",  # Firebrick
  "ID=SynBlk_3" = "#A9B2C3",  # MediumSeaGreen
  "ID=SynBlk_4" = "#F8D568",  # MediumSeaGreen
  "ID=SynBlk_5" = "#9999CC",  
  "ID=SynBlk_6" = "#00CC99",  
  "ID=SynBlk_7" ="#8B8680", # Gray
  "ID=NonSynBlk_M_1" = "#5ab4ac",  # M
  "ID=NonSynBlk_M_2" = "#5ab4ac",  # M
  "ID=NonSynBlk_M_3" = "#5ab4ac",  # M
  "ID=NonSynBlk_M_4" = "#5ab4ac",  # M
  "ID=NonSynBlk_M_5" = "#5ab4ac",  # M
  "ID=NonSynBlk_M_6" = "#5ab4ac",  # M
  "ID=NonSynBlk_m_1" = "#d8b365",  # m
  "ID=NonSynBlk_m_2" = "#d8b365",  # m
  "ID=NonSynBlk_m_3" = "#d8b365",  # m
  "ID=NonSynBlk_m_4" = "#d8b365",  # m
  "ID=NonSynBlk_m_5" = "#d8b365",  # m
  "ID=NonSynBlk_m_6" = "#d8b365"  # m
)

colors_M <- c(
  "ID=SynBlk_1" ="#8B8680", # Gray
  "ID=NonSynBlk_M_1" = "#5ab4ac",  # M
  "ID=SynBlk_2" = "#9999CC",  # Firebrick
  "ID=NonSynBlk_M_2" = "#5ab4ac",  # M
  "ID=SynBlk_3" ="#2EB4E6",  # MediumSeaGreen
  "ID=NonSynBlk_M_3" = "#5ab4ac",  # M
  "ID=SynBlk_4" = "#9999CC",  # MediumSeaGreen
  "ID=NonSynBlk_M_4" = "#5ab4ac",  # M
  "ID=SynBlk_5" = "#2EB4E6",  # Firebrick
  "ID=NonSynBlk_M_5" = "#5ab4ac",  # M
  "ID=SynBlk_6" = "#9999CC",  # Firebrick
  "ID=NonSynBlk_M_6" = "#5ab4ac",  # M
  "ID=SynBlk_7" ="#8B8680" # Gray
)
# Negativo "#9999CC" 
# Positivo "#2EB4E6"
colors_m <- c(
  "ID=SynBlk_1" ="#8B8680", # Gray
  "ID=NonSynBlk_m_1" = "#d8b365",  # m
  "ID=SynBlk_5" = "#2EB4E6",  # MediumSeaGreen
  "ID=NonSynBlk_m_2" = "#d8b365",  # m
  "ID=SynBlk_4" = "#9999CC",  # MediumSeaGreen
  "ID=NonSynBlk_m_3" = "#d8b365",  # m
  "ID=SynBlk_3" = "#2EB4E6",  # Firebrick
  "ID=NonSynBlk_m_4" = "#d8b365",  # m
  "ID=SynBlk_2" = "#9999CC",  # Firebrick
  "ID=NonSynBlk_m_5" = "#d8b365",  # m
  "ID=SynBlk_6" = "#9999CC",  # Firebrick
  "ID=NonSynBlk_m_6" = "#d8b365",  # m
  "ID=SynBlk_7" ="#8B8680" # Gray
)
# Negativo "#9999CC" 
# Positivo "#2EB4E6"

# Plot
ggplot(chromosomes, aes(xmin = X4/10^6, xmax = X5/10^6, ymin = ymin, ymax = ymax, fill = X9)) +
  geom_rect(color = "black") +
  scale_fill_manual(values = colors_synblk) + 
  theme_minimal() +
  labs(title = "Chromosome Segments", x = "Genomic Position", y = "Chromosome") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


# Define plot limits
plot(c(0, 40), c(-1, 3), type = "n", xlab = "Genomic Position", ylab = "", axes = FALSE, main = "Chromosomes with Segments")
# Draw Chromosomes as Background Bars
rect(0, 1.5,37903724/10^6, 2.5, col = "gray90", border = FALSE)  # Chromosome 1 (Upper)
rect(0, -0.5, 36710837/10^6, 0.5, col = "gray90", border =FALSE) # Chromosome 2 (Lower)

# Add minor ticks at every 1 Mb interval
axis(1, at = seq(0, 38, by = 1), las = 1, cex.axis = 0.001, labels = FALSE, tck = -0.02)
# Add major ticks with labels every 5 Mb interval
axis(1, at = seq(0, 38, by = 5), las = 1, cex.axis = 0.7, tck = -0.04)
# Add Y-axis (optional, to show chromosome positions)
#axis(2, at = c(-0.5, 2.5), labels = c("Male Hap 2", "Male Hap 1"), las = 2)

# Function to Draw Chromosome Blocks
draw_blocks <- function(data, y_min, y_max, colors) {
  for (i in 1:nrow(data)) {
    rect(data$X4[i]/10^6, y_min, data$X5[i]/10^6, y_max, col = colors[i], border = "grey40")
  }
}
draw_blocks(chr_12_M, 1.5, 2.5, colors_M)
draw_blocks(chr_12_m, -0.5, 0.5, colors_m)


# Plot the polygons representing the alignments.
chr_12_M_syn <- chr_12_M %>%  filter(X3 != "NonSynBlk") %>% arrange(X9)
chr_12_m_syn <- chr_12_m %>%  filter(X3 != "NonSynBlk") %>% arrange(X9)

for (x in 1:7) {
  if (chr_12_M_syn[x, 7] == "+") {
    x_positions <- c(chr_12_m_syn[x, 5] / 10^6, chr_12_m_syn[x, 4] / 10^6, 
                     chr_12_M_syn[x, 4] / 10^6, chr_12_M_syn[x, 5] / 10^6)
  } else {
    x_positions <- c(chr_12_m_syn[x, 5] / 10^6, chr_12_m_syn[x, 4] / 10^6, 
                     chr_12_M_syn[x, 5] / 10^6, chr_12_M_syn[x, 4] / 10^6)
  }
  
  polygon(x = x_positions, 
          y = c(0.5, 0.5, 1.5, 1.5),
          col = adjustcolor('black', alpha.f = 0.1), 
          border = FALSE)
}