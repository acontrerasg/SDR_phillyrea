#!/usr/bin/env Rscript 
library(vroom)
library(tidyverse)
library(gtools)

# Functions
# Aggregate bed wins
aggregate_bed_windows <- function(bed, win_size = 100000) {
  library(dplyr)
  
  new_bed <- bed %>%
    mutate(window = floor(start / win_size)) %>%
    group_by(chr, window) %>%
    summarise(
      start = min(start),
      end = max(end),
      kmer_count = sum(kmer_count),
      .groups = "drop"
    )
  
  return(new_bed)
}

# Compute cumulative positions
compute_cumulative_positions <- function(bed) {
  library(dplyr)
  
  bed_cum <- bed %>%
    # Compute chromosome sizes
    group_by(chr) %>%
    summarise(chr_len = max(start), .groups = "drop") %>%
    # Compute cumulative offsets
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(chr, tot) %>%
    # Join back to original data
    left_join(bed, by = "chr") %>%
    # Order and compute cumulative position
    arrange(chr, start) %>%
    mutate(startcum = start + tot)
  
  return(bed_cum)
}

# Plot  Manhattan  of Kmer counts 
plot_manhattan_kmers <- function(don,
                                 title = "K-mer Manhattan plot",
                                 colors = c("black", "grey40"),
                                 shade_colors = c("grey80", "white"),
                                 threshold = 500,
                                 highlight_col = "#EB3742")
{
  # Chromosome axis centers + boundaries
  axisdf <- don %>%
    group_by(chr) %>%
    summarize(
      center = (max(startcum) + min(startcum)) / 2,
      xmin = min(startcum),
      xmax = max(startcum),
      .groups = "drop"
    ) %>%
    arrange(center)
  
  # alternating shading colors (user-defined)
  axisdf$fill <- rep(shade_colors, length.out = nrow(axisdf))
  
  # Chromosome colors for points
  chr_levels <- axisdf$chr
  n_chr <- length(chr_levels)
  
  col_vec <- rep(colors, length.out = n_chr)
  names(col_vec) <- chr_levels
  
  # Split data
  don_low  <- subset(don, kmer_count <= threshold)
  don_high <- subset(don, kmer_count > threshold)
  
  # Plot
  p <- ggplot() +
    
    # Background chromosome shading
    geom_rect(
      data = axisdf,
      aes(xmin = xmin, xmax = xmax,
          ymin = -Inf, ymax = Inf,
          fill = fill),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    scale_fill_identity() +
    
    # Base points
    geom_point(
      data = don_low,
      aes(x = startcum, y = kmer_count, color = as.factor(chr)),
      alpha = 1, size = 0.8
    ) +
    
    # Highlight points
    geom_point(
      data = don_high,
      aes(x = startcum, y = kmer_count),
      color = highlight_col, size = 1
    ) +
    
    scale_color_manual(values = col_vec) +
    
    scale_x_continuous(
      breaks = axisdf$center,
      labels = axisdf$chr
    ) +
    
    scale_y_continuous(expand = c(0, 0)) +
    
    ggthemes::theme_few(base_size = 20) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text = element_text(colour = "black")
    ) +
    
    ylab("K-mer count") +
    xlab("Chromosome") +
    ggtitle(title)
  
  return(p)
}


# MAIN 
#bed_path="PATH/TO/FILE"
# MAIN 
remote="PATH/TO/FOLDER"
path="PATH/TO/FOLDER/"
bed_path=paste0(remote,path,"Male_unique_hap1.sorted_intersection_counts.bed")
bed <- vroom(bed_path, col_names = c("chr", "start", "end", "kmer_count"))
bed <-  bed %>% mutate(chr = as.numeric(str_remove(string = chr, pattern = "Pan_040.h1-"))) %>% 
  mutate(chr = factor(chr, levels = 1:23, ordered = TRUE)) 

win_size=100000
bed_agg <- aggregate_bed_windows(bed, win_size = win_size)
com_pos_bed <- compute_cumulative_positions(bed_agg)
plot_manhattan_kmers(com_pos_bed)

col = c("grey60", "grey90")

plot_manhattan_kmers(com_pos_bed, colors = col, highlight_col = "#CC474B")

# Zoom in in chrtomosome 12
bed_agg %>%
  filter(chr == 12) %>%
  mutate(start_mbp = start / 1e6) %>%
  ggplot(aes(x = start_mbp, y = kmer_count)) +
  geom_line(linewidth = 1, color = "#5fc78eff") +
  scale_x_continuous(
    breaks = seq(
      floor(min(bed_agg$start[bed_agg$chr == 12]) / 1e6),
      ceiling(max(bed_agg$start[bed_agg$chr == 12]) / 1e6),
      by = 5)
  ) +
  xlab("Mbp") +
  ylab("K-mer count") +
  ggthemes::theme_few(base_size = 20) +
  theme(
    # Text
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    # Axis lines & ticks
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    # Panel border (optional but nice)
    panel.border = element_rect(color = "black", fill = NA),
  )

