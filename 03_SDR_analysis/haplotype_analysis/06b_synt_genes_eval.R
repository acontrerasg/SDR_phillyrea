#!/usr/bin/env Rscript
library(tidyverse)
library(patchwork)
remote="PATH/TO/FOLDER"
path="PATH/TO/FOLDER"

# Merging function
get_rbh <- function(df_avb, df_bva) {
  # Column name constants matching DIAMOND verbose TSV headers
  COL_QUERY    <- "Query ID"
  COL_SUBJECT  <- "Subject ID"
  COL_PIDENT   <- "Percentage of identical matches"
  COL_BITSCORE <- "Bit score"
  COL_EVALUE   <- "Expected value"
  # ---- Slim down and rename for merging ------------------------------------
  avb_slim <- data.frame(
    protein_A     = df_avb[[COL_QUERY]],
    protein_B     = df_avb[[COL_SUBJECT]],
    pident_AvsB   = df_avb[[COL_PIDENT]],
    bitscore_AvsB = df_avb[[COL_BITSCORE]],
    evalue_AvsB   = df_avb[[COL_EVALUE]],
    stringsAsFactors = FALSE
  )
  
  bva_slim <- data.frame(
    protein_B     = df_bva[[COL_QUERY]],
    protein_A     = df_bva[[COL_SUBJECT]],
    pident_BvsA   = df_bva[[COL_PIDENT]],
    bitscore_BvsA = df_bva[[COL_BITSCORE]],
    evalue_BvsA   = df_bva[[COL_EVALUE]],
    stringsAsFactors = FALSE
  )
  
  # Inner join on both IDs — pair must agree in both directions ----------
  merge(avb_slim, bva_slim, by = c("protein_A", "protein_B")) %>%  as_tibble()
}
# Function for plotting SDR gene counts
# Function to plot dodged counts per SDR block with custom order
plot_dodged_counts <- function(data, sdr_col = "SDR_block_A", counting_col = "Identity", order_SDR_block = NULL) {
  
  # Convert column names to symbols for tidy evaluation
  sdr_sym <- sym(sdr_col)
  cnt_sym <- sym(counting_col)
  
  plot_data <- data %>%
    filter(!is.na(!!sdr_sym)) %>%
    group_by(!!sdr_sym, !!cnt_sym) %>%
    tally() %>%
    ungroup()
  
  # If order provided, convert SDR column to factor with that order
  if (!is.null(order_SDR_block)) {
    plot_data <- plot_data %>%
      mutate(!!sdr_sym := factor(!!sdr_sym, levels = order_SDR_block))
  }
  
  ggplot(plot_data, aes(x = !!sdr_sym, y = n, fill = !!cnt_sym)) +
    geom_col(position = position_dodge(width = 0.8), color="grey40") +
    scale_fill_grey() +
    geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.5) +
    labs(
      y = "Count",
      x = sdr_col,
      fill = counting_col
    ) +
    theme_minimal()
}
# PLot rearragmenets
plot_synteny_bump <- function(df, genome_labels = c(A = "Genome A", B = "Genome B")) {
  
  plot_df <- df %>%
    mutate(
      rank_A = rank(start_A, ties.method = "first", na.last = "keep"),
      rank_B = rank(start_B, ties.method = "first", na.last = "keep")
    ) %>% 
    pivot_longer(
      cols = c(rank_A, rank_B, SDR_block_A, SDR_block_B),
      names_to = c(".value", "genome"),
      names_pattern = "(rank|SDR_block)_(A|B)"
    )
  # create numeric x positions
  plot_df$genome_x <- as.numeric(factor(plot_df$genome))
  p <- plot_df %>%  
    ggplot(aes(x = genome, y = rank, group = protein_A)) +
    # Syntenic lines (dodge)
    geom_line(
      data = subset(plot_df,Identity == "Homolog"),
      aes(color = Synteny, group = protein_A, linetype = Synteny),
      alpha = 1,
      linewidth = 0.5
    ) +
    scale_color_manual(values = c(
      "Syntenic" = "grey20",
      "Rearranged" = "grey80",
      "Non-syntenic" = "white"
    )) +
    scale_linetype_manual(values = c(
      "Syntenic" = "solid",
      "Rearranged" = "dashed",
      "Non-syntenic" = "solid"
    )) + 
    ggnewscale::new_scale_color() +
    geom_point(
      aes(
        x = genome_x + ifelse(genome_x %% 2 == 0,
                              ifelse(Identity == "Homolog", 0.05, 0.1),
                              ifelse(Identity == "Homolog", -0.05,  -0.1)),
        color = SDR_block,
        shape = Identity,
        size = Identity
      )) +
    scale_shape_manual(values = c(
      "Homolog" = 95,
      "Private" = 95
    )) + 
    scale_size_manual(values = c(
      "Homolog" = 6,
      "Private" = 2.5
    )) +
    scale_x_discrete(labels = genome_labels, expand = expansion(mult = c(1, 1))) +
    scale_y_continuous(expand = expansion(mult = c(0.5, 0.5))) +
    theme_minimal() +
    guides(size = "none") +
    labs(
      x = "",
      y = "Relative gene position",
      color = "SDR_block",
      shape = "Identity"
    )
  
  return(p)
}

# Get diamond output
protein_h1vsh2 <- vroom::vroom(file = paste0(remote,path,"protein_h1vsh2.tsv"), col_names = T, delim = "\t")
# Get diamond output
protein_h2vsh1 <- vroom::vroom(file = paste0(remote,path,"protein_h2vsh1.tsv"), col_names = T, delim = "\t")

# Merged dmd results
protein_hits <- get_rbh(protein_h1vsh2, protein_h2vsh1) %>%  
  mutate(protein_A = str_remove(protein_A, "-mRNA-1")) %>%
  mutate(protein_B = str_remove(protein_B, "-mRNA-1"))


# Get bed files:
bed_cols_h1 <- c("chr_A", "start_A", "end_A", "protein_A", "SDR_block_A")
bed_cols_h2 <- c("chr_B", "start_B", "end_B", "protein_B", "SDR_block_B")
protein_h1_bed <- vroom::vroom(file = paste0(remote,path,"gene_h1.bed"), col_names = bed_cols_h1)
protein_h2_bed <- vroom::vroom(file = paste0(remote,path,"gene_h2.bed"), col_names = bed_cols_h2)


coordinates_protein_hits <- full_join(protein_h1_bed, protein_hits, by="protein_A") %>% 
  full_join(., protein_h2_bed, by="protein_B")

# Get homology  by identiity and Synteny by block procedence
coordinates_protein_hits <- coordinates_protein_hits %>% 
  mutate(
    Identity = if_else(
      !is.na(SDR_block_A) & !is.na(SDR_block_B) &
        pident_AvsB > 80 & pident_BvsA > 80,
      "Homolog", "Private"
    ),
    Synteny = case_when(
      Identity =="Homolog"& SDR_block_A == SDR_block_B ~ "Syntenic",
      Identity == "Homolog"& SDR_block_A != SDR_block_B ~ "Rearranged",
      TRUE ~ "Non-syntenic"
    )
  )

# Remove "ID=" prefix from SDR_block_A and SDR_block_B
coordinates_protein_hits <- coordinates_protein_hits %>%
  mutate(
    SDR_block_A = sub("^ID=", "", SDR_block_A),
    SDR_block_B = sub("^ID=", "", SDR_block_B)
  )

# Get Rearrgameed but syntenic proteins
coordinates_protein_hits %>%   
  group_by(Identity, Synteny) %>%  
  tally()


#  Do plots
a <- plot_synteny_bump(coordinates_protein_hits, genome_labels = c(A = "SDR-Y", B = "SDR-X"))

order_SDR_block <- c("M1","m1","TIB1","m2","M2","B2","M3","m3","TIB3","m4", 
                     "M4", "TB4", "M5","m5","I5","m6", "M6")

b <- plot_dodged_counts(coordinates_protein_hits, sdr_col = "SDR_block_A", 
                        order_SDR_block = order_SDR_block, counting_col = "Synteny")
c <- plot_dodged_counts(coordinates_protein_hits, sdr_col = "SDR_block_B",
                        order_SDR_block = order_SDR_block, counting_col = "Synteny")
# Patchwork
(b / c)

# Final classification file:
Final_classification_table <- coordinates_protein_hits %>% 
  select(protein_A, protein_B, Identity, Synteny, SDR_block_A, SDR_block_B) %>% 
  dplyr::rename(
    hap1 = protein_A,
    hap2 = protein_B
  )

write_delim(Final_classification_table,
            file=paste0(remote,path,"syntenic_classification_SDR.tsv"),
            delim = "\t", col_names = T,quote = "none")