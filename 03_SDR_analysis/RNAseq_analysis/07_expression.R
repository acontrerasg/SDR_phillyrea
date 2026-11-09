#!/usr/bin Rscript
library(tidyverse)
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
library(vroom)
library(patchwork)
library(gtools)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

# Color schemes:
sex_colors <- c("Male" = "#5ab4ac", "Herm" = "#d8b365", "Hermaphrodite" = "#d8b365")

## FUNCTIONS
# Define column names for GFF format
gff_cols <- c(
  "seqid",      # chromosome or scaffold
  "source",     # name of the program or database
  "type",       # feature type
  "start",      # start position
  "end",        # end position
  "score",      # numeric value
  "strand",     # + (forward) or - (reverse)
  "phase",      # frame 0,1,2
  "attributes"  # semicolon-separated key-value pairs
)


#############
# FUNCTIONS #
#############
# Generate degs

generate_deg_phenotype <- function(
    count_path,
    phenotype_df,
    remove_samples = NULL,
    phenotype_levels = NULL,
    dds_name = "dds",
    res_name = "res",
    df_name = "df"
) {
  
  # Load Count Matrix
  cts <- as.matrix(
    read.csv(count_path, sep = "\t", row.names = "GENEID")
  )
  
  
  cts_removed <- cts[startsWith(rownames(cts), "__"), ]
  cts <- cts[!startsWith(rownames(cts), "__"), ]
  
  # Optionally remove samples
  if (!is.null(remove_samples)) {
    cts <- cts[, !colnames(cts) %in% remove_samples, drop = FALSE]
  }
  
  # Prepare stats DF used as colData
  df_total <- data.frame(
    sample = colnames(cts),
    total_count = colSums(cts)
  )
  
  df_removed <- data.frame(
    feature = rownames(cts_removed),
    cts_removed,
    check.names = FALSE
  )
  rownames(df_removed) <- NULL
  
  df_removed_t <- data.frame(t(df_removed[, -1, drop = FALSE]))
  colnames(df_removed_t) <- df_removed$feature
  df_removed_t <- tibble::rownames_to_column(df_removed_t, "sample")
  
  df_stats <- dplyr::full_join(df_removed_t, df_total, by = "sample") %>%
    mutate(sample = stringr::str_remove(sample, "_counts"))
  
  
  df_stats <- dplyr::left_join(df_stats, phenotype_df, by = "sample")
  
  # Load phenotype file
  if (!is.null(remove_samples)) {
    df_stats <- df_stats %>% dplyr::filter(!sample %in% remove_samples)
  }
  
  rownames(df_stats) <- df_stats$sample
  
  # Clean column names in count matrix ----
  colnames(cts) <- gsub("_counts", "", colnames(cts))
  
  stopifnot(all(rownames(df_stats) == colnames(cts)))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData   = df_stats,
    design    = ~ phenotype
  )
  
  #  User-defined factor levels
  if (!is.null(phenotype_levels)) {
    dds$phenotype <- factor(dds$phenotype, levels = phenotype_levels)
  } else {
    dds$phenotype <- factor(dds$phenotype)
  }
  
  # Pre-filter low count genes:
  min_group_size <- min(table(dds$phenotype))
  keep <- rowSums(counts(dds) >= 10) >= min_group_size
  dds <- dds[keep, ]
  #keep <- rowSums(counts(dds) >= 10) >= 6
  #dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Return objects under user-provided names ----
  df <- as.data.frame(res)
  df$gene <- rownames(res)
  
  out <- list()
  out[[dds_name]] <- dds
  out[[res_name]] <- res
  out[[df_name]] <- df
  
  return(out)
}

# Generate customized summaries for DEGs
summarize_degs <- function(res,
                           padj_cutoff = 0.01,
                           lfc_cutoff = 1,
                           lfc_col = "log2FoldChange",
                           padj_col = "padj") {
  
  # Convert to dataframe if needed
  res_df <- as.data.frame(res)
  
  # Remove NA padj
  res_df <- res_df[!is.na(res_df[[padj_col]]), ]
  
  # Define significance
  up <- sum(res_df[[padj_col]] < padj_cutoff &
              res_df[[lfc_col]] > lfc_cutoff,
            na.rm = TRUE)
  
  down <- sum(res_df[[padj_col]] < padj_cutoff &
                res_df[[lfc_col]] < -lfc_cutoff,
              na.rm = TRUE)
  
  total_sig <- up + down
  
  summary_table <- data.frame(
    Category = c("Upregulated", "Downregulated", "Total Significant"),
    Count = c(up, down, total_sig)
  )
  
  return(summary_table)
}

# Generate Volcano plot
plot_volcano <- function(df,
                         top_n = 10,
                         padj_cutoff = 0.01,
                         lfc_cutoff = 1,
                         lfc_col = "log2FoldChange",
                         padj_col = "padj",
                         gene_col = "gene",
                         title = "Volcano Plot",
                         alpha = 0.7,
                         point_size = 2,
                         colors = c("Upregulated" = "#5ab4ac",
                                    "Downregulated" = "#d8b365",
                                    "noDEG" = "grey80")) {
  
  # Remove NA padj
  df <- df %>%
    filter(!is.na(.data[[padj_col]]))
  
  # Create DEG classification
  df <- df %>%
    mutate(
      DEG_class = case_when(
        .data[[padj_col]] < padj_cutoff &
          .data[[lfc_col]] > lfc_cutoff  ~ "Upregulated",
        
        .data[[padj_col]] < padj_cutoff &
          .data[[lfc_col]] < -lfc_cutoff ~ "Downregulated",
        
        TRUE ~ "noDEG"
      ),
      neglog10 = -log10(.data[[padj_col]])
    )
  
  # Select top genes ONLY among DEGs
  top_genes <- df %>%
    filter(DEG_class != "noDEG") %>%
    arrange(.data[[padj_col]]) %>%
    slice_head(n = top_n)
  
  p <- ggplot(df,
              aes(x = .data[[lfc_col]],
                  y = neglog10,
                  color = DEG_class)) +
    
    geom_point(alpha = alpha, size = point_size) +
    
    scale_color_manual(values = colors) +
    
    # Threshold lines
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               color = "grey40") +
    
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               color = "grey40") +
    
    # Label top DEGs
    geom_label_repel(
      data = top_genes,
      aes(
        x = .data[[lfc_col]],
        y = neglog10,
        label = .data[[gene_col]]
      ),
      inherit.aes = FALSE,
      color = "black",
      fill = "white",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      size = 3.5,
      max.overlaps = Inf
    ) +
    
    labs(title = title,
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value",
         color = NULL) +
    
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

plot_deg_lfc <- function(df,
                         title = "",
                         legend_title = "Expression",
                         colors = c(
                           noDEG = "grey60",
                           Up = "#5ab4ac",
                           Down = "#d8b365"
                         ),
                         jitter_width = 0.1) {
  
  df %>%
    dplyr::mutate(
      point_class = dplyr::case_when(
        DEG_class == "noDEG" ~ "noDEG",
        DEG_class == "DEG" & log2FoldChange > 0 ~ "Up",
        DEG_class == "DEG" & log2FoldChange < 0 ~ "Down",
        TRUE ~ "noDEG"
      )
    ) %>%
    ggplot(aes(x = Identity,
               y = log2FoldChange,
               color = point_class)) +
    geom_jitter(width = jitter_width) +
    scale_color_manual(
      values = colors,
      name = legend_title
    ) +
    theme_light() +
    ggtitle(title)
}


########
# PATHS
remote='PATH/TO/FOLDER'
path_counts_1='PATH/TO/FILE'
path_counts_2='PATH/TO/FILE'
# Fet phenotype
path="PATH/TO/FOLDER/"
file="phenotype.tsv"
phenotype=vroom::vroom(paste0(remote,path,file), delim = "\t", col_names = T)
phenotype <- phenotype %>%
  dplyr::mutate(
    dplyr::across(where(is.character), as.factor)
  )
############################################

# Call DESEQ2 objects with contrasts using the sex phenotype

out_h1_sex <- generate_deg_phenotype(
  count_path = paste0(remote, path_counts_1),
  phenotype_df = phenotype, 
  remove_samples = c("Pa03", "Pa63", "Pa64", "Pa65"),
  phenotype_levels = c("Hermaphrodite", "Male")
)

out_h2_sex <- generate_deg_phenotype(
  count_path = paste0(remote, path_counts_2),
  phenotype_df = phenotype, 
  remove_samples = c("Pa03", "Pa63", "Pa64", "Pa65"),
  phenotype_levels = c("Hermaphrodite", "Male")
)

# Apply custom tresholds of padj 0.01 and lfc of 1.
summarize_degs(out_h1_sex$res, padj_cutoff = 0.01 , lfc_cutoff = 1)
summarize_degs(out_h2_sex$res, padj_cutoff = 0.01 , lfc_cutoff = 1)
