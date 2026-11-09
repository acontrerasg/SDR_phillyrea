#!/usr/bin/env Rscript 
# Function
make_count_stats <- function(remote, path){
  
  # Count matrix
  cts <- as.matrix(
    read.csv(
      paste0(remote, path),
      sep = "\t",
      row.names = "GENEID"
    )
  )
  
  # Remove extra rows
  cts_removed <- cts[startsWith(rownames(cts), "__"), , drop = FALSE]
  cts <- cts[!startsWith(rownames(cts), "__"), , drop = FALSE]
  
  # Total counts
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
  
  df_removed_t <- df_removed_t %>%
    tibble::rownames_to_column("sample")
  
  df_stats <- dplyr::full_join(df_removed_t, df_total, by = "sample") %>%
    dplyr::mutate(sample = stringr::str_remove(sample, "_counts"))
  
  return(df_stats)
}
# Plot Function
make_count_plots <- function(df_stats, color_scheme,plot_title = NULL){
  
  # Total count plot
  a <- df_stats %>%
    ggplot(aes(
      y = sample,
      x = total_count / 10^6,
      fill = phenotype
    )) +
    geom_col() +
    scale_fill_manual(values = color_scheme) +
    theme_light() +
    xlab("Total count in millions") +
    ggtitle("Total Count Plot")
  
  # Initialize plot list
  plot_list <- list(
    Plot_total = a
  )
  
  # Generate scatterplots
  for (col in colnames(df_stats)) {
    
    if (col %in% c(
      "total_count",
      "SI",
      "phenotype",
      "sample",
      "SI_geno"
    )) {
      next
    }
    
    p <- df_stats %>%
      ggplot(aes(
        x = total_count / 10^6,
        y = .data[[col]] / 10^6,
        color = phenotype,
        label = sample
      )) +
      geom_point() +
      geom_text(vjust = -0.5, size = 3) +
      scale_color_manual(values = color_scheme) +
      theme_light() +
      xlab("Total Count in Millions") +
      ylab(paste(col, "in Millions")) +
      ggtitle(paste("Plot of", col, "vs Total Count"))
    
    plot_list[[col]] <- p
  }
  
  # Combine plots
  combined_plot <- patchwork::wrap_plots(plot_list) +
    patchwork::plot_layout(guides = "collect")
  
  # Add overall title if provided
  if (!is.null(plot_title)) {
    combined_plot <- combined_plot +
      patchwork::plot_annotation(
        title = plot_title
      )
  }
  
  
  return(combined_plot)
}

######################

# Color scheme
color_scheme <- c(Herm="#d8b365", Male="#5ab4ac")
# Phenotype Samples
remote='PATH/TO/FOLDER'
path="PATH/TO/FOLDER/"
file="phenotype.tsv"
phenotype=vroom::vroom(paste0(remote,path,file), delim = "\t", col_names = T)
phenotype <- phenotype %>%
  mutate(phenotype = ifelse(phenotype == "Male", "Male", "Herm"))

# COUNT H1
# Path to matrix counts
remote='PATH/TO/FOLDER'
path="PATH/TO/FILE"

df_stats_h1 <- make_count_stats(remote,path)
df_stats_h1 <- dplyr::full_join(df_stats_h1, phenotype, by="sample") 

combined_plot <- make_count_plots(df_stats_h1, color_scheme, plot_title = "Counts H1")
combined_plot

### Count H2
# Path to matrix counts
remote='PATH/TO/FOLDER'
path="PATH/TO/FILE"
#
df_stats_h2 <- make_count_stats(remote,path)
df_stats_h2 <- dplyr::full_join(df_stats_h2, phenotype, by="sample") 

combined_plot <- make_count_plots(df_stats_h2, color_scheme, plot_title = "Counts H2")
combined_plot
