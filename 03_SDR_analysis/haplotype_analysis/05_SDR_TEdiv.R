#!/usr/bin Rscript
library(dplyr)
library(stringr)
library(ggplot2)
library(vroom)
library(patchwork)
library(ggstatsplot)

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
# FUNCTIONS
# Get SDR annotation
process_TE_gff <- function(remote, path, file, gff_cols) {
  ## read + annotate
  TE_anno <- vroom::vroom(file = paste0(remote, path, file),
                           delim = "\t",
                           col_names = gff_cols,
                          comment = "#") %>%
    mutate(
      TEID = str_remove(str_remove(attributes, "ID="), ";.*$"),
      Perc_Div = as.numeric(str_remove(str_remove(attributes, ".*Perc_Div="), ";.*")))
  ## return final tibble
  return(TE_anno)
}
# Wilcox test 
wilcox_test_by_group <- function(df, group_col, value_col, condition_col) {
  
  # Validate that columns exist in dataframe
  if (!group_col %in% names(df)) {
    stop(paste("Column", group_col, "not found in dataframe"))
  }
  if (!value_col %in% names(df)) {
    stop(paste("Column", value_col, "not found in dataframe"))
  }
  if (!condition_col %in% names(df)) {
    stop(paste("Column", condition_col, "not found in dataframe"))
  }
  
  # Check that condition column has exactly 2 unique values
  unique_conditions <- unique(df[[condition_col]])
  unique_conditions <- unique_conditions[!is.na(unique_conditions)]
  
  if (length(unique_conditions) != 2) {
    stop(paste("Condition column", condition_col, "must have exactly 2 unique values.",
               "Found:", length(unique_conditions), 
               paste0("(", paste(unique_conditions, collapse = ", "), ")")))
  }
  
  condition1 <- unique_conditions[1]
  condition2 <- unique_conditions[2]
  
  # Get all unique groups
  group_vector <- unique(df[[group_col]])
  group_vector <- group_vector[!is.na(group_vector)]
  
  # Print header
  cat("Wilcoxon Test Results by", group_col, "\n")
  cat("Comparing:", condition1, "vs", condition2, "\n")
  cat("=================================\n\n")
  
  # Initialize results dataframe
  results <- data.frame()
  
  # Loop through each group
  for (group_name in group_vector) {
    # Subset data for this group
    group_data <- df[df[[group_col]] == group_name & !is.na(df[[group_col]]), ]
    
    # Extract values for each condition
    group1_values <- group_data[[value_col]][group_data[[condition_col]] == condition1 & 
                                               !is.na(group_data[[condition_col]])]
    group2_values <- group_data[[value_col]][group_data[[condition_col]] == condition2 & 
                                               !is.na(group_data[[condition_col]])]
    
    # Remove NAs from values
    group1_values <- group1_values[!is.na(group1_values)]
    group2_values <- group2_values[!is.na(group2_values)]
    
    # Only test if both groups have data
    if (length(group1_values) > 0 & length(group2_values) > 0) {
      # Perform Wilcoxon rank-sum test
      test_result <- wilcox.test(group1_values, group2_values)
      
      # Print results
      cat(group_name, "\n")
      cat("  ", condition1, ": n =", length(group1_values), 
          ", median =", round(median(group1_values), 2), "\n")
      cat("  ", condition2, ": n =", length(group2_values), 
          ", median =", round(median(group2_values), 2), "\n")
      cat("  W =", test_result$statistic, 
          ", p-value =", format(test_result$p.value, scientific = TRUE))
      
      if (test_result$p.value < 0.05) {
        cat(" *\n")
      } else {
        cat("\n")
      }
      cat("\n")
      
      # Store results
      result_row <- data.frame(
        Group = group_name,
        n_cond1 = length(group1_values),
        n_cond2 = length(group2_values),
        median_cond1 = median(group1_values),
        median_cond2 = median(group2_values),
        W_statistic = as.numeric(test_result$statistic),
        p_value = test_result$p.value,
        significant = ifelse(test_result$p.value < 0.05, "Yes", "No"),
        stringsAsFactors = FALSE
      )
      
      # Rename columns to include actual condition names
      names(result_row)[2] <- paste0("n_", condition1)
      names(result_row)[3] <- paste0("n_", condition2)
      names(result_row)[4] <- paste0("median_", condition1)
      names(result_row)[5] <- paste0("median_", condition2)
      
      results <- rbind(results, result_row)
      
    } else {
      cat(group_name, "\n")
      cat("  Warning: Insufficient data for comparison\n")
      cat("    ", condition1, ": n =", length(group1_values), "\n")
      cat("    ", condition2, ": n =", length(group2_values), "\n\n")
    }
  }
  
  # Print summary table
  if (nrow(results) > 0) {
    cat("\n=== Summary Table ===\n")
    print(results, row.names = FALSE)
  } else {
    cat("\nNo valid comparisons could be performed.\n")
  }
  
  # Return results invisibly
  invisible(results)
}

# Get SDR annotation:
remote='PATH/TO/FOLDER'
path="PATH/TO/FOLDER"
# Hap1
file_h1="Pan_040.hap1-rt_SDR_TEs.gff3"
SDR_TEs_hap1 <- process_TE_gff(remote, path, file_h1, gff_cols) %>% filter(!grepl("CA3p|CA5p", attributes))
# Hap2
file_h2="Pan_040.hap2-rt_SDR_TEs.gff3"
SDR_TEs_hap2 <- process_TE_gff(remote, path, file_h2, gff_cols)  %>% filter(!grepl("CA3p|CA5p", attributes))

# Hap1 TEs
path="PATH/TO/FOLDER"
All_TEs_hap1 <- process_TE_gff(remote, path, "Pan_040.hap1-rt_Final_TEanno_Chrom.gff3", gff_cols)
# Hap2 TEs
All_TEs_hap2 <- process_TE_gff(remote, path, "Pan_040.hap2-rt_Final_TEanno_Chrom.gff3", gff_cols)

# Add SDR info:
All_TEs_hap1 <- All_TEs_hap1 %>%
  mutate(SDR = ifelse(TEID %in% SDR_TEs_hap1$TEID, "SDR", "background"))

All_TEs_hap2 <- All_TEs_hap2 %>%
  mutate(SDR = ifelse(TEID %in% SDR_TEs_hap2$TEID, "SDR", "background"))

# Prepare the data
test_data1 <- All_TEs_hap1 %>%
  filter(type %in% c("DNA/DTM", "LTR/Gypsy", "LTR/Copia"))
test_data2 <- All_TEs_hap2 %>%
  filter(type %in% c("DNA/DTM", "LTR/Gypsy", "LTR/Copia"))

# Calculate statistics
wilcox_test_by_group(
  df = test_data1,
  group_col = "type",
  value_col = "Perc_Div",
  condition_col = "SDR"
)

wilcox_test_by_group(
  df = test_data2,
  group_col = "type",
  value_col = "Perc_Div",
  condition_col = "SDR"
)

# Plot Divergence
plot1 <- grouped_ggbetweenstats(
  data = test_data1,
  x = SDR,
  y = Perc_Div,
  grouping.var = type,
  xlab = "Region",
  ylab = "Percent Divergence",
  type = "nonparametric",  # Wilcoxon test
  pairwise.display = "all",
  centrality.type="nonparametric", 
  violin.args = list(width = 0),
  boxplot.args = list(width = 0.2),
  point.args = list(alpha = 0),  # Remove dots completely
  package = "ggsci",
  ggtheme = ggplot2::theme_bw(), 
  palette = "nrc_npg"
)
# Create grouped plot without dots
plot2 <- grouped_ggbetweenstats(
  data = test_data2,
  x = SDR,
  y = Perc_Div,
  grouping.var = type,
  xlab = "Region",
  ylab = "Percent Divergence",
  type = "nonparametric",  # Wilcoxon test
  pairwise.display = "all",
  centrality.type="nonparametric", 
  violin.args = list(width = 0),
  boxplot.args = list(width = 0.2),
  point.args = list(alpha = 0),  # Remove dots completely
  package = "ggsci",
  ggtheme = ggplot2::theme_bw(),
  palette = "nrc_npg"
)
plot1/plot2
