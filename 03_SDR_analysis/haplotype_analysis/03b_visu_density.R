#!/usr/bin/env Rscript
remote_path="PATH/TO/FOLDER"
path="PATH/TO/FOLDER/"



# Define helper functions
process_density_data <- function(file_name, col_names) {
  vroom(paste0(file_name), col_names = col_names) %>%
    mutate(
      winID = str_remove(winID, pattern = "ID=.*_"),
      winID = case_when(
        winID == 1 ~ "CA3p",
        winID == 2 ~ "TIB1",
        winID == 3 ~ "B2",
        winID == 4 ~ "TIB3",
        winID == 5 ~ "TB4",
        winID == 6 ~ "I5",
        winID == 7 ~ "CA5p",
        TRUE ~ "background" # Renaming all other values as "background"
      ),
      # Set factor with custom order
      winID = factor(winID, levels = c("background", "CA3p", "TIB1", "B2", "TIB3", "TB4", "I5", "CA5p"))
    )
}



create_density_plot <- function(data, title) {
  # Define a custom color palette for winID levels
  winid_colors <- c(
    "background" = "black",  
    "CA3p" = "#8B8680",        
    "TIB1" = "#9999CC",        
    "B2" = "#2EB4E6",         
    "TIB3" ="#9999CC",        
    "TB4" = "#2EB4E6",         
    "I5" = "#9999CC",         
    "CA5p" = "#8B8680"         
  )
  
  data %>%
    ggplot(aes(x = winID, y = density, fill = winID)) +
    geom_boxplot() +
    facet_wrap(~type) +
    coord_flip() + 
    scale_fill_manual(values = winid_colors) +
    labs(
      title = title,
      x = NULL,
      y = "Density per 100 kb window"
    ) +
    theme_light(base_size = 16) +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black"),
      panel.border = element_rect(color = "black"),
      legend.position = "none"
    )
}

# Column names for both datasets
col_names <- c("chr", "start", "end", "winID", "density", "type")

# Process and plot haplotype 1 data
h1_density <- process_density_data(
  paste0(remote_path,path,"Pan_040_h1_100000_SDR_density_results.tsv"), 
  col_names
)
h1_plot <- create_density_plot(h1_density, "Density Haplotype 1")
print(h1_plot)

# Process and plot haplotype 2 data
h2_density <- process_density_data(
  paste0(remote_path,path,"Pan_040_h2_100000_SDR_density_results.tsv"), 
  col_names
)
h2_plot <- create_density_plot(h2_density, "Density Haplotype 2")
print(h2_plot)