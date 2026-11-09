library(pafr)
library(vroom)
library(dplyr)
library(patchwork)
library(stringr)
library(gtools)
# Remote path

path="PATH/TO/FOLDER"

# List all .paf files in the folder
paf_files <- list.files(path, pattern = "_asm5.paf$", full.names = TRUE)

# Apply the read_paf function to each file and store the results in a list
paf_data_list <- lapply(paf_files, read_paf)

combined_paf_asm5 <- bind_rows(paf_data_list, .id = "source_file")

# Filter secondary  and low quality mapping rows
combined_paf_asm5 <- combined_paf_asm5 %>%  filter(tp == "P") %>%  filter(mapq == 60)
View(combined_paf_asm5)

# Dotplots all
a <- combined_paf_asm5 %>%
  filter(source_file == 1) %>%
  mutate(
    # Remove prefix "Pla.h1-" and suffix "_RagTag"
    clean_qname = str_remove_all(qname, "^Pla\\.h1-|_RagTag$"),
    # Extract numeric part (after cleaning)
    order_num = as.numeric(clean_qname)
  ) %>%
  arrange(order_num) %>%
  dotplot(label_seqs = TRUE, 
          order_by = "qstart", 
          line_size = 1, 
          alignment_colour = "dodgerblue") +
  theme_light(base_size = 16)

b <- combined_paf_asm5 %>%
  filter(source_file == 2) %>%
  mutate(
    # Remove prefix "Pla.h1-" and suffix "_RagTag"
    clean_qname = str_remove_all(qname, "^Pla\\.h1-|_RagTag$"),
    # Extract numeric part (after cleaning)
    order_num = as.numeric(clean_qname)
  ) %>%
  arrange(order_num) %>%
  dotplot(label_seqs = TRUE, 
          order_by = "qstart", 
          line_size = 1, 
          alignment_colour = "dodgerblue") +
  theme_light(base_size = 16)

c <- combined_paf_asm5 %>%
  filter(source_file == 3) %>%
  mutate(
    # Remove prefix "Pla.h1-" and suffix "_RagTag"
    clean_qname = str_remove_all(qname, "^Pla\\.h2-|_RagTag$"),
    # Extract numeric part (after cleaning)
    order_num = as.numeric(clean_qname)
  ) %>%
  arrange(order_num) %>%
  dotplot(label_seqs = TRUE, 
          order_by = "qstart", 
          line_size = 1, 
          alignment_colour = "dodgerblue") +
  theme_light(base_size = 16)

d <- combined_paf_asm5 %>%
  filter(source_file == 1) %>%
  mutate(
    # Remove prefix "Pla.h1-" and suffix "_RagTag"
    clean_qname = str_remove_all(qname, "^Pla\\.h2-|_RagTag$"),
    # Extract numeric part (after cleaning)
    order_num = as.numeric(clean_qname)
  ) %>%
  arrange(order_num) %>%
  dotplot(label_seqs = TRUE, 
          order_by = "qstart", 
          line_size = 1, 
          alignment_colour = "dodgerblue") +
  theme_light(base_size = 16)

# Merge plots
(a | b) / (c | d)

# Dotplots of chrom 12 
a <- combined_paf_asm5 %>% 
  filter(source_file == 1) %>% filter(alen > 10^1) %>%   
  filter(tname=="Pan_040.h1-12") %>% 
  filter(qname=="Pla.h1-12_RagTag") %>%  
  dotplot(label_seqs = TRUE)

b <- combined_paf_asm5 %>% 
  filter(source_file == 2) %>% filter(alen > 10^1) %>%   
  filter(tname=="Pan_040.h2-12") %>% 
  filter(qname=="Pla.h1-12_RagTag") %>%  
  dotplot(label_seqs = TRUE)

c <- combined_paf_asm5 %>% 
  filter(source_file == 3) %>% filter(alen > 10^1) %>%   
  filter(tname=="Pan_040.h1-12") %>% 
  filter(qname=="Pla.h2-12_RagTag") %>%  
  dotplot(label_seqs = TRUE)

d <- combined_paf_asm5 %>% 
  filter(source_file == 4) %>% filter(alen > 10^1) %>%   
  filter(tname=="Pan_040.h2-12") %>% 
  filter(qname=="Pla.h2-12_RagTag") %>%  
  dotplot(label_seqs = TRUE)
# Merge plots
(a | b) / ( c | d)