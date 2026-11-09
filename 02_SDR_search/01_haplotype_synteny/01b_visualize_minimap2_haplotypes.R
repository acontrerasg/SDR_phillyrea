#!/usr/bin/env Rscript
library(pafr)
library(vroom)
library(dplyr)
library(patchwork)

remote_path="PATH/TO/FOLDER"
minimap_path="PATH/TO/FOLDER/" # nolint

#load up all paf files merged together
all_paf <- pafr::read_paf(paste0(remote_path,minimap_path,"Pang_haplotypes_asm5.paf"))

# Create function for synteny plot in pafr
create_synteny_plot <- function(q_chrom, t_chrom, title) {
  all_paf %>%
    filter(tp == "P") %>%
    filter(mapq == 60) %>% 
    filter(alen > 10^5) %>% 
    plot_synteny(
      q_chrom = q_chrom, 
      t_chrom = t_chrom, 
      centre = FALSE
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(
      title = title,
      x = "Position in Mb"
    )
}

plots <- list(
  create_synteny_plot("Pan_040.h2-12", "Pan_040.h1-12", "Synteny Alignment between \n Male hap1 and hap2"),
  create_synteny_plot("Pan_01N25.h2-12", "Pan_01N25.h1-12", "Synteny Alignment between \n Hermaphrodite hap1 and hap2"),
  
  create_synteny_plot("Pan_01N25.h1-12", "Pan_040.h1-12", "Synteny Alignment between \n hermaphrodite hap1 and Male hap1"),
  create_synteny_plot("Pan_01N25.h2-12", "Pan_040.h1-12", "Synteny Alignment between \n hermaphrodite hap2 and Male hap1"),
  
  create_synteny_plot("Pan_01N25.h1-12", "Pan_040.h2-12", "Synteny Alignment between \n hermaphrodite hap1 and Male hap2"),
  create_synteny_plot("Pan_01N25.h2-12", "Pan_040.h2-12", "Synteny Alignment between \n hermaphrodite hap2 and Male hap2")
  
)


(plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plots[[6]]) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  )

# Combine plots
# Reduce(`+`, plots)
combined_plot_2x3 <- (plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plots[[6]]) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  )

combined_plot_3x2 <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]]) /
  (plots[[5]] | plots[[6]]) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
  )
# 2x3 layout - wider than tall
ggsave("synteny_plots_2x3.pdf", 
       combined_plot_2x3, 
       width = 24,      
       height = 16,     
       units = "in",
       limitsize = FALSE)
# 3x2 layout - taller than wide
ggsave("synteny_plots_3x2.pdf", 
       combined_plot_3x2, 
       width = 16,      
       height = 24,     
       units = "in",
       limitsize = FALSE)


# Narrow down the SDR:
all_paf %>%
  filter(tp == "P") %>% 
  filter(alen > 10^4) %>% 
  filter(mapq == 60) %>% 
  filter(qname == "Pan_040.h2-12") %>%
  filter(tname == "Pan_040.h1-12") %>%
  filter(qstart > 12.5*10^6 | qstart < 7*10^6) %>%
  plot_synteny(q_chrom = "Pan_040.h2-12", t_chrom = "Pan_040.h1-12")


# Get the  syntenic regions  when comparing with Pan H2
data <- all_paf %>%
  filter(tp == "P") %>% 
  filter(alen > 10^4) %>% 
  filter(mapq == 60) %>% 
  filter(qname == "Pan_040.h2-12") %>%
  filter(tname == "Pan_040.h1-12") %>%
  filter(qstart > 12.5*10^6 | qstart < 7*10^6)  %>%  as_tibble()
# Transform to bed:
syntenyc_region_040_h2 <- data %>%
  select(
    chrom = qname,
    start = qstart,
    end = qend,
    name = tname,      
    score = mapq,
    strand = strand
  )
syntenyc_region_040_h1 <- data %>%
  select(
    chrom = tname,
    start = tstart,
    end = tend,
    name = qname,      # using query name as the name field
    score = mapq,
    strand = strand
  )
################################################################
# FOR HAP1 01N25 
data <- all_paf %>%
  filter(tp == "P") %>% 
  filter(alen > 10^5) %>% 
  filter(mapq == 60) %>% 
  filter(qname == "Pan_01N25.h1-12") %>%
  filter(tname == "Pan_040.h1-12") %>%
  filter(qstart > 13*10^6 | qstart < 7*10^6) %>%
  as_tibble()
syntenyc_region_01N25_h1 <- data %>%
  select(
    chrom = qname,
    start = qstart,
    end = qend,
    name = tname,      
    score = mapq,
    strand = strand
  )

data <- all_paf %>%
  filter(tp == "P") %>% 
  filter(alen > 10^5) %>% 
  filter(mapq == 60) %>% 
  filter(qname == "Pan_01N25.h2-12") %>%
  filter(tname == "Pan_040.h1-12") %>%
  filter(qstart > 12.5*10^6 | qstart < 7*10^6) %>%
  as_tibble()
syntenyc_region_01N25_h2 <- data %>%
  select(
    chrom = qname,
    start = qstart,
    end = qend,
    name = tname,      
    score = mapq,
    strand = strand
  )

#### Write syntenic regions to bed:
readr::write_tsv(syntenyc_region_040_h1,paste0(remote_path,minimap_path,"Pang_040_h1_syntenic.coord.bed"), col_names = FALSE)
readr::write_tsv(syntenyc_region_040_h2,paste0(remote_path,minimap_path,"Pang_040_h2_syntenic.coord.bed"), col_names = FALSE)
readr::write_tsv(syntenyc_region_01N25_h1,paste0(remote_path,minimap_path,"Pang_01N25_h1_syntenic.coord.bed"), col_names = FALSE)
readr::write_tsv(syntenyc_region_01N25_h2,paste0(remote_path,minimap_path,"Pang_01N25_h2_syntenic.coord.bed"), col_names = FALSE)




# Analize Male hap1 vs hap2 rearragments


all_paf %>%
  filter(tp == "P") %>% 
  filter(alen > 10^4) %>% 
  filter(mapq == 60) %>% 
  filter(qname == "Pan_040.h2-12") %>%
  filter(tname == "Pan_040.h1-12") %>%
  filter(qstart < 12.5*10^6 ) %>% 
  filter(qstart > 7*10^6) %>% 
  ggplot(aes(x=alen, y=dv, color=strand)) + 
  geom_point() + 
  scale_color_manual(
    values = c("+" = "blue", "-" = "red"),
    labels = c("+" = "Collinear", "-" = "Inverted"),
    name = NULL) + 
  scale_x_log10(breaks = c(1e4, 1e5, 1e6), 
                labels = c("10Kbp", "100Kbp", "1Mbp")) +
  xlab("SDR Rearrangement Length") +
  ylab("Approximate per-base sequence divergence") +
  theme_light(base_size = 12)

syntenic_bits <- all_paf %>%
  filter(tp == "P") %>% 
  filter(alen > 10^5) %>% 
  filter(mapq == 60) %>% 
  filter(qname == "Pan_040.h2-12") %>%
  filter(tname == "Pan_040.h1-12") %>%
  filter(qstart > 12.5*10^6 | qstart < 7*10^6) 

# Save file with rearramgents

rearragment_data <- all_paf %>%
  filter(tp == "P") %>% 
  filter(alen > 10^4) %>% 
  filter(mapq == 60) %>% 
  filter(qname == "Pan_040.h2-12") %>%
  filter(tname == "Pan_040.h1-12") %>%
  filter(qstart < 12.5*10^6 ) %>% 
  filter(qstart > 7*10^6)  %>%  as_tibble()


readr::write_tsv(rearragment_data,paste0(remote_path,minimap_path,"rearragments_Male_SDR.tsv"), col_names = TRUE)



all_paf %>%
  filter(tp == "P") %>% 
  filter(alen > 10^5) %>% 
  filter(mapq == 60) %>% 
  filter(qname == "Pan_040.h2-12") %>%
  filter(tname == "Pan_040.h1-12") %>%
  filter(qstart < 12.5*10^6 ) %>% 
  filter(qstart > 7*10^6) %>% 
  ggplot(aes(x=alen, y=dv, color=strand)) + 
  geom_point()  +  geom_point(data = syntenic_bits, aes(x=alen, y=dv, color="green")) + 
  scale_x_log10(breaks = c(1e4, 1e5, 1e6,1e7), 
                labels = c("10Kbp", "100Kbp", "1Mbp","10Mbp")) +
  xlab("SDR Rearrangement Length") +
  ylab("Approximate per-base \n sequence divergence") +
  theme_light(base_size = 12) +    scale_color_manual(
    values = c("+" = "blue", "-" = "red"),
    labels = c("+" = "Collinear", "-" = "Inverted"),
    name = NULL)
