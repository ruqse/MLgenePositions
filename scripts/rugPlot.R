
# The biomaRt package requires BiocMananger for installation
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")


# Load necessary libraries
library(tidyverse)
library(data.table)
library(biomaRt)
library(glue)
library(ggplot2)

# Set working directory to the parent directory of the active RStudio document
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/../"))


# Configure biomaRt for the wormbase database
mart <- useMart("parasite_mart", host = "https://parasite.wormbase.org", 
                port = 443, dataset = "wbps_gene")



# Define attributes to retrieve
attributes <- c('wbps_gene_id', 'external_gene_id', 'chromosome_name', 'start_position', 'end_position')


# # Define gene names and remove duplicates
gene_names <- c("glc-1", "glc-2", "glc-3","glc-5", "avr-14", "avr-15", 
                "glc-4",  "unc-7", "unc-9", "osm-1", "osm-5", "dyf-11",
                "che-1", "che-2", "che-3", "che-6", "che-10", "che-11", 
                "che-12", "che-13", "che-14", "osm-1", "osm-3", "osm-5", 
                "osm-6", "osm-12", "daf-10", "daf-19", "dyf-1", "dyf-2", 
                "dyf-3", "dyf-4", "dyf-5", "dyf-6", "dyf-7", "dyf-9", "dyf-10", 
                "dyf-11", "dyf-13", "mec-1", "mec-8", "bbs-1", "bbs-8", "dhc-3", 
                "che-2", "che-3", "che-10", "che-11", "che-12", "che-13", "daf-6", 
                "daf-10", "daf-19", "mec-8", "osm-1", "osm-3", "osm-5", "osm-6", "dyf-1", 
                "dyf-2", "dyf-3", "dyf-4", "dyf-5", "dyf-6", "dyf-7", "dyf-8", "dyf-9", "dyf-10", 
                "dyf-11", "dyf-12", "dyf-13", "unc-33", "unc-44", "unc-13", "che-3", "daf-10", 
                "osm-3", "che-3", "daf-10", "osm-3", "pgp-1", "pgp-2", "pgp-3", "pgp-5", "pgp-6",
                "pgp-9", "pgp-10", "pgp-11", "pgp-12","pgp-14","mrp-1", "mrp-3", "mrp-6", "mrp-8", 
                "haf-4", "haf-9","pmp-4", "pmp-5","lgc-37") # Replace as desired

gene_names <- unique(gene_names)


# s# Define filters and corresponding values
filters <- c("gene_name", "species_id_1010") 
values <- list(gene_names, "caelegprjna13758")

# # Fetch gene data
gene_data <- getBM(attributes = attributes, 
                   filters = filters, 
                   values = values, 
                   mart = mart)

# Process the fetched gene data for better readability and clarity
candidates <- gene_data %>% 
  dplyr::rename(
    GeneID = wbps_gene_id,
    Gene = external_gene_id,
    CHROM = chromosome_name,
    startPOS = start_position,
    endPOS = end_position
  ) %>% 
  dplyr::select(-GeneID) %>% 
  dplyr::mutate(
    Gene_class = case_when(
      Gene %in% c("glc-1", "glc-2", "glc-3", "glc-4","avr-14", "avr-15") ~ "GluCl",
      startsWith(Gene, "pgp") | startsWith(Gene, "pmp") | startsWith(Gene, "haf") ~ "Efflux pump",
      startsWith(Gene,"lgc")~ "GABA", TRUE ~ "Dyf"
    )
  )

  # factorize levels of Gene_class

candidates$Gene_class <- factor(candidates$Gene_class, levels = c("GluCl", "GABA",
                                                                  "Efflux pump", "Dyf"))

# Define the Gene_class colors
rug_colors <- c("GluCl" = "red", 
                "Efflux pump" = "blue", 
                "GABA" = "#FF8C00", 
                "Dyf" = "grey")


# get the chromosome lengths for c.elegans
chr_lens <- data.table::fread("data/c_elegans_chr_lengths.tsv") %>% 
  dplyr::select(CHROM, startPOS = start, endPOS = stop)

# shape chr_lens for plotting
chr_lens2 <- rbind(chr_lens) %>%
  tidyr::pivot_longer(cols = -CHROM, names_to = "poop", values_to = "POS") %>%
  dplyr::select(-poop) 

# plot the gene locations and the chromosome lengths
rugPlot <- ggplot() + 
  theme_bw() + 
  geom_bar(data = chr_lens2, mapping = aes(x = POS/1000000), alpha = 0.0) +
  facet_grid(.~CHROM, scales = "free_x", space = "free") +
  geom_segment(data = candidates, aes(x = startPOS/1000000, xend = endPOS/1000000, y = -Inf, yend = Inf,
                                      color = Gene_class), linetype=1, size = 0.4)+
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  scale_color_manual(values = rug_colors)+
  #facet_grid(Gene_class ~ CHROM, scales = "free_x", space = "free")+
  labs(x =  "Genomic position (Mb)", y = NULL)+
  theme(legend.position = "none", 
        panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_text(size = 12, face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),# remove y axis labels
        axis.ticks.y = element_blank(),
        # axis.ticks.x = element_blank(),# remove y axis ticks
        strip.text.x = element_blank(),
        strip.text.y = element_blank())

