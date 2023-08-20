
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


# load mappings data
load("data/filtered.mappings.results.rda")

# load the independent test cuttoffs
load("data/independent_test_cutoff.rda")


# initialize a list to store the plots
plot_list <- list()

# loop over each trait
for (ivm in unique(filtered.mappings.results$trait)){
  IVM.mappings <- filtered.mappings.results %>% 
    dplyr::filter(trait %in% ivm) 
  
  
  # 
  # do we have mito mapping?
  IVM.mito_check <- IVM.mappings %>%
    na.omit()
  
  ## MANHATTAN PLOTS ##
  for.plot.IVM <- IVM.mappings %>%
    dplyr::mutate(CHROM = as.factor(CHROM)) %>%
    {
      if(!("MtDNA" %in% IVM.mito_check$CHROM)) dplyr::filter(., CHROM != "MtDNA") else .
    }
  BF.IVM <- IVM.mappings %>% 
    dplyr::group_by(trait, algorithm) %>% 
    dplyr::filter(log10p != 0) %>% 
    dplyr::distinct(marker, log10p) %>%
    dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
    dplyr::ungroup() %>%
    dplyr::select(BF) %>%
    unique(.) %>%
    dplyr::slice(1) %>% # BF can be slightly different between loco and inbred... but just plot one (5.46 v 5.47...)
    as.numeric()
  
  EIGEN <- independent_test_cutoff
  BF.frame.IVM <- IVM.mappings %>%
    dplyr::select(trait) %>%
    dplyr::filter(!duplicated(trait)) %>%
    dplyr::mutate(BF = BF.IVM, EIGEN  = EIGEN, user = unique(IVM.mappings$BF)[1])
  
  # if user selected a different threshold, use that, otherwise plot BF and EIGEN
  if(BF.frame.IVM$user %in% c(BF.frame.IVM$BF, BF.frame.IVM$EIGEN)) {
    for.plot.ann.IVM <- for.plot.IVM %>%
      dplyr::mutate(sig = case_when(log10p > BF.frame.IVM$BF ~ "BF",
                                    log10p > BF.frame.IVM$EIGEN ~ "EIGEN",
                                    TRUE ~ "NONSIG"))
    
    sig.colors <- c("#EE4266","#EE4266", "black") # changed here to differentiate BF (red) and EIGEN (#EE4266)
    names(sig.colors) <- c("BF","EIGEN", "NONSIG")
  } else {
    for.plot.ann.IVM <- for.plot.IVM %>%
      dplyr::mutate(sig = case_when(log10p > BF.frame.IVM$user ~ "user",
                                    TRUE ~ "NONSIG"))
    
    sig.colors <- c("#EE4266", "black")
    names(sig.colors) <- c("user", "NONSIG")
  }
  
  test.IVM <- BF.frame.IVM %>%
    tidyr::pivot_longer(BF:user) %>%
    dplyr::distinct() %>%
    dplyr::filter(name %in% names(sig.colors)) %>% 
    dplyr::filter(!name %in% "BF")
  
  # are we plotting mito or no?
  if("MtDNA" %in% unique(for.plot.ann.IVM$CHROM)) {
    facet_scales <- "fixed"
  } else {
    facet_scales <- "free"
  }
  
  manhattan <- ggplot() + 
    theme_bw() + 
    geom_bar(data = chr_lens2, mapping = aes(x = POS/1000000), alpha = 0.0) +
    facet_grid( .~ CHROM, scales = "free_x", space = "free") +
    geom_point(data = for.plot.ann.IVM, 
               mapping = aes(x = POS/1000000, 
                             y = log10p,
                             colour = sig,
                             alpha = sig)) +
    scale_alpha_manual(values = c("BF" = 1, "EIGEN" = 1, "user" = 1, "NONSIG" = 0.25), guide = "none") +
    scale_colour_manual(values = sig.colors, guide = "none") +  # no legend for points color
    scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
    scale_y_continuous(breaks = seq(0, round(max(for.plot.ann.IVM$log10p)), by = 2))+
    geom_hline(data = test.IVM, aes(yintercept = value, linetype = name)) + 
    scale_linetype_manual(values = c("BF" = 1, "EIGEN" = 3, "user" = 2), guide = "none") +  # no legend for linetypes
    labs(x =  NULL,
         y = expression(bold(-log[10](italic(p))))) +
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          plot.title = element_text(size = 10),
          axis.title.x = element_blank(),  # remove x axis title
          axis.text.x = element_blank(),  # remove x axis text
          axis.ticks.x = element_blank(),  # remove x axis ticks
          axis.title.y = element_text(size = 10),
          strip.text = element_text(size = 10, face="bold")
          #               strip.text.y = element_text(angle = 0, size = 12),
          #                strip.background = element_rect(fill = "grey", colour = "black")
    )+
    facet_grid(algorithm ~ CHROM, scales = "free_x", space = "free")
  
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
  
  
  
  # adjust mergins of the manhattan plot and the rug plot
  manhattan <- manhattan + theme(plot.margin = margin(5.5, 5.5, 0, 8, "pt"))
  rugPlot <- rugPlot + theme(plot.margin = margin(0, 5.5, 5.5, 5.5, "pt"))
  combined_plot_sized <- manhattan + rugPlot + patchwork::plot_layout(ncol = 1, heights = c(7, 1))
  plot_list[[ivm]] <- combined_plot_sized
  
 ggsave(paste0("Figures/singlePlots/",ivm,".png"), combined_plot_sized , width=7.5, height=3.5, units = "in", dpi = 600)
  
}    


# combine all plots into one figure
cowplots <- cowplot::plot_grid(plot_list$length_Ivermectin_12,plot_list$length_Ivermectin_08, align="V", ncol = 1, labels = c("A","B"))

ggsave(paste0("Figures/Ivermectin2021GWA.png"), cowplots , width=7.5, height=5, units = "in", dpi = 600)


