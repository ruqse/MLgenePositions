# Vsualisation of C. elegans mappings after Ivermectin drug exposure

## Overview
This script provides a visualization of the mappings for C. elegans genes implicated in Ivermectin drug exposure. The data is sourced from the [**AndersenLab**](https://andersenlab.org/).

The *Figures.R* script fetches genomic coordinates of genes implicated in Ivermectin

## Requirements
The following R packages are required to run the script:

tidyverse
data.table
biomaRt
glue
ggplot2

## Instructions
Make sure all the required libraries are installed. The script contains commands for installing the required libraries.
Update the gene_names variable with the genes you are interested in visualizing.
Ensure you have the chromosome lengths file (data/c_elegans_chr_lengths.tsv), the mappings data (data/filtered.mappings.results.rda), and the independent test cutoffs (data/independent_test_cutoff.rda).
Run the Figures.R script.

## Output
The script will produce a series of Manhattan plots representing the genomic mappings of genes of interest. Genes are categorized into classes such as "GluCl", "GABA", "Efflux pump", and "Dyf". The plots also highlight significant genomic positions based on different significance thresholds.

## Notes
The script sets the working directory to the parent directory of the active RStudio document. Make sure the script is placed appropriately with respect to the data files.
Genomic data is fetched from the wormbase database using the biomaRt package.

Hope this helps! If you have additional features or functionalities that you'd like to highlight in the README, let me know!
