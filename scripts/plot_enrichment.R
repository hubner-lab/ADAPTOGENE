#!/usr/bin/env Rscript
# Generate enrichment visualization plots (emapplot, cnetplot, dotplot)

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(stringr)
    library(clusterProfiler)
    library(enrichplot)
    library(ggplot2)
    library(qs)
})

source("/pipeline/scripts/R/lib/enrichment_plots.R")

args = commandArgs(trailingOnly = TRUE)
################
INTERMEDIATE_DIR = args[1]
PLOTS_DIR        = args[2]
TOP_TERMS        = as.numeric(args[3])
PLOT_WIDTH       = as.numeric(args[4])
PLOT_HEIGHT      = as.numeric(args[5])
CNET_LABEL       = if (length(args) >= 6) args[6] else "gene_id"
GENES_FILE       = if (length(args) >= 7) args[7] else "NULL"
TOP_PLOT_REGIONS = if (length(args) >= 8) as.numeric(args[8]) else 0
REGIONS_FILE     = if (length(args) >= 9) args[9] else "NULL"
################

message('INFO: Generating enrichment plots')
message(paste0('INFO: .qs directory: ', INTERMEDIATE_DIR))
message(paste0('INFO: Plots directory: ', PLOTS_DIR))

label_map      <- build_label_map(GENES_FILE, CNET_LABEL)
top_region_ids <- top_regions_for_plotting(REGIONS_FILE, TOP_PLOT_REGIONS)

total <- process_all_enrichment(INTERMEDIATE_DIR, PLOTS_DIR,
                                 TOP_TERMS, PLOT_WIDTH, PLOT_HEIGHT,
                                 label_map      = label_map,
                                 top_region_ids = top_region_ids)

message(if (total > 0) paste0('INFO: Created ', total, ' plots')
        else 'WARNING: No plots created')
message('INFO: Complete')
