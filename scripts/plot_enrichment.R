#!/usr/bin/env Rscript
# Generate enrichment visualization plots (emapplot, cnetplot, dotplot)
# Creates plots for each region's enrichment results

library(data.table)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(qs)

args = commandArgs(trailingOnly=TRUE)
################
INTERMEDIATE_DIR = args[1]   # Base intermediate directory with .qs files: intermediate/enrichment/
PLOTS_DIR = args[2]          # Base plots directory: plots/association/enrichment/
TOP_TERMS = args[3] %>% as.numeric   # Number of top terms to show
PLOT_WIDTH = args[4] %>% as.numeric  # Plot width in inches
PLOT_HEIGHT = args[5] %>% as.numeric # Plot height in inches
################

message('INFO: Generating enrichment visualization plots')
message(paste0('INFO: Enrichment .qs directory: ', INTERMEDIATE_DIR))
message(paste0('INFO: Plots directory: ', PLOTS_DIR))
message(paste0('INFO: Top terms to show: ', TOP_TERMS))

# Check if enrichment directory exists
if (!dir.exists(INTERMEDIATE_DIR)) {
    message('WARNING: Enrichment .qs directory does not exist, skipping plotting')
    quit(save = "no", status = 0)
}

######################## Functions

# Create emapplot (GO term network)
create_emapplot <- function(enrich_obj, n_terms, output_base) {
    tryCatch({
        # Calculate pairwise term similarity FIRST (required for emapplot)
        # Use enrichplot::pairwise_termsim to add similarity matrix
        enrich_obj <- enrichplot::pairwise_termsim(enrich_obj)

        # emapplot with showCategory to limit number of terms
        p <- enrichplot::emapplot(enrich_obj, showCategory = n_terms) +
            ggtitle("GO Term Similarity Network") +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

        # Save PNG
        ggsave(paste0(output_base, '.png'), plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)
        message(paste0('INFO: Saved emapplot PNG to ', output_base, '.png'))

        # Save SVG
        ggsave(paste0(output_base, '.svg'), plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
        message(paste0('INFO: Saved emapplot SVG to ', output_base, '.svg'))

    }, error = function(e) {
        message(paste0('WARNING: Could not create emapplot: ', e$message))
        message('INFO: This may occur if there are too few terms or terms are not related')
    })
}

# Create cnetplot (gene-concept network)
create_cnetplot <- function(enrich_obj, n_terms, output_base) {
    tryCatch({
        # cnetplot without color.params (simplified)
        # Just show top N categories by default
        p <- enrichplot::cnetplot(enrich_obj, showCategory = n_terms) +
            ggtitle("Gene-Concept Network") +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

        # Save PNG
        ggsave(paste0(output_base, '.png'), plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)
        message(paste0('INFO: Saved cnetplot PNG to ', output_base, '.png'))

        # Save SVG
        ggsave(paste0(output_base, '.svg'), plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
        message(paste0('INFO: Saved cnetplot SVG to ', output_base, '.svg'))

    }, error = function(e) {
        message(paste0('WARNING: Could not create cnetplot: ', e$message))
    })
}

# Create dotplot (enrichment dotplot)
create_dotplot <- function(enrich_obj, n_terms, output_base) {
    tryCatch({
        p <- enrichplot::dotplot(enrich_obj, showCategory = n_terms) +
            ggtitle("GO Enrichment") +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

        # Save PNG
        ggsave(paste0(output_base, '.png'), plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 300)
        message(paste0('INFO: Saved dotplot PNG to ', output_base, '.png'))

        # Save SVG
        ggsave(paste0(output_base, '.svg'), plot = p, width = PLOT_WIDTH, height = PLOT_HEIGHT)
        message(paste0('INFO: Saved dotplot SVG to ', output_base, '.svg'))

    }, error = function(e) {
        message(paste0('WARNING: Could not create dotplot: ', e$message))
    })
}

# Process all enrichment files in a directory
process_enrichment_files <- function(enrichment_dir, plots_dir) {
    # Find all per-region enrichment .qs files
    region_files <- list.files(enrichment_dir,
                                pattern = "^Region_.*_enrichment\\.qs$",
                                full.names = TRUE,
                                recursive = TRUE)

    if (length(region_files) == 0) {
        message('WARNING: No per-region enrichment .qs files found')
        return(0)
    }

    message(paste0('INFO: Found ', length(region_files), ' region enrichment files'))

    plots_created <- 0

    for (file in region_files) {
        message(paste0('INFO: Processing ', basename(file)))

        # Extract trait and region_id from file path
        # Path format: enrichment/bio_1/Region_3:15669175-15669175_bio_1_enrichment.qs
        rel_path <- sub(paste0(enrichment_dir, '/?'), '', file)
        trait <- dirname(rel_path)
        filename <- basename(file)
        region_id <- sub('Region_', '', sub('_enrichment\\.qs$', '', filename))

        # Create output directory for this trait
        trait_plot_dir <- file.path(plots_dir, trait)
        if (!dir.exists(trait_plot_dir)) {
            dir.create(trait_plot_dir, recursive = TRUE)
            message(paste0('INFO: Created directory: ', trait_plot_dir))
        }

        # Load enrichResult object from .qs file
        enrich_obj <- qs::qread(file)

        if (is.null(enrich_obj) || !inherits(enrich_obj, "enrichResult")) {
            message(paste0('WARNING: Invalid enrichment object in ', basename(file)))
            next
        }

        n_terms <- nrow(enrich_obj@result)
        message(paste0('INFO: Region ', region_id, ': ', n_terms, ' enriched terms'))

        # Use minimum of TOP_TERMS and actual number of terms
        n_plot <- min(TOP_TERMS, n_terms)

        # Generate plots (save both PNG and SVG)
        base_name <- paste0('Region_', region_id)

        # Dotplot (always works)
        dotplot_base <- file.path(trait_plot_dir, paste0(base_name, '_dotplot'))
        create_dotplot(enrich_obj, n_plot, dotplot_base)
        plots_created <- plots_created + 2  # PNG + SVG

        # Cnetplot (needs gene information)
        cnetplot_base <- file.path(trait_plot_dir, paste0(base_name, '_cnetplot'))
        create_cnetplot(enrich_obj, n_plot, cnetplot_base)
        plots_created <- plots_created + 2  # PNG + SVG

        # Emapplot (needs multiple terms and similarity calculation)
        if (n_terms >= 2) {
            emapplot_base <- file.path(trait_plot_dir, paste0(base_name, '_emapplot'))
            create_emapplot(enrich_obj, n_plot, emapplot_base)
            plots_created <- plots_created + 2  # PNG + SVG
        } else {
            message('INFO: Skipping emapplot (need at least 2 terms)')
        }
    }

    return(plots_created)
}

#################### Main

# Process all enrichment results
total_plots <- process_enrichment_files(INTERMEDIATE_DIR, PLOTS_DIR)

if (total_plots > 0) {
    message(paste0('INFO: Created ', total_plots, ' plots'))
} else {
    message('WARNING: No plots were created')
}

message('INFO: Complete')
