#!/usr/bin/env Rscript
# Write pipeline summary table for the current mode
# Appends/updates a running Pipeline_summary.tsv with metrics from the current step

library(data.table)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
################
MODE = args[1]           # Pipeline mode: processing, structure, structure_K, association, maladaptation
OUTPUT = args[2]         # Pipeline_summary.tsv path
# Remaining args are mode-specific input files
################

message(paste0('INFO: Writing summary for mode: ', MODE))

######################## Functions

# Read existing summary or create empty one
read_or_create_summary <- function(output_path) {
    if (file.exists(output_path) && file.size(output_path) > 0) {
        existing <- fread(output_path)
        message(paste0('INFO: Loaded existing summary with ', nrow(existing), ' rows'))
        return(existing)
    }
    return(data.table(step = character(), metric = character(), value = character()))
}

# Add rows to summary, replacing any existing rows for the same step
update_summary <- function(existing, new_rows) {
    # Remove old rows for this step
    step_name <- unique(new_rows$step)
    existing <- existing[!step %in% step_name]
    # Append new rows
    rbind(existing, new_rows)
}

# Create a row
row <- function(step, metric, value) {
    data.table(step = step, metric = metric, value = as.character(value))
}

#################### Main

# Read existing summary
summary_dt <- read_or_create_summary(OUTPUT)

if (MODE == 'processing') {
    # args: MODE OUTPUT vcf_filt vcf_ld samples_list samples_filtered samples_removed
    VCF_FILT = args[3]
    VCF_LD = args[4]
    SAMPLES_LIST = args[5]
    SAMPLES_FILTERED = args[6]
    SAMPLES_REMOVED = args[7]

    # Count samples
    n_samples_total <- length(readLines(SAMPLES_LIST))
    n_samples_filtered <- length(readLines(SAMPLES_FILTERED))
    n_samples_removed <- length(readLines(SAMPLES_REMOVED))

    # Count SNPs from VCF files (non-header lines)
    count_vcf_snps <- function(vcf_path) {
        con <- file(vcf_path, 'r')
        n <- 0L
        while (length(line <- readLines(con, n = 10000)) > 0) {
            n <- n + sum(!startsWith(line, '#'))
        }
        close(con)
        return(n)
    }

    n_snps_filtered <- count_vcf_snps(VCF_FILT)
    n_snps_ld <- count_vcf_snps(VCF_LD)

    new_rows <- rbind(
        row('processing', 'samples_total', n_samples_total),
        row('processing', 'samples_after_filtering', n_samples_filtered),
        row('processing', 'samples_removed', n_samples_removed),
        row('processing', 'snps_after_filtering', n_snps_filtered),
        row('processing', 'snps_after_ld_pruning', n_snps_ld)
    )

} else if (MODE == 'structure') {
    # args: MODE OUTPUT cross_entropy_plot K_START K_END
    K_START = args[3] %>% as.numeric
    K_END = args[4] %>% as.numeric

    new_rows <- rbind(
        row('structure', 'K_range', paste0(K_START, '-', K_END))
    )

} else if (MODE == 'structure_K') {
    # args: MODE OUTPUT K_BEST climate_site predictors
    K_BEST = args[3]
    CLIMATE_SITE = args[4]
    PREDICTORS = args[5]

    predictors_list <- str_split(PREDICTORS, ',')[[1]]

    # Count climate variables
    climate <- fread(CLIMATE_SITE)
    n_climate_vars <- ncol(climate) - 4  # Subtract site, sample, lat, lon columns

    new_rows <- rbind(
        row('structure_K', 'K_best', K_BEST),
        row('structure_K', 'climate_predictors', PREDICTORS),
        row('structure_K', 'n_climate_variables', n_climate_vars)
    )

} else if (MODE == 'association') {
    # args: MODE OUTPUT selected_snps regions_per_trait regions_combined genes_per_region
    #       enrichment_done [optional]
    SELECTED_SNPS = args[3]
    REGIONS_PER_TRAIT = args[4]
    REGIONS_COMBINED = args[5]
    GENES_PER_REGION = args[6]
    ENRICHMENT_DONE = if (length(args) >= 7) args[7] else "NULL"

    # Count selected SNPs
    snps <- fread(SELECTED_SNPS)
    n_selected_snps <- nrow(snps)

    # Count per-method SNPs
    method_cols <- setdiff(colnames(snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
    method_rows <- list()
    for (m in method_cols) {
        n_method <- sum(!is.na(snps[[m]]) & snps[[m]] != '' & snps[[m]] != '""')
        method_rows <- c(method_rows, list(row('association', paste0('sig_snps_', m), n_method)))
    }

    # Count regions
    regions_trait <- fread(REGIONS_PER_TRAIT)
    regions_combined <- fread(REGIONS_COMBINED)

    # Count genes
    genes <- fread(GENES_PER_REGION)
    n_genes <- if (nrow(genes) > 0 && 'gene_id' %in% colnames(genes)) n_distinct(genes$gene_id) else 0

    # Count enrichment
    enrichment_summary <- ''
    if (ENRICHMENT_DONE != "NULL" && file.exists(ENRICHMENT_DONE)) {
        enrichment_summary <- 'completed'
    }

    new_rows <- rbind(
        row('association', 'selected_snps_total', n_selected_snps),
        do.call(rbind, method_rows),
        row('association', 'regions_per_trait', nrow(regions_trait)),
        row('association', 'regions_combined', nrow(regions_combined)),
        row('association', 'genes_found', n_genes),
        row('association', 'enrichment_status', enrichment_summary)
    )

    # Add per-trait region counts
    if (nrow(regions_trait) > 0 && 'trait' %in% colnames(regions_trait)) {
        trait_counts <- regions_trait %>%
            dplyr::group_by(trait) %>%
            dplyr::summarise(n = dplyr::n(), .groups = 'drop')
        for (i in seq_len(nrow(trait_counts))) {
            new_rows <- rbind(new_rows,
                row('association', paste0('regions_', trait_counts$trait[i]), trait_counts$n[i]))
        }
    }

} else if (MODE == 'maladaptation') {
    # args: MODE OUTPUT gf_adaptive offset_site_values
    GF_ADAPTIVE = args[3]
    OFFSET_SITE = args[4]

    # Load offset values for range
    offset <- fread(OFFSET_SITE)
    offset_col <- setdiff(colnames(offset), c('site', 'sample', 'latitude', 'longitude'))
    if (length(offset_col) > 0) {
        offset_vals <- offset[[offset_col[1]]]
        offset_min <- round(min(offset_vals, na.rm = TRUE), 4)
        offset_max <- round(max(offset_vals, na.rm = TRUE), 4)
        offset_mean <- round(mean(offset_vals, na.rm = TRUE), 4)
    } else {
        offset_min <- offset_max <- offset_mean <- NA
    }

    new_rows <- rbind(
        row('maladaptation', 'gf_model', 'completed'),
        row('maladaptation', 'offset_min', offset_min),
        row('maladaptation', 'offset_max', offset_max),
        row('maladaptation', 'offset_mean', offset_mean)
    )

} else if (MODE == 'association_phenotypes') {
    # args: MODE OUTPUT missing_summary selected_snps regions_per_trait regions_combined genes_per_region [enrichment_done]
    MISSING_SUMMARY = args[3]
    SELECTED_SNPS = args[4]
    REGIONS_PER_TRAIT = args[5]
    REGIONS_COMBINED = args[6]
    GENES_PER_REGION = args[7]
    ENRICHMENT_DONE = if (length(args) >= 8) args[8] else "NULL"

    # Read missing value summary
    missing <- fread(MISSING_SUMMARY)
    n_traits <- nrow(missing)

    # Count selected SNPs
    snps <- fread(SELECTED_SNPS)
    n_selected_snps <- nrow(snps)

    # Count regions
    regions_trait <- fread(REGIONS_PER_TRAIT)
    regions_combined <- fread(REGIONS_COMBINED)

    # Count genes
    genes <- fread(GENES_PER_REGION)
    n_genes <- if (nrow(genes) > 0 && 'gene_id' %in% colnames(genes)) n_distinct(genes$gene_id) else 0

    # Enrichment status
    enrichment_summary <- ''
    if (ENRICHMENT_DONE != "NULL" && file.exists(ENRICHMENT_DONE)) {
        enrichment_summary <- 'completed'
    }

    new_rows <- rbind(
        row('association_phenotypes', 'n_phenotype_traits', n_traits),
        row('association_phenotypes', 'missing_strategy', missing$strategy[1]),
        row('association_phenotypes', 'selected_snps_total', n_selected_snps),
        row('association_phenotypes', 'regions_per_trait', nrow(regions_trait)),
        row('association_phenotypes', 'regions_combined', nrow(regions_combined)),
        row('association_phenotypes', 'genes_found', n_genes),
        row('association_phenotypes', 'enrichment_status', enrichment_summary)
    )

    # Per-trait missing info
    for (i in seq_len(nrow(missing))) {
        t <- missing$trait[i]
        new_rows <- rbind(new_rows,
            row('association_phenotypes', paste0(t, '_n_available'), missing$n_available[i]),
            row('association_phenotypes', paste0(t, '_n_missing'), missing$n_missing[i]))
        if (missing$n_missing[i] > 0 && nchar(missing$missing_samples[i]) > 0) {
            new_rows <- rbind(new_rows,
                row('association_phenotypes', paste0(t, '_dropped_samples'), missing$missing_samples[i]))
        }
    }

    # Per-trait region counts
    if (nrow(regions_trait) > 0 && 'trait' %in% colnames(regions_trait)) {
        trait_counts <- regions_trait %>%
            dplyr::group_by(trait) %>%
            dplyr::summarise(n = dplyr::n(), .groups = 'drop')
        for (i in seq_len(nrow(trait_counts))) {
            new_rows <- rbind(new_rows,
                row('association_phenotypes', paste0('regions_', trait_counts$trait[i]), trait_counts$n[i]))
        }
    }

} else if (MODE == 'overlapping') {
    # args: MODE OUTPUT overlap_summary selected_snps regions_per_trait regions_combined genes_per_region [enrichment_done]
    OVERLAP_SUMMARY = args[3]
    SELECTED_SNPS = args[4]
    REGIONS_PER_TRAIT = args[5]
    REGIONS_COMBINED = args[6]
    GENES_PER_REGION = args[7]
    ENRICHMENT_DONE = if (length(args) >= 8) args[8] else "NULL"

    # Read overlap summary (last rows are SUMMARY rows)
    overlap <- fread(OVERLAP_SUMMARY)
    summary_rows <- overlap[gea_region_id == 'SUMMARY']
    overlap_pairs <- overlap[gea_region_id != 'SUMMARY']

    n_gea <- as.integer(summary_rows[gwas_region_id == 'gea_total', overlap_length])
    n_gwas <- as.integer(summary_rows[gwas_region_id == 'gwas_total', overlap_length])
    n_gea_overlap <- as.integer(summary_rows[gwas_region_id == 'gea_overlapping', overlap_length])
    n_gwas_overlap <- as.integer(summary_rows[gwas_region_id == 'gwas_overlapping', overlap_length])
    n_pairs <- as.integer(summary_rows[gwas_region_id == 'overlap_pairs', overlap_length])

    gea_pct <- summary_rows[gwas_region_id == 'gea_total', overlap_pct]
    gwas_pct <- summary_rows[gwas_region_id == 'gwas_total', overlap_pct]

    # Count merged SNPs and regions
    snps <- fread(SELECTED_SNPS)
    regions_trait <- fread(REGIONS_PER_TRAIT)
    regions_combined <- fread(REGIONS_COMBINED)
    genes <- fread(GENES_PER_REGION)
    n_genes <- if (nrow(genes) > 0 && 'gene_id' %in% colnames(genes)) n_distinct(genes$gene_id) else 0

    enrichment_summary <- ''
    if (ENRICHMENT_DONE != "NULL" && file.exists(ENRICHMENT_DONE)) {
        enrichment_summary <- 'completed'
    }

    new_rows <- rbind(
        row('overlapping', 'gea_regions', n_gea),
        row('overlapping', 'gwas_regions', n_gwas),
        row('overlapping', 'overlap_pairs', n_pairs),
        row('overlapping', 'gea_overlapping', n_gea_overlap),
        row('overlapping', 'gwas_overlapping', n_gwas_overlap),
        row('overlapping', 'gea_overlap_pct', gea_pct),
        row('overlapping', 'gwas_overlap_pct', gwas_pct),
        row('overlapping', 'merged_snps_total', nrow(snps)),
        row('overlapping', 'new_regions_per_trait', nrow(regions_trait)),
        row('overlapping', 'new_regions_combined', nrow(regions_combined)),
        row('overlapping', 'genes_found', n_genes),
        row('overlapping', 'enrichment_status', enrichment_summary)
    )

} else {
    message(paste0('WARNING: Unknown mode "', MODE, '", skipping summary'))
    quit(save = "no", status = 0)
}

# Update summary
summary_dt <- update_summary(summary_dt, new_rows)

# Save
fwrite(summary_dt, OUTPUT, sep = '\t')
message(paste0('INFO: Saved summary with ', nrow(summary_dt), ' rows to ', OUTPUT))
message('INFO: Complete')
