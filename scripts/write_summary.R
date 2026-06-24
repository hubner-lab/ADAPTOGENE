#!/usr/bin/env Rscript
# Write pipeline summary table for the current mode
# Appends/updates a running Pipeline_summary.tsv with metrics from the current step

library(data.table)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
################
MODE = args[1]           # Pipeline mode: processing, prestructure, structure, gea, gwas, gea_x_gwas, maladaptation
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
    #       qc_raw_summary filtering_summary het_outlier_sd has_dp [depth_sample depth_site]
    VCF_FILT          = args[3]
    VCF_LD            = args[4]
    SAMPLES_LIST      = args[5]
    SAMPLES_FILTERED  = args[6]
    SAMPLES_REMOVED   = args[7]
    QC_RAW_SUMMARY    = args[8]
    FILTERING_SUMMARY = args[9]
    HET_OUTLIER_SD    = args[10]
    HAS_DP            = args[11]
    DEPTH_SUMMARY     = if (length(args) >= 12) args[12] else 'NULL'

    # Count samples
    n_samples_total    <- length(readLines(SAMPLES_LIST))
    n_samples_filtered <- length(readLines(SAMPLES_FILTERED))
    n_samples_removed  <- length(readLines(SAMPLES_REMOVED))

    # Count SNPs from VCF files (non-header lines)
    count_vcf_snps <- function(vcf_path) {
        con <- file(vcf_path, 'r')
        n <- 0L
        while (length(line <- readLines(con, n = 10000)) > 0)
            n <- n + sum(!startsWith(line, '#'))
        close(con)
        return(n)
    }
    n_snps_filtered <- count_vcf_snps(VCF_FILT)
    n_snps_ld       <- count_vcf_snps(VCF_LD)

    # Parse raw bcftools stats output (SN/TSTV lines)
    parse_bcftools_stats <- function(path) {
        if (!file.exists(path) || file.size(path) == 0) return(list(snps=NA, titv=NA))
        lines <- readLines(path)
        snps <- NA; titv <- NA
        for (ln in lines) {
            if (startsWith(ln, 'SN')) {
                parts <- strsplit(ln, '\t')[[1]]
                if (length(parts) >= 4 && grepl('number of SNPs', parts[3]))
                    snps <- trimws(parts[4])
            } else if (startsWith(ln, 'TSTV')) {
                parts <- strsplit(ln, '\t')[[1]]
                if (length(parts) >= 5) titv <- trimws(parts[5])
            }
        }
        list(snps=snps, titv=titv)
    }
    bcf <- parse_bcftools_stats(QC_RAW_SUMMARY)
    n_raw_snps <- bcf$snps
    titv       <- bcf$titv

    # Read filtering summary for detailed attrition
    filt_sum <- fread(FILTERING_SUMMARY)

    # Het outlier count from filtering summary
    n_het_removed <- 0L
    if (HET_OUTLIER_SD != 'NULL' && 'After het outlier removal' %in% filt_sum$stage) {
        n_before <- filt_sum[stage == 'After sample missingness filter', n_samples]
        n_after  <- filt_sum[stage == 'After het outlier removal', n_samples]
        if (length(n_before) > 0 && length(n_after) > 0)
            n_het_removed <- as.integer(n_before) - as.integer(n_after)
    }

    # Depth info
    depth_qc_status <- if (HAS_DP == 'TRUE') 'available' else 'not_provided'

    new_rows <- rbind(
        row('processing', 'samples_total',              n_samples_total),
        row('processing', 'samples_after_filtering',    n_samples_filtered),
        row('processing', 'samples_removed',            n_samples_removed),
        row('processing', 'samples_het_outliers_removed', n_het_removed),
        row('processing', 'snps_raw',                   n_raw_snps),
        row('processing', 'snps_after_filtering',       n_snps_filtered),
        row('processing', 'snps_after_ld_pruning',      n_snps_ld),
        row('processing', 'titv_ratio',                 titv),
        row('processing', 'depth_qc',                   depth_qc_status)
    )

    # Append depth stats if available
    if (HAS_DP == 'TRUE' && DEPTH_SUMMARY != 'NULL' && file.exists(DEPTH_SUMMARY)) {
        depth_dt <- fread(DEPTH_SUMMARY)
        for (i in seq_len(nrow(depth_dt))) {
            new_rows <- rbind(new_rows, row('processing', depth_dt$metric[i], depth_dt$value[i]))
        }
    }

} else if (MODE == 'prestructure') {
    # args: MODE OUTPUT cross_entropy_plot K_START K_END
    K_START = args[3] %>% as.numeric
    K_END = args[4] %>% as.numeric

    new_rows <- rbind(
        row('prestructure',  'K_range', paste0(K_START, '-', K_END))
    )

} else if (MODE == 'structure') {
    # args: MODE OUTPUT K_BEST climate_site predictors ld_decay_path ld_decay_group_by ld_decay_scope
    K_BEST = args[3]
    CLIMATE_SITE = args[4]
    PREDICTORS = args[5]
    LD_DECAY_PATH = if (length(args) >= 6) args[6] else "NULL"
    LD_DECAY_GROUP_BY = if (length(args) >= 7) args[7] else "NULL"
    LD_DECAY_SCOPE_VAL = if (length(args) >= 8) args[8] else "NULL"

    new_rows <- rbind(row('structure', 'K_best', K_BEST))

    if (CLIMATE_SITE != 'NULL') {
        predictors_list <- str_split(PREDICTORS, ',')[[1]]
        climate <- fread(CLIMATE_SITE)
        n_climate_vars <- ncol(climate) - 4  # Subtract site, sample, lat, lon columns

        new_rows <- rbind(new_rows,
            row('structure', 'climate_predictors', PREDICTORS),
            row('structure', 'n_climate_variables', n_climate_vars)
        )
    }

    # LD decay summary
    if (LD_DECAY_PATH != 'NULL' && file.exists(LD_DECAY_PATH)) {
        ld_table <- fread(LD_DECAY_PATH)
        gw_all <- ld_table[group == 'All' & scope == 'genome_wide']
        if (nrow(gw_all) > 0) {
            new_rows <- rbind(new_rows,
                row('structure', 'ld_decay_half_decay_genome_wide', gw_all$half_decay_bp[1]),
                row('structure', 'ld_decay_r2_02_genome_wide', gw_all$r2_02_bp[1])
            )
        }
        new_rows <- rbind(new_rows,
            row('structure', 'ld_decay_group_by', LD_DECAY_GROUP_BY),
            row('structure', 'ld_decay_groups_analyzed', nrow(ld_table[scope == 'genome_wide'])),
            row('structure', 'ld_decay_scope', LD_DECAY_SCOPE_VAL)
        )
    }

} else if (MODE %in% c('gea', 'gwas')) {
    # args: MODE OUTPUT selected_snps regions_per_trait regions_combined genes_per_region
    #        [region_dist_spec r2_threshold ld_decay_group ld_decay_path]
    is_gwas <- MODE == 'gwas'
    arg_offset <- if (is_gwas) 1L else 0L  # gwas has an extra missing_summary arg before selected_snps
    MISSING_SUMMARY   <- if (is_gwas) args[3] else NULL
    SELECTED_SNPS     <- args[3 + arg_offset]
    REGIONS_PER_TRAIT <- args[4 + arg_offset]
    REGIONS_COMBINED  <- args[5 + arg_offset]
    GENES_PER_REGION  <- args[6 + arg_offset]
    REGION_DIST_SPEC  <- if (length(args) >= 7 + arg_offset) args[7 + arg_offset] else 'unknown'
    R2_THRESHOLD_VAL  <- if (length(args) >= 8 + arg_offset) as.numeric(args[8 + arg_offset]) else 0.2
    LD_DECAY_GROUP_VAL <- if (length(args) >= 9 + arg_offset) args[9 + arg_offset] else 'All'
    LD_DECAY_PATH_VAL <- if (length(args) >= 10 + arg_offset) args[10 + arg_offset] else 'NULL'

    # Count selected SNPs
    snps <- fread(SELECTED_SNPS)
    n_selected_snps <- nrow(snps)

    # Count per-method SNPs
    method_cols <- setdiff(colnames(snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
    method_rows <- list()
    for (m in method_cols) {
        n_method <- sum(!is.na(snps[[m]]) & snps[[m]] != '' & snps[[m]] != '""')
        method_rows <- c(method_rows, list(row(MODE, paste0('sig_snps_', m), n_method)))
    }

    # Count regions
    regions_trait    <- fread(REGIONS_PER_TRAIT)
    regions_combined <- fread(REGIONS_COMBINED)

    # Count genes
    genes   <- fread(GENES_PER_REGION)
    n_genes <- if (nrow(genes) > 0 && 'gene_id' %in% colnames(genes)) n_distinct(genes$gene_id) else 0

    new_rows <- rbind(
        row(MODE, 'selected_snps_total', n_selected_snps),
        do.call(rbind, method_rows),
        row(MODE, 'regions_per_trait', nrow(regions_trait)),
        row(MODE, 'regions_combined',  nrow(regions_combined)),
        row(MODE, 'genes_found',       n_genes),
        row(MODE, 'region_distance_spec',    REGION_DIST_SPEC),
        row(MODE, 'region_r2_threshold',     R2_THRESHOLD_VAL),
        row(MODE, 'region_ld_decay_group',   LD_DECAY_GROUP_VAL)
    )

    # Per-chromosome distances from LD decay table (when auto mode)
    if (LD_DECAY_PATH_VAL != 'NULL' && file.exists(LD_DECAY_PATH_VAL) &&
        REGION_DIST_SPEC %in% c('auto_per_chromosome', 'auto_genome_wide', 'auto')) {
        ld_table <- fread(LD_DECAY_PATH_VAL)
        if (REGION_DIST_SPEC == 'auto_per_chromosome') {
            chr_rows <- ld_table[group == LD_DECAY_GROUP_VAL & !scope %in% c('genome_wide')]
            if (nrow(chr_rows) > 0) {
                new_rows <- rbind(new_rows,
                    row(MODE, 'region_distance_chr_min_bp',
                        round(min(chr_rows$r2_02_bp, na.rm = TRUE))),
                    row(MODE, 'region_distance_chr_median_bp',
                        round(median(chr_rows$r2_02_bp, na.rm = TRUE))),
                    row(MODE, 'region_distance_chr_max_bp',
                        round(max(chr_rows$r2_02_bp, na.rm = TRUE)))
                )
            }
        } else {
            gw_row <- ld_table[group == LD_DECAY_GROUP_VAL & scope == 'genome_wide']
            if (nrow(gw_row) > 0) {
                new_rows <- rbind(new_rows,
                    row(MODE, 'region_distance_genome_wide_bp', round(gw_row$r2_02_bp[1])))
            }
        }
    }

    # Add per-trait region counts
    if (nrow(regions_trait) > 0 && 'trait' %in% colnames(regions_trait)) {
        trait_counts <- regions_trait %>%
            dplyr::group_by(trait) %>%
            dplyr::summarise(n = dplyr::n(), .groups = 'drop')
        for (i in seq_len(nrow(trait_counts))) {
            new_rows <- rbind(new_rows,
                row(MODE, paste0('regions_', trait_counts$trait[i]), trait_counts$n[i]))
        }
    }

    # GWAS-specific: missing value summary
    if (is_gwas && !is.null(MISSING_SUMMARY) && file.exists(MISSING_SUMMARY)) {
        missing <- fread(MISSING_SUMMARY)
        new_rows <- rbind(new_rows,
            row('gwas', 'n_phenotype_traits', nrow(missing)),
            row('gwas', 'missing_strategy',   missing$strategy[1])
        )
        for (i in seq_len(nrow(missing))) {
            t <- missing$trait[i]
            new_rows <- rbind(new_rows,
                row('gwas', paste0(t, '_n_available'), missing$n_available[i]),
                row('gwas', paste0(t, '_n_missing'),   missing$n_missing[i]))
            if (missing$n_missing[i] > 0 && nchar(missing$missing_samples[i]) > 0) {
                new_rows <- rbind(new_rows,
                    row('gwas', paste0(t, '_dropped_samples'), missing$missing_samples[i]))
            }
        }
    }

} else if (MODE == 'maladaptation') {
    # args: MODE OUTPUT offset_site_values (space-joined paths, one per snp_set x spatial_tag)
    OFFSET_SITE = args[3]

    # Split space-joined paths into individual paths
    offset_paths <- strsplit(trimws(OFFSET_SITE), "\\s+")[[1]]
    offset_paths <- offset_paths[nchar(offset_paths) > 0]

    # Emit per-tag offset stats (one row per SNP set x spatial tag)
    tag_rows <- lapply(offset_paths, function(p) {
        tag <- basename(dirname(p))
        if (!file.exists(p)) {
            return(list(
                row('maladaptation', sprintf('offset_min_%s', tag), NA),
                row('maladaptation', sprintf('offset_max_%s', tag), NA),
                row('maladaptation', sprintf('offset_mean_%s', tag), NA)
            ))
        }
        offset <- fread(p, colClasses = c(site = 'character', sample = 'character'))
        offset_col <- setdiff(colnames(offset), c('site', 'sample', 'latitude', 'longitude'))
        if (length(offset_col) > 0) {
            vals <- offset[[offset_col[1]]]
            list(
                row('maladaptation', sprintf('offset_min_%s', tag),  round(min(vals,  na.rm = TRUE), 4)),
                row('maladaptation', sprintf('offset_max_%s', tag),  round(max(vals,  na.rm = TRUE), 4)),
                row('maladaptation', sprintf('offset_mean_%s', tag), round(mean(vals, na.rm = TRUE), 4))
            )
        } else {
            list(
                row('maladaptation', sprintf('offset_min_%s', tag),  NA),
                row('maladaptation', sprintf('offset_max_%s', tag),  NA),
                row('maladaptation', sprintf('offset_mean_%s', tag), NA)
            )
        }
    })

    new_rows <- rbind(
        row('maladaptation', 'gf_model', 'completed'),
        rbindlist(lapply(unlist(tag_rows, recursive = FALSE), as.data.table))
    )

} else if (MODE == 'gea_x_gwas') {
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
        row('gea_x_gwas', 'gea_regions', n_gea),
        row('gea_x_gwas', 'gwas_regions', n_gwas),
        row('gea_x_gwas', 'overlap_pairs', n_pairs),
        row('gea_x_gwas', 'gea_overlapping', n_gea_overlap),
        row('gea_x_gwas', 'gwas_overlapping', n_gwas_overlap),
        row('gea_x_gwas', 'gea_overlap_pct', gea_pct),
        row('gea_x_gwas', 'gwas_overlap_pct', gwas_pct),
        row('gea_x_gwas', 'merged_snps_total', nrow(snps)),
        row('gea_x_gwas', 'new_regions_per_trait', nrow(regions_trait)),
        row('gea_x_gwas', 'new_regions_combined', nrow(regions_combined)),
        row('gea_x_gwas', 'genes_found', n_genes),
        row('gea_x_gwas', 'enrichment_status', enrichment_summary)
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
