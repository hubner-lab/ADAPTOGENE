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

} else if (MODE == 'structure') {
    # args: MODE OUTPUT cross_entropy_plot K_START K_END
    K_START = args[3] %>% as.numeric
    K_END = args[4] %>% as.numeric

    new_rows <- rbind(
        row('structure', 'K_range', paste0(K_START, '-', K_END))
    )

} else if (MODE == 'structure_K') {
    # args: MODE OUTPUT K_BEST climate_site predictors ld_decay_path ld_decay_group_by ld_decay_scope
    K_BEST = args[3]
    CLIMATE_SITE = args[4]
    PREDICTORS = args[5]
    LD_DECAY_PATH = if (length(args) >= 6) args[6] else "NULL"
    LD_DECAY_GROUP_BY = if (length(args) >= 7) args[7] else "NULL"
    LD_DECAY_SCOPE_VAL = if (length(args) >= 8) args[8] else "NULL"

    new_rows <- rbind(row('structure_K', 'K_best', K_BEST))

    if (CLIMATE_SITE != 'NULL') {
        predictors_list <- str_split(PREDICTORS, ',')[[1]]
        climate <- fread(CLIMATE_SITE)
        n_climate_vars <- ncol(climate) - 4  # Subtract site, sample, lat, lon columns

        new_rows <- rbind(new_rows,
            row('structure_K', 'climate_predictors', PREDICTORS),
            row('structure_K', 'n_climate_variables', n_climate_vars)
        )
    }

    # LD decay summary
    if (LD_DECAY_PATH != 'NULL' && file.exists(LD_DECAY_PATH)) {
        ld_table <- fread(LD_DECAY_PATH)
        gw_all <- ld_table[group == 'All' & scope == 'genome_wide']
        if (nrow(gw_all) > 0) {
            new_rows <- rbind(new_rows,
                row('structure_K', 'ld_decay_half_decay_genome_wide', gw_all$half_decay_bp[1]),
                row('structure_K', 'ld_decay_r2_02_genome_wide', gw_all$r2_02_bp[1])
            )
        }
        new_rows <- rbind(new_rows,
            row('structure_K', 'ld_decay_group_by', LD_DECAY_GROUP_BY),
            row('structure_K', 'ld_decay_groups_analyzed', nrow(ld_table[scope == 'genome_wide'])),
            row('structure_K', 'ld_decay_scope', LD_DECAY_SCOPE_VAL)
        )
    }

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

} else if (MODE == 'haplotype_scan') {
    TAGS <- args[3]
    tags_list <- strsplit(TAGS, ",")[[1]]

    new_rows <- rbind(
        row('haplotype_scan', 'tags_analyzed', paste(tags_list, collapse = "; ")),
        row('haplotype_scan', 'n_tag_combos', length(tags_list))
    )

    # Read scan_status.tsv for each tag and extract comprehensive stats
    for (tag in tags_list) {
        # Try to read scan status file
        inter_base <- dirname(dirname(OUTPUT))  # Go up from tables/ to results/
        status_file <- file.path(inter_base, 'intermediate', paste0('haplotype_', tag), 'scan_status.tsv')

        if (file.exists(status_file)) {
            status <- fread(status_file)
            summary_row <- status[region_id == 'SUMMARY']
            detail <- status[region_id != 'SUMMARY']

            if (nrow(summary_row) > 0) {
                n_total <- summary_row$start[1]
                n_success <- summary_row$end[1]
                n_skipped <- summary_row$snp_count[1]
                n_failed <- as.integer(summary_row$status[1])
            } else {
                n_total <- nrow(detail)
                n_success <- sum(detail$status == 'success')
                n_skipped <- sum(detail$status == 'skipped_low_snps')
                n_failed <- sum(grepl('failed', detail$status))
            }

            # Average SNPs per successful region
            avg_snps <- if (n_success > 0) {
                mean(detail[status == 'success']$snp_count, na.rm = TRUE)
            } else { 0 }

            # Most common failure reason
            failure_counts <- detail[status != 'success', .N, by = failure_reason][order(-N)]
            top_failure <- if (nrow(failure_counts) > 0) {
                substr(failure_counts$failure_reason[1], 1, 100)
            } else { '' }

            new_rows <- rbind(new_rows,
                row(paste0('haplotype_scan_', tag), 'regions_attempted', n_total),
                row(paste0('haplotype_scan_', tag), 'regions_succeeded', n_success),
                row(paste0('haplotype_scan_', tag), 'regions_skipped_low_snps', n_skipped),
                row(paste0('haplotype_scan_', tag), 'regions_failed', n_failed),
                row(paste0('haplotype_scan_', tag), 'avg_snps_per_success_region', sprintf('%.1f', avg_snps))
            )

            if (top_failure != '') {
                new_rows <- rbind(new_rows,
                    row(paste0('haplotype_scan_', tag), 'most_common_failure', top_failure)
                )
            }
        } else {
            # Fallback: just count regions file if scan_status.tsv doesn't exist
            regions_file <- file.path(dirname(OUTPUT), paste0("haplotype_", tag, "/Selected_regions.tsv"))
            if (file.exists(regions_file)) {
                reg <- fread(regions_file)
                new_rows <- rbind(new_rows,
                    row(paste0('haplotype_scan_', tag), 'regions_selected', nrow(reg))
                )
            }
        }
    }

} else if (MODE == 'haplotype') {
    TAGS <- args[3]
    EPSILON <- args[4]
    tags_list <- strsplit(TAGS, ",")[[1]]

    new_rows <- rbind(
        row('haplotype', 'epsilon_selected', EPSILON),
        row('haplotype', 'tags_analyzed', paste(tags_list, collapse = "; "))
    )

    for (tag in tags_list) {
        assign_file <- file.path(dirname(OUTPUT), paste0("haplotype_", tag, "/Haplotype_assignments.tsv"))
        if (file.exists(assign_file)) {
            assigns <- fread(assign_file)
            n_haps <- if ('haplotype' %in% names(assigns)) n_distinct(assigns$haplotype) else 0
            n_regions <- if ('region_id' %in% names(assigns)) n_distinct(assigns$region_id) else 0
            new_rows <- rbind(new_rows,
                row('haplotype', paste0('haplotype_groups_', tag), n_haps),
                row('haplotype', paste0('regions_with_results_', tag), n_regions)
            )
        }
    }

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
