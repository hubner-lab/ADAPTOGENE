#!/usr/bin/env Rscript
# Custom Manhattan plot for association analysis
# Generates clean, customizable Manhattan plots using ggplot2
# Optionally highlights significant regions with method-specific shapes

library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(qvalue)
library(scales)

args = commandArgs(trailingOnly=TRUE)
################################
ASSOC_TABLE = args[1]   # SNPID chr pos TRAITS...
ADJUST = args[2]        # e.g., "bonf_0.05" or "qval_0.05"
Kbest = args[3] %>% as.numeric
METHOD = args[4]        # EMMAX or LFMM (current plot method)
TRAIT = args[5]         # single trait name, e.g., "bio_1"
PLOT_DIR = args[6]      # output directory
REGIONS_FILE = args[7]  # optional: Regions.tsv for highlighting (can be "NULL")
SIGSNPS_FILES = args[8] # comma-separated list of sigSNPs files from all methods
################################

# Handle optional files
if (is.na(REGIONS_FILE) || REGIONS_FILE == "NULL" || REGIONS_FILE == "") {
    REGIONS_FILE <- NULL
}
if (is.na(SIGSNPS_FILES) || SIGSNPS_FILES == "NULL" || SIGSNPS_FILES == "") {
    SIGSNPS_FILES <- NULL
}

message(paste0('INFO: Manhattan plot for ', METHOD, ' - ', TRAIT))
message(paste0('INFO: K = ', Kbest))
message(paste0('INFO: Adjustment: ', ADJUST))

################################ Functions

FUN_max_pvalue_fdr <- function(pvalues, pval_threshold) {
    qvalues_result <- qvalue(pvalues)
    significant_pvalues <- pvalues[qvalues_result$qvalues < pval_threshold]
    if (length(significant_pvalues) > 0) {
        return(max(significant_pvalues))
    } else {
        return(NULL)
    }
}

FUN_max_pvalue_top <- function(pvalues, topN) {
    if (length(pvalues) < topN) {
        stop("topN is larger than the number of available p-values")
    }
    sorted_pvalues <- sort(pvalues, decreasing = FALSE)
    return(max(sorted_pvalues[1:topN]))
}

# Create cumulative position for Manhattan plot
prepare_manhattan_data <- function(df, chr_col = "chr", pos_col = "pos", pval_col) {
    # Get chromosome order (natural sort)
    chr_order <- df[[chr_col]] %>% unique() %>%
        .[order(as.numeric(str_extract(., "\\d+")))]

    df <- df %>%
        dplyr::mutate(chr_f = factor(.data[[chr_col]], levels = chr_order))

    # Calculate cumulative positions
    chr_lengths <- df %>%
        group_by(chr_f) %>%
        summarise(chr_len = max(.data[[pos_col]]), .groups = 'drop')

    chr_lengths <- chr_lengths %>%
        dplyr::mutate(
            tot = cumsum(as.numeric(chr_len)) - chr_len,
            center = tot + chr_len / 2
        )

    df <- df %>%
        left_join(chr_lengths %>% dplyr::select(chr_f, tot), by = "chr_f") %>%
        dplyr::mutate(
            pos_cum = .data[[pos_col]] + tot,
            log10p = -log10(.data[[pval_col]])
        )

    list(data = df, chr_info = chr_lengths)
}

# Manhattan plot theme
theme_manhattan <- function() {
    theme_minimal() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(color = "grey30", linewidth = 0.3),
        axis.ticks.y = element_line(color = "grey30", linewidth = 0.3),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        legend.position = "none",
        plot.margin = margin(10, 15, 10, 10)
    )
}

# Generate color palette for chromosomes
get_chr_colors <- function(n_chr) {
    # Alternating colors: dark blue and light blue-grey
    rep(c("#2166AC", "#92C5DE"), length.out = n_chr)
}

# Generate color palette for highlighted regions
get_region_colors <- function(n_regions) {
    if (n_regions <= 10) {
        # Distinct colors for up to 10 regions
        colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                   "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5")
        return(colors[1:n_regions])
    } else {
        return(scales::hue_pal()(n_regions))
    }
}

################################ Main

# Load association data
message('INFO: Loading association data')
snps_assoc <- fread(ASSOC_TABLE, sep = '\t', header = T)
snps_assoc$chr <- as.character(snps_assoc$chr)
message(paste0('INFO: Loaded ', nrow(snps_assoc), ' SNPs'))

# Check if trait exists
if (!(TRAIT %in% colnames(snps_assoc))) {
    stop(paste0("Trait '", TRAIT, "' not found in association table"))
}

# Parse adjustment parameters
adjustment <- ADJUST %>% str_split('_') %>% unlist %>% .[1]
pval_threshold <- ADJUST %>% str_split('_') %>% unlist %>% .[2] %>% as.numeric

message(paste0('INFO: Adjustment method: ', adjustment))
message(paste0('INFO: Threshold: ', pval_threshold))

# Calculate threshold based on adjustment method
original_threshold <- pval_threshold
if (adjustment == 'bonf') {
    pval_threshold <- original_threshold / nrow(snps_assoc)
}
if (adjustment == 'qval') {
    pval_threshold <- FUN_max_pvalue_fdr(snps_assoc[[TRAIT]], original_threshold)
    if (is.null(pval_threshold)) {
        message('WARNING: No significant SNPs at q-value threshold')
        pval_threshold <- original_threshold / nrow(snps_assoc)  # fallback to Bonferroni
    }
}
if (adjustment == 'top') {
    pval_threshold <- FUN_max_pvalue_top(snps_assoc[[TRAIT]], original_threshold)
}

message(paste0('INFO: Final p-value threshold: ', pval_threshold))
threshold_log10 <- -log10(pval_threshold)

# Prepare data for Manhattan plot
message('INFO: Preparing Manhattan data')
plot_data <- snps_assoc %>%
    dplyr::select(SNPID, chr, pos, !!sym(TRAIT)) %>%
    dplyr::rename(pvalue = !!sym(TRAIT))

manhattan_prep <- prepare_manhattan_data(plot_data, "chr", "pos", "pvalue")
plot_df <- manhattan_prep$data
chr_info <- manhattan_prep$chr_info

# Get chromosome colors
n_chr <- length(unique(plot_df$chr_f))
chr_colors <- get_chr_colors(n_chr)
names(chr_colors) <- levels(plot_df$chr_f)

# Calculate y-axis limit (add some padding)
y_max <- max(plot_df$log10p, na.rm = TRUE) * 1.1
y_max <- max(y_max, threshold_log10 * 1.2)  # Ensure threshold is visible

################################ Plot 1: Simple Manhattan

message('INFO: Generating simple Manhattan plot')

# Mark significant SNPs for current method
plot_df <- plot_df %>%
    dplyr::mutate(is_significant = log10p >= threshold_log10)

n_sig <- sum(plot_df$is_significant)
message(paste0('INFO: Found ', n_sig, ' significant SNPs above threshold'))

p_simple <- ggplot() +
    # Non-significant SNPs (chromosome colors)
    geom_point(data = plot_df %>% filter(!is_significant),
               aes(x = pos_cum, y = log10p, color = chr_f),
               alpha = 0.7, size = 0.8) +
    scale_color_manual(values = chr_colors, guide = "none") +
    # Significant SNPs (red, larger)
    geom_point(data = plot_df %>% filter(is_significant),
               aes(x = pos_cum, y = log10p),
               color = "red", alpha = 0.9, size = 1.8) +
    scale_x_continuous(
        labels = chr_info$chr_f,
        breaks = chr_info$center,
        expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.05)),
        limits = c(0, y_max)
    ) +
    geom_hline(yintercept = threshold_log10, linetype = "dashed",
               color = "red", linewidth = 0.5) +
    labs(
        title = paste0(METHOD, " - ", TRAIT),
        subtitle = paste0("K = ", Kbest, ", threshold: ", adjustment, " ", original_threshold,
                         " (", n_sig, " significant SNPs)"),
        x = "Chromosome",
        y = expression(-log[10](p-value))
    ) +
    theme_manhattan()

# Save simple plot in PNG and SVG formats
simple_base <- paste0("Manhattan_", TRAIT, "_K", Kbest, "_", ADJUST)

ggsave(file.path(PLOT_DIR, paste0(simple_base, ".png")), p_simple,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".svg")), p_simple,
       width = 10, height = 4)
message(paste0('INFO: Saved simple Manhattan plot: ', simple_base, ' (.png, .svg)'))

################################ Plot 2: Manhattan with Regions (if regions file provided)

regions_generated <- FALSE

if (!is.null(REGIONS_FILE) && file.exists(REGIONS_FILE)) {
    message('INFO: Loading regions for highlighting')
    regions <- fread(REGIONS_FILE)
    if ('chr' %in% colnames(regions)) regions$chr <- as.character(regions$chr)

    if (nrow(regions) > 0) {
        # Filter regions that contain the current trait (exact match, not substring)
        # Use word boundary to avoid "bio_1" matching "bio_12"
        regions_trait <- regions %>%
            dplyr::filter(str_detect(traits, paste0("(^|,)", TRAIT, "(,|$)")))

        message(paste0('INFO: Found ', nrow(regions_trait), ' regions for trait ', TRAIT))

        if (nrow(regions_trait) > 0) {
            # Select top 10 regions by SNP count (most significant)
            MAX_REGIONS <- 10
            if (nrow(regions_trait) > MAX_REGIONS) {
                message(paste0('INFO: Limiting to top ', MAX_REGIONS, ' regions by SNP count'))
                regions_trait <- regions_trait %>%
                    dplyr::arrange(desc(snp_count), min_pvalue) %>%
                    head(MAX_REGIONS)
            }

            # Add cumulative positions to regions
            # Ensure chr types match (fread may read chr as integer)
            regions_trait <- regions_trait %>%
                dplyr::mutate(chr = as.character(chr)) %>%
                left_join(chr_info %>%
                         dplyr::mutate(chr = as.character(chr_f)) %>%
                         dplyr::select(chr, tot),
                         by = "chr") %>%
                dplyr::mutate(
                    start_cum = start + tot,
                    end_cum = end + tot
                )

            # Ensure all region boxes have a visible minimum width
            # Use 1.5% of total genome span as minimum half-width
            total_span <- max(plot_df$pos_cum) - min(plot_df$pos_cum)
            min_half_width <- total_span * 0.015

            regions_trait <- regions_trait %>%
                dplyr::mutate(
                    region_center = (start_cum + end_cum) / 2,
                    region_half_width = pmax((end_cum - start_cum) / 2, min_half_width),
                    start_cum = region_center - region_half_width,
                    end_cum = region_center + region_half_width
                ) %>%
                dplyr::select(-region_center, -region_half_width)

            # Load sigSNPs from all methods for per-trait method attribution
            method_snps <- list()
            if (!is.null(SIGSNPS_FILES)) {
                message('INFO: Loading sigSNPs files for per-trait method attribution')
                sigsnps_vec <- str_split(SIGSNPS_FILES, ',')[[1]]
                sigsnps_vec <- sigsnps_vec[sigsnps_vec != ""]

                for (sigfile in sigsnps_vec) {
                    if (file.exists(sigfile)) {
                        dt <- fread(sigfile)
                        if ('chr' %in% colnames(dt)) dt$chr <- as.character(dt$chr)
                        if (nrow(dt) > 0) {
                            # Extract method from filename (e.g., "EMMAX" from "EMMAX_pvalues_K5_sigSNPs_bonf_0.05.tsv")
                            file_method <- basename(sigfile) %>% str_extract("^[A-Z]+")
                            # Filter for current trait only
                            dt_trait <- dt %>% dplyr::filter(trait == TRAIT)
                            if (nrow(dt_trait) > 0) {
                                method_snps[[file_method]] <- dt_trait$SNPID
                                message(paste0('INFO: ', file_method, ' has ', nrow(dt_trait), ' sigSNPs for ', TRAIT))
                            }
                        }
                    }
                }
            }

            # Build per-trait method attribution lookup
            snp_method_lookup <- NULL
            if (length(method_snps) > 0) {
                all_sig_snps <- unique(unlist(method_snps))
                snp_method_lookup <- data.frame(
                    SNPID = all_sig_snps,
                    stringsAsFactors = FALSE
                )
                # Determine which methods detected each SNP for this trait
                snp_method_lookup$source_method <- sapply(all_sig_snps, function(snp) {
                    methods_with_snp <- names(method_snps)[sapply(method_snps, function(m) snp %in% m)]
                    paste(methods_with_snp, collapse = ',')
                })
                message(paste0('INFO: Built method lookup for ', nrow(snp_method_lookup), ' sigSNPs for ', TRAIT))
            }

            # Initialize region marking
            plot_df_regions <- plot_df %>%
                dplyr::mutate(
                    in_region = FALSE,
                    region_id = NA_character_,
                    source_method = NA_character_  # Which method(s) detected this SNP for THIS trait
                )

            # Mark SNPs that are in regions
            for (i in 1:nrow(regions_trait)) {
                r <- regions_trait[i, ]

                # Find SNPs within the region boundaries (with some padding)
                in_this_region <- as.character(plot_df_regions$chr) == as.character(r$chr) &
                                  plot_df_regions$pos >= (r$start - 500000) &
                                  plot_df_regions$pos <= (r$end + 500000)

                # If we have per-trait method lookup, use it
                if (!is.null(snp_method_lookup) && nrow(snp_method_lookup) > 0) {
                    # Get SNPs that are significant (in any method) for this trait AND in this region
                    region_sig_snps <- snp_method_lookup %>%
                        dplyr::filter(SNPID %in% plot_df_regions$SNPID[in_this_region])

                    # Mark SNPs with their trait-specific method attribution
                    for (j in seq_len(nrow(region_sig_snps))) {
                        snp_row <- which(plot_df_regions$SNPID == region_sig_snps$SNPID[j])
                        if (length(snp_row) > 0) {
                            plot_df_regions$in_region[snp_row] <- TRUE
                            plot_df_regions$region_id[snp_row] <- r$region_id
                            plot_df_regions$source_method[snp_row] <- region_sig_snps$source_method[j]
                        }
                    }
                } else {
                    # No method lookup - use only significant SNPs for current method
                    snps_to_mark <- in_this_region & plot_df_regions$is_significant
                    plot_df_regions$in_region[snps_to_mark] <- TRUE
                    plot_df_regions$region_id[snps_to_mark] <- r$region_id
                    plot_df_regions$source_method[snps_to_mark] <- METHOD
                }
            }

            # Determine shape based on source method
            # Circle (16) = current method only
            # Triangle (17) = other method only
            # Diamond (18) = both methods
            plot_df_regions <- plot_df_regions %>%
                dplyr::mutate(
                    point_shape = case_when(
                        !in_region ~ NA_integer_,
                        is.na(source_method) ~ 16L,  # Default circle
                        source_method == METHOD ~ 16L,  # Current method only = circle
                        str_detect(source_method, METHOD) & str_detect(source_method, ",") ~ 18L,  # Both methods = diamond
                        str_detect(source_method, ",") ~ 18L,  # Multiple methods = diamond
                        TRUE ~ 17L  # Other method only = triangle
                    )
                )

            # Create region colors
            region_ids <- unique(na.omit(plot_df_regions$region_id))
            region_colors <- get_region_colors(length(region_ids))
            names(region_colors) <- region_ids

            n_in_regions <- sum(plot_df_regions$in_region)
            message(paste0('INFO: ', n_in_regions, ' SNPs to highlight in regions'))

            # Count by source (per-trait method attribution)
            if (!is.null(snp_method_lookup) && n_in_regions > 0) {
                n_current <- sum(plot_df_regions$in_region & plot_df_regions$source_method == METHOD, na.rm = TRUE)
                n_both <- sum(plot_df_regions$in_region & str_detect(plot_df_regions$source_method, ","), na.rm = TRUE)
                n_other <- n_in_regions - n_current - n_both
                message(paste0('INFO: Per-trait sources for ', TRAIT, ': ', METHOD, '=', n_current, ', other=', n_other, ', both=', n_both))
            }

            # Create the plot with region highlighting
            message('INFO: Generating Manhattan plot with regions')

            # Prepare data subsets
            df_background <- plot_df_regions %>% filter(!in_region & !is_significant)
            df_sig_not_region <- plot_df_regions %>% filter(is_significant & !in_region)
            df_in_region <- plot_df_regions %>% filter(in_region)

            p_regions <- ggplot() +
                # Region rectangles (background highlight) - draw first
                geom_rect(data = regions_trait,
                         aes(xmin = start_cum, xmax = end_cum,
                             ymin = 0, ymax = y_max, fill = region_id),
                         alpha = 0.1) +
                scale_fill_manual(values = region_colors, guide = "none") +
                # Background SNPs (non-significant, not in region, chromosome colors)
                geom_point(data = df_background,
                          aes(x = pos_cum, y = log10p, color = chr_f),
                          alpha = 0.5, size = 0.6) +
                scale_color_manual(values = chr_colors, guide = "none") +
                # Significant SNPs not in any region (red)
                geom_point(data = df_sig_not_region,
                          aes(x = pos_cum, y = log10p),
                          color = "red", alpha = 0.8, size = 1.5) +
                # SNPs in regions - colored by region, shaped by source method
                ggnewscale::new_scale_color()

            # Add region SNPs with different shapes
            if (nrow(df_in_region) > 0) {
                # Current method (circle)
                df_current <- df_in_region %>% filter(point_shape == 16)
                if (nrow(df_current) > 0) {
                    p_regions <- p_regions +
                        geom_point(data = df_current,
                                  aes(x = pos_cum, y = log10p, color = region_id),
                                  shape = 16, alpha = 0.9, size = 2.5)
                }

                # Other method (triangle)
                df_other <- df_in_region %>% filter(point_shape == 17)
                if (nrow(df_other) > 0) {
                    p_regions <- p_regions +
                        geom_point(data = df_other,
                                  aes(x = pos_cum, y = log10p, color = region_id),
                                  shape = 17, alpha = 0.9, size = 2.5)
                }

                # Both methods (diamond)
                df_both <- df_in_region %>% filter(point_shape == 18)
                if (nrow(df_both) > 0) {
                    p_regions <- p_regions +
                        geom_point(data = df_both,
                                  aes(x = pos_cum, y = log10p, color = region_id),
                                  shape = 18, alpha = 0.9, size = 3.0)
                }
            }

            p_regions <- p_regions +
                scale_color_manual(values = region_colors, name = "Region") +
                # Threshold line
                geom_hline(yintercept = threshold_log10, linetype = "dashed",
                          color = "red", linewidth = 0.5) +
                # Axis settings
                scale_x_continuous(
                    labels = chr_info$chr_f,
                    breaks = chr_info$center,
                    expand = c(0.01, 0.01)
                ) +
                scale_y_continuous(
                    expand = expansion(mult = c(0, 0.05)),
                    limits = c(0, y_max)
                ) +
                labs(
                    title = paste0(METHOD, " - ", TRAIT, " (with regions)"),
                    subtitle = paste0("K = ", Kbest, ", showing ", nrow(regions_trait), " regions",
                                     " | shapes: \u25CF=", METHOD, ", \u25B2=other, \u25C6=both"),
                    x = "Chromosome",
                    y = expression(-log[10](p-value))
                ) +
                theme_manhattan() +
                theme(
                    legend.position = "right",
                    legend.text = element_text(size = 6),
                    legend.title = element_text(size = 8),
                    legend.key.size = unit(0.4, "cm")
                ) +
                guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))

            # Save regions plot in PNG and SVG formats
            regions_base <- paste0("Manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, "_regions")
            ggsave(file.path(PLOT_DIR, paste0(regions_base, ".png")), p_regions,
                   width = 12, height = 4.5, dpi = 300)
            ggsave(file.path(PLOT_DIR, paste0(regions_base, ".svg")), p_regions,
                   width = 12, height = 4.5)
            message(paste0('INFO: Saved regions Manhattan plot: ', regions_base, ' (.png, .svg)'))
            regions_generated <- TRUE

        } else {
            message('INFO: No regions found for this trait')
        }
    } else {
        message('INFO: Regions file is empty')
    }
} else {
    message('INFO: No regions file provided or file not found')
}

# If no regions plot was generated, save the simple plot as the regions version
# This ensures the output files exist for snakemake
if (!regions_generated) {
    message('INFO: Saving simple plot as regions version (no regions to highlight)')
    p_no_regions <- p_simple +
        labs(subtitle = paste0("K = ", Kbest, ", threshold: ", adjustment, " ", original_threshold, " (no significant regions)"))
    regions_base <- paste0("Manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, "_regions")
    ggsave(file.path(PLOT_DIR, paste0(regions_base, ".png")), p_no_regions,
           width = 10, height = 4, dpi = 300)
    ggsave(file.path(PLOT_DIR, paste0(regions_base, ".svg")), p_no_regions,
           width = 10, height = 4)
}

message('INFO: Manhattan plot complete')
