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
library(scattermore)
library(jsonlite)

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
SCATTERMORE_THRESHOLD = args[9] %>% as.numeric  # Use scattermore when SNPs > this
ALL_PREDICTORS = args[10]  # Comma-separated list of all trait names (for consistent Okabe-Ito colors)
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

# Okabe-Ito colorblind-safe palette for traits
get_trait_colors <- function(traits) {
    okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999999")
    if (length(traits) <= 8) {
        setNames(okabe_ito[1:length(traits)], traits)
    } else {
        setNames(viridis::turbo(length(traits)), traits)
    }
}

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

    # Add gaps between chromosomes (2% of mean chr length)
    chr_gap <- mean(chr_lengths$chr_len) * 0.02

    chr_lengths <- chr_lengths %>%
        dplyr::mutate(
            chr_idx = dplyr::row_number(),
            tot = cumsum(as.numeric(chr_len)) - chr_len + (chr_idx - 1) * chr_gap,
            center = tot + chr_len / 2
        ) %>%
        dplyr::select(-chr_idx)

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

# Smart geom layer: uses scattermore for large datasets, geom_point for small
# scattermore doesn't work with ggplot2 color scales, so colors must be pre-computed
# Parameters:
#   data: data frame with pos_cum, log10p, and optionally chr_f for coloring
#   chr_colors: named vector of colors (names = chromosome levels)
#   use_scattermore: whether to use scattermore (TRUE for large datasets)
#   pointsize: for scattermore, use 3-4 for visibility (non-integers like 3.2 look better)
#   pixels: raster resolution, should match output size (width*dpi, height*dpi)
add_scatter_layer <- function(data, chr_colors, alpha = 0.5, size = 0.8,
                              use_scattermore = FALSE, pointsize = 10) {
    if (nrow(data) == 0) {
        return(geom_blank())
    }

    if (use_scattermore) {
        # scattermore needs pre-computed colors (doesn't work with scale_color_manual)
        data$point_color <- chr_colors[as.character(data$chr_f)]
        # pixels should match output: 10in x 4in @ 300dpi = 3000x1200
        geom_scattermore(data = data,
                         aes(x = pos_cum, y = log10p, color = point_color),
                         alpha = alpha, pointsize = pointsize,
                         pixels = c(3000, 1200), interpolate = FALSE)
    } else {
        geom_point(data = data,
                   aes(x = pos_cum, y = log10p, color = chr_f),
                   alpha = alpha, size = size)
    }
}

################################ Compute trait color

# Parse all predictors and compute trait-specific Okabe-Ito color
all_traits <- str_split(ALL_PREDICTORS, ',')[[1]] %>% str_trim()
all_trait_colors <- get_trait_colors(all_traits)
sig_color <- unname(all_trait_colors[TRAIT])
if (is.na(sig_color)) sig_color <- "#E69F00"  # fallback
message(paste0('INFO: Trait color: ', sig_color))

################################ Main

# Load association data
message('INFO: Loading association data')
snps_assoc <- fread(ASSOC_TABLE, sep = '\t', header = T)
snps_assoc$chr <- as.character(snps_assoc$chr)
n_snps <- nrow(snps_assoc)
message(paste0('INFO: Loaded ', n_snps, ' SNPs'))

# Decide rendering method based on SNP count
use_scattermore <- n_snps > SCATTERMORE_THRESHOLD
if (use_scattermore) {
    message(paste0('INFO: Using scattermore for fast rendering (', n_snps, ' > ', SCATTERMORE_THRESHOLD, ')'))
} else {
    message(paste0('INFO: Using standard geom_point (', n_snps, ' <= ', SCATTERMORE_THRESHOLD, ')'))}


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
    # Non-significant SNPs (chromosome colors) - use scattermore if large dataset
    add_scatter_layer(data = plot_df %>% filter(!is_significant),
                      chr_colors = chr_colors,
                      alpha = 0.5, size = 0.8,
                      use_scattermore = use_scattermore)

# Add appropriate color scale based on rendering method
if (use_scattermore) {
    p_simple <- p_simple + scale_color_identity()
} else {
    p_simple <- p_simple + scale_color_manual(values = chr_colors, guide = "none")
}

p_simple <- p_simple +
    # Significant SNPs (trait-colored, larger) - always use geom_point for proper styling
    geom_point(data = plot_df %>% filter(is_significant),
               aes(x = pos_cum, y = log10p),
               color = sig_color, alpha = 0.9, size = 2.5) +
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
simple_base <- paste0("manhattan_", TRAIT, "_K", Kbest, "_", ADJUST)

ggsave(file.path(PLOT_DIR, paste0(simple_base, ".png")), p_simple,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".svg")), p_simple,
       width = 10, height = 4)
message(paste0('INFO: Saved simple Manhattan plot: ', simple_base, ' (.png, .svg)'))

################################ QQ Plot (per-method, single trait)

message('INFO: Generating QQ plot')

qq_df <- plot_df %>%
    arrange(pvalue) %>%
    mutate(
        expected = -log10(ppoints(n())),
        observed = log10p
    )

# 95% confidence interval band
n_qq <- nrow(qq_df)
ci_qq <- data.frame(
    expected = -log10(ppoints(n_qq)),
    ci_lower = -log10(qbeta(0.975, 1:n_qq, n_qq:1)),
    ci_upper = -log10(qbeta(0.025, 1:n_qq, n_qq:1))
)

p_qq <- ggplot() +
    geom_ribbon(data = ci_qq,
                aes(x = expected, ymin = ci_lower, ymax = ci_upper),
                fill = "grey80", alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    # Non-significant points (alternating chr colors)
    geom_point(data = qq_df %>% filter(!is_significant),
               aes(x = expected, y = observed, color = chr_f),
               alpha = 0.5, size = 1.5) +
    scale_color_manual(values = chr_colors, guide = "none") +
    # Significant points (trait color)
    geom_point(data = qq_df %>% filter(is_significant),
               aes(x = expected, y = observed),
               color = sig_color, size = 2.5) +
    coord_fixed() +
    labs(
        title = paste0(METHOD, " - ", TRAIT, " QQ Plot"),
        subtitle = paste0("K = ", Kbest, ", ", n_sig, " significant SNPs"),
        x = expression(Expected ~ -log[10](p-value)),
        y = expression(Observed ~ -log[10](p-value))
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    )

qq_base <- paste0("qq_", TRAIT, "_K", Kbest, "_", ADJUST)
ggsave(file.path(PLOT_DIR, paste0(qq_base, ".png")), p_qq,
       width = 6, height = 6, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(qq_base, ".svg")), p_qq,
       width = 6, height = 6)
message(paste0('INFO: Saved QQ plot: ', qq_base, ' (.png, .svg)'))

################################ Background PNG + Coords JSON (for Shiny plotly overlay)

message('INFO: Generating background Manhattan for Shiny overlay')

# Compute the exact axis ranges the ggplot renders (matching expand parameters)
x_min_data <- min(plot_df$pos_cum)
x_max_data <- max(plot_df$pos_cum)
x_span     <- x_max_data - x_min_data
bg_x_lo    <- x_min_data - 0.01 * x_span   # expand = c(0.01, 0.01)
bg_x_hi    <- x_max_data + 0.01 * x_span

# y: limits = c(0, y_max), expand = mult c(0, 0.05) — use 0 bottom so both simple and
# regions versions share one coords.json (regions uses mult c(0, 0.05) anyway)
bg_y_lo <- 0
bg_y_hi <- y_max * 1.05

# Background-only plot: non-sig SNPs + threshold line, no axes (Shiny draws its own)
p_bg <- ggplot() +
    add_scatter_layer(data = plot_df %>% dplyr::filter(!is_significant),
                      chr_colors = chr_colors,
                      alpha = 0.5, size = 0.8,
                      use_scattermore = use_scattermore) +
    geom_hline(yintercept = threshold_log10, linetype = "dashed",
               color = "red", linewidth = 0.5) +
    scale_x_continuous(limits = c(bg_x_lo, bg_x_hi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(bg_y_lo, bg_y_hi), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

if (use_scattermore) {
    p_bg <- p_bg + scale_color_identity()
} else {
    p_bg <- p_bg + scale_color_manual(values = chr_colors, guide = "none")
}

bg_base <- paste0("manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, "_background")
ggsave(file.path(PLOT_DIR, paste0(bg_base, ".png")), p_bg,
       width = 10, height = 4, dpi = 300)
message(paste0('INFO: Saved background Manhattan: ', bg_base, '.png'))

# Coordinate mapping JSON for plotly axis alignment in Shiny
coords_list <- list(
    chr_offsets    = setNames(as.list(chr_info$tot),     as.character(chr_info$chr_f)),
    chr_lengths    = setNames(as.list(chr_info$chr_len), as.character(chr_info$chr_f)),
    gap_fraction   = 0.02,
    x_range        = c(bg_x_lo, bg_x_hi),
    y_range        = c(bg_y_lo, bg_y_hi),
    bonferroni_y   = threshold_log10,
    plot_width_px  = 3000L,
    plot_height_px = 1200L
)
coords_file <- file.path(PLOT_DIR,
    paste0("manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, "_coords.json"))
jsonlite::write_json(coords_list, coords_file, auto_unbox = TRUE, digits = 6)
message(paste0('INFO: Saved coords JSON: ', basename(coords_file)))

################################ Plot 2: Manhattan with Regions (if regions file provided)

regions_generated <- FALSE

if (!is.null(REGIONS_FILE) && file.exists(REGIONS_FILE)) {
    message('INFO: Loading regions for highlighting')
    regions <- fread(REGIONS_FILE)
    if ('chr' %in% colnames(regions)) regions$chr <- as.character(regions$chr)

    if (nrow(regions) > 0) {
        # Filter regions for the current trait
        # Per-trait regions have a "trait" column with exact matches
        regions_trait <- regions %>%
            dplyr::filter(trait == TRAIT)

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

            # Region colors with chr:start-end labels for legend
            region_ids <- unique(regions_trait$region_id)
            region_colors <- get_region_colors(length(region_ids))
            names(region_colors) <- region_ids

            # Create display labels: chr:start-end
            region_labels <- regions_trait %>%
                dplyr::distinct(region_id, .keep_all = TRUE) %>%
                dplyr::mutate(label = paste0(chr, ":", scales::comma(start), "-", scales::comma(end))) %>%
                dplyr::select(region_id, label)
            region_label_map <- setNames(region_labels$label, region_labels$region_id)

            # Remap region_id to display label in data and colors
            regions_trait$region_id <- region_label_map[regions_trait$region_id]
            region_colors <- setNames(region_colors, region_label_map[names(region_colors)])

            # Create the plot with region highlighting
            message('INFO: Generating Manhattan plot with regions')

            # Prepare data subsets
            df_background <- plot_df %>% filter(!is_significant)
            df_sig_reg <- plot_df %>% filter(is_significant)

            p_regions <- ggplot() +
                # Region rectangles (background highlight) - draw first
                geom_rect(data = regions_trait,
                         aes(xmin = start_cum, xmax = end_cum,
                             ymin = 0, ymax = y_max, fill = region_id),
                         alpha = 0.15) +
                scale_fill_manual(values = region_colors, name = "Region") +
                # Background SNPs (non-significant) - use scattermore if large
                add_scatter_layer(data = df_background,
                                  chr_colors = chr_colors,
                                  alpha = 0.5, size = 0.8,
                                  use_scattermore = use_scattermore)

            # Add appropriate color scale based on rendering method
            if (use_scattermore) {
                p_regions <- p_regions + scale_color_identity()
            } else {
                p_regions <- p_regions + scale_color_manual(values = chr_colors, guide = "none")
            }

            p_regions <- p_regions +
                # All significant SNPs - trait-colored
                geom_point(data = df_sig_reg,
                          aes(x = pos_cum, y = log10p),
                          color = sig_color, alpha = 0.9, size = 2.5) +
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
                    subtitle = paste0("K = ", Kbest, ", showing ", nrow(regions_trait), " regions"),
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
                guides(
                    fill = guide_legend(ncol = 1, override.aes = list(alpha = 0.3))
                )

            # Save regions plot in PNG and SVG formats
            regions_base <- paste0("manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, "_regions")
            ggsave(file.path(PLOT_DIR, paste0(regions_base, ".png")), p_regions,
                   width = 12, height = 4.5, dpi = 300)
            ggsave(file.path(PLOT_DIR, paste0(regions_base, ".svg")), p_regions,
                   width = 12, height = 4.5)
            message(paste0('INFO: Saved regions Manhattan plot: ', regions_base, ' (.png, .svg)'))

            # Background regions version for Shiny (no sig SNPs, no axes)
            p_bg_regions <- ggplot() +
                geom_rect(data = regions_trait,
                         aes(xmin = start_cum, xmax = end_cum,
                             ymin = bg_y_lo, ymax = bg_y_hi, fill = region_id),
                         alpha = 0.15) +
                scale_fill_manual(values = region_colors, guide = "none") +
                add_scatter_layer(data = df_background,
                                  chr_colors = chr_colors,
                                  alpha = 0.5, size = 0.8,
                                  use_scattermore = use_scattermore) +
                geom_hline(yintercept = threshold_log10, linetype = "dashed",
                           color = "red", linewidth = 0.5) +
                scale_x_continuous(limits = c(bg_x_lo, bg_x_hi), expand = c(0, 0)) +
                scale_y_continuous(limits = c(bg_y_lo, bg_y_hi), expand = c(0, 0)) +
                theme_void() +
                theme(plot.background = element_rect(fill = "white", color = NA))

            if (use_scattermore) {
                p_bg_regions <- p_bg_regions + scale_color_identity()
            } else {
                p_bg_regions <- p_bg_regions + scale_color_manual(values = chr_colors, guide = "none")
            }

            bg_reg_base <- paste0("manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, "_regions_background")
            ggsave(file.path(PLOT_DIR, paste0(bg_reg_base, ".png")), p_bg_regions,
                   width = 12, height = 4.5, dpi = 300)
            message(paste0('INFO: Saved background regions Manhattan: ', bg_reg_base, '.png'))

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
    regions_base <- paste0("manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, "_regions")
    ggsave(file.path(PLOT_DIR, paste0(regions_base, ".png")), p_no_regions,
           width = 10, height = 4, dpi = 300)
    ggsave(file.path(PLOT_DIR, paste0(regions_base, ".svg")), p_no_regions,
           width = 10, height = 4)
}

# Ensure _regions_background.png exists when a regions file was provided
# (may not be created if no regions match the current trait)
if (!is.null(REGIONS_FILE)) {
    expected_bg_reg <- file.path(PLOT_DIR,
        paste0("manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, "_regions_background.png"))
    if (!file.exists(expected_bg_reg)) {
        file.copy(file.path(PLOT_DIR, paste0(bg_base, ".png")), expected_bg_reg)
        message('INFO: No regions for this trait - copied simple background as regions background')
    }
}

message('INFO: Manhattan plot complete')
