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

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/utils/manhattan_utils.R")

args = commandArgs(trailingOnly=TRUE)
################################
ASSOC_TABLE    = args[1]   # SNPID chr pos TRAITS...
ADJUST         = args[2]   # e.g., "bonf_0.05" or "qval_0.05"
Kbest          = args[3] %>% as.numeric
METHOD         = args[4]   # EMMAX or LFMM (current plot method)
TRAIT          = args[5]   # single trait name, e.g., "bio_1"
PLOT_DIR       = args[6]   # output directory
ALL_PREDICTORS = args[7]   # Comma-separated list of all trait names
REGIME         = if (length(args) >= 8) tolower(args[8]) else "snp"  # "snp" or "wza"
################################

# WZA mode uses same output format but with "wza_" filename prefix
is_wza    <- REGIME == "wza"
name_pfx  <- if (is_wza) "wza_" else ""

message(paste0('INFO: Manhattan plot for ', METHOD, ' - ', TRAIT, ' [regime=', REGIME, ']'))
message(paste0('INFO: K = ', Kbest))
message(paste0('INFO: Adjustment: ', ADJUST))
if (is_wza) message('INFO: WZA regime — Bonferroni denominator = number of windows')

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

message(paste0('INFO: Using scattermore for rendering (', n_snps, ' SNPs)'))

# Check if trait exists
if (!(TRAIT %in% colnames(snps_assoc))) {
    stop(paste0("Trait '", TRAIT, "' not found in association table"))
}

# Parse adjustment parameters
adjustment <- ADJUST %>% str_split('_') %>% unlist %>% .[1]
original_threshold <- ADJUST %>% str_split('_') %>% unlist %>% .[2] %>% as.numeric

message(paste0('INFO: Adjustment method: ', adjustment))
message(paste0('INFO: Threshold: ', original_threshold))

# WZA: filter to rows that have the trait column filled (windows with ≥1 SNP)
if (is_wza && "n_snps" %in% colnames(snps_assoc)) {
    snps_assoc <- snps_assoc[!is.na(snps_assoc[[TRAIT]]), ]
    n_snps <- nrow(snps_assoc)
    message(paste0('INFO: After NA filter: ', n_snps, ' windows'))
}

# Calculate threshold based on adjustment method
pval_threshold <- compute_pval_threshold(snps_assoc[[TRAIT]], adjustment, original_threshold)

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
    # Non-significant SNPs (chromosome colors) - scattermore for performance
    add_scatter_layer(data = plot_df %>% filter(!is_significant),
                      chr_colors = chr_colors,
                      alpha = 0.5) +
    scale_color_identity()

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
simple_base <- paste0("manhattan_", name_pfx, TRAIT, "_K", Kbest, "_", ADJUST)

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
    labs(
        title = paste0(METHOD, " - ", TRAIT, " QQ Plot"),
        subtitle = paste0("K = ", Kbest, ", ", n_sig, " significant SNPs"),
        x = expression(Expected ~ -log[10](p-value)),
        y = expression(Observed ~ -log[10](p-value))
    ) +
    theme_minimal() +
    theme(
        plot.title    = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.margin   = margin(4, 8, 4, 4)
    )

qq_base <- paste0("qq_", name_pfx, TRAIT, "_K", Kbest, "_", ADJUST)
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
                      alpha = 0.5) +
    scale_color_identity() +
    geom_hline(yintercept = threshold_log10, linetype = "dashed",
               color = "red", linewidth = 0.5) +
    scale_x_continuous(limits = c(bg_x_lo, bg_x_hi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(bg_y_lo, bg_y_hi), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

bg_base <- paste0("manhattan_", name_pfx, TRAIT, "_K", Kbest, "_", ADJUST, "_background")
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
    paste0("manhattan_", name_pfx, TRAIT, "_K", Kbest, "_", ADJUST, "_coords.json"))
jsonlite::write_json(coords_list, coords_file, auto_unbox = TRUE, digits = 6)
message(paste0('INFO: Saved coords JSON: ', basename(coords_file)))

message('INFO: Manhattan plot complete')
