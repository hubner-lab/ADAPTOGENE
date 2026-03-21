#!/usr/bin/env Rscript
# Combined Manhattan plot for all traits and methods
# Shows all GWAS results on a single plot with:
# - Color = trait (Okabe-Ito colorblind-safe palette)
# - Shape = method (up to 10 distinct shapes)
# - Every (SNP, trait, method) detection shown individually (no "Multi" collapse)
# - Regions marked with rectangles (from Regions_climate_combined.tsv)
# - Always uses scattermore for background SNPs

library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(qvalue)
library(scales)
library(scattermore)
library(ggnewscale)

args = commandArgs(trailingOnly=TRUE)
################################
ASSOC_FILES_STR = args[1]   # Comma-separated "METHOD:ADJUST:FILEPATH"
PREDICTORS = args[2]        # Comma-separated trait names
Kbest = args[3] %>% as.numeric
REGIONS_FILE = args[4]      # Regions_climate_combined.tsv
PLOT_DIR = args[5]          # Output directory
################################

message('INFO: Combined Manhattan plot for all traits and methods')
message(paste0('INFO: K = ', Kbest))

# Parse ASSOC_FILES_STR into list
assoc_list <- str_split(ASSOC_FILES_STR, ',')[[1]]
assoc_info <- list()
for (item in assoc_list) {
    parts <- str_split(item, ':')[[1]]
    if (length(parts) == 3) {
        method <- parts[1]
        adjust <- parts[2]
        filepath <- parts[3]
        assoc_info[[method]] <- list(method = method, adjust = adjust, filepath = filepath)
    }
}

if (length(assoc_info) == 0) {
    stop("No valid association files provided")
}

message(paste0('INFO: Found ', length(assoc_info), ' methods: ', paste(names(assoc_info), collapse = ', ')))

# Parse predictors
traits <- str_split(PREDICTORS, ',')[[1]] %>% str_trim()
message(paste0('INFO: Traits: ', paste(traits, collapse = ', ')))

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

# Distinct shapes for methods (up to 10)
get_method_shapes <- function(methods) {
    shapes <- c(16, 17, 15, 18, 3, 4, 8, 6, 9, 14)
    # circle, triangle, square, diamond, cross, X, star, inv-triangle, diamond+cross, square+triangle
    setNames(shapes[1:length(methods)], methods)
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
        plot.margin = margin(10, 15, 10, 10)
    )
}

# Generate color palette for regions
get_region_colors <- function(n_regions) {
    if (n_regions <= 10) {
        colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                   "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5")
        return(colors[1:n_regions])
    } else {
        return(scales::hue_pal()(n_regions))
    }
}

################################ Main

# Load and combine all association data
message('INFO: Loading association data for all methods and traits')

all_data <- list()
method_thresholds <- list()

for (method_name in names(assoc_info)) {
    info <- assoc_info[[method_name]]
    message(paste0('INFO: Loading ', method_name, ' from ', info$filepath))

    snps_assoc <- fread(info$filepath, sep = '\t', header = TRUE)
    snps_assoc$chr <- as.character(snps_assoc$chr)

    # Parse adjustment parameters
    adjustment <- str_split(info$adjust, '_')[[1]][1]
    pval_threshold <- str_split(info$adjust, '_')[[1]][2] %>% as.numeric()

    message(paste0('INFO: ', method_name, ' adjustment: ', adjustment, ' ', pval_threshold))

    # For each trait, create long-format data
    for (trait in traits) {
        if (!(trait %in% colnames(snps_assoc))) {
            message(paste0('WARNING: Trait ', trait, ' not found in ', method_name, ' data'))
            next
        }

        # Calculate threshold
        original_threshold <- pval_threshold
        threshold_value <- pval_threshold

        if (adjustment == 'bonf') {
            threshold_value <- original_threshold / nrow(snps_assoc)
        }
        if (adjustment == 'qval') {
            threshold_value <- FUN_max_pvalue_fdr(snps_assoc[[trait]], original_threshold)
            if (is.null(threshold_value)) {
                message(paste0('WARNING: No significant SNPs at q-value threshold for ', method_name, ' ', trait))
                threshold_value <- original_threshold / nrow(snps_assoc)
            }
        }
        if (adjustment == 'top') {
            threshold_value <- FUN_max_pvalue_top(snps_assoc[[trait]], original_threshold)
        }

        threshold_log10 <- -log10(threshold_value)
        method_thresholds[[paste0(method_name, "_", trait)]] <- threshold_log10

        # Create long-format data
        trait_data <- snps_assoc %>%
            dplyr::select(SNPID, chr, pos, !!sym(trait)) %>%
            dplyr::rename(pvalue = !!sym(trait)) %>%
            dplyr::mutate(
                method = method_name,
                trait = trait,
                log10p = -log10(pvalue),
                is_significant = log10p >= threshold_log10
            )

        all_data[[paste0(method_name, "_", trait)]] <- trait_data

        n_sig <- sum(trait_data$is_significant)
        message(paste0('INFO: ', method_name, ' - ', trait, ': ', n_sig, ' significant SNPs'))
    }
}

# Combine all data
plot_data_all <- bind_rows(all_data)
message(paste0('INFO: Combined data: ', nrow(plot_data_all), ' total data points'))

# Calculate minimum threshold (most stringent)
min_threshold <- min(unlist(method_thresholds))
message(paste0('INFO: Minimum threshold (most stringent): -log10(p) = ', round(min_threshold, 2)))

# Prepare manhattan data (using first trait for chromosome info)
reference_data <- plot_data_all %>%
    dplyr::filter(trait == traits[1], method == names(assoc_info)[1]) %>%
    dplyr::select(SNPID, chr, pos, pvalue)

manhattan_prep <- prepare_manhattan_data(reference_data, "chr", "pos", "pvalue")
chr_info <- manhattan_prep$chr_info

# Add cumulative positions to all data
plot_data_all <- plot_data_all %>%
    left_join(chr_info %>%
              dplyr::mutate(chr = as.character(chr_f)) %>%
              dplyr::select(chr, tot),
              by = "chr") %>%
    dplyr::mutate(pos_cum = pos + tot)

# Get palettes
trait_colors <- get_trait_colors(traits)
method_shapes <- get_method_shapes(names(assoc_info))

# Calculate y-axis limit
y_max <- max(plot_data_all$log10p, na.rm = TRUE) * 1.1

message('INFO: Generating combined Manhattan plot (simple version)')

# Background SNPs: non-significant, from first method only (avoid overplotting)
df_background <- plot_data_all %>%
    dplyr::filter(!is_significant) %>%
    dplyr::filter(method == names(assoc_info)[1])

# Significant SNPs: every (SNP, trait, method) combination shown individually
df_significant <- plot_data_all %>%
    dplyr::filter(is_significant)

n_sig_points <- nrow(df_significant)
n_sig_unique <- length(unique(df_significant$SNPID))
message(paste0('INFO: Total significant points (SNP x trait x method): ', n_sig_points))
message(paste0('INFO: Unique significant SNPs: ', n_sig_unique))

# Simple plot — per-trait background scattermore (chr gaps for separation)
p_simple <- ggplot()

# Per-trait background SNPs (scattermore, first method only)
first_method <- names(assoc_info)[1]
for (t in traits) {
    df_t <- df_background %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_simple <- p_simple +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = log10p, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 1500), interpolate = FALSE)
}

p_simple <- p_simple +
    scale_color_identity() +
    # 3. Significant SNPs - color by trait, shape by method
    ggnewscale::new_scale_color() +
    geom_point(data = df_significant,
               aes(x = pos_cum, y = log10p, color = trait, shape = method),
               alpha = 0.9, size = 2.5, stroke = 0.5) +
    scale_color_manual(values = trait_colors, name = "Trait") +
    scale_shape_manual(values = method_shapes, name = "Method") +
    scale_x_continuous(
        labels = chr_info$chr_f,
        breaks = chr_info$center,
        expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.05)),
        limits = c(0, y_max)
    ) +
    geom_hline(yintercept = min_threshold, linetype = "dashed",
               color = "red", linewidth = 0.5, alpha = 0.5) +
    labs(
        title = "Combined Association Analysis (All Traits & Methods)",
        subtitle = paste0("K = ", Kbest, " | ", n_sig_unique, " significant SNPs, ",
                         n_sig_points, " detections across ",
                         length(names(assoc_info)), " methods"),
        x = "Chromosome",
        y = expression(-log[10](p-value))
    ) +
    theme_manhattan() +
    theme(
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.5, "cm")
    ) +
    guides(
        color = guide_legend(order = 1),
        shape = guide_legend(order = 2)
    )

# Save simple plot
simple_base <- paste0("manhattan_combined_K", Kbest)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".png")), p_simple,
       width = 14, height = 5, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".svg")), p_simple,
       width = 14, height = 5)
message(paste0('INFO: Saved combined simple plot: ', simple_base))

################################ Plot with Regions

if (file.exists(REGIONS_FILE)) {
    message('INFO: Loading combined climate regions')
    regions <- fread(REGIONS_FILE)
    if ('chr' %in% colnames(regions)) regions$chr <- as.character(regions$chr)

    if (nrow(regions) > 0) {
        message(paste0('INFO: Found ', nrow(regions), ' combined regions'))

        # Limit to top 10 regions by SNP count
        MAX_REGIONS <- 10
        if (nrow(regions) > MAX_REGIONS) {
            message(paste0('INFO: Limiting to top ', MAX_REGIONS, ' regions'))
            regions <- regions %>%
                dplyr::arrange(desc(snp_count), min_pvalue) %>%
                head(MAX_REGIONS)
        }

        # Add cumulative positions
        regions <- regions %>%
            dplyr::mutate(chr = as.character(chr)) %>%
            left_join(chr_info %>%
                     dplyr::mutate(chr = as.character(chr_f)) %>%
                     dplyr::select(chr, tot),
                     by = "chr") %>%
            dplyr::mutate(
                start_cum = start + tot,
                end_cum = end + tot
            )

        # Ensure minimum visible width for regions
        total_span <- max(plot_data_all$pos_cum) - min(plot_data_all$pos_cum)
        min_half_width <- total_span * 0.015

        regions <- regions %>%
            dplyr::mutate(
                region_center = (start_cum + end_cum) / 2,
                region_half_width = pmax((end_cum - start_cum) / 2, min_half_width),
                start_cum = region_center - region_half_width,
                end_cum = region_center + region_half_width
            ) %>%
            dplyr::select(-region_center, -region_half_width)

        # Region colors with chr:start-end labels for legend
        region_ids <- unique(regions$region_id)
        region_colors <- get_region_colors(length(region_ids))
        names(region_colors) <- region_ids

        # Create display labels: chr:start-end
        region_labels <- regions %>%
            dplyr::distinct(region_id, .keep_all = TRUE) %>%
            dplyr::mutate(label = paste0(chr, ":", scales::comma(start), "-", scales::comma(end))) %>%
            dplyr::select(region_id, label)
        region_label_map <- setNames(region_labels$label, region_labels$region_id)

        # Remap region_id to display label in data and colors
        regions$region_id <- region_label_map[regions$region_id]
        region_colors <- setNames(region_colors, region_label_map[names(region_colors)])

        # Prepare data subsets (same as simple version)
        df_background_reg <- plot_data_all %>%
            dplyr::filter(!is_significant) %>%
            dplyr::filter(method == names(assoc_info)[1])

        df_sig_reg <- plot_data_all %>%
            dplyr::filter(is_significant)

        # Plot with regions — per-trait background + region rectangles (chr gaps for separation)
        message('INFO: Generating combined Manhattan plot with regions')

        p_regions <- ggplot() +
            # 1. Region rectangles (with legend)
            geom_rect(data = regions,
                     aes(xmin = start_cum, xmax = end_cum,
                         ymin = 0, ymax = y_max, fill = region_id),
                     alpha = 0.15) +
            scale_fill_manual(values = region_colors, name = "Region")

        # 3. Per-trait background SNPs
        for (t in traits) {
            df_t <- df_background_reg %>% dplyr::filter(trait == t)
            df_t$point_color <- unname(trait_colors[t])
            p_regions <- p_regions +
                geom_scattermore(data = df_t,
                                 aes(x = pos_cum, y = log10p, color = point_color),
                                 alpha = 0.5, pointsize = 10,
                                 pixels = c(4000, 1500), interpolate = FALSE)
        }

        p_regions <- p_regions +
            scale_color_identity() +
            # 4. All significant SNPs - color by trait, shape by method
            ggnewscale::new_scale_color() +
            geom_point(data = df_sig_reg,
                      aes(x = pos_cum, y = log10p, color = trait, shape = method),
                      alpha = 0.9, size = 2.5, stroke = 0.5) +
            scale_color_manual(values = trait_colors, name = "Trait") +
            scale_shape_manual(values = method_shapes, name = "Method") +
            # Threshold line
            geom_hline(yintercept = min_threshold, linetype = "dashed",
                      color = "red", linewidth = 0.5, alpha = 0.5) +
            # Axes
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
                title = "Combined Association Analysis with Regions (All Traits & Methods)",
                subtitle = paste0("K = ", Kbest, " | ", length(region_ids), " regions | color=trait, shape=method"),
                x = "Chromosome",
                y = expression(-log[10](p-value))
            ) +
            theme_manhattan() +
            theme(
                legend.position = "right",
                legend.text = element_text(size = 7),
                legend.title = element_text(size = 8, face = "bold"),
                legend.key.size = unit(0.4, "cm"),
                legend.box = "vertical"
            ) +
            guides(
                fill = guide_legend(ncol = 1, order = 3, override.aes = list(alpha = 0.3)),
                color = guide_legend(ncol = 1, order = 1, override.aes = list(size = 2)),
                shape = guide_legend(ncol = 1, order = 2, override.aes = list(size = 2))
            )

        # Save regions plot
        regions_base <- paste0("manhattan_combined_K", Kbest, "_regions")
        ggsave(file.path(PLOT_DIR, paste0(regions_base, ".png")), p_regions,
               width = 15, height = 5.5, dpi = 300)
        ggsave(file.path(PLOT_DIR, paste0(regions_base, ".svg")), p_regions,
               width = 15, height = 5.5)
        message(paste0('INFO: Saved combined regions plot: ', regions_base))

    } else {
        message('INFO: Regions file is empty')
    }
} else {
    message('INFO: Regions file not found')
}

################################ QQ Plot (Combined)

message('INFO: Generating combined QQ plot')

# For each (trait, method) combo: sort p-values, compute expected uniform quantiles
qq_data <- plot_data_all %>%
    group_by(trait, method) %>%
    arrange(pvalue) %>%
    mutate(
        expected = -log10(ppoints(n())),
        observed = -log10(pvalue)
    ) %>%
    ungroup()

# 95% confidence interval band (under uniform null)
n_snps_ci <- max(table(paste(qq_data$trait, qq_data$method)))
ci_data <- data.frame(
    expected = -log10(ppoints(n_snps_ci)),
    ci_lower = -log10(qbeta(0.975, 1:n_snps_ci, n_snps_ci:1)),
    ci_upper = -log10(qbeta(0.025, 1:n_snps_ci, n_snps_ci:1))
)

p_qq <- ggplot() +
    geom_ribbon(data = ci_data,
                aes(x = expected, ymin = ci_lower, ymax = ci_upper),
                fill = "grey80", alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    # Non-significant points
    geom_point(data = qq_data %>% filter(!is_significant),
               aes(x = expected, y = observed, color = trait, shape = method),
               alpha = 0.5, size = 1.5) +
    # Significant points (larger)
    geom_point(data = qq_data %>% filter(is_significant),
               aes(x = expected, y = observed, color = trait, shape = method),
               size = 2.5) +
    scale_color_manual(values = trait_colors, name = "Trait") +
    scale_shape_manual(values = method_shapes, name = "Method") +
    coord_fixed() +
    labs(
        title = "Combined QQ Plot (All Traits & Methods)",
        subtitle = paste0("K = ", Kbest, " | ", length(names(assoc_info)), " methods, ",
                         length(traits), " traits"),
        x = expression(Expected ~ -log[10](p-value)),
        y = expression(Observed ~ -log[10](p-value))
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.5, "cm")
    ) +
    guides(
        color = guide_legend(order = 1),
        shape = guide_legend(order = 2)
    )

qq_base <- paste0("qq_combined_K", Kbest)
ggsave(file.path(PLOT_DIR, paste0(qq_base, ".png")), p_qq,
       width = 7, height = 7, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(qq_base, ".svg")), p_qq,
       width = 7, height = 7)
message(paste0('INFO: Saved combined QQ plot: ', qq_base))

message('INFO: Combined Manhattan plot complete')
