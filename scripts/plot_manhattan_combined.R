#!/usr/bin/env Rscript
# Combined Manhattan plot for all traits and methods
# Shows all GWAS results on a single plot with:
# - Different colors for different traits
# - Different shapes for different methods
# - Regions marked with rectangles (from Regions_climate_combined.tsv)
# - Always uses scattermore for performance

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
        plot.margin = margin(10, 15, 10, 10)
    )
}

# Generate color palette for trait-method combinations
get_trait_method_colors <- function(traits, methods) {
    # Define distinct colors for each trait-method combination
    # Use colorblind-friendly palette
    all_combos <- expand.grid(trait = traits, method = methods, stringsAsFactors = FALSE)
    all_combos$combo <- paste0(all_combos$trait, "_", all_combos$method)

    # Assign colors - using a mix of warm and cool colors
    if (nrow(all_combos) <= 4) {
        combo_colors <- c(
            "#E41A1C",  # Red
            "#FF7F00",  # Orange
            "#377EB8",  # Blue
            "#4DAF4A"   # Green
        )
    } else {
        combo_colors <- scales::hue_pal()(nrow(all_combos))
    }

    names(combo_colors) <- all_combos$combo[1:length(combo_colors)]

    # Add special color for multi-association
    combo_colors["Multi"] <- "#000000"  # Black for multiple associations

    return(combo_colors)
}

# Generate color palette for chromosomes (background) - light pastel colors
get_chr_colors <- function(n_chr) {
    # Light pastel alternating colors for subtle but visible background
    rep(c("#B3D9FF", "#E6F2FF"), length.out = n_chr)
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

# Get colors
trait_method_colors <- get_trait_method_colors(traits, names(assoc_info))
chr_colors <- get_chr_colors(length(unique(chr_info$chr_f)))
names(chr_colors) <- levels(chr_info$chr_f)

# Calculate y-axis limit
y_max <- max(plot_data_all$log10p, na.rm = TRUE) * 1.1

message('INFO: Generating combined Manhattan plot (simple version)')

# Add trait-method combination column
plot_data_all <- plot_data_all %>%
    dplyr::mutate(trait_method = paste0(trait, "_", method))

# Simple plot - background SNPs only (non-significant)
df_background <- plot_data_all %>%
    dplyr::filter(!is_significant) %>%
    # Take only one method per trait to avoid overplotting background
    dplyr::filter(method == names(assoc_info)[1])

# Create chromosome background colors for df_background
df_background$point_color <- chr_colors[as.character(df_background$chr)]

# Significant SNPs - determine if multi-association
# For each unique SNP position, count how many trait-method combinations it's significant in
sig_snp_summary <- plot_data_all %>%
    dplyr::filter(is_significant) %>%
    group_by(SNPID, chr, pos, pos_cum) %>%
    summarise(
        n_associations = n(),
        trait_methods = paste(trait_method, collapse = ","),
        max_log10p = max(log10p),
        .groups = 'drop'
    ) %>%
    dplyr::mutate(
        display_color = ifelse(n_associations > 1, "Multi",
                              strsplit(trait_methods, ",")[[1]][1])
    )

# Join back to get colors
df_significant <- sig_snp_summary %>%
    dplyr::mutate(
        color_group = display_color,
        log10p = max_log10p
    )

n_sig_snps <- nrow(df_significant)
n_multi <- sum(df_significant$color_group == "Multi")
message(paste0('INFO: Total unique significant SNPs: ', n_sig_snps))
message(paste0('INFO: Multi-association SNPs: ', n_multi))

# Simple plot
p_simple <- ggplot() +
    # Background SNPs (chromosome colors) - always use scattermore
    geom_scattermore(data = df_background,
                     aes(x = pos_cum, y = log10p, color = point_color),
                     alpha = 0.9, pointsize = 5,
                     pixels = c(4000, 1500), interpolate = FALSE) +
    scale_color_identity() +
    # Significant SNPs - colored by trait-method combination
    ggnewscale::new_scale_color() +
    geom_point(data = df_significant,
               aes(x = pos_cum, y = log10p, color = color_group),
               alpha = 0.9, size = 2.5, shape = 16) +
    scale_color_manual(values = trait_method_colors,
                      name = "Trait-Method",
                      breaks = names(trait_method_colors)) +
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
        subtitle = paste0("K = ", Kbest, " | ", n_sig_snps, " significant SNPs (",
                         n_multi, " multi-association)"),
        x = "Chromosome",
        y = expression(-log[10](p-value))
    ) +
    theme_manhattan() +
    theme(
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"),
        legend.key.size = unit(0.5, "cm")
    )

# Save simple plot
simple_base <- paste0("Manhattan_all_traits_combined_K", Kbest)
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

        # Mark SNPs in regions
        plot_data_regions <- plot_data_all %>%
            dplyr::mutate(
                in_region = FALSE,
                region_id = NA_character_
            )

        for (i in 1:nrow(regions)) {
            r <- regions[i, ]
            # Find significant SNPs in this region
            in_this_region <- as.character(plot_data_regions$chr) == as.character(r$chr) &
                              plot_data_regions$pos >= (r$start - 500000) &
                              plot_data_regions$pos <= (r$end + 500000) &
                              plot_data_regions$is_significant

            plot_data_regions$in_region[in_this_region] <- TRUE
            plot_data_regions$region_id[in_this_region] <- r$region_id
        }

        # Region colors
        region_ids <- unique(regions$region_id)
        region_colors <- get_region_colors(length(region_ids))
        names(region_colors) <- region_ids

        n_in_regions <- sum(plot_data_regions$in_region)
        message(paste0('INFO: ', n_in_regions, ' significant SNPs in regions'))

        # Prepare data subsets
        df_background_reg <- plot_data_regions %>%
            dplyr::filter(!is_significant) %>%
            dplyr::filter(method == names(assoc_info)[1])
        df_background_reg$point_color <- chr_colors[as.character(df_background_reg$chr)]

        # For significant SNPs, aggregate by position to handle multi-associations
        sig_snps_regions <- plot_data_regions %>%
            dplyr::filter(is_significant) %>%
            group_by(SNPID, chr, pos, pos_cum, in_region, region_id) %>%
            summarise(
                n_associations = n(),
                trait_methods = paste(trait_method, collapse = ","),
                max_log10p = max(log10p),
                .groups = 'drop'
            ) %>%
            dplyr::mutate(
                display_color = ifelse(n_associations > 1, "Multi",
                                      strsplit(trait_methods, ",")[[1]][1]),
                log10p = max_log10p
            )

        df_sig_not_region <- sig_snps_regions %>%
            dplyr::filter(!in_region)

        df_in_region <- sig_snps_regions %>%
            dplyr::filter(in_region)

        # Calculate statistics
        n_sig_unique <- nrow(sig_snps_regions)
        n_multi_regions <- sum(sig_snps_regions$display_color == "Multi")

        # Plot with regions
        message('INFO: Generating combined Manhattan plot with regions')

        p_regions <- ggplot() +
            # Region rectangles (background)
            geom_rect(data = regions,
                     aes(xmin = start_cum, xmax = end_cum,
                         ymin = 0, ymax = y_max, fill = region_id),
                     alpha = 0.15) +
            scale_fill_manual(values = region_colors, guide = "none") +
            # Background SNPs - always scattermore
            geom_scattermore(data = df_background_reg,
                             aes(x = pos_cum, y = log10p, color = point_color),
                             alpha = 0.9, pointsize = 5,
                             pixels = c(4000, 1500), interpolate = FALSE) +
            scale_color_identity() +
            # Significant SNPs NOT in regions - colored by trait-method
            ggnewscale::new_scale_color() +
            geom_point(data = df_sig_not_region,
                      aes(x = pos_cum, y = log10p, color = display_color),
                      alpha = 0.8, size = 2.0, shape = 16) +
            scale_color_manual(values = trait_method_colors, name = "Trait-Method") +
            # SNPs IN regions - colored by region (larger)
            ggnewscale::new_scale_color() +
            geom_point(data = df_in_region,
                      aes(x = pos_cum, y = log10p, color = region_id),
                      alpha = 0.95, size = 3.0, shape = 18) +
            scale_color_manual(values = region_colors, name = "Region") +
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
                subtitle = paste0("K = ", Kbest, " | ", n_sig_unique, " significant SNPs (",
                                 n_multi_regions, " multi-association) | ",
                                 nrow(regions), " regions highlighted"),
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
                fill = "none",
                color = guide_legend(ncol = 1, override.aes = list(size = 2))
            )

        # Save regions plot
        regions_base <- paste0("Manhattan_all_traits_combined_K", Kbest, "_regions")
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

message('INFO: Combined Manhattan plot complete')
