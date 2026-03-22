#!/usr/bin/env Rscript
# Static Miami plot combining GEA (top) and GWAS (bottom)
# GEA data plotted with positive -log10(p), GWAS with negative -log10(p)
# Color = trait (Okabe-Ito), Shape = method, scattermore mandatory for background

library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(qvalue)
library(scales)
library(scattermore)
library(ggnewscale)
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)
################################
GEA_FILES_STR = args[1]       # Comma-sep "METHOD:ADJUST:FILEPATH"
GWAS_FILES_STR = args[2]      # Comma-sep "METHOD:ADJUST:FILEPATH"
PREDICTORS_GEA = args[3]      # Comma-sep trait names (bio_1,bio_12)
PREDICTORS_GWAS = args[4]     # Comma-sep trait names (height,flowering_time)
Kbest = args[5] %>% as.numeric
REGIONS_FILE = args[6]        # overlap Regions_combined.tsv or "NULL"
SIGSNPS_FILE = args[7]        # overlap selected_snps_all.tsv or "NULL"
PLOT_DIR = args[8]
SCATTERMORE_THRESHOLD = args[9] %>% as.numeric
################################

# Handle optional files
if (is.na(REGIONS_FILE) || REGIONS_FILE == "NULL" || REGIONS_FILE == "") {
    REGIONS_FILE <- NULL
}
if (is.na(SIGSNPS_FILE) || SIGSNPS_FILE == "NULL" || SIGSNPS_FILE == "") {
    SIGSNPS_FILE <- NULL
}

message('INFO: Static Miami plot (GEA top, GWAS bottom)')
message(paste0('INFO: K = ', Kbest))

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
    setNames(shapes[1:length(methods)], methods)
}

# Create cumulative position for Manhattan plot with inter-chromosome gaps
prepare_manhattan_data <- function(df, chr_col = "chr", pos_col = "pos", pval_col) {
    chr_order <- df[[chr_col]] %>% unique() %>%
        .[order(as.numeric(str_extract(., "\\d+")))]

    df <- df %>%
        dplyr::mutate(chr_f = factor(.data[[chr_col]], levels = chr_order))

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

# Manhattan plot theme adapted for Miami
theme_miami <- function() {
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

################################ Parse input files

parse_assoc_str <- function(files_str) {
    assoc_list <- str_split(files_str, ',')[[1]]
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
    return(assoc_info)
}

gea_info <- parse_assoc_str(GEA_FILES_STR)
gwas_info <- parse_assoc_str(GWAS_FILES_STR)

gea_traits <- str_split(PREDICTORS_GEA, ',')[[1]] %>% str_trim()
gwas_traits <- str_split(PREDICTORS_GWAS, ',')[[1]] %>% str_trim()

message(paste0('INFO: GEA methods: ', paste(names(gea_info), collapse = ', ')))
message(paste0('INFO: GWAS methods: ', paste(names(gwas_info), collapse = ', ')))
message(paste0('INFO: GEA traits: ', paste(gea_traits, collapse = ', ')))
message(paste0('INFO: GWAS traits: ', paste(gwas_traits, collapse = ', ')))

################################ Load data

load_assoc_data <- function(assoc_info, traits, panel_label) {
    all_data <- list()
    method_thresholds <- list()

    for (method_name in names(assoc_info)) {
        info <- assoc_info[[method_name]]
        message(paste0('INFO: Loading ', panel_label, ' ', method_name, ' from ', info$filepath))

        snps_assoc <- fread(info$filepath, sep = '\t', header = TRUE)
        snps_assoc$chr <- as.character(snps_assoc$chr)

        adjustment <- str_split(info$adjust, '_')[[1]][1]
        pval_threshold <- str_split(info$adjust, '_')[[1]][2] %>% as.numeric()

        for (trait in traits) {
            if (!(trait %in% colnames(snps_assoc))) {
                message(paste0('WARNING: Trait ', trait, ' not found in ', method_name))
                next
            }

            original_threshold <- pval_threshold
            threshold_value <- pval_threshold

            if (adjustment == 'bonf') {
                threshold_value <- original_threshold / nrow(snps_assoc)
            }
            if (adjustment == 'qval') {
                threshold_value <- FUN_max_pvalue_fdr(snps_assoc[[trait]], original_threshold)
                if (is.null(threshold_value)) {
                    threshold_value <- original_threshold / nrow(snps_assoc)
                }
            }
            if (adjustment == 'top') {
                threshold_value <- FUN_max_pvalue_top(snps_assoc[[trait]], original_threshold)
            }

            threshold_log10 <- -log10(threshold_value)
            method_thresholds[[paste0(method_name, "_", trait)]] <- threshold_log10

            trait_data <- snps_assoc %>%
                dplyr::select(SNPID, chr, pos, !!sym(trait)) %>%
                dplyr::rename(pvalue = !!sym(trait)) %>%
                dplyr::mutate(
                    method = method_name,
                    trait = trait,
                    log10p = -log10(pvalue),
                    is_significant = log10p >= threshold_log10,
                    panel = panel_label
                )

            all_data[[paste0(method_name, "_", trait)]] <- trait_data

            n_sig <- sum(trait_data$is_significant)
            message(paste0('INFO: ', panel_label, ' ', method_name, ' - ', trait, ': ', n_sig, ' significant SNPs'))
        }
    }

    list(data = bind_rows(all_data), thresholds = method_thresholds)
}

gea_result <- load_assoc_data(gea_info, gea_traits, "GEA")
gwas_result <- load_assoc_data(gwas_info, gwas_traits, "GWAS")

gea_data <- gea_result$data
gwas_data <- gwas_result$data

# Calculate thresholds (min per panel)
gea_threshold <- min(unlist(gea_result$thresholds))
gwas_threshold <- min(unlist(gwas_result$thresholds))
message(paste0('INFO: GEA threshold: -log10(p) = ', round(gea_threshold, 2)))
message(paste0('INFO: GWAS threshold: -log10(p) = ', round(gwas_threshold, 2)))

################################ Prepare positions

# Use GEA data for chromosome coordinate system
reference_data <- gea_data %>%
    dplyr::filter(trait == gea_traits[1], method == names(gea_info)[1]) %>%
    dplyr::select(SNPID, chr, pos, pvalue)

manhattan_prep <- prepare_manhattan_data(reference_data, "chr", "pos", "pvalue")
chr_info <- manhattan_prep$chr_info

# Add cumulative positions to both panels
add_cum_pos <- function(df) {
    df %>%
        left_join(chr_info %>%
                  dplyr::mutate(chr = as.character(chr_f)) %>%
                  dplyr::select(chr, tot),
                  by = "chr") %>%
        dplyr::mutate(pos_cum = pos + tot)
}

gea_data <- add_cum_pos(gea_data)
gwas_data <- add_cum_pos(gwas_data)

# For Miami: GEA positive, GWAS negative
gea_data <- gea_data %>% dplyr::mutate(miami_y = log10p)
gwas_data <- gwas_data %>% dplyr::mutate(miami_y = -log10p)

################################ Palettes

# All traits (union of GEA + GWAS) for consistent colors
all_traits <- unique(c(gea_traits, gwas_traits))
trait_colors <- get_trait_colors(all_traits)

# All methods (union) for consistent shapes
all_methods <- unique(c(names(gea_info), names(gwas_info)))
method_shapes <- get_method_shapes(all_methods)

################################ Y-axis limits

y_max_gea <- max(gea_data$log10p, na.rm = TRUE) * 1.1
y_max_gwas <- max(gwas_data$log10p, na.rm = TRUE) * 1.1
y_limit <- max(y_max_gea, y_max_gwas)

################################ Build simple Miami plot

message('INFO: Generating simple Miami plot')

# Background SNPs (non-sig, first method only per panel)
gea_bg <- gea_data %>%
    dplyr::filter(!is_significant) %>%
    dplyr::filter(method == names(gea_info)[1])

gwas_bg <- gwas_data %>%
    dplyr::filter(!is_significant) %>%
    dplyr::filter(method == names(gwas_info)[1])

# Sig SNPs
gea_sig <- gea_data %>% dplyr::filter(is_significant)
gwas_sig <- gwas_data %>% dplyr::filter(is_significant)

n_gea_sig <- length(unique(gea_sig$SNPID))
n_gwas_sig <- length(unique(gwas_sig$SNPID))
message(paste0('INFO: GEA significant: ', n_gea_sig, ' unique SNPs, ', nrow(gea_sig), ' detections'))
message(paste0('INFO: GWAS significant: ', n_gwas_sig, ' unique SNPs, ', nrow(gwas_sig), ' detections'))

# Simple Miami — per-trait background scattermore (chr gaps for separation)
p_miami <- ggplot()

# 2. GEA background per trait
first_gea_method <- names(gea_info)[1]
for (t in gea_traits) {
    df_t <- gea_bg %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_miami <- p_miami +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = miami_y, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 2400), interpolate = FALSE)
}

# 3. GWAS background per trait
first_gwas_method <- names(gwas_info)[1]
for (t in gwas_traits) {
    df_t <- gwas_bg %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_miami <- p_miami +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = miami_y, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 2400), interpolate = FALSE)
}

p_miami <- p_miami +
    scale_color_identity() +
    # 4. GEA sig SNPs
    ggnewscale::new_scale_color() +
    geom_point(data = gea_sig,
               aes(x = pos_cum, y = miami_y, color = trait, shape = method),
               alpha = 0.9, size = 2.5, stroke = 0.5) +
    # GWAS sig SNPs
    geom_point(data = gwas_sig,
               aes(x = pos_cum, y = miami_y, color = trait, shape = method),
               alpha = 0.9, size = 2.5, stroke = 0.5) +
    scale_color_manual(values = trait_colors, name = "Trait") +
    scale_shape_manual(values = method_shapes, name = "Method") +
    # Threshold lines
    geom_hline(yintercept = gea_threshold, linetype = "dashed",
               color = "red", linewidth = 0.4, alpha = 0.5) +
    geom_hline(yintercept = -gwas_threshold, linetype = "dashed",
               color = "red", linewidth = 0.4, alpha = 0.5) +
    # Bold zero line separating GEA and GWAS
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
    # Axes
    scale_x_continuous(
        labels = chr_info$chr_f,
        breaks = chr_info$center,
        expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
        limits = c(-y_limit, y_limit),
        labels = function(x) abs(x),
        expand = expansion(mult = c(0.02, 0.02))
    ) +
    # Labels
    labs(
        title = "Miami Plot: GEA (top) vs GWAS (bottom)",
        subtitle = paste0("K = ", Kbest, " | GEA: ", n_gea_sig, " sig SNPs | GWAS: ", n_gwas_sig, " sig SNPs"),
        x = "Chromosome",
        y = expression(-log[10](p-value))
    ) +
    # Annotate panel labels
    annotate("text", x = -Inf, y = y_limit * 0.95, label = "GEA",
             hjust = -0.3, fontface = "bold", size = 4, color = "grey30") +
    annotate("text", x = -Inf, y = -y_limit * 0.95, label = "GWAS",
             hjust = -0.3, fontface = "bold", size = 4, color = "grey30") +
    theme_miami() +
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

# Save simple Miami
simple_base <- paste0("miami_combined_K", Kbest)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".png")), p_miami,
       width = 14, height = 8, dpi = 300)
ggsave(file.path(PLOT_DIR, paste0(simple_base, ".svg")), p_miami,
       width = 14, height = 8)
message(paste0('INFO: Saved simple Miami plot: ', simple_base))

################################ Background PNG + Coords JSON (for Shiny plotly overlay)

message('INFO: Generating Miami background for Shiny overlay')

all_pos_cum <- c(gea_data$pos_cum, gwas_data$pos_cum)
x_min_data  <- min(all_pos_cum)
x_max_data  <- max(all_pos_cum)
x_span      <- x_max_data - x_min_data
bg_x_lo     <- x_min_data - 0.01 * x_span
bg_x_hi     <- x_max_data + 0.01 * x_span

# y: limits = c(-y_limit, y_limit), expand = mult c(0.02, 0.02)
bg_y_lo <- -y_limit * 1.04
bg_y_hi <-  y_limit * 1.04

p_bg_miami <- ggplot()

for (t in gea_traits) {
    df_t <- gea_bg %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_bg_miami <- p_bg_miami +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = miami_y, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 2400), interpolate = FALSE)
}

for (t in gwas_traits) {
    df_t <- gwas_bg %>% dplyr::filter(trait == t)
    df_t$point_color <- unname(trait_colors[t])
    p_bg_miami <- p_bg_miami +
        geom_scattermore(data = df_t,
                         aes(x = pos_cum, y = miami_y, color = point_color),
                         alpha = 0.5, pointsize = 10,
                         pixels = c(4000, 2400), interpolate = FALSE)
}

p_bg_miami <- p_bg_miami +
    scale_color_identity() +
    geom_hline(yintercept =  gea_threshold,  linetype = "dashed",
               color = "red", linewidth = 0.4, alpha = 0.5) +
    geom_hline(yintercept = -gwas_threshold, linetype = "dashed",
               color = "red", linewidth = 0.4, alpha = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
    scale_x_continuous(limits = c(bg_x_lo, bg_x_hi), expand = c(0, 0)) +
    scale_y_continuous(limits = c(bg_y_lo, bg_y_hi), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

bg_miami_base <- paste0("miami_combined_K", Kbest, "_background")
ggsave(file.path(PLOT_DIR, paste0(bg_miami_base, ".png")), p_bg_miami,
       width = 14, height = 8, dpi = 300)
message(paste0('INFO: Saved Miami background: ', bg_miami_base, '.png'))

coords_miami <- list(
    chr_offsets      = setNames(as.list(chr_info$tot),     as.character(chr_info$chr_f)),
    chr_lengths      = setNames(as.list(chr_info$chr_len), as.character(chr_info$chr_f)),
    gap_fraction     = 0.02,
    x_range          = c(bg_x_lo, bg_x_hi),
    y_range          = c(bg_y_lo, bg_y_hi),
    gea_threshold_y  = gea_threshold,
    gwas_threshold_y = -gwas_threshold,
    is_miami         = TRUE,
    plot_width_px    = 4200L,
    plot_height_px   = 2400L
)
coords_miami_file <- file.path(PLOT_DIR, paste0("miami_combined_K", Kbest, "_coords.json"))
jsonlite::write_json(coords_miami, coords_miami_file, auto_unbox = TRUE, digits = 6)
message(paste0('INFO: Saved Miami coords JSON: ', basename(coords_miami_file)))

################################ Miami plot with Regions

if (!is.null(REGIONS_FILE) && file.exists(REGIONS_FILE)) {
    message('INFO: Loading overlap regions')
    regions <- fread(REGIONS_FILE)
    if ('chr' %in% colnames(regions)) regions$chr <- as.character(regions$chr)

    if (nrow(regions) > 0) {
        message(paste0('INFO: Found ', nrow(regions), ' combined regions'))

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

        # Minimum visible width
        total_span <- max(c(gea_data$pos_cum, gwas_data$pos_cum)) - min(c(gea_data$pos_cum, gwas_data$pos_cum))
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

        # Sig SNPs (all, colored by trait — regions shown as rectangles only)
        gea_sig_r <- gea_data %>% dplyr::filter(is_significant)
        gwas_sig_r <- gwas_data %>% dplyr::filter(is_significant)

        # Regions Miami — per-trait background + region rectangles (chr gaps for separation)
        p_miami_r <- ggplot() +
            # 1. Region rectangles (with legend)
            geom_rect(data = regions,
                     aes(xmin = start_cum, xmax = end_cum,
                         ymin = -y_limit, ymax = y_limit, fill = region_id),
                     alpha = 0.15) +
            scale_fill_manual(values = region_colors, name = "Region")

        # 3. GEA background per trait
        for (t in gea_traits) {
            df_t <- gea_bg %>% dplyr::filter(trait == t)
            df_t$point_color <- unname(trait_colors[t])
            p_miami_r <- p_miami_r +
                geom_scattermore(data = df_t,
                                 aes(x = pos_cum, y = miami_y, color = point_color),
                                 alpha = 0.5, pointsize = 10,
                                 pixels = c(4000, 2400), interpolate = FALSE)
        }

        # 4. GWAS background per trait
        for (t in gwas_traits) {
            df_t <- gwas_bg %>% dplyr::filter(trait == t)
            df_t$point_color <- unname(trait_colors[t])
            p_miami_r <- p_miami_r +
                geom_scattermore(data = df_t,
                                 aes(x = pos_cum, y = miami_y, color = point_color),
                                 alpha = 0.5, pointsize = 10,
                                 pixels = c(4000, 2400), interpolate = FALSE)
        }

        p_miami_r <- p_miami_r +
            scale_color_identity() +
            # 5. All sig SNPs - color by trait, shape by method
            ggnewscale::new_scale_color() +
            geom_point(data = bind_rows(gea_sig_r, gwas_sig_r),
                      aes(x = pos_cum, y = miami_y, color = trait, shape = method),
                      alpha = 0.9, size = 2.5, stroke = 0.5) +
            scale_color_manual(values = trait_colors, name = "Trait") +
            scale_shape_manual(values = method_shapes, name = "Method") +
            # Threshold lines
            geom_hline(yintercept = gea_threshold, linetype = "dashed",
                       color = "red", linewidth = 0.4, alpha = 0.5) +
            geom_hline(yintercept = -gwas_threshold, linetype = "dashed",
                       color = "red", linewidth = 0.4, alpha = 0.5) +
            # Bold zero line
            geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
            # Axes
            scale_x_continuous(
                labels = chr_info$chr_f,
                breaks = chr_info$center,
                expand = c(0.01, 0.01)
            ) +
            scale_y_continuous(
                limits = c(-y_limit, y_limit),
                labels = function(x) abs(x),
                expand = expansion(mult = c(0.02, 0.02))
            ) +
            labs(
                title = "Miami Plot with Regions: GEA (top) vs GWAS (bottom)",
                subtitle = paste0("K = ", Kbest, " | ", length(region_ids), " regions | color=trait, shape=method"),
                x = "Chromosome",
                y = expression(-log[10](p-value))
            ) +
            annotate("text", x = -Inf, y = y_limit * 0.95, label = "GEA",
                     hjust = -0.3, fontface = "bold", size = 4, color = "grey30") +
            annotate("text", x = -Inf, y = -y_limit * 0.95, label = "GWAS",
                     hjust = -0.3, fontface = "bold", size = 4, color = "grey30") +
            theme_miami() +
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

        # Save regions Miami
        regions_base <- paste0("miami_combined_K", Kbest, "_regions")
        ggsave(file.path(PLOT_DIR, paste0(regions_base, ".png")), p_miami_r,
               width = 14, height = 8, dpi = 300)
        ggsave(file.path(PLOT_DIR, paste0(regions_base, ".svg")), p_miami_r,
               width = 14, height = 8)
        message(paste0('INFO: Saved Miami regions plot: ', regions_base))

        # Background regions version for Shiny (no sig SNPs, no axes)
        p_bg_miami_r <- ggplot() +
            geom_rect(data = regions,
                     aes(xmin = start_cum, xmax = end_cum,
                         ymin = bg_y_lo, ymax = bg_y_hi, fill = region_id),
                     alpha = 0.15) +
            scale_fill_manual(values = region_colors, guide = "none")

        for (t in gea_traits) {
            df_t <- gea_bg %>% dplyr::filter(trait == t)
            df_t$point_color <- unname(trait_colors[t])
            p_bg_miami_r <- p_bg_miami_r +
                geom_scattermore(data = df_t,
                                 aes(x = pos_cum, y = miami_y, color = point_color),
                                 alpha = 0.5, pointsize = 10,
                                 pixels = c(4000, 2400), interpolate = FALSE)
        }

        for (t in gwas_traits) {
            df_t <- gwas_bg %>% dplyr::filter(trait == t)
            df_t$point_color <- unname(trait_colors[t])
            p_bg_miami_r <- p_bg_miami_r +
                geom_scattermore(data = df_t,
                                 aes(x = pos_cum, y = miami_y, color = point_color),
                                 alpha = 0.5, pointsize = 10,
                                 pixels = c(4000, 2400), interpolate = FALSE)
        }

        p_bg_miami_r <- p_bg_miami_r +
            scale_color_identity() +
            geom_hline(yintercept =  gea_threshold,  linetype = "dashed",
                       color = "red", linewidth = 0.4, alpha = 0.5) +
            geom_hline(yintercept = -gwas_threshold, linetype = "dashed",
                       color = "red", linewidth = 0.4, alpha = 0.5) +
            geom_hline(yintercept = 0, linewidth = 0.8, color = "black") +
            scale_x_continuous(limits = c(bg_x_lo, bg_x_hi), expand = c(0, 0)) +
            scale_y_continuous(limits = c(bg_y_lo, bg_y_hi), expand = c(0, 0)) +
            theme_void() +
            theme(plot.background = element_rect(fill = "white", color = NA))

        bg_miami_r_base <- paste0("miami_combined_K", Kbest, "_regions_background")
        ggsave(file.path(PLOT_DIR, paste0(bg_miami_r_base, ".png")), p_bg_miami_r,
               width = 14, height = 8, dpi = 300)
        message(paste0('INFO: Saved Miami background regions: ', bg_miami_r_base, '.png'))

    } else {
        message('INFO: Regions file is empty')
    }
} else {
    message('INFO: No regions file provided or not found')
}

# Ensure _regions_background.png exists (may be skipped if regions file is NULL/empty)
if (!is.null(REGIONS_FILE)) {
    expected_bg_miami_r <- file.path(PLOT_DIR,
        paste0("miami_combined_K", Kbest, "_regions_background.png"))
    if (!file.exists(expected_bg_miami_r)) {
        file.copy(file.path(PLOT_DIR, paste0(bg_miami_base, ".png")), expected_bg_miami_r)
        message('INFO: No regions - copied simple Miami background as regions background')
    }
}

message('INFO: Miami plot complete')
