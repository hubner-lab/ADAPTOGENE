#!/usr/bin/env Rscript
# Publication-ready GEA Manhattan plot — all bioclimatic traits on one panel
# Temperature traits in warm palette, precipitation traits in cool palette
# Non-significant SNPs rendered with scattermore; significant SNPs with geom_point
# Highlighted article regions as colored boxes

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(scattermore)
library(ggnewscale)

################################ Configurable parameters

PVALUES_FILE <- "WBDC2_results/tables/association/EMMAX/EMMAX_pvalues_K5.tsv"
OUTPUT_FILE  <- "WBDC2_results/plots/Manhattan_GEA_all_bioclim.svg"

# Climate traits only (12 bioclim variables)
CLIMATE_TRAITS <- c("bio_1","bio_2","bio_3","bio_4","bio_5","bio_6",
                     "bio_8","bio_9","bio_10","bio_13","bio_14","bio_15")

# Human-readable labels for legend
TRAIT_LABELS <- c(
  bio_1  = "BIO1 Annual Mean Temp",
  bio_2  = "BIO2 Mean Diurnal Range",
  bio_3  = "BIO3 Isothermality",
  bio_4  = "BIO4 Temp Seasonality",
  bio_5  = "BIO5 Max Temp Warmest Mo.",
  bio_6  = "BIO6 Min Temp Coldest Mo.",
  bio_8  = "BIO8 Mean Temp Wettest Q.",
  bio_9  = "BIO9 Mean Temp Driest Q.",
  bio_10 = "BIO10 Mean Temp Warmest Q.",
  bio_13 = "BIO13 Precip Wettest Mo.",
  bio_14 = "BIO14 Precip Driest Mo.",
  bio_15 = "BIO15 Precip Seasonality"
)

# Temperature traits — warm palette (9 hand-picked distinguishable warm tones)
TEMP_TRAITS  <- c("bio_1","bio_2","bio_3","bio_4","bio_5","bio_6","bio_8","bio_9","bio_10")
TEMP_COLORS  <- c(
  bio_1  = "#B71C1C",  # deep red
  bio_2  = "#E53935",  # red
  bio_3  = "#FF6F00",  # amber
  bio_4  = "#FF8F00",  # dark amber
  bio_5  = "#F57C00",  # orange
  bio_6  = "#EF6C00",  # deep orange
  bio_8  = "#FFB300",  # gold
  bio_9  = "#FDD835",  # yellow
  bio_10 = "#D84315"   # burnt orange
)

# Precipitation traits — cool palette (3 blues/teals)
PRECIP_TRAITS  <- c("bio_13","bio_14","bio_15")
PRECIP_COLORS  <- c(
  bio_13 = "#1565C0",  # dark blue
  bio_14 = "#0097A7",  # teal
  bio_15 = "#42A5F5"   # light blue
)

ALL_COLORS <- c(TEMP_COLORS, PRECIP_COLORS)

# Plot dimensions (single-column figure)
PLOT_WIDTH_MM  <- 180
PLOT_HEIGHT_MM <- 80

# Significance: Bonferroni
BONF_ALPHA <- 0.05

################################ Article regions to highlight

# Define as data.table: chr, start, end, short label for legend
REGIONS <- data.table(
  chr   = c("3H",         "7H",         "1H",         "4H",         "2H",         "6H",         "7H"),
  start = c(536121931,    418230754,    50809274,     567646849,    163488648,    349559726,    390121100),
  end   = c(538306268,    425997112,    55247125,     569646859,    169726331,    350668888,    395267626),
  label = c("3H:536-538", "7H:418-426", "1H:51-55",  "4H:568-570", "2H:163-170", "6H:349-351", "7H:390-395")
)

REGION_COLORS <- c(
  "3H:536-538" = "#E41A1C",  # red
  "7H:418-426" = "#377EB8",  # blue
  "1H:51-55"   = "#4DAF4A",  # green
  "4H:568-570" = "#984EA3",  # purple
  "2H:163-170" = "#FF7F00",  # orange
  "6H:349-351" = "#A6CEE3",  # light blue
  "7H:390-395" = "#B2DF8A"   # light green
)

################################ Load data

message("INFO: Loading p-values from ", PVALUES_FILE)
cols_needed <- c("SNPID", "chr", "pos", CLIMATE_TRAITS)
pval <- fread(PVALUES_FILE, select = cols_needed, sep = "\t", header = TRUE)
pval[, chr := as.character(chr)]
n_snps <- nrow(pval)
message("INFO: Loaded ", n_snps, " SNPs x ", length(CLIMATE_TRAITS), " traits")

################################ Cumulative chromosome positions

chr_order <- unique(pval$chr)[order(as.numeric(str_extract(unique(pval$chr), "\\d+")))]
pval[, chr_f := factor(chr, levels = chr_order)]

chr_info <- pval[, .(chr_len = max(pos)), by = chr_f]
chr_info[, tot := cumsum(as.numeric(chr_len)) - chr_len]
chr_info[, center := tot + chr_len / 2]

pval <- merge(pval, chr_info[, .(chr_f, tot)], by = "chr_f", all.x = TRUE)
pval[, pos_cum := pos + tot]

################################ Melt to long format

message("INFO: Melting to long format")
dt_long <- melt(pval,
                id.vars = c("SNPID", "chr", "chr_f", "pos", "pos_cum", "tot"),
                measure.vars = CLIMATE_TRAITS,
                variable.name = "trait",
                value.name = "pvalue")
dt_long[, trait := as.character(trait)]

# -log10 transform
dt_long[, log10p := -log10(pvalue)]

# Bonferroni threshold
bonf_threshold <- BONF_ALPHA / n_snps
threshold_log10 <- -log10(bonf_threshold)
message("INFO: Bonferroni threshold: ", signif(bonf_threshold, 4),
        " (-log10 = ", round(threshold_log10, 2), ")")

# Flag significant
dt_long[, is_sig := log10p >= threshold_log10]
n_sig <- sum(dt_long$is_sig)
message("INFO: ", n_sig, " significant trait-SNP pairs across all traits")

################################ Split sig / non-sig

df_nonsig <- dt_long[is_sig == FALSE]
df_sig    <- dt_long[is_sig == TRUE]

message("INFO: Non-significant points: ", nrow(df_nonsig))
message("INFO: Significant points: ", nrow(df_sig))

# Pre-compute hex colors for non-sig (scattermore needs pre-assigned colors)
# Vectorized: compute alpha-adjusted color per trait (12 values), then lookup
NONSIG_ALPHA_COLORS <- vapply(ALL_COLORS, function(col) {
  adjustcolor(col, alpha.f = 0.35)
}, character(1))
df_nonsig[, hex_color := NONSIG_ALPHA_COLORS[trait]]

# Split non-sig into temp and precip for layering order
# (temp drawn first = bottom, precip drawn second = on top)
df_nonsig_temp   <- df_nonsig[trait %in% TEMP_TRAITS]
df_nonsig_precip <- df_nonsig[trait %in% PRECIP_TRAITS]
message("INFO: Non-sig temp: ", nrow(df_nonsig_temp),
        ", precip: ", nrow(df_nonsig_precip))

################################ Prepare regions with cumulative positions

chr_tot_lookup <- chr_info[, .(chr = as.character(chr_f), tot)]
REGIONS[, chr := as.character(chr)]
REGIONS <- merge(REGIONS, chr_tot_lookup, by = "chr", all.x = TRUE)
REGIONS[, start_cum := start + tot]
REGIONS[, end_cum   := end   + tot]

# Ensure minimum visible width (1.5% of genome span)
total_span <- max(pval$pos_cum) - min(pval$pos_cum)
min_half_width <- total_span * 0.008
REGIONS[, `:=`(
  region_center     = (start_cum + end_cum) / 2,
  region_half_width = pmax((end_cum - start_cum) / 2, min_half_width)
)]
REGIONS[, start_cum := region_center - region_half_width]
REGIONS[, end_cum   := region_center + region_half_width]

message("INFO: ", nrow(REGIONS), " article regions prepared")

################################ Build plot

message("INFO: Building publication Manhattan plot")

# Y-axis range
y_max <- max(dt_long$log10p, na.rm = TRUE) * 1.08
y_max <- max(y_max, threshold_log10 * 1.15)

# Trait factor for legend ordering (temp first, then precip)
trait_order <- c(TEMP_TRAITS, PRECIP_TRAITS)
df_sig[, trait_f := factor(trait, levels = trait_order, labels = TRAIT_LABELS[trait_order])]

# Region factor for legend ordering
REGIONS[, label_f := factor(label, levels = REGIONS$label)]

p <- ggplot() +

  # Layer 0: Region highlight boxes (behind everything)
  geom_rect(
    data = REGIONS,
    aes(xmin = start_cum, xmax = end_cum, ymin = 0, ymax = y_max, fill = label_f),
    alpha = 0.15
  ) +
  scale_fill_manual(
    values = REGION_COLORS,
    name = "Region"
  ) +

  # Layer 1a: non-significant TEMPERATURE SNPs (bottom raster)
  geom_scattermore(
    data = df_nonsig_temp,
    aes(x = pos_cum, y = log10p, color = hex_color),
    pointsize = 3,
    pixels = c(3200, 1200),
    interpolate = FALSE
  ) +
  scale_color_identity() +

  # Layer 1b: non-significant PRECIPITATION SNPs (on top of temp)
  ggnewscale::new_scale_color() +
  geom_scattermore(
    data = df_nonsig_precip,
    aes(x = pos_cum, y = log10p, color = hex_color),
    pointsize = 3,
    pixels = c(3200, 1200),
    interpolate = FALSE
  ) +
  scale_color_identity() +

  # Layer 2: significant SNPs — standard geom_point with black outline
  ggnewscale::new_scale_color() +
  geom_point(
    data = df_sig,
    aes(x = pos_cum, y = log10p, color = trait_f),
    shape = 21, fill = NA, size = 1.0, stroke = 0.15, alpha = 0.9
  ) +
  # Filled layer on top
  geom_point(
    data = df_sig,
    aes(x = pos_cum, y = log10p, color = trait_f),
    shape = 16, size = 0.8, alpha = 0.9
  ) +
  scale_color_manual(
    values = setNames(ALL_COLORS[trait_order], TRAIT_LABELS[trait_order]),
    name = NULL
  ) +

  # Bonferroni threshold line
  geom_hline(yintercept = threshold_log10, linetype = "dashed",
             color = "#D32F2F", linewidth = 0.35) +

  # X-axis: chromosome labels
  scale_x_continuous(
    labels = chr_info$chr_f,
    breaks = chr_info$center,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05)),
    limits = c(0, y_max)
  ) +

  # Labels (no title — publication figures use captions)
  labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +

  # Publication theme
  theme_minimal(base_family = "sans") +
  theme(
    # Grid
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
    # Axes
    axis.line.x  = element_line(color = "grey30", linewidth = 0.3),
    axis.line.y  = element_line(color = "grey30", linewidth = 0.3),
    axis.ticks.x = element_line(color = "grey30", linewidth = 0.3),
    axis.ticks.y = element_line(color = "grey30", linewidth = 0.3),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title.x = element_text(size = 10, margin = margin(t = 4)),
    axis.title.y = element_text(size = 10, margin = margin(r = 4)),
    # Legend
    legend.position   = "right",
    legend.text        = element_text(size = 6.5),
    legend.title       = element_text(size = 7.5, face = "bold"),
    legend.key.size    = unit(3, "mm"),
    legend.spacing.y   = unit(0.5, "mm"),
    legend.margin      = margin(0, 0, 0, 2),
    legend.box         = "vertical",
    legend.box.spacing = unit(2, "mm"),
    # Margins
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(
    fill  = guide_legend(order = 1, ncol = 1, override.aes = list(alpha = 0.3)),
    color = guide_legend(order = 2, ncol = 1, override.aes = list(size = 2.5))
  )

################################ Save

message("INFO: Saving SVG to ", OUTPUT_FILE)
dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)
ggsave(OUTPUT_FILE, p,
       width = PLOT_WIDTH_MM, height = PLOT_HEIGHT_MM, units = "mm",
       device = "svg")

# Also save PNG for quick preview
png_file <- sub("\\.svg$", ".png", OUTPUT_FILE)
message("INFO: Saving PNG to ", png_file)
ggsave(png_file, p,
       width = PLOT_WIDTH_MM, height = PLOT_HEIGHT_MM, units = "mm",
       dpi = 300, device = "png")

message("INFO: Done. ", n_sig, " significant points plotted across ",
        length(unique(df_sig$trait)), " traits.")
