#!/usr/bin/env Rscript
# Unified PieMap plotting script
# Handles both climate-based piemaps and genetic offset piemaps
library(qs)
library(data.table)
library(dplyr)
library(see)
library(ggplot2)
library(ggrepel)
library(raster)
library(ggspatial)
library(scatterpie)
library(ggnewscale)
library(svglite)

args = commandArgs(trailingOnly=TRUE)
#################################
RASTER_PATH = args[1]       # path to .grd file
RASTER_LAYER = args[2]      # layer name (e.g., "bio_1") or "1" for single-layer rasters
RASTER_LEGEND = args[3]     # legend title for raster (e.g., "bio_1" or "Genetic Offset")
SAMPLES = args[4]           # metadata file
CLUSTERS = args[5]          # clusters file
SIZE_TRAIT = args[6]        # trait file for pie sizing, or "NULL" for uniform
SIZE_TRAIT_NAME = args[7]   # name for size legend, or "NULL"
PALETTE = args[8] %>% as.numeric
PALETTE_REV = args[9] %>% as.logical
PIE_ALPHA = args[10] %>% as.numeric
POP_LABEL = args[11]
POP_LABEL_SIZE = args[12] %>% as.numeric
PLOT_DIR = args[13]
INTER_DIR = args[14]
OUTPUT_PREFIX = args[15]
#################################

set.seed(42)

message(paste0('INFO: Creating PieMap: ', OUTPUT_PREFIX))

# Color palette for pie charts
my_colors <- c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222",
               "#b3446c", "#2b3d26", "#dcd300", "#875692", "#848482", "#c2b280", "#882d17", "#654522",
               "#e25822", "#008856", "#f99379", "#604e97", "#f6a600")

############################# Functions

#' Calculate autoscaled pie size range based on map extent
#' @param raster_layer Raster layer to get extent from
#' @return List with min_radius, max_radius, mid_radius (in map units)
calculate_pie_size_range <- function(raster_layer) {
  map_extent <- extent(raster_layer)
  x_range <- map_extent@xmax - map_extent@xmin
  y_range <- map_extent@ymax - map_extent@ymin

  # Use the smaller dimension to ensure pies fit well
  map_size <- min(c(x_range, y_range))

  # Define min/max as proportion of map size
  min_radius <- map_size * 0.02   # 2% of map size
  max_radius <- map_size * 0.06   # 6% of map size
  mid_radius <- (min_radius + max_radius) / 2

  message(paste0('INFO: Map extent - X: ', round(x_range, 2), ', Y: ', round(y_range, 2)))
  message(paste0('INFO: Autoscaled pie sizes - min: ', round(min_radius, 4),
                 ', mid: ', round(mid_radius, 4), ', max: ', round(max_radius, 4)))

  return(list(
    min_radius = min_radius,
    max_radius = max_radius,
    mid_radius = mid_radius
  ))
}

#' Scale trait values to pie radius using min-max normalization
#' @param trait_values Numeric vector of trait values
#' @param min_radius Minimum pie radius
#' @param max_radius Maximum pie radius
#' @return Numeric vector of scaled radii
scale_trait_to_radius <- function(trait_values, min_radius, max_radius) {
  trait_min <- min(trait_values, na.rm = TRUE)
  trait_max <- max(trait_values, na.rm = TRUE)

  # Handle case where all values are the same
  if (trait_max == trait_min) {
    return(rep((min_radius + max_radius) / 2, length(trait_values)))
  }

  # Min-max scaling to [min_radius, max_radius]
  scaled <- min_radius +
    (trait_values - trait_min) / (trait_max - trait_min) *
    (max_radius - min_radius)

  return(scaled)
}

#' Create a PieMap plot
#' @param samples Sample metadata with site, sample, latitude, longitude
#' @param clusters Cluster assignments per sample
#' @param raster_layer Raster layer for background
#' @param trait Optional data.frame with site and trait value columns for pie sizing
#' @param trait_name Name for size legend (e.g., "Tajima's D")
#' @param pie_size_range List with min_radius, max_radius, mid_radius
#' @param pie_alpha Alpha for pie charts
#' @param palette_num Colorhex palette number
#' @param palette_reverse Reverse palette direction
#' @param pop_labels Show population labels
#' @param pop_label_size Size of population labels
#' @param legend_map_title Title for raster legend
#' @return ggplot object
create_piemap <- function(samples,
                          clusters,
                          raster_layer,
                          trait = NULL,
                          trait_name = NULL,
                          pie_size_range = NULL,
                          pie_alpha = 0.6,
                          palette_num = 1022614,
                          palette_reverse = FALSE,
                          pop_labels = FALSE,
                          pop_label_size = 5,
                          legend_map_title = 'Value') {

  # Calculate pie size range if not provided
  if (is.null(pie_size_range)) {
    pie_size_range <- calculate_pie_size_range(raster_layer)
  }

  # Get cluster columns
  cluster_cols <- grep('C[0-9]+', colnames(clusters), value = TRUE)
  n_clust <- length(cluster_cols)

  # Process cluster info by population
  clusters_by_pop <- clusters %>%
    dplyr::select(-sample) %>%
    dplyr::group_by(site) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::left_join(samples[, c('site', 'latitude', 'longitude')] %>% unique, by = 'site')

  # Apply pie sizing based on trait or use midpoint for uniform size
  if (is.null(trait) || trait_name == 'NULL' || is.null(trait_name)) {
    # No trait: use midpoint size for all pies
    clusters_by_pop$pie_radius <- pie_size_range$mid_radius
    message('INFO: Using uniform midpoint pie size (no trait provided)')
    trait_name <- NULL
  } else {
    # Trait provided: min-max scale trait values to pie radius
    clusters_by_pop <- clusters_by_pop %>%
      dplyr::left_join(trait, by = 'site') %>%
      dplyr::rename(trait_value = colnames(trait)[2])

    clusters_by_pop$pie_radius <- scale_trait_to_radius(
      clusters_by_pop$trait_value,
      pie_size_range$min_radius,
      pie_size_range$max_radius
    )
    message(paste0('INFO: Scaled ', trait_name, ' values to pie radius [',
                   round(pie_size_range$min_radius, 4), ', ',
                   round(pie_size_range$max_radius, 4), ']'))
  }

  # Build plot
  gPlot <- ggplot() +
    layer_spatial(raster_layer, aes(fill = after_stat(band1))) +
    scale_fill_colorhex_c(na.value = NA, palette = palette_num, reverse = palette_reverse) +
    labs(fill = legend_map_title) +
    ggnewscale::new_scale_fill() +
    geom_scatterpie(data = clusters_by_pop,
                    aes(x = longitude, y = latitude, group = site, r = pie_radius),
                    cols = cluster_cols,
                    color = 'black',
                    alpha = pie_alpha) +
    scale_fill_manual(values = my_colors[1:n_clust]) +
    theme_bw() +
    labs(fill = 'Clusters')

  # Add population labels
  if (pop_labels) {
    label_data <- clusters_by_pop %>%
      dplyr::distinct(site, .keep_all = TRUE)
    gPlot <- gPlot +
      geom_text_repel(data = label_data,
                      aes(x = longitude, y = latitude, label = site),
                      alpha = 0.7,
                      color = 'black',
                      size = pop_label_size,
                      box.padding = 0.5,
                      point.padding = 0.3,
                      segment.color = 'grey50',
                      segment.alpha = 0.5,
                      max.overlaps = Inf)
  }

  # Add size legend showing actual trait values
  if (!is.null(trait_name)) {
    trait_vals <- clusters_by_pop$trait_value
    trait_min <- min(trait_vals, na.rm = TRUE)
    trait_max <- max(trait_vals, na.rm = TRUE)
    trait_mid <- (trait_min + trait_max) / 2

    # Create break values and corresponding point sizes (in mm for geom_point)
    size_min <- 2
    size_max <- 6
    size_mid <- (size_min + size_max) / 2

    legend_data <- data.frame(
      x = rep(mean(clusters_by_pop$longitude), 3),
      y = rep(mean(clusters_by_pop$latitude), 3),
      trait_value = c(trait_min, trait_mid, trait_max),
      point_size = c(size_min, size_mid, size_max)
    )

    # Format legend values based on trait type
    break_vals <- c(trait_min, trait_mid, trait_max)
    if (grepl("Pi|pi|diversity", trait_name, ignore.case = TRUE)) {
      break_labels <- formatC(break_vals, format = "e", digits = 2)
    } else {
      break_labels <- formatC(break_vals, format = "f", digits = 2)
    }

    gPlot <- gPlot +
      geom_point(data = legend_data,
                 aes(x = x, y = y, size = trait_value),
                 alpha = 0) +
      scale_size_continuous(
        name = trait_name,
        range = c(size_min, size_max),
        breaks = break_vals,
        labels = break_labels
      ) +
      guides(size = guide_legend(override.aes = list(alpha = 1, color = 'grey30')))

    message(paste0('INFO: Added size legend for ', trait_name,
                   ' [', round(trait_min, 3), ', ', round(trait_max, 3), ']'))
  }

  return(gPlot)
}

################################# MAIN

# Load raster
raster_data <- tryCatch({
  brick(RASTER_PATH)
}, error = function(e) {
  raster(RASTER_PATH)
})

# Extract layer
if (nlayers(raster_data) > 1) {
  # Multi-layer raster (e.g., climate stack)
  if (!RASTER_LAYER %in% names(raster_data)) {
    stop(paste0("ERROR: Layer '", RASTER_LAYER, "' not found in raster. Available: ",
                paste(names(raster_data), collapse = ", ")))
  }
  rlayer <- raster_data[[RASTER_LAYER]]
  message(paste0('INFO: Using layer "', RASTER_LAYER, '" from multi-layer raster'))
} else {
  # Single-layer raster (e.g., genetic offset)
  rlayer <- raster_data
  message('INFO: Using single-layer raster')
}

# Load data
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

clusters <- fread(CLUSTERS, colClasses = c("site" = "character"))

# Load size trait if provided
if (SIZE_TRAIT != "NULL" && file.exists(SIZE_TRAIT)) {
  size_trait <- fread(SIZE_TRAIT, colClasses = c("site" = "character"))
  message(paste0('INFO: Loaded size trait from ', SIZE_TRAIT))
} else {
  size_trait <- NULL
  message('INFO: No size trait provided, using uniform pie size')
}

# Parse population label setting
pop_labels <- toupper(substr(POP_LABEL, 1, 1)) == 'T'

# Calculate pie size range based on map extent
pie_size_range <- calculate_pie_size_range(rlayer)

# Create plot
gPlot <- create_piemap(
  samples = samples,
  clusters = clusters,
  raster_layer = rlayer,
  trait = size_trait,
  trait_name = if (SIZE_TRAIT_NAME != "NULL") SIZE_TRAIT_NAME else NULL,
  pie_size_range = pie_size_range,
  pie_alpha = PIE_ALPHA,
  palette_num = PALETTE,
  palette_reverse = PALETTE_REV,
  pop_labels = pop_labels,
  pop_label_size = POP_LABEL_SIZE,
  legend_map_title = RASTER_LEGEND
)

# Save outputs
ggsave(paste0(PLOT_DIR, OUTPUT_PREFIX, '.png'), gPlot)
ggsave(paste0(PLOT_DIR, OUTPUT_PREFIX, '.svg'), gPlot,
       device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
qsave(gPlot, paste0(INTER_DIR, OUTPUT_PREFIX, '.qs'))

message(paste0('INFO: PieMap complete: ', OUTPUT_PREFIX))
