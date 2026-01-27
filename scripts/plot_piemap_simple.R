library(qs)
library(data.table)
library(dplyr)
library(see)
library(ggplot2)
library(raster)
library(ggspatial)
library(scatterpie)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
#################################
SAMPLES = args[1]
CLUSTERS = args[2]
RASTER_LAYER = args[3]    # .grd file
BIO = args[4]             # single bio variable name (e.g., "bio_1")
PALLETE_MAP = args[5] %>% as.numeric
PALLETE_MAP_reverse = args[6] %>% as.logical
PIE_ALPHA = args[7] %>% as.numeric
POP_LABELS = args[8] %>% as.logical
POP_LABEL_SIZE = args[9] %>% as.numeric
PLOT_DIR = args[10]
INTER_DIR = args[11]
# Note: PIE_SIZE and PIE_RESCALE removed - now autoscaled based on map extent
#################################

set.seed(42)

# Color palette
my.colors = c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222",
              "#b3446c", "#2b3d26", "#dcd300", "#875692", "#848482", "#c2b280", "#882d17", "#654522", "#e25822",
              "#008856", "#f99379", "#604e97", "#f6a600")

############################# Function

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
  # These proportions work well for typical population genomics maps
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

FUN_pie_map_simple <- function(samples,
                               clusters,
                               raster_layer,
                               pie_size_range = NULL,  # list with min_radius, max_radius, mid_radius
                               colors = my.colors[1:5],
                               pallete_num = 1022614,
                               pallete_reverse = FALSE,
                               mapcolor_lim = NULL,
                               pie_alpha = 0.6,
                               pop.labels = FALSE,
                               pop.label.size = 5,
                               legend_map_title = 'Climate',
                               legend_pie_title = 'Clusters',
                               legend_position = 'right') {

  # Calculate pie size range if not provided
  if (is.null(pie_size_range)) {
    pie_size_range <- calculate_pie_size_range(raster_layer)
  }

  # Process cluster info by population
  clusters.by.pop <- clusters %>%
    dplyr::select(-sample) %>%
    dplyr::group_by(site) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::left_join(samples[, c('site', 'latitude', 'longitude')] %>% unique, by = 'site')

  # Use midpoint size for uniform pies (no trait scaling)
  clusters.by.pop$pie_radius <- pie_size_range$mid_radius
  message(paste0('INFO: Using uniform midpoint pie size: ', round(pie_size_range$mid_radius, 4)))

  # Plot
  gPlot <- ggplot() +
    layer_spatial(raster_layer, aes(fill = after_stat(band1))) +
    scale_fill_colorhex_c(na.value = NA, palette = pallete_num, reverse = pallete_reverse, limits = mapcolor_lim) +
    labs(fill = legend_map_title) +
    ggnewscale::new_scale_fill() +
    geom_scatterpie(data = clusters.by.pop,
                    aes(x = longitude, y = latitude, group = site, r = pie_radius),
                    cols = grep('C[0-9]+', colnames(clusters), value = TRUE),
                    color = 'black',
                    alpha = pie_alpha) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    labs(fill = legend_pie_title) +
    theme(legend.position = legend_position)

  if (pop.labels) {
    gPlot <- gPlot +
      geom_text(data = clusters.by.pop,
                aes(x = longitude, y = latitude, label = site),
                alpha = pie_alpha,
                color = 'white',
                size = pop.label.size)
  }

  return(gPlot)
}

################################# MAIN

# Load raster stack
raster_stack <- brick(RASTER_LAYER)

# Check if bio variable exists in raster
if (!BIO %in% names(raster_stack)) {
  stop(paste0("ERROR: Variable '", BIO, "' not found in raster. Available: ",
              paste(names(raster_stack), collapse = ", ")))
}

# Extract single layer
rlayer <- raster_stack[[BIO]]

# Load data
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

clusters <- fread(CLUSTERS,
                  colClasses = c("site" = "character"))
Nclust <- ncol(clusters) - 2  # minus site and sample columns

# Calculate autoscaled pie size range based on map extent
message('INFO: Calculating autoscaled pie sizes based on map extent')
pie_size_range <- calculate_pie_size_range(rlayer)

# Generate simple piemap with uniform midpoint pie size
message('INFO: Creating simple PieMap with uniform midpoint pie size')
gPieMap <- FUN_pie_map_simple(samples,
                              clusters,
                              rlayer,
                              pie_size_range = pie_size_range,
                              pie_alpha = PIE_ALPHA,
                              colors = my.colors[1:Nclust],
                              pallete_num = PALLETE_MAP,
                              pallete_reverse = PALLETE_MAP_reverse,
                              pop.labels = POP_LABELS,
                              pop.label.size = POP_LABEL_SIZE,
                              legend_map_title = BIO)

ggsave(paste0(PLOT_DIR, "PieMap_", BIO, ".png"), gPieMap)
ggsave(paste0(PLOT_DIR, "PieMap_", BIO, ".svg"), gPieMap,
       device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
qsave(gPieMap, paste0(INTER_DIR, "PieMap_", BIO, ".qs"))

message('INFO: Simple PieMap complete')
