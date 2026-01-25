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
PIE_SIZE = args[7] %>% as.numeric
PIE_ALPHA = args[8] %>% as.numeric
POP_LABELS = args[9] %>% as.logical
PIE_RESCALE = args[10] %>% as.numeric
POP_LABEL_SIZE = args[11] %>% as.numeric
PLOT_DIR = args[12]
INTER_DIR = args[13]
#################################

set.seed(42)

# Color palette
my.colors = c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222",
              "#b3446c", "#2b3d26", "#dcd300", "#875692", "#848482", "#c2b280", "#882d17", "#654522", "#e25822",
              "#008856", "#f99379", "#604e97", "#f6a600")

############################# Function
FUN_pie_map_simple <- function(samples,
                               clusters,
                               raster_layer,
                               colors = my.colors[1:5],
                               pallete_num = 1022614,
                               pallete_reverse = FALSE,
                               mapcolor_lim = NULL,
                               pie.size = 10,
                               pie_alpha = 0.6,
                               pop.labels = FALSE,
                               pop.label.size = 5,
                               legend_map_title = 'Climate',
                               legend_pie_title = 'Clusters',
                               legend_position = 'right') {

  # Process cluster info by population
  clusters.by.pop <- clusters %>%
    dplyr::select(-sample) %>%
    dplyr::group_by(site) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::left_join(samples[, c('site', 'latitude', 'longitude')] %>% unique, by = 'site')

  # Calculate uniform pie size based on map extent
  map_extent <- extent(raster_layer)
  x_range <- map_extent@xmax - map_extent@xmin
  y_range <- map_extent@ymax - map_extent@ymin
  map_size <- mean(c(x_range, y_range))
  scaling_factor <- map_size * 0.05  # uniform size for all pies
  clusters.by.pop$pie_radius <- scaling_factor

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

# Generate simple piemap
gPieMap <- FUN_pie_map_simple(samples,
                              clusters,
                              rlayer,
                              pie.size = PIE_SIZE,
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
