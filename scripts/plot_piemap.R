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
IBD = args[4]
TAJIMA = args[5]
PI_DIVERSITY = args[6]
CUSTOM_TRAIT = args[7]    # path or 'NULL'
BIO = args[8]             # single bio variable name (e.g., "bio_1")
PALLETE_MAP = args[9] %>% as.numeric
PALLETE_MAP_reverse = args[10] %>% as.logical
PIE_ALPHA = args[11] %>% as.numeric
IBD_ALPHA = args[12] %>% as.numeric
IBD_COLOR = args[13]
POP_LABELS = args[14] %>% as.logical
POP_LABEL_SIZE = args[15] %>% as.numeric
PLOT_DIR = args[16]
INTER_DIR = args[17]
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

FUN_pie_map <- function(samples,
                        clusters,
                        raster_layer,
                        trait = NULL,
                        pie_size_range = NULL,  # list with min_radius, max_radius, mid_radius
                        ibd.dt = NULL,
                        ibd.line_alpha = 0.6,
                        ibd.line_color = 'black',
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

  # Apply pie sizing based on trait or use midpoint for uniform size
  if (is.null(trait)) {
    # No trait: use midpoint size for all pies
    clusters.by.pop$pie_radius <- pie_size_range$mid_radius
    message('INFO: Using uniform midpoint pie size (no trait provided)')
  } else {
    # Trait provided: min-max scale trait values to pie radius
    clusters.by.pop <- clusters.by.pop %>%
      dplyr::left_join(trait, by = 'site') %>%
      dplyr::rename(trait_value = colnames(trait)[2])

    clusters.by.pop$pie_radius <- scale_trait_to_radius(
      clusters.by.pop$trait_value,
      pie_size_range$min_radius,
      pie_size_range$max_radius
    )
    message(paste0('INFO: Scaled trait values to pie radius [',
                   round(pie_size_range$min_radius, 4), ', ',
                   round(pie_size_range$max_radius, 4), ']'))
  }

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

  if (!is.null(ibd.dt) && nrow(ibd.dt) > 0) {
    gPlot <- gPlot +
      lapply(1:nrow(ibd.dt), function(i) {
        geom_line(data = clusters.by.pop %>% dplyr::filter(site %in% c(ibd.dt[i, ]$pop1, ibd.dt[i, ]$pop2)),
                  aes(x = longitude, y = latitude),
                  color = ibd.line_color,
                  alpha = ibd.line_alpha)
      })
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

ibd.dt <- fread(IBD,
                colClasses = c("pop1" = "character",
                               "pop2" = "character"))

tajima <- fread(TAJIMA,
                colClasses = c("site" = "character"))

pi_diversity <- fread(PI_DIVERSITY,
                      colClasses = c("site" = "character"))

# Calculate autoscaled pie size range based on map extent (once for all plots)
message('INFO: Calculating autoscaled pie sizes based on map extent')
pie_size_range <- calculate_pie_size_range(rlayer)

# TajimaD plot - pie size scaled by Tajima's D values
message('INFO: Creating TajimaD PieMap with trait-scaled pie sizes')
gTajima <- FUN_pie_map(samples,
                       clusters,
                       rlayer,
                       trait = tajima,
                       pie_size_range = pie_size_range,
                       ibd.dt = ibd.dt,
                       ibd.line_alpha = IBD_ALPHA,
                       ibd.line_color = IBD_COLOR,
                       pie_alpha = PIE_ALPHA,
                       colors = my.colors[1:Nclust],
                       pallete_num = PALLETE_MAP,
                       pallete_reverse = PALLETE_MAP_reverse,
                       pop.labels = POP_LABELS,
                       pop.label.size = POP_LABEL_SIZE,
                       legend_map_title = BIO)

ggsave(paste0(PLOT_DIR, "PieMap_", BIO, "_TajimaD.png"), gTajima)
ggsave(paste0(PLOT_DIR, "PieMap_", BIO, "_TajimaD.svg"), gTajima,
       device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
qsave(gTajima, paste0(INTER_DIR, "PieMap_", BIO, "_TajimaD.qs"))

# PiDiversity plot - pie size scaled by Pi diversity values
message('INFO: Creating PiDiversity PieMap with trait-scaled pie sizes')
gDiversity <- FUN_pie_map(samples,
                          clusters,
                          rlayer,
                          trait = pi_diversity,
                          pie_size_range = pie_size_range,
                          ibd.dt = ibd.dt,
                          ibd.line_alpha = IBD_ALPHA,
                          ibd.line_color = IBD_COLOR,
                          pie_alpha = PIE_ALPHA,
                          colors = my.colors[1:Nclust],
                          pallete_num = PALLETE_MAP,
                          pallete_reverse = PALLETE_MAP_reverse,
                          pop.labels = POP_LABELS,
                          pop.label.size = POP_LABEL_SIZE,
                          legend_map_title = BIO)

ggsave(paste0(PLOT_DIR, "PieMap_", BIO, "_PiDiversity.png"), gDiversity)
ggsave(paste0(PLOT_DIR, "PieMap_", BIO, "_PiDiversity.svg"), gDiversity,
       device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
qsave(gDiversity, paste0(INTER_DIR, "PieMap_", BIO, "_PiDiversity.qs"))

# Custom trait plot (if file exists)
if (file.exists(CUSTOM_TRAIT) && CUSTOM_TRAIT != "NULL") {
  custom_trait <- fread(CUSTOM_TRAIT,
                        colClasses = c("site" = "character"))

  trait_name <- custom_trait %>% dplyr::select(-site) %>% colnames

  message(paste0('INFO: Creating Custom trait (', trait_name, ') PieMap with trait-scaled pie sizes'))
  gCustom <- FUN_pie_map(samples,
                         clusters,
                         rlayer,
                         trait = custom_trait,
                         pie_size_range = pie_size_range,
                         ibd.dt = ibd.dt,
                         ibd.line_alpha = IBD_ALPHA,
                         ibd.line_color = IBD_COLOR,
                         pie_alpha = PIE_ALPHA,
                         colors = my.colors[1:Nclust],
                         pallete_num = PALLETE_MAP,
                         pallete_reverse = PALLETE_MAP_reverse,
                         pop.labels = POP_LABELS,
                         pop.label.size = POP_LABEL_SIZE,
                         legend_map_title = BIO)

  ggsave(paste0(PLOT_DIR, "PieMap_", BIO, "_", trait_name, ".png"), gCustom)
  ggsave(paste0(PLOT_DIR, "PieMap_", BIO, "_", trait_name, ".svg"), gCustom,
         device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
  qsave(gCustom, paste0(INTER_DIR, "PieMap_", BIO, "_", trait_name, ".qs"))
}

message('INFO: PieMap plots complete')
