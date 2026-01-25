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
PIE_SIZE = args[11] %>% as.numeric
PIE_ALPHA = args[12] %>% as.numeric
IBD_ALPHA = args[13] %>% as.numeric
IBD_COLOR = args[14]
POP_LABELS = args[15] %>% as.logical
PIE_RESCALE = args[16] %>% as.numeric
POP_LABEL_SIZE = args[17] %>% as.numeric
PLOT_DIR = args[18]
INTER_DIR = args[19]
#################################

set.seed(42)

# Color palette
my.colors = c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222",
              "#b3446c", "#2b3d26", "#dcd300", "#875692", "#848482", "#c2b280", "#882d17", "#654522", "#e25822",
              "#008856", "#f99379", "#604e97", "#f6a600")

############################# Function
FUN_pie_map <- function(samples,
                        clusters,
                        raster_layer,
                        trait = NULL,
                        ibd.dt = NULL,
                        ibd.line_alpha = 0.6,
                        ibd.line_color = 'black',
                        colors = my.colors[1:5],
                        pallete_num = 1022614,
                        pallete_reverse = FALSE,
                        mapcolor_lim = NULL,
                        pie.size = 10,
                        pie_alpha = 0.6,
                        pie_rescale = 1,
                        pop.labels = FALSE,
                        pop.label.size = 5,
                        legend_map_title = 'Climate',
                        legend_pie_title = 'Clusters',
                        legend_position = 'right') {

  if (is.null(trait)) {
    trait = data.frame(site = clusters$site,
                       trait = pie.size)
  }

  # Process cluster info by population
  clusters.by.pop <- clusters %>%
    dplyr::select(-sample) %>%
    dplyr::group_by(site) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::left_join(samples[, c('site', 'latitude', 'longitude')] %>% unique, by = 'site') %>%
    dplyr::left_join(trait, by = 'site') %>%
    dplyr::rename(trait = colnames(trait)[2])

  # Pie size scaling
  map_extent <- extent(raster_layer)
  x_range <- map_extent@xmax - map_extent@xmin
  y_range <- map_extent@ymax - map_extent@ymin
  map_size <- mean(c(x_range, y_range))

  min_size <- 0.2
  max_size <- 1

  clusters.by.pop$trait_rescaled <- min_size +
    (clusters.by.pop$trait - min(clusters.by.pop$trait, na.rm = TRUE)) /
    (max(clusters.by.pop$trait, na.rm = TRUE) - min(clusters.by.pop$trait, na.rm = TRUE)) *
    (max_size - min_size)

  scaling_factor <- map_size * 0.1
  clusters.by.pop$trait_scaled <- clusters.by.pop$trait_rescaled * scaling_factor

  # Plot

  gPlot <- ggplot() +
    layer_spatial(raster_layer, aes(fill = after_stat(band1))) +
    scale_fill_colorhex_c(na.value = NA, palette = pallete_num, reverse = pallete_reverse, limits = mapcolor_lim) +
    labs(fill = legend_map_title) +
    ggnewscale::new_scale_fill() +
    geom_scatterpie(data = clusters.by.pop,
                    aes(x = longitude, y = latitude, group = site, r = trait_scaled),
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

# TajimaD plot
gTajima <- FUN_pie_map(samples,
                       clusters,
                       rlayer,
                       ibd.dt = ibd.dt,
                       ibd.line_alpha = IBD_ALPHA,
                       ibd.line_color = IBD_COLOR,
                       trait = tajima,
                       pie.size = PIE_SIZE,
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

# PiDiversity plot
gDiversity <- FUN_pie_map(samples,
                          clusters,
                          rlayer,
                          ibd.dt = ibd.dt,
                          ibd.line_alpha = IBD_ALPHA,
                          ibd.line_color = IBD_COLOR,
                          trait = pi_diversity,
                          pie.size = PIE_SIZE,
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
  
  gCustom <- FUN_pie_map(samples,
                         clusters,
                         rlayer,
                         ibd.dt = ibd.dt,
                         ibd.line_alpha = IBD_ALPHA,
                         ibd.line_color = IBD_COLOR,
                         trait = custom_trait,
                         pie.size = PIE_SIZE,
                         pie_alpha = PIE_ALPHA,
                         colors = my.colors[1:Nclust],
                         pallete_num = PALLETE_MAP,
                         pallete_reverse = PALLETE_MAP_reverse,
                         pie_rescale = PIE_RESCALE,
                         pop.labels = POP_LABELS,
                         pop.label.size = POP_LABEL_SIZE,
                         legend_map_title = BIO)
  
  ggsave(paste0(PLOT_DIR, "PieMap_", BIO, "_", trait_name, ".png"), gCustom)
  ggsave(paste0(PLOT_DIR, "PieMap_", BIO, "_", trait_name, ".svg"), gCustom,
         device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
  qsave(gCustom, paste0(INTER_DIR, "PieMap_", BIO, "_", trait_name, ".qs"))
}
