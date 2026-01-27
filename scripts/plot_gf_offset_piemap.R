library(dplyr)
library(data.table)
library(ggplot2)
library(raster)
library(ggspatial)
library(scatterpie)
library(see)
library(ggnewscale)
library(ggrepel)
library(svglite)
library(qs)
library(stringr)
args = commandArgs(trailingOnly=TRUE)
###############
OFFSET_RASTER = args[1]
SAMPLES = args[2]
CLUSTERS = args[3]
IBD = args[4]
TAJIMA = args[5]
PALETTE = args[6] %>% as.numeric
PALETTE_REV = args[7]
PIE_ALPHA = args[8] %>% as.numeric
IBD_ALPHA = args[9] %>% as.numeric
IBD_COLOR = args[10]
POP_LABEL = args[11]
POP_LABEL_SIZE = args[12] %>% as.numeric
PLOT_DIR = args[13]
INTER_DIR = args[14]
SUFFIX = args[15]
###############

message('INFO: Plotting Genetic Offset PieMap')

PALETTE_REV <- as.logical(PALETTE_REV)

# Colors from Polychrome palette
my_colors <- c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222",
               "#b3446c", "#2b3d26", "#dcd300", "#875692", "#848482", "#c2b280", "#882d17", "#654522",
               "#e25822", "#008856", "#f99379", "#604e97", "#f6a600")

# Load data
rast_offset <- raster(OFFSET_RASTER)

samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

clusters <- fread(CLUSTERS, colClasses = c("site" = "character"))
cluster_cols <- grep('C[0-9]+', colnames(clusters), value = TRUE)
n_clust <- length(cluster_cols)

ibd_dt <- fread(IBD, colClasses = c("pop1" = "character", "pop2" = "character"))

tajima <- fread(TAJIMA, colClasses = c("site" = "character"))

# Process clusters by population
trait <- tajima
trait[[2]] <- scales::rescale(trait[[2]], c(0.3, 1)) / 10

clusters_by_pop <- clusters %>%
  dplyr::select(-sample) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::left_join(samples[, c('site', 'latitude', 'longitude')] %>% unique, by = 'site') %>%
  dplyr::left_join(trait, by = 'site') %>%
  dplyr::rename(trait = colnames(trait)[2])

# Build plot
gPlot <- ggplot() +
  layer_spatial(rast_offset, aes(fill = after_stat(band1))) +
  scale_fill_colorhex_c(na.value = NA, palette = PALETTE, reverse = TRUE) +
  labs(fill = 'Genetic Offset') +
  ggnewscale::new_scale_fill() +
  geom_scatterpie(data = clusters_by_pop,
                  aes(x = longitude, y = latitude, group = site, r = trait),
                  cols = cluster_cols,
                  color = 'black',
                  alpha = PIE_ALPHA) +
  scale_fill_manual(values = my_colors[1:n_clust]) +
  theme_bw() +
  labs(fill = 'Clusters')

# Add IBD lines
if (nrow(ibd_dt) > 0) {
  gPlot <- gPlot +
    lapply(1:nrow(ibd_dt), function(i) {
      geom_line(data = clusters_by_pop %>%
                  dplyr::filter(site %in% c(ibd_dt[i, ]$pop1, ibd_dt[i, ]$pop2)),
                aes(x = longitude, y = latitude),
                color = IBD_COLOR,
                alpha = IBD_ALPHA)
    })
}

# Add population labels
if (toupper(substr(POP_LABEL, 1, 1)) == 'T') {
  label_data <- clusters_by_pop %>%
    dplyr::distinct(site, .keep_all = TRUE)
  gPlot <- gPlot +
    geom_text_repel(data = label_data,
                    aes(x = longitude, y = latitude, label = site),
                    alpha = 0.7, color = 'black', size = POP_LABEL_SIZE,
                    box.padding = 0.5, point.padding = 0.3,
                    segment.color = 'grey50', segment.alpha = 0.5,
                    max.overlaps = Inf)
}

# Save
ggsave(paste0(PLOT_DIR, 'GeneticOffsetPieMap_', SUFFIX, '.png'), gPlot)
ggsave(paste0(PLOT_DIR, 'GeneticOffsetPieMap_', SUFFIX, '.svg'), gPlot,
       device = svglite::svglite, bg = "transparent", fix_text_size = FALSE)
qsave(gPlot, paste0(INTER_DIR, 'GeneticOffsetPieMap_', SUFFIX, '.qs'))

message('INFO: Genetic Offset PieMap complete')
