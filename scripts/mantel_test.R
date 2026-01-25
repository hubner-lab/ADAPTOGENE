library(dplyr)
library(data.table)
library(parallel)
library(vegan)
library(stringr)
library(geosphere)
library(ggplot2)
library(qs)

args = commandArgs(trailingOnly=TRUE)
####################################
SAMPLES = args[1]
CLUSTERS = args[2]
ENV = args[3]
PREDICTORS_SELECTED = args[4] %>% str_split(',') %>% unlist
PLOT_DIR = args[5]
INTER_DIR = args[6]
####################################

#################################### Functions
plot_mantel_comparison <- function(ibe, ibd, ibe_partial, ibd_partial) {
  # Prepare data frame for plotting
  data <- data.frame(
    test_type = c("Simple Mantel", "Partial Mantel", "Simple Mantel", "Partial Mantel"),
    comparison = c("GEN ~ GEO", "GEN ~ GEO | ENV", "GEN ~ ENV", "GEN ~ ENV | GEO"),
    correlation = c(ibd$statistic, ibd_partial$statistic,
                    ibe$statistic, ibe_partial$statistic),
    p_value = c(ibd$signif, ibd_partial$signif,
                ibe$signif, ibe_partial$signif),
    x_pos = 1:4,
    color_group = c("GEO", "GEO", "ENV", "ENV")
  )

  # Create connecting segments data
  segments_data <- data.frame(
    x = c(1, 3),
    xend = c(2, 4),
    y = c(ibd$statistic, ibe$statistic),
    yend = c(ibd_partial$statistic, ibe_partial$statistic)
  )

  # Create the plot
  p <- ggplot() +
    geom_segment(data = segments_data,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "gray50",
                 linewidth = 0.5,
                 linetype = "dashed") +
    geom_point(data = data,
               aes(x = x_pos, y = correlation, color = color_group, shape = test_type),
               size = 4) +
    geom_text(data = data,
              aes(x = x_pos, y = correlation,
                  label = sprintf("r = %.3f\np = %.3f", correlation, p_value)),
              vjust = -1.2,
              size = 3) +
    scale_color_manual(values = c("GEO" = "#56B4E9", "ENV" = "#E69F00")) +
    scale_shape_manual(values = c("Simple Mantel" = 16, "Partial Mantel" = 17)) +
    scale_x_continuous(breaks = 1:4,
                       labels = c("GEN ~ GEO", "GEN ~ GEO | ENV",
                                  "GEN ~ ENV", "GEN ~ ENV | GEO"),
                       limits = c(0.5, 4.5)) +
    scale_y_continuous(limits = c(
      min(c(data$correlation, 0)) - 0.1,
      max(data$correlation) + 0.2
    )) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      y = "Mantel Correlation Coefficient (r)",
      title = "Simple and Partial Mantel Test Correlations"
    )

  return(p)
}

########################################## Main

# Load + process
geo = fread(SAMPLES) %>%
  dplyr::select(latitude, longitude)

env = fread(ENV) %>%
  dplyr::select(all_of(PREDICTORS_SELECTED)) %>%
  scale

clust = fread(CLUSTERS) %>%
  dplyr::select(-sample, -site)

# Calculate distances
env_dist = vegdist(env,
                   method = "euclidean",
                   binary = FALSE,
                   diag = FALSE,
                   upper = FALSE,
                   na.rm = FALSE)

clust_dist = vegdist(clust,
                     method = 'hellinger')

geo_dist = distm(geo,
                 fun = distVincentyEllipsoid)

# Calculate IBD, IBE and partial versions
ibe = mantel.partial(clust_dist,
                     env_dist,
                     geo_dist,
                     method = 'pearson',
                     permutations = 999)

ibd = mantel.partial(clust_dist,
                     geo_dist,
                     env_dist,
                     method = 'pearson',
                     permutations = 999)

ibd_full = mantel(clust_dist,
                  geo_dist,
                  method = 'pearson',
                  permutations = 999)

ibe_full = mantel(clust_dist,
                  env_dist,
                  method = 'pearson',
                  permutations = 999)

# Plot
gMantel = plot_mantel_comparison(ibe_full, ibd_full, ibe, ibd)

# Save
ggsave(paste0(PLOT_DIR, 'MantelTest.png'), gMantel)
ggsave(paste0(PLOT_DIR, 'MantelTest.svg'), gMantel,
       device = svglite::svglite, bg = 'transparent')
qsave(gMantel, paste0(INTER_DIR, 'MantelTest.qs'))
