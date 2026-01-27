library(dplyr)
library(data.table)
library(vegan)
library(stringr)
library(geosphere)
library(ggplot2)
library(ggpubr)
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

analyze_mantel <- function(geo, env, clust) {
  # Calculate distance matrices
  env_dist <- vegdist(env,
                      method = "euclidean",
                      binary = FALSE,
                      diag = FALSE,
                      upper = FALSE,
                      na.rm = FALSE)

  clust_dist <- vegdist(clust,
                        method = 'hellinger')

  geo_dist <- distm(geo, fun = distVincentyEllipsoid)

  # Simple Mantel tests
  ibd_full <- mantel(clust_dist,
                     geo_dist,
                     method = 'pearson',
                     permutations = 999)

  ibe_full <- mantel(clust_dist,
                     env_dist,
                     method = 'pearson',
                     permutations = 999)

  # Partial Mantel tests
  ibd_partial <- mantel.partial(clust_dist,
                                geo_dist,
                                env_dist,
                                method = 'pearson',
                                permutations = 999)

  ibe_partial <- mantel.partial(clust_dist,
                                env_dist,
                                geo_dist,
                                method = 'pearson',
                                permutations = 999)

  # Calculate variance components (R-squared)
  r2_geo_total <- max(0, ibd_full$statistic^2)
  r2_env_total <- max(0, ibe_full$statistic^2)
  r2_geo_pure <- max(0, ibd_partial$statistic^2)
  r2_env_pure <- max(0, ibe_partial$statistic^2)

  # Shared variance (intersection)
  r2_shared <- max(0, r2_geo_total + r2_env_total - r2_geo_pure - r2_env_pure)

  # Unexplained variance
  r2_total_explained <- r2_geo_pure + r2_env_pure + r2_shared
  r2_unexplained <- max(0, 1 - r2_total_explained)

  # Normalize if components don't sum to 1
  total_check <- r2_geo_pure + r2_env_pure + r2_shared + r2_unexplained
  if (abs(total_check - 1) > 0.01) {
    message('WARNING: Variance components do not sum to 1. Normalizing...')
    r2_geo_pure <- r2_geo_pure / total_check
    r2_env_pure <- r2_env_pure / total_check
    r2_shared <- r2_shared / total_check
    r2_unexplained <- r2_unexplained / total_check
  }

  message(sprintf('INFO: Variance components - Geography: %.3f, Environment: %.3f, Shared: %.3f, Unexplained: %.3f',
                  r2_geo_pure, r2_env_pure, r2_shared, r2_unexplained))

  results <- list(
    mantel_tests = list(
      ibd_full = ibd_full,
      ibe_full = ibe_full,
      ibd_partial = ibd_partial,
      ibe_partial = ibe_partial
    ),
    variance_components = list(
      geo_pure = r2_geo_pure,
      env_pure = r2_env_pure,
      shared = r2_shared,
      unexplained = r2_unexplained
    )
  )

  return(results)
}

create_variance_pie <- function(results) {
  geo_pure <- results$variance_components$geo_pure
  env_pure <- results$variance_components$env_pure
  shared <- results$variance_components$shared
  unexplained <- results$variance_components$unexplained

  # Define colors as named vector for consistent mapping
  component_colors <- c(
    "Geography Only" = "#56B4E9",
    "Environment Only" = "#E69F00",
    "Geography \u00d7 Environment" = "#9B59B6",
    "Unexplained" = "#999999"
  )

  pie_data <- data.frame(
    component = names(component_colors),
    value = c(geo_pure, env_pure, shared, unexplained),
    stringsAsFactors = FALSE
  )

  # Filter negligible components
  pie_data <- pie_data[pie_data$value > 0.001, ]
  pie_data$percentage <- pie_data$value / sum(pie_data$value) * 100
  pie_data$percent_label <- sprintf("%.1f%%", pie_data$percentage)

  colors <- component_colors[pie_data$component]

  p <- ggpubr::ggpie(
    data = pie_data,
    x = "value",
    label = "percent_label",
    fill = "component",
    color = "white",
    palette = unname(colors),
    lab.pos = "in",
    lab.font = c(4, "bold", "black")
  ) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      panel.background = element_blank(),
      plot.background = element_blank()
    ) +
    guides(fill = guide_legend(title = NULL))

  # Log breakdown
  for (i in 1:nrow(pie_data)) {
    message(sprintf('INFO: %s: %.3f (%.1f%%)',
                    pie_data$component[i],
                    pie_data$value[i],
                    pie_data$percentage[i]))
  }

  return(p)
}

########################################## Main

# Load + process
geo_raw <- fread(SAMPLES)
geo_cols <- colnames(geo_raw)

lat_col <- geo_cols[grepl("^lat", geo_cols, ignore.case = TRUE)]
lon_col <- geo_cols[grepl("^lon", geo_cols, ignore.case = TRUE)]

if (length(lat_col) == 0) stop("No latitude column found (expected column starting with 'lat')")
if (length(lon_col) == 0) stop("No longitude column found (expected column starting with 'lon')")

if (length(lat_col) > 1) { message('WARNING: Multiple latitude columns, using: ', lat_col[1]); lat_col <- lat_col[1] }
if (length(lon_col) > 1) { message('WARNING: Multiple longitude columns, using: ', lon_col[1]); lon_col <- lon_col[1] }

geo <- geo_raw %>%
  dplyr::select(longitude = all_of(lon_col), latitude = all_of(lat_col))

if (any(is.na(geo))) {
  message('WARNING: Missing values in geographic coordinates, removing affected rows')
  geo <- na.omit(geo)
}

message(sprintf('INFO: Geographic data - Longitude: [%.2f, %.2f], Latitude: [%.2f, %.2f]',
                min(geo$longitude), max(geo$longitude), min(geo$latitude), max(geo$latitude)))

env <- fread(ENV) %>%
  dplyr::select(all_of(PREDICTORS_SELECTED)) %>%
  scale()

clust <- fread(CLUSTERS) %>%
  dplyr::select(-sample, -site)

# Run Mantel analysis
results <- analyze_mantel(geo, env, clust)

# Create variance pie chart
gMantel <- create_variance_pie(results)

# Save
ggsave(paste0(PLOT_DIR, 'MantelTest.png'), gMantel, width = 10, height = 8, dpi = 300)
ggsave(paste0(PLOT_DIR, 'MantelTest.svg'), gMantel,
       device = svglite::svglite, bg = 'transparent', fix_text_size = FALSE)
qsave(results, paste0(INTER_DIR, 'MantelTest.qs'))

message('INFO: Mantel test complete')
