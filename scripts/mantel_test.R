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
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
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
geo_raw <- fread(SAMPLES, colClasses = c("site" = "character", "sample" = "character"))
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

env_raw <- fread(ENV) %>%
  dplyr::select(all_of(PREDICTORS_SELECTED))

# Align the ancestry matrix to the metadata BY SAMPLE ID, not by position.
# clusters_K{k}.tsv covers every sample that reached sNMF, while SAMPLES here is
# metadata_climate_valid.tsv — the climate-valid subset — so the two tables differ
# in length whenever any sample lacks coordinates or climate values (46 vs 47 on
# the shipped test dataset). The previous positional read subset `clust` with a
# logical vector one element shorter than its own row count, which R recycles, so
# every ancestry row after the dropped sample was paired with the wrong site's
# geography and climate — silently, since nothing compared the lengths.
clust_raw <- fread(CLUSTERS, colClasses = c("sample" = "character", "site" = "character"))
clust_idx <- match(geo_raw$sample, clust_raw$sample)
if (anyNA(clust_idx)) {
  stop('ERROR: ', sum(is.na(clust_idx)), ' sample(s) in ', basename(SAMPLES),
       ' have no row in ', basename(CLUSTERS),
       ' — cannot align the ancestry matrix to the climate/geography tables.')
}
if (nrow(clust_raw) != length(clust_idx)) {
  message(sprintf('INFO: ancestry table has %d samples, %d are climate-valid — subsetting by sample ID',
                  nrow(clust_raw), length(clust_idx)))
}
clust <- clust_raw[clust_idx, ] %>%
  dplyr::select(-sample, -site)

# Remove rows with NA climate values (affects geo, env, clust equally)
complete_rows <- complete.cases(env_raw)
if (any(!complete_rows)) {
  n_removed <- sum(!complete_rows)
  message(sprintf('WARNING: Removing %d samples with missing climate values (%d remain)',
                  n_removed, sum(complete_rows)))
  env_raw <- env_raw[complete_rows, ]
  geo <- geo[complete_rows, ]
  clust <- clust[complete_rows, ]
}

# SITE-LEVEL AGGREGATION. Every one of the three matrices below is per SAMPLE,
# but IBD/IBE is a question about SITES: geographic coordinates and climate values
# are identical for every sample sequenced at a site, so a 30-sample site
# contributes 435 within-site pairs at distance ~0 to both the geographic and the
# environmental matrix, and Mantel's permutation test counts those as independent
# units — the p-value is inflated by sampling effort, not by geography. The same
# imbalance also weights scale()'s centre and sd by sample size, and unlike the
# other consumers of the climate table this one turns the scaled columns into a
# EUCLIDEAN DISTANCE, where column scale is meaningful and the weighting is
# therefore NOT absorbed (cf. the note in download_climate_present.R).
# So: collapse to one row per site first, then scale. Ancestry is averaged per
# site (Q rows sum to 1, so the mean is still a valid ancestry composition —
# the standard population-level unit for a Hellinger distance); coordinates are
# averaged too, in case a site's samples carry slightly jittered coordinates.
site_vec <- geo_raw$site[complete_rows]
stopifnot(length(site_vec) == nrow(env_raw), nrow(geo) == nrow(env_raw),
          nrow(clust) == nrow(env_raw))

site_levels <- unique(site_vec)
n_sites     <- length(site_levels)
agg_mean <- function(df) {
    as.data.frame(do.call(rbind, lapply(site_levels, function(s)
        colMeans(as.matrix(df[site_vec == s, , drop = FALSE]), na.rm = TRUE))))
}
message(sprintf('INFO: collapsing %d samples to %d sites for the Mantel tests',
                length(site_vec), n_sites))
env_raw <- agg_mean(env_raw)
geo     <- agg_mean(geo)
clust   <- agg_mean(clust)

if (n_sites < 4) {
    stop('ERROR: Mantel test needs at least 4 sites (', n_sites, ' present) — ',
         'a distance matrix over 3 sites has 3 pairs and no permutation test is ',
         'meaningful. Disable Population.calc_stats or add sites.')
}
if (n_sites < 8) {
    message(sprintf(paste0('WARNING: only %d sites — %d pairwise distances, and the ',
                           'permutation test cannot resolve p below 1/%d. Read the ',
                           'variance components as descriptive, not as a test.'),
                    n_sites, n_sites * (n_sites - 1) / 2, factorial(n_sites)))
}

env <- scale(env_raw)

# Drop zero-variance predictors (scale() produces NaN when sd=0)
zero_var <- apply(env, 2, function(x) all(is.nan(x)))
if (any(zero_var)) {
  dropped <- colnames(env)[zero_var]
  message(sprintf('WARNING: Dropping %d zero-variance predictor(s): %s',
                  length(dropped), paste(dropped, collapse = ', ')))
  env <- env[, !zero_var, drop = FALSE]
}

if (ncol(env) == 0) {
  stop('ERROR: All predictors have zero variance — cannot compute environmental distances. ',
       'Check that PREDICTORS_SELECTED contains variables with variation across sites.')
}

# Run Mantel analysis
results <- analyze_mantel(geo, env, clust)

# Create variance pie chart
gMantel <- create_variance_pie(results)

# Save
ggsave(paste0(PLOT_DIR, 'mantel_test.png'), gMantel, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(paste0(PLOT_DIR, 'mantel_test.svg'), gMantel,
       device = svglite::svglite, bg = 'white', fix_text_size = FALSE)
qsave(results, paste0(INTER_DIR, 'mantel_test.qs'))

message('INFO: Mantel test complete')
