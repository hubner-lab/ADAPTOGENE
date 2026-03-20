library(terra)
library(data.table)
library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)
#################################
SAMPLES = args[1]
MODEL_RASTERS = args[2]  # comma-separated list of per-model .tif files
N_MODELS = args[3] %>% as.numeric
PRESENT_RASTER = args[4]  # present climate raster for cell count validation
PRESENT_ALL = args[5]     # present climate all values for cell count validation
OUTPUT_RASTER = args[6]
OUTPUT_ALL = args[7]
OUTPUT_SITE = args[8]
#################################

message('INFO: Merging future climate rasters from individual models')

# Load samples for site extraction
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# Load all per-model rasters
model_files <- MODEL_RASTERS %>% str_split(',') %>% unlist
message(paste0('INFO: Loading ', length(model_files), ' model rasters'))

clim_list <- lapply(model_files, function(f) {
  message(paste0('INFO:   Loading: ', f))
  rast(f)
})

# Get layer names from first model
layer_names <- names(clim_list[[1]])
n_layers <- length(layer_names)
message(paste0('INFO: ', n_layers, ' bioclimatic layers per model'))

# Average across models
if (length(clim_list) == 1) {
  clim_future <- clim_list[[1]]
} else {
  # Stack all models and average per layer
  all_layers <- rast(clim_list)
  clim_future <- tapp(all_layers, rep(1:n_layers, N_MODELS), fun = "mean", na.rm = TRUE)
  names(clim_future) <- layer_names
}

message('INFO: Model averaging complete')

# Extract all values with cell IDs for raster remapping
# terra::extract with cell indices doesn't return ID column (unlike raster package)
clim_future_all <- terra::extract(clim_future, 1:ncell(clim_future))
clim_future_all$ID <- 1:ncell(clim_future)
clim_future_all <- na.omit(clim_future_all) %>%
  dplyr::select(ID, everything())

# Validate cell count matches present climate
message('INFO: Validating cell count matches present climate...')
present_all <- fread(PRESENT_ALL)

future_cells <- nrow(clim_future_all)
present_cells <- nrow(present_all)

message(paste0('INFO: Present climate cells: ', present_cells))
message(paste0('INFO: Future climate cells: ', future_cells))

if (future_cells != present_cells) {
  stop(paste0('ERROR: Cell count mismatch - present: ', present_cells, ', future: ', future_cells))
}
message('INFO: Cell count validation passed')

# Save merged raster
writeRaster(clim_future,
            filename = OUTPUT_RASTER,
            overwrite = TRUE,
            gdal = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
message(paste0('INFO: Saved merged raster: ', OUTPUT_RASTER))

# Save all values
clim_future_all %>%
  fwrite(OUTPUT_ALL, sep = '\t')
message(paste0('INFO: Saved all values: ', OUTPUT_ALL))

# Extract site values
coords <- vect(data.frame(long = samples$longitude, lat = samples$latitude),
               geom = c("long", "lat"), crs = "EPSG:4326")

clim_future_site <- terra::extract(clim_future, coords)[, -1] %>% as.data.frame

# Validate: check for samples with NA future climate values
na_rows <- which(rowSums(is.na(clim_future_site)) > 0)
if (length(na_rows) > 0) {
    bad_samples <- samples[na_rows, ]
    message("ERROR: Future climate extraction returned NA for the following samples:")
    for (i in seq_len(nrow(bad_samples))) {
        msg <- if (bad_samples$latitude[i] == 0 && bad_samples$longitude[i] == 0) {
            "likely placeholder (0,0)"
        } else {
            "falls on ocean/NoData pixel"
        }
        message(paste0("  - ", bad_samples$sample[i], " (", bad_samples$site[i],
                       ") at (", bad_samples$latitude[i], ", ", bad_samples$longitude[i], "): ", msg))
    }
    stop(paste0("Future climate extraction failed for ", length(na_rows), " samples. ",
                "Fix coordinates in input metadata or remove these samples."))
}

# Save site values (with sample column for traceability)
cbind(sample = samples$sample, clim_future_site) %>%
  fwrite(OUTPUT_SITE, sep = '\t')
message(paste0('INFO: Saved site values: ', OUTPUT_SITE))
message('INFO: Future climate merge complete')
