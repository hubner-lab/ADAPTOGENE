library(raster)
library(data.table)
library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)
#################################
SAMPLES = args[1]
MODEL_RASTERS = args[2]  # comma-separated list of per-model .grd files
N_MODELS = args[3] %>% as.numeric
OUTPUT_RASTER = args[4]
OUTPUT_ALL = args[5]
OUTPUT_SITE = args[6]
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
  stack(f)
})

# Get layer names from first model
layer_names <- names(clim_list[[1]])
n_layers <- length(layer_names)
message(paste0('INFO: ', n_layers, ' bioclimatic layers per model'))

# Average across models
if (length(clim_list) == 1) {
  clim_future <- clim_list[[1]]
} else {
  clim_future <- clim_list %>%
    stack %>%
    stackApply(rep(1:n_layers, N_MODELS), mean, na.rm = TRUE) %>%
    setNames(layer_names)
}

message('INFO: Model averaging complete')

# Save merged raster
writeRaster(clim_future,
            filename = OUTPUT_RASTER,
            format = "raster",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
message(paste0('INFO: Saved merged raster: ', OUTPUT_RASTER))

# Extract all values (na.omit to match present climate format)
clim_future_all <- raster::extract(clim_future, 1:ncell(clim_future), df = TRUE) %>%
  na.omit
clim_future_all %>%
  fwrite(OUTPUT_ALL, sep = '\t')
message(paste0('INFO: Saved all values: ', OUTPUT_ALL))

# Extract site values
coords <- data.frame(long = samples$longitude,
                     lat = samples$latitude)
coordinates(coords) <- ~long + lat

clim_future_site <- raster::extract(clim_future, coords) %>% as.data.frame
clim_future_site %>%
  fwrite(OUTPUT_SITE, sep = '\t')
message(paste0('INFO: Saved site values: ', OUTPUT_SITE))
message('INFO: Future climate merge complete')
