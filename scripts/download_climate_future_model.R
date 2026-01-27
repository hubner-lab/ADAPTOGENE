library(geodata)
library(raster)
library(data.table)
library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)
#################################
SAMPLES = args[1]
CROP_REGION = args[2]
GAP = args[3] %>% as.numeric
INDIR = args[4]
SSP = args[5]
YEAR = args[6]
MODEL = args[7]  # single model name
RESOLUTION = args[8] %>% as.numeric
OUTPUT = args[9]  # output .grd file path
#################################

# Load samples for crop region
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# Define crop region
if (CROP_REGION == "auto") {
  min_lat <- min(samples$latitude) - GAP
  max_lat <- max(samples$latitude) + GAP
  min_long <- min(samples$longitude) - GAP
  max_long <- max(samples$longitude) + GAP
  crop_region <- c(min_long, max_long, min_lat, max_lat)
} else {
  crop_region <- CROP_REGION %>% str_split(',') %>% unlist %>% as.numeric
}

message(paste0('INFO: Downloading future climate for model: ', MODEL))
message(paste0('INFO: SSP=', SSP, ', YEAR=', YEAR, ', Resolution=', RESOLUTION))
message(paste0('INFO: Crop region: minLon=', crop_region[1], ', maxLon=', crop_region[2],
               ', minLat=', crop_region[3], ', maxLat=', crop_region[4]))

# Download CMIP6 bioclimatic data for this model
clim_model <- cmip6_world(
  lon = crop_region[1:2],
  lat = crop_region[3:4],
  model = MODEL,
  ssp = SSP,
  time = YEAR,
  var = 'bioc',
  res = RESOLUTION,
  path = INDIR
) %>%
  raster::crop(., crop_region) %>%
  setNames(
    names(.) %>%
      str_split('_') %>%
      lapply(function(x) x %>% rev %>% .[1]) %>%
      unlist %>%
      paste0('bio_', .)
  ) %>%
  stack %>%
  raster::crop(., crop_region)

message(paste0('INFO: Download complete for model: ', MODEL))
message(paste0('INFO: Layers: ', paste(names(clim_model), collapse = ', ')))

# Save per-model raster
writeRaster(clim_model,
            filename = OUTPUT,
            format = "raster",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

message(paste0('INFO: Saved raster to: ', OUTPUT))
