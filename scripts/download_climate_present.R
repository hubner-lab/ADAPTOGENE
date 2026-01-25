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
DATA_DIR = args[4]      # directory for storing downloaded climate data
RESOLUTION = args[5] %>% as.numeric
INTER_DIR = args[6]     # for raster output
TABLES_DIR = args[7]    # for tsv outputs
#################################

#################################### Functions
FUN_define_crop_region <- function(samples, gap) {
  min_lat <- min(samples$latitude) - gap
  max_lat <- max(samples$latitude) + gap
  min_long <- min(samples$longitude) - gap
  max_long <- max(samples$longitude) + gap
  return(c(min_long, max_long, min_lat, max_lat))
}

FUN_present_climate <- function(samples,
                                data_dir,
                                crop_region,
                                resolution = 0.5) {
  
  # Download world map
  worldclim_global('bio', resolution, data_dir)
  
  # Construct folder path based on resolution
  resolution_format <- switch(as.character(resolution),
                              '0.5' = "30s",
                              '2.5' = "2.5m",
                              '5'  = "5m",
                              '10' = "10m")
  
  clim.list <- dir(paste0(data_dir, '/wc2.1_', resolution_format),
  #clim.list <- dir(paste0(data_dir, '/climate/wc2.1_', resolution_format), # depending on the version of geodata!
                   full.names = TRUE,
                   pattern = paste0("^wc2\\.1_", resolution_format, "_bio_[0-9]+\\.tif$"))
  message(clim.list %>% str)
  
  # Stack from tif
  clim.layer <- stack(clim.list)
  message('INFO: Present climate stack of tif files loaded')
  message(samples %>% str)
  
  # Extract coordinates
  coords <- data.frame(long = samples$longitude,
                       lat = samples$latitude)
  coordinates(coords) <- ~ long + lat
  
  # Crop
  if (length(crop_region) == 1 && crop_region == 'world') {
    clim.present <- clim.layer
  } else {
    clim.present <- raster::crop(clim.layer, crop_region)
  }
  
  names(clim.present) <- gsub(paste0("wc2\\.1_", resolution_format, "_"), "", names(clim.present))
  
  # Extract data for coordinates
  clim.present.coords <- raster::extract(clim.present, coords) %>%
    as.data.frame
  
  # Extract each pixel values
  clim.land <- raster::extract(clim.present,
                               1:ncell(clim.present),
                               df = TRUE) %>%
    na.omit %>%
    setNames(c('ID', names(clim.present)))
  
  return(list(RasterStack = clim.present,
              SiteValues = clim.present.coords,
              AllValues = clim.land))
}

################################### MAIN
# Load Sample info
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# Define crop regions
if (CROP_REGION == "auto") {
  crop_region <- FUN_define_crop_region(samples, gap = GAP)
  message(paste0('INFO: Extract climate data in region minLongitude=', crop_region[1],
                 '; maxLongitude=', crop_region[2],
                 '; minLatitude=', crop_region[3],
                 '; maxLatitude=', crop_region[4]))
} else if (CROP_REGION == "world") {
  crop_region <- 'world'
} else {
  crop_region <- CROP_REGION %>% str_split(',') %>% unlist %>% as.numeric
}

message(crop_region %>% str)

# Download and process bioclimatic raster and variables
clim_present <- FUN_present_climate(samples,
                                    DATA_DIR,
                                    crop_region,
                                    RESOLUTION)

message('INFO: Loading climate present data complete, saving to the disk')
message(clim_present %>% str)

# Save raster
writeRaster(clim_present$RasterStack,
            filename = paste0(INTER_DIR, 'climate_present.grd'),
            format = "raster",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

# Save all values
clim_present$AllValues %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_all.tsv'), sep = '\t')

# Save site values
clim_present$SiteValues %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_site.tsv'), sep = '\t')

# Save site values scaled
clim_present$SiteValues %>%
  dplyr::mutate_all(scale) %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_site_scaled.tsv'), sep = '\t')
