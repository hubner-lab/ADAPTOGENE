library(terra)
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
RASTER_DIR = args[6]    # for raster output
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

# Download WorldClim data directly from geodata.ucdavis.edu
FUN_download_worldclim <- function(data_dir, resolution = 0.5) {
  # Convert resolution to WorldClim format
  resolution_format <- switch(as.character(resolution),
                              '0.5' = "30s",
                              '2.5' = "2.5m",
                              '5'  = "5m",
                              '10' = "10m")

  tif_dir <- file.path(data_dir, paste0("wc2.1_", resolution_format))
  zip_file <- file.path(data_dir, paste0("wc2.1_", resolution_format, "_bio.zip"))

  # Check if tif files already exist
  expected_tifs <- file.path(tif_dir, paste0("wc2.1_", resolution_format, "_bio_", 1:19, ".tif"))

  if (all(file.exists(expected_tifs))) {
    message("INFO: All bioclimatic tif files already exist, skipping download")
    return(tif_dir)
  }

  # Create directory if needed
  if (!dir.exists(tif_dir)) {
    dir.create(tif_dir, recursive = TRUE)
  }

  # Download zip if not exists
  if (!file.exists(zip_file)) {
    url <- paste0("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_",
                  resolution_format, "_bio.zip")
    message(paste0("INFO: Downloading WorldClim bioclimatic data from: ", url))
    old_timeout <- getOption("timeout")
    options(timeout = 600)
    on.exit(options(timeout = old_timeout), add = TRUE)
    download.file(url, zip_file, mode = "wb", quiet = FALSE)
  } else {
    message("INFO: Zip file already exists, skipping download")
  }

  # Unzip using system unzip (R's unzip can fail on large files)
  message("INFO: Extracting tif files...")
  exit_code <- system2("unzip", args = c("-o", "-j", shQuote(zip_file), "-d", shQuote(tif_dir)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(exit_code, "status")) && attr(exit_code, "status") != 0) {
    # Fallback to R unzip
    warning("System unzip failed, trying R unzip...")
    unzip(zip_file, exdir = tif_dir)
  }

  return(tif_dir)
}

FUN_present_climate <- function(samples,
                                data_dir,
                                crop_region,
                                resolution = 0.5) {

  # Download/check WorldClim data
  tif_dir <- FUN_download_worldclim(data_dir, resolution)

  # Construct folder path based on resolution
  resolution_format <- switch(as.character(resolution),
                              '0.5' = "30s",
                              '2.5' = "2.5m",
                              '5'  = "5m",
                              '10' = "10m")

  clim.list <- dir(tif_dir,
                   full.names = TRUE,
                   pattern = paste0("^wc2\\.1_", resolution_format, "_bio_[0-9]+\\.tif$"))
  message(clim.list %>% str)

  if (length(clim.list) == 0) {
    stop(paste0("ERROR: No bioclimatic tif files found in ", tif_dir))
  }

  # Load as SpatRaster
  clim.layer <- rast(clim.list)
  message('INFO: Present climate SpatRaster loaded')
  message(samples %>% str)

  # Create SpatVector for site coordinates
  coords <- vect(data.frame(long = samples$longitude, lat = samples$latitude),
                 geom = c("long", "lat"), crs = "EPSG:4326")

  # Crop
  if (length(crop_region) == 1 && crop_region == 'world') {
    clim.present <- clim.layer
  } else {
    clim.present <- terra::crop(clim.layer, ext(crop_region))
  }

  names(clim.present) <- gsub(paste0("wc2\\.1_", resolution_format, "_"), "", names(clim.present))

  # Extract data for coordinates
  clim.present.coords <- terra::extract(clim.present, coords)[, -1] %>%
    as.data.frame

  # Extract all pixel values
  clim.land <- terra::extract(clim.present, 1:ncell(clim.present)) %>%
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

# Validate: check for samples with NA climate values
site_values <- clim_present$SiteValues
na_rows <- which(rowSums(is.na(site_values)) > 0)
if (length(na_rows) > 0) {
    bad_samples <- samples[na_rows, ]
    message("ERROR: Climate extraction returned NA for the following samples:")
    for (i in seq_len(nrow(bad_samples))) {
        msg <- if (bad_samples$latitude[i] == 0 && bad_samples$longitude[i] == 0) {
            "likely placeholder (0,0)"
        } else {
            "falls on ocean/NoData pixel"
        }
        message(paste0("  - ", bad_samples$sample[i], " (", bad_samples$site[i],
                       ") at (", bad_samples$latitude[i], ", ", bad_samples$longitude[i], "): ", msg))
    }
    stop(paste0("Climate extraction failed for ", length(na_rows), " samples. ",
                "Fix coordinates in input metadata or remove these samples."))
}

# Save raster
writeRaster(clim_present$RasterStack,
            filename = paste0(RASTER_DIR, 'climate_present_rasterstack.tif'),
            overwrite = TRUE,
            gdal = c("INTERLEAVE=BAND", "COMPRESS=LZW"))

# Save all values
clim_present$AllValues %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_all.tsv'), sep = '\t')

# Save site values (with sample column for traceability)
cbind(sample = samples$sample, site_values) %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_site.tsv'), sep = '\t')

# Save site values scaled (with sample column)
cbind(sample = samples$sample, site_values %>% dplyr::mutate_all(scale)) %>%
  fwrite(paste0(TABLES_DIR, 'climate_present_site_scaled.tsv'), sep = '\t')
