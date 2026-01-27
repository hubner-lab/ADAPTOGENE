library(qs)
library(data.table)
library(dplyr)
library(raster)
library(stringr)
library(gradientForest)
args = commandArgs(trailingOnly=TRUE)
###############
GF_PATH = args[1]
PREDICTORS_SELECTED = args[2] %>% str_split(',') %>% unlist
FUTURE_ALL = args[3]
PRESENT_ALL = args[4]
PRESENT_RASTER = args[5]
SAMPLES = args[6]
OUTPUT_RASTER = args[7]
OUTPUT_MAP_VALUES = args[8]
OUTPUT_SITE_VALUES = args[9]
###############

message('INFO: Calculating Genetic Offset')

# Load GF model
gf <- qread(GF_PATH)

# Load climate data (select only predictors + ID)
future_climate_all <- fread(FUTURE_ALL) %>% dplyr::select(ID, !!PREDICTORS_SELECTED)
present_climate_all <- fread(PRESENT_ALL) %>% dplyr::select(ID, !!PREDICTORS_SELECTED)

message(paste0('INFO: Present climate cells: ', nrow(present_climate_all)))
message(paste0('INFO: Future climate cells: ', nrow(future_climate_all)))

# Predict cumulative importance for present and future
pred <- predict(gf, present_climate_all[, -1])
pred.future <- predict(gf, future_climate_all[, -1])

# Calculate genetic offset: GO = sqrt(sum((future - present)^2))
genetic_offset <- sapply(1:ncol(pred.future), function(i) {
  (pred.future[, i] - pred[, i])^2
}) %>%
  rowSums(na.rm = TRUE) %>%
  sqrt()

message(paste0('INFO: Genetic offset range: ', round(min(genetic_offset, na.rm = TRUE), 4),
               ' - ', round(max(genetic_offset, na.rm = TRUE), 4)))

# Create genetic offset raster
clim_present <- brick(PRESENT_RASTER)
rast_offset <- clim_present[[PREDICTORS_SELECTED[1]]]
rast_offset[present_climate_all$ID] <- genetic_offset

# Save offset raster
writeRaster(rast_offset,
            filename = OUTPUT_RASTER,
            format = "raster",
            overwrite = TRUE,
            options = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
message(paste0('INFO: Saved offset raster: ', OUTPUT_RASTER))

# Save GO values for whole map (matrix format)
GO_values <- matrix(values(rast_offset), ncol = ncol(rast_offset), byrow = TRUE)
GO_values %>%
  as.data.table %>%
  fwrite(OUTPUT_MAP_VALUES, sep = '\t', col.names = FALSE)
message(paste0('INFO: Saved map GO values: ', OUTPUT_MAP_VALUES))

# Save GO values per sampling site
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))
coords <- data.frame(long = samples$longitude, lat = samples$latitude)
coordinates(coords) <- ~long + lat

go_site <- raster::extract(rast_offset, coords) %>%
  as.data.frame %>%
  cbind(samples[, c('site', 'sample')] %>% unique, .) %>%
  setNames(c('site', 'sample', 'genetic_offset')) %>%
  dplyr::arrange(desc(genetic_offset))

go_site %>%
  fwrite(OUTPUT_SITE_VALUES, sep = '\t')
message(paste0('INFO: Saved site GO values: ', OUTPUT_SITE_VALUES))
message('INFO: Genetic offset calculation complete')
