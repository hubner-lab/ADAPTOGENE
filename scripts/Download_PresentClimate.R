library(geodata)
library(raster)
library(ggplot2)
library(data.table)
library(dplyr)
library(qs)
library(stringr)
args = commandArgs(trailingOnly=TRUE)
#################################
SAMPLES = args[1]
CROP_REGION = args[2]
GAP = args[3] %>% as.numeric
INTERMEDIATE = args[4] # data dir for storing downloaded climate data
RESOLUTION = args[5] %>% as.numeric
#################################### Functions
# Plots settings
plot_theme <-
        theme_classic(base_size=12, base_family = 'Helvetica') +
             theme(plot.title = element_text(size = 12, face = 'bold'),
                   axis.text = element_text(face = 'bold', size = 12),
                   axis.title = element_text(face = "bold", size = 16),
                   legend.text = element_text(size = 15),
                   legend.title = element_text(face = "bold", size = 15)
                   )

FUN_define_crop_region <- function(samples, gap = GAP) {

  # Calculate min/max lat and long with the gap
  min_lat <- min(samples$latitude) - gap
  max_lat <- max(samples$latitude) + gap
  min_long <- min(samples$longitude) - gap
  max_long <- max(samples$longitude) + gap

  # Return the crop region as a list or named vector
  return(c( min_long, max_long, min_lat, max_lat))
}

FUN_present_climate <- function(samples, # SpatialPoints
				INTERMEDIATE,
                                crop_region, # long long lat lat
				resolution = 0.5, # Valid resolutions are 10, 5, 2.5, and 0.5 (minutes of a degree)
                                PCs = NULL){ #  number of main PCs and NULL if raw climate needed

  # browser()
  # Download world map
  worldclim_global('bio', resolution, INTERMEDIATE) 
  # specify path
   # Construct folder path based on resolution
  resolution_format <- switch(resolution %>% as.character,
                       '0.5' = "30s",
                       '2.5' = "2.5m",
                       '5'  = "5m",
                       '10' = "10m")
  #TODO change paths for different resolutions
  clim.list <- dir(paste0(INTERMEDIATE,'/climate/wc2.1_', resolution_format), 
		   full.names=T, 
		   pattern = paste0("^wc2\\.1_", resolution_format, "_bio_[0-9]+\\.tif$"))
  											message(clim.list %>% str)
  # Stack from tif
  clim.layer <-  stack(clim.list)
  								message('INFO: Present climate stack of tif files loaded')
  								message(samples %>% str)
  # crop
  ## Extract coordinates
  coords <-
	  data.frame(long = samples$longitude,
        	     lat = samples$latitude)
  coordinates(coords)= ~ long + lat

  if(length(crop_region) == 1 && crop_region == 'world'){ # very RAM consuming!!!
	clim.present = clim.layer
  } else {
  	clim.present <- raster::crop(clim.layer, crop_region)
  }

  names(clim.present) <- gsub(paste0("wc2\\.1_", resolution_format, "_"), "", names(clim.present))
  # Extract data for coordinates
  clim.present.coords <-
    raster::extract(clim.present, coords) %>%
    as.data.frame
  
  # Extract each pixels values
  clim.land <- 
    raster::extract(clim.present,
                    1:ncell(clim.present), 
                    df = TRUE) %>% 
    na.omit %>%
    setNames(c('ID', clim.present %>% names))

  
  #TODO Second part for PCA 
  # Return bio layers if PCA is NULL
          if(is.null(PCs)){
            
	       	     return(list(clim.present,
        	                clim.present.coords,
        	                clim.land) %>%
                	     setNames(c('RasterStack', 'SiteValues', 'AllValues')))
                   }
          # If PCA is any number return rasterPCA layers
          else{
          
          # rasterPCA
          clim.present.pca <- rasterPCA(clim.present %>% raster::scale(center = T, scale = T) # PCA is affected by scale  #TODO MAKE PCA ON SELECTED VARS ONLY
                              , spca=TRUE)
          # Print summary
          print(summary(clim.present.pca$model))
          # Print loadings of choosed PCs
          print(knitr::kable(round(clim.present.pca$model$loadings[,1:PCs] %>% as.data.frame %>% setNames(paste0('CLIM', 1:PCs)),3)))
          
          # Get first N PCs
          clim.present.pca.coords <-
            raster::extract(clim.present.pca$map, coords)[,1:PCs] %>%  # extract first 3 PCs
            as.data.frame %>%
            setNames(paste0('CLIM', 1:PCs))
          
          clim.land.pca <- 
            raster::extract(clim.present.pca[1:PCs]$map,
                            1:ncell(clim.present.pca[1:PCs]$map), 
                            df = TRUE) %>% 
            na.omit %>%
            setNames(c('ID', clim.present.pca[1:PCs]$map %>% names))
          
          
            return(list(clim.present.pca[1:PCs]$map, 
                        clim.present.pca.coords,
                        clim.land.pca) %>%
                     setNames(c('RasterStack', 'SiteValues', 'AllValues'))) # return list(raster, coords_values)
          }
}

################################### MAIN
# Load Sample info
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character", 
				'sample' = 'character',
				'latitude' = 'numeric',
				'longitude' = 'numeric')) 

# Define crop regions
if(CROP_REGION == "auto"){
	crop_region = FUN_define_crop_region(samples, gap = GAP)
							message(paste0('INFO: Extract climate data in region minLongitude=', crop_region[1], '; maxLongitude=', crop_region[2], '; minLatitude=', crop_region[3], '; maxLatitude=', crop_region[4]))
} else if (CROP_REGION == "world") {
    crop_region = 'world'
} else {
	crop_region = CROP_REGION %>% str_split(',') %>% unlist %>% as.numeric
}
							message(crop_region %>% str) # tmp
# Download and process bioclimatic raster and variables
clim_present <- FUN_present_climate(samples,
				    INTERMEDIATE,
				    crop_region,
				    RESOLUTION) #TODO add PCs option
							message('INFO: Loading climate present data complete, saving to the disk')
							message(clim_present %>% str)
# Save raster
writeRaster(clim_present$RasterStack, 
	    filename=paste0('intermediate/Climate_present_RasterStack.grd'), 
	    format="raster", 
	    overwrite=TRUE,
	    options=c("INTERLEAVE=BAND","COMPRESS=LZW"))
#qsave(clim_present$RasterStack, paste0(INTERMEDIATE, 'Climate_present_RasterStack.qs'))

# Save all values
clim_present$AllValues %>%
	fwrite('tables/Climate_present_all.tsv', sep = '\t')
# Save region values
clim_present$SiteValues %>%
	fwrite('tables/Climate_present_region.tsv', sep = '\t')
# Save region values scaled
clim_present$SiteValues %>%
	dplyr::mutate_all(scale) %>%
	fwrite('tables/Climate_present_region_scaled.tsv', sep = '\t')
# Plot clim factors distributions
for (bio in colnames(clim_present$SiteValues)){
	gDens <-
		clim_present$SiteValues %>%
				ggplot(aes(x = .data[[bio]])) +
					geom_density(fill = 'blue') + 
						labs(title = bio, x = "Values", y = "Density") +
						plot_theme
	ggsave(paste0('plots/DensityPlot_', bio, '.png'), gDens)
	qsave(gDens, paste0('intermediate/DensityPlot_', bio, '.qs'))
}

