library(geodata)
library(raster)
library(data.table)
library(dplyr)
library(stringr)
library(parallel)
args = commandArgs(trailingOnly=TRUE)
#################################
SAMPLES = args[1]
CROP_REGION = args[2]
GAP = args[3] %>% as.numeric
INDIR = args[4] # data dir for storing downloaded climate data
SSP = args[5] # 126, 245, 370, 585
YEAR = args[6] # year range 2061-2080
MODELS = args[7] %>% str_split(',') %>% unlist # c('MPI-ESM1-2-HR', 'UKESM1-0-LL', 'IPSL-CM6A-LR')
CPU = args[8] # threads of downloading from geodata
RESOLUTION = args[9] %>% as.numeric
#################################
#################################### Functions
FUN_define_crop_region <- function(samples, gap = GAP) {

  # Calculate min/max lat and long with the gap
  min_lat <- min(samples$latitude) - gap
  max_lat <- max(samples$latitude) + gap
  min_long <- min(samples$longitude) - gap
  max_long <- max(samples$longitude) + gap

  # Return the crop region as a list or named vector
  return(c( min_long, max_long, min_lat, max_lat))
}

#TODO make for sophisticated, now just existed
FUN_future_climate <- function(samples, # SpatialPoints 
				INDIR,
                                crop_region, # long long lat lat
                                PCs = NULL,
				ssp = 585, # 126, 245, 370, 585 + ??? 119, 434, 340, 460 (according to https://cds.climate.copernicus.eu/datasets/projections-cmip6?tab=download) !!! TODO check it, check for bioclim which models exists + year ranges (wait while geodata will become to work) !!!
				year = "2061-2080", # 
				models_cmip6 = c('MPI-ESM1-2-HR', 'UKESM1-0-LL', 'IPSL-CM6A-LR'),
				CPU = 1
				,resolution = 0.5 # SKIPPED
				){ #  number of main PCs and NULL if raw climate needed

	#TODO no needed for renaming
	#resolution_format <- switch(resolution %>% as.character,
        #               '0.5' = "30s",
        #               '2.5' = "2.5m",
        #               '5'  = "5m",
        #               '10' = "10m")	
								message(resolution)	
	clim.future_list <-
		mclapply(models_cmip6, function(model){
	
		   cmip6_world(lon = crop_region[1:2],
		              lat = crop_region[3:4],
		            model = model,
		            ssp = ssp,
		            time = year,
		            var = 'bioc',
			    res = resolution, #TODO later add step with aggregating to desired resolution
			    path = INDIR, 
		            quite = F) %>%
		     raster::crop(., crop_region) %>%
		   setNames(names(.) %>% str_split('_') %>% lapply(function(x) x %>% rev %>% .[1]) %>% unlist %>% paste0('bio_', .)) %>%
		   stack %>%
		   raster::crop(., crop_region) # don't needed?
		}, mc.cores = length(models_cmip6)) #TODO add mclapply
										message('INFO: Loading future climate data complete, now processing')

	clim.future <-
	  clim.future_list %>%
	    stack %>%
	    stackApply(rep(1:19, 3), mean, na.rm = T) %>% # should add names?
	    setNames(names(clim.future_list[[1]]))

    	clim.future.values = raster::extract(clim.future, 1:ncell(clim.future), df = T)

	coords <- data.frame(long = samples$longitude,
	                     lat = samples$latitude)
	coordinates(coords)= ~ long + lat

	clim.future.sitevalues = raster::extract(clim.future, coords)  %>% as.data.frame

  	return(list(RasterStack = clim.future,
		    AllValues = clim.future.values,
		    SiteValues = clim.future.sitevalues
		))


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
} else {
		
		crop_region = CROP_REGION %>% str_split(',') %>% unlist %>% as.numeric
}
							message(crop_region %>% str)
							message(paste0('INFO: Extract climate data in region minLongitude=', crop_region[1], '; maxLongitude=', crop_region[2], '; minLatitude=', crop_region[3], '; maxLatitude=', crop_region[4]))

# Download and process bioclimatic raster and variables
clim_future <- FUN_future_climate(samples, 
				  INDIR, 
				  crop_region,
				  PCs = NULL, #TODO
				  ssp = SSP,
				  year = YEAR,
				  models_cmip6 = MODELS,
				  CPU = CPU
				  ,resolution = RESOLUTION
					)
				message('INFO: Loading climate future data complete, saving to the disk')
								message(clim_future %>% str)
# Save raster
writeRaster(clim_future$RasterStack, 
	    filename=paste0('intermediate/Climate_future_year', YEAR, '_ssp', SSP, '_RasterStack.grd'), 
	    format="raster", 
	    overwrite=TRUE,
	    options=c("INTERLEAVE=BAND","COMPRESS=LZW"))

# Save all values
clim_future$AllValues %>%
	fwrite(paste0('tables/Climate_future_year', YEAR, '_ssp', SSP, '_all.tsv' ), sep = '\t')
# Save region values
clim_future$SiteValues %>%
	fwrite(paste0('tables/Climate_future_year', YEAR, '_ssp', SSP, '_site.tsv' ), sep = '\t') #TODO rename on site?
