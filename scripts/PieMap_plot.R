library(qs)
library(data.table)
library(dplyr)
library(see)
library(ggplot2)
library(raster)
library(ggspatial)
library(scatterpie)
library(magrittr)
args = commandArgs(trailingOnly=TRUE)
# Print the command that called the R script
					message(paste("Rscript called with:", commandArgs()[1], paste(args, collapse = " ")))
#################################
SAMPLES = args[1]
CLUSTERS = args[2]
RASTER_LAYER = args[3] # grd file
IBD = args[4]
TAJIMA = args[5]
PI_DIVERSITY = args[6]
CUSTOM_TRAIT = args[7]

PALLETE_MAP = args[8] %>% as.numeric# num of pallete
PALLETE_MAP_reverse = args[9] %>% as.logical # reverse pallete colors or not

PIE_SIZE = args[10] %>% as.numeric
PIE_ALPHA = args[11] %>% as.numeric

IBD_ALPHA = args[12] %>% as.numeric
IBD_COLOR = args[13]

POP_LABELS = args[14] %>% as.logical

PIE_RESCALE = args[15] %>% as.numeric
POP_LABEL_SIZE = args[16] %>% as.numeric
#################################
set.seed(42)
# Choose colors
my.colors = c("#f3c300", "#a1caf1", "#be0032", "#8db600", "#e68fac", "#0067a5", "#f38400", "#222222", 
  "#b3446c", "#2b3d26", "#dcd300",  "#875692", "#848482", "#c2b280", "#882d17", "#654522", "#e25822",
  "#008856", "#f99379", "#604e97", "#f6a600") # from Polychrome pallete
############################# Function

FUN_pie_map <- function(samples,
                        clusters, 
                        raster_layer, 
                        # layer_name = 'bio_12',
                        trait = NULL, # 2 column, site and trait (in this order)
  
  			ibd.dt = NULL, # 2 column dataframe, pop1 and pop2 columns
                        ibd.line_alpha = 0.6,
                        ibd.line_color = 'black',
  
  			colors = wes_palettes[["FantasticFox1"]][2:5], # should match the number of clusters
                        pallete_num = 1022614,
                        pallete_reverse = F,
                        mapcolor_lim = NULL, # c(MIN, MAX)
                        
                        pie.size = 10,
                        pie_alpha = 0.6,
                                            
                        
                        pie_rescale = 1,
                        pop.labels = F,
			pop.label.size = 5,
                        legend_map_title = 'Precipitation',
                        legend_pie_title = 'Clusters',
                        legend_position_pie.x = 0,
                        legend_position_pie.y = 0,
                        legend_position = 'right'
                        ) {
  # browser()
  
  if(is.null(trait)){ # create blank trait if there is no trait
    trait = data.frame(site = clusters$site,
                       trait = pie.size)
  } 
  
  # Process cluster info by population
  clusters.by.pop <-
    clusters %>%
    dplyr::select(-sample) %>% 
    dplyr::group_by(site) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::left_join(samples[,c('site', 'latitude', 'longitude')] %>% unique, by = 'site') %>%
    dplyr::left_join(trait, by = 'site') %>% # Join with trait
    dplyr::rename(trait = colnames(trait)[2])
  
  # Pie size scaling adjusting
  # Get map extent from raster
  map_extent <- extent(raster_layer)
  x_range <- map_extent@xmax - map_extent@xmin
  y_range <- map_extent@ymax - map_extent@ymin
  
  # Use mean range as "map size"
  map_size <- mean(c(x_range, y_range))
  
  # Custom range for pie sizes
  min_size <- 0.2
  max_size <- 1
  
  # Rescale trait to [min_size, max_size]
  clusters.by.pop$trait_rescaled <- min_size +
    (clusters.by.pop$trait - min(clusters.by.pop$trait, na.rm = TRUE)) /
    (max(clusters.by.pop$trait, na.rm = TRUE) - min(clusters.by.pop$trait, na.rm = TRUE)) *
    (max_size - min_size)
  # Scale pies relative to map size
  scaling_factor <- map_size * 0.1   # adjust 0.02 for relative size
  clusters.by.pop$trait_scaled <- clusters.by.pop$trait_rescaled * scaling_factor
	#  message(trait %>% str) 
 # message(trait[[2]] %>% range)
  # return(clusters.by.pop)

  # Plot
  gPlot <-
    ggplot() +
    # Mao
    layer_spatial(raster_layer, aes(fill = after_stat(band1))) +
    scale_fill_colorhex_c(na.value = NA, palette = pallete_num, reverse = pallete_reverse, limits = mapcolor_lim) +
    labs(fill = legend_map_title) +
    
    ggnewscale::new_scale_fill() + # refresh color scale
    
    # Pop structyre
      geom_scatterpie(data = clusters.by.pop, 
                      aes(x = longitude, y = latitude, group = site
                          , r = trait_scaled # trait_scaled if draw equal to smth (once make for comparison between species, will not be basic functionality)
                          ),
                          cols = grep('C[0-9]+', colnames(clusters), value = T),
                          color = 'black',
                          alpha = pie_alpha
                          #,pie_scale = pie.size
                          ) + 
      scale_fill_manual(values = colors) +
      theme_bw() +
      labs(fill  = legend_pie_title) +
      theme(legend.position = legend_position) #+
      # geom_scatterpie_legend(clusters.by.pop$trait, x = legend_position_pie.x, y = legend_position_pie.y)
  
  # Add labels
  if(pop.labels){ 
    
    gPlot <-
      gPlot + 
        geom_text(data = clusters.by.pop, 
                  aes(x = longitude, y = latitude, label = site), 
                  alpha = pie_alpha,
		  color = 'white',
		  size = pop.label.size
	) 
    
    }
  
  if(!is.null(ibd.dt)){
    # IBD/IBE
    gPlot <-
      gPlot +
      lapply(1:nrow(ibd.dt), function(i){
          # geom_encircle(data = clusters.by.pop %>% dplyr::filter(site %in% pair), 
                        # aes(x = longitude, y = latitude),
                        # color = 'black', spread = 0.0000001, alpha = 0.2, expand = 0.001)
          geom_line(data = clusters.by.pop %>% dplyr::filter(site %in% c(ibd.dt[i,]$pop1, ibd.dt[i,]$pop2) ),
                    aes(x = longitude, y = latitude),
                    color = ibd.line_color,
                    alpha = ibd.line_alpha
                    )
      })
  }
  
  return(gPlot)
}

################################# MAIN
# Load raster stack from qs object
raster_layer <- brick(RASTER_LAYER)
					message(raster_layer %>% str)
# Load SAMPLES
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character", 
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))
# Load CLUSTERS
clusters <- fread(CLUSTERS,
                  colClasses = c("site" = "character"))
Nclust = nrow(clusters) - 2
# Load IBD
ibd.dt <- fread(IBD,
		colClasses = c("pop1" = "character",
			       "pop2" = "character"))
# Load TajimaD
tajima <- fread(TAJIMA,
		colClasses = c("site" = "character"))

# Load Pi diversity 
pi_diversity <- fread(PI_DIVERSITY,
		      colClasses = c("site" = "character"))

if(file.exists(CUSTOM_TRAIT)){
	custom_trait <- fread(CUSTOM_TRAIT,
		      	      colClasses = c("site" = "character"))
}
# Plot PieMap
# Iterate through all layer names and perform an operation
for (bio in names(raster_layer)) {
  rlayer <- raster_layer[[bio]]  # Extract the layer by name
   
   # No trait
   
   gNotrait <-
	 FUN_pie_map(samples, 
	 	    clusters, 
		    rlayer,
		    
		    ibd.dt = ibd.dt,
		    ibd.line_alpha = IBD_ALPHA,
		    ibd.line_color = IBD_COLOR,
		    
		    trait = NULL,
		    pie.size = PIE_SIZE,
		    pie_alpha = PIE_ALPHA,
		    colors = my.colors[1:Nclust],
	
		    pallete_num = PALLETE_MAP,
		    pallete_reverse = PALLETE_MAP_reverse,

		    pop.labels = POP_LABELS 
	#TODO make legend properly, maybe on top
	)
   ggsave(paste0("plots/PieMap_", bio, '.png'), gNotrait)
   ggsave(paste0("plots/PieMap_", bio, '.svg'), gNotrait,
   	device = svglite::svglite,
   	bg = 'transparent',
   	fix_text_size = F)
   qsave(gNotrait, paste0("intermediate/PieMap_", bio, '_Notrait.qs'))
   
   # Tajima
   gTajima <-
	 FUN_pie_map(samples, 
	 	    clusters, 
		    rlayer,
		    
		    ibd.dt = ibd.dt,
		    ibd.line_alpha = IBD_ALPHA,
		    ibd.line_color = IBD_COLOR,
		    
		    trait = tajima,
		    pie.size = PIE_SIZE,
		    pie_alpha = PIE_ALPHA,
		    colors = my.colors[1:Nclust],
	
		    pallete_num = PALLETE_MAP,
		    pallete_reverse = PALLETE_MAP_reverse,

		    pop.labels = POP_LABELS 
	#TODO make legend properly, maybe on top
	)
   ggsave(paste0("plots/PieMap_", bio, '_TajimaD.png'), gTajima)
   ggsave(paste0("plots/PieMap_", bio, '_TajimaD.svg'), gTajima,
   	device = svglite::svglite,
   	bg = 'transparent',
   	fix_text_size = F)
   qsave(gTajima, paste0("intermediate/PieMap_", bio, '_TajimaD.qs'))

   # Diversity
   gDiversity <-
         FUN_pie_map(samples,
                    clusters,
                    rlayer,

                    ibd.dt = ibd.dt,
                    ibd.line_alpha = IBD_ALPHA,
                    ibd.line_color = IBD_COLOR,

                    trait = pi_diversity,
                    pie.size = PIE_SIZE,
                    pie_alpha = PIE_ALPHA,
                    colors = my.colors[1:Nclust],

                    pallete_num = PALLETE_MAP,
                    pallete_reverse = PALLETE_MAP_reverse,

		    pop.labels = POP_LABELS 
        #TODO make legend properly, maybe on top
        )

   ggsave(paste0("plots/PieMap_", bio, '_PiDiversity.png'), gDiversity)
    
   ggsave(paste0("plots/PieMap_", bio, '_PiDiversity.svg'), gDiversity,
	       device = svglite::svglite,
	       bg = "transparent",       # прозрачный фон (по желанию))
	       fix_text_size = F
		)

   qsave(gDiversity, paste0("intermediate/PieMap_", bio, '_PiDiversity.qs'))

   if(file.exists(CUSTOM_TRAIT)){
   
	gCustom <-
		FUN_pie_map(samples,
			    clusters,
			    rlayer,
			    ibd.dt = ibd.dt,
			    ibd.line_alpha = IBD_ALPHA,
			    ibd.line_color = IBD_COLOR,
			    trait = custom_trait,
			    pie.size = PIE_SIZE,
			    pie_alpha = PIE_ALPHA,
			    colors = my.colors[1:Nclust],
			    pallete_num = PALLETE_MAP,
			    pallete_reverse = PALLETE_MAP_reverse,
			    pie_rescale = PIE_RESCALE,
                            pop.labels = POP_LABELS,
			    pop.label.size = POP_LABEL_SIZE
		)

	ggsave(paste0("plots/PieMap_", bio, '_', custom_trait %>% dplyr::select(-site) %>% colnames,'.png'), gCustom)
	ggsave(paste0("plots/PieMap_", bio, '_', custom_trait %>% dplyr::select(-site) %>% colnames,'.svg'), gCustom, 
	       device = svglite::svglite,
	       bg = "transparent",       # прозрачный фон (по желанию))
	       fix_text_size = F
		)
	qsave(gCustom, paste0("intermediate/PieMap_", bio, '_', custom_trait %>% dplyr::select(-site) %>% colnames,'.qs'))
   }
}
