library(stringr)
library(tibble)
library(qs)
library(data.table)
library(dplyr)
library(see)
library(ggplot2)
library(raster)
library(ggspatial)
library(scatterpie)
library(gradientForest)
library(ggpubr)
library(forcats)
library(svglite)
args = commandArgs(trailingOnly=TRUE)
###############
gf = qread(args[1])
gf_random = qread(args[2])
PREDICTORS_SELECTED = args[3] %>% str_split(',') %>% unlist
future_climate_all = fread(args[4]) %>% dplyr::select(ID, !!PREDICTORS_SELECTED) 
present_climate_all = fread(args[5]) %>% dplyr::select(ID, !!PREDICTORS_SELECTED)
CLIM_LAYER = args[6] ; clim.present <- brick(CLIM_LAYER) # grd file
SAMPLES = args[7]
CLUSTERS = args[8]
IBD = args[9]
TAJIMA = args[10]
PALLETE_MAP = args[11] %>% as.numeric # num of pallete
PIE_SIZE = args[12] %>% as.numeric
PIE_ALPHA = args[13] %>% as.numeric
IBD_ALPHA = args[14] %>% as.numeric
IBD_COLOR = args[15]
# For out file labelling
PCNM = args[16]
GF_SUFFIX = args[17]
###############
plot_theme = theme_classic(base_size=12, base_family = 'Helvetica') + 
             theme(plot.title = element_text(size = 17, face = 'bold'),
                   axis.text = element_text(face = "bold", size = 18),
                   axis.title = element_text(face = "bold", size = 22),
                   legend.text = element_text(size = 18),
                   legend.title = element_text(face = "bold", size = 22)
                   )

SUFFIX = paste0(GF_SUFFIX, '_', PCNM, 'PCNM')
###################### Function
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

                        rescale.pie = c(0.3, 1, 10),
                        pop.labels = F,
                        legend_map_title = 'DefaultTitle',
                        legend_pie_title = 'Clusters',
                        legend_position_pie.x = 0,
                        legend_position_pie.y = 0,
                        legend_position = 'right'
                        ) {
  # browser()

  if(is.null(trait)){ # create blank trait if there is no trait
    trait = data.frame(site = clusters$site,
                       trait = pie.size)
  } else {trait[[2]] <- scales::rescale(trait[[2]], c(rescale.pie[1], rescale.pie[2])) / rescale.pie[3] } #  GS / TAJIMA

  # Process cluster info by population
  clusters.by.pop <-
    clusters %>%
    dplyr::select(-sample) %>%
    dplyr::group_by(site) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::left_join(samples[,c('site', 'latitude', 'longitude')] %>% unique, by = 'site') %>%
    dplyr::left_join(trait, by = 'site') %>% # Join with trait
    dplyr::rename(trait = colnames(trait)[2])

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
                          , r = trait
                          ),
                          cols = grep('C[0-9]+', colnames(clusters), value = T),
                          color = 'black',
                          alpha = pie_alpha
                          #,pie_scale = pie.scale
                          ) +
      scale_fill_manual(values = colors) +
      # ID
      # geom_text(data = clusters.by.pop,
                # aes(x = longitude, y = latitude, label = site),
                # color = 'white',
                # size = 2) +
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
                  alpha = pie_alpha)

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

###################### MAIN
# CumImp
											message("INFO: Plotting CumImp plot")
gList <-
  lapply(PREDICTORS_SELECTED, function(bio){
    gf_random.cumimp <- cumimp(gf_random, predictor=bio, type="Overall")

    gf.cumimp <- cumimp(gf, predictor=bio, type="Overall" )



    df = data.frame(x = c(gf_random.cumimp$x,
                          gf.cumimp$x ),
                  y = c(gf_random.cumimp$y
                        ,gf.cumimp$y ),
                  Model = c(rep('Neutral', length(gf_random.cumimp$x)),
                            rep('Adaptive', length(gf.cumimp$x)) ) )
                # Plot
    gPlot <-
      df %>%
        ggplot(aes(x = x, y = y, color = Model)) +
          geom_line() +
          labs(x = bio, y = 'Cumulative importance') +
          scale_color_manual(values = c('red', 'blue', 'green' #TODO remember why this colors
                                             , 'black'
                                             #, 'brown', 'orange'
                                             # ,'pink', 'orange', 'yellow', 'purple'
                                             )) +
          plot_theme +
      theme(axis.text = element_text(face = NULL, size = 12),
            axis.title.y = element_text(face = NULL, size = 12),
            axis.title.x = element_text(face = NULL, size = 16),
            legend.title = element_blank())

    return(gPlot)
	})
						#	message(gList %>% str)
gCumImp <- ggarrange(plotlist = gList, nrow = 3, ncol = 3, common.legend = T) #TODO make automatic calculation of plot size

ggsave(paste0('plots/gradientForest_CumulativeImportance_', SUFFIX, '.png'), 
       gCumImp)
ggsave(paste0('plots/gradientForest_CumulativeImportance_', SUFFIX, '.svg'), 
       gCumImp,
       device = svglite::svglite,
       bg = "transparent",       # прозрачный фон (по желанию))
       fix_text_size = F
	)
qsave(gCumImp, paste0('intermediate/gradientForest_CumulativeImportance_', SUFFIX, '.qs'))

# Importance
									message("INFO: Plotting Overall Importance plot")
imp = importance(gf, type = 'Weighted') %>%
  enframe() %>%
  dplyr::mutate(name = as.factor(name))

gAdapt <-
  ggplot(imp, aes(y = fct_reorder(name, value), x = value)) +
    geom_bar(stat = 'identity', fill = 'red', color = 'black') +
    plot_theme +
    labs(x = expression(paste("R"^2, " weighted importance ")),
         title = 'Adaptive') +
    xlim(c(0, imp$value %>% max + 0.005)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(colour = 'black', size = 12),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 14))

imp_random = importance(gf_random, type = 'Weighted') %>%
  enframe() %>%
  dplyr::mutate(name = as.factor(name))

gAll <-
  ggplot(imp_random, aes(y = fct_reorder(name, value), x = value)) +
    geom_bar(stat = 'identity', fill = 'blue', color = 'black') +
    plot_theme +
    labs(x = expression(paste("R"^2, " weighted importance ")),
         title = 'Neutral') +
    xlim(c(0, imp$value %>% max + 0.005))  +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(colour = 'black', size = 12),
          axis.text.x = element_text(size = 10),
          axis.title.x = element_text(size = 14))

gImp <- ggpubr::ggarrange(gAdapt, gAll, ncol = 2)
ggsave(paste0('plots/gradientForest_OverallImportance_', SUFFIX, '.png'), 
       gImp)
ggsave(paste0('plots/gradientForest_OverallImportance_', SUFFIX, '.svg'), 
       gImp,
       device = svglite::svglite,
       bg = "transparent",       # прозрачный фон (по желанию))
       fix_text_size = F
	)
qsave(gImp, paste0('intermediate/gradientForest_OverallImportance_', SUFFIX, '.qs'))

# GeneticOffset
										message("INFO: Calculating Genetic Offset")
										message(present_climate_all %>% str)
										message(future_climate_all %>% str)
# Predict
pred <- predict(gf, present_climate_all[,-1])
pred.future <- predict(gf, future_climate_all[,-1])
										message(pred %>% str)
										message(pred.future %>% str)
# Genetic offset
genetic.offset.adaptive <- 
	sapply(1:ncol(pred.future), function(i) { (pred.future[,i] - pred[,i])^2 }) %>% 
	rowSums(na.rm = T) %>% 
	sqrt()
# Genetic offset
#genetic.offset.adaptive <- sqrt((pred.future[,1]-pred[,1])^2 + 
#                                 (pred.future[,2]-pred[,2])^2 +
#                                  (pred.future[,3]-pred[,3])^2 + 
#                                  (pred.future[,4]-pred[,4])^2 +
#                                  (pred.future[,5]-pred[,5])^2 +
#                                  (pred.future[,6]-pred[,6])^2 +
#                                  (pred.future[,7]-pred[,7])^2
#                                )

#Define raster properties
rast.offset <- clim.present[[PREDICTORS_SELECTED[1]]]
										message(clim.present %>% str)
										message(clim.present %>% names)
										message('RAST')
										message(rast.offset %>% str)
										message(rast.offset %>% names)
										message(genetic.offset.adaptive %>% str)
#Assign genetic offset values (difference between future and present predictions) to raster
rast.offset[present_climate_all$ID] <- genetic.offset.adaptive

# SAVE GO values for whole map
GO_values <- matrix(values(rast.offset), ncol = ncol(rast.offset), byrow = TRUE)
										message(GO_values %>% str)
GO_values %>%
	as.data.table %>%
	fwrite(paste0('tables/gradientForest_GeneticOffesetValues_', SUFFIX, '.tsv'),
		sep = '\t',
		col.names = F)

# SAVE GO values
#raster::extract(rast.offset.585, coords) %>%  # extract first 3 PCs
#    as.data.frame %>%
#    cbind(geo.df$site, .) %>%
#    unique %>%
#    setNames(c('site', 'GO.585')) %>%
#    dplyr::arrange(GO.585) #%T>%
#  # fwrite('../results/GO.585_values.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
										message("INFO: Plotting Genetic Offset")
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
# Load SAMPLES
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character", 
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))
# Plot
gPieMapGO <-
       FUN_pie_map(samples, 
                   clusters, 
                   rast.offset, 
                   
		   ibd.dt = ibd.dt,
                   ibd.line_alpha = IBD_ALPHA,
                   ibd.line_color = IBD_COLOR,
                    
                   trait = tajima,
                   pie.size = PIE_SIZE,
                   pie_alpha = PIE_ALPHA,
                   colors = my.colors[1:Nclust],
		   
		   legend_map_title = 'GeneticOffset',

                   pallete_num = PALLETE_MAP,
                   pallete_reverse = T, #TODO

		   #TEMP
		   #mapcolor_lim = c(0, 0.07)
			)
ggsave(paste0('plots/gradientForest_GeneticOffesetPieMap_', SUFFIX, '.png'), 
       gPieMapGO)
ggsave(paste0('plots/gradientForest_GeneticOffesetPieMap_', SUFFIX, '.svg'), 
       gPieMapGO,
       device = svglite::svglite,
       bg = "transparent",       # прозрачный фон (по желанию))
       fix_text_size = F
	)
qsave(gPieMapGO, paste0('intermediate/gradientForest_GeneticOffesetPieMap_', SUFFIX, '.qs'))
