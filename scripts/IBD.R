library(dplyr)
library(data.table)
library(parallel)
library(vegan) # for Mantel test
library(stringr)
library(geosphere) # distm function
args = commandArgs(trailingOnly=TRUE)
set.seed(42)
#######################
CLUSTERS = args[1]
SAMPLES = args[2]
CPU = args[3]
#######################


############################## Function
FUN_ibd <- function(clusters, samples, CPU){

  # browser()
  # IBD refers to the phenomenon where genetic similarity decreases as geographic distance between populations increases.

  # To test for IBD/IBE, we can use the Mantel test to compare a matrix of genetic distances among populations with a matrix of geographic/environment distances among populations. If the two matrices are positively correlated, this would suggest that IBD/IBE is playing a role in shaping the genetic differentiation among populations.
                                                                        message(paste0('INFO: RUN IBD analysis with: n = ', 
										       samples$site %>% unique %>% length, 
										       '; GeoMatrix = ', nrow(samples)))

  ibd <-
    mclapply(samples$site %>% unique, function(i){
      mclapply(samples$site %>% unique, function(i2){

          # Select populations
          clust_df = clusters[clusters$site %in% c(i,i2),] %>% dplyr::select(starts_with('C'))
          geo_df = samples[samples$site %in% c(i,i2),] %>% dplyr::select(latitude, longitude)
#							message(paste0('INFO: clust: ', clust_df %>% str))
#	  						message(paste0('INFO: geo: ', geo_df %>% str))
          clust_dist = vegdist(clust_df, method = 'jaccard')
#							message(paste0('INFO: clust_dist: ', clust_dist %>% str))
          geo_dist = distm(geo_df, fun=distVincentyEllipsoid) # In this case, it uses the Vincenty distance, which is a method for calculating the distance between two points on the surface of an ellipsoid. It is considered to be the most accurate method for calculating distances on the Earth's surface, and it's recommended when the accuracy is important.
#	  						message(paste0('INFO: geo_dist: ',geo_dist %>% str))


          # Perform Mantel test for cheking IBD
          ibd = mantel(clust_dist, geo_dist,
                 method = 'pearson',
                 permutations=999)

          res_dt <- data.table(Comparison = paste0(i,' vs ', i2),
                     IBD_pval = ibd$signif,
                     IBD_stat = ibd$statistic)

	  return(res_dt)
      }, mc.cores = CPU, mc.allow.recursive = T) %>% do.call(rbind, .)
    }, mc.cores = CPU, mc.allow.recursive = T) %>% do.call(rbind, . ) %>%
    na.omit								
								message("INFO: IBD calculation completed, process the result")
								message(ibd %>% str)
  # Look on not isolated:
  ibd.list <-
    ibd[ibd$IBD_pval > 0.05,] %>%
    .$Comparison %>%
    str_split(' vs ') %>%
    lapply(sort) %>%
    unique

  return(list(raw = ibd,
              not_ibd = ibd.list)
         )
}
############################ MAIN
# Load tables
#TODO here we filtered out all populations which contain only 1 individual, due to impossibility of calculate IBD metric for it
samples <- fread(SAMPLES, 
		colClasses = c("site" = "character", 'sample' = 'character')) %>%
	dplyr::group_by(site) %>%
	dplyr::filter(n() > 1) %>%
	dplyr::ungroup()
clusters <- fread(CLUSTERS,
		  colClasses = c("site" = "character")) %>%
	dplyr::group_by(site) %>%
	dplyr::filter(n() > 1) %>%
	dplyr::ungroup()
# Calculate IBD
ibd <- FUN_ibd(clusters, samples, CPU)

# Save
ibd$raw %>%
	fwrite('tables/IBD_raw.tsv', sep = '\t')
				#			message(ibd$not_ibd %>% str)
ibd$not_ibd %>%	
	do.call(rbind, . ) %>%
	as.data.table %>%
	dplyr::rename(pop1 = V1, pop2 = V2) %>%
	fwrite('tables/IBD_notIsolated.tsv', sep = '\t') # not isolated 

#TODO plot heatmap
