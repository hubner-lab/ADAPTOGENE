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
TABLES_DIR = args[4]
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
        tryCatch({
          # Select populations
          clust_df = clusters[clusters$site %in% c(i,i2),] %>% dplyr::select(starts_with('C'))
          geo_df = samples[samples$site %in% c(i,i2),] %>% dplyr::select(longitude, latitude)

          clust_dist = vegdist(clust_df, method = 'jaccard')
          geo_dist = distm(geo_df, fun=distVincentyEllipsoid)

          # Perform Mantel test for cheking IBD
          ibd = mantel(clust_dist, geo_dist,
                 method = 'pearson',
                 permutations=999)

          res_dt <- data.table(Comparison = paste0(i,' vs ', i2),
                     IBD_pval = ibd$signif,
                     IBD_stat = ibd$statistic)

	  return(res_dt)
        }, error = function(e) {
          warning(paste("IBD error for", i, "vs", i2, ":", e$message))
          return(NULL)
        })
      }, mc.cores = CPU, mc.allow.recursive = T) %>% Filter(Negate(is.null), .) %>% do.call(rbind, .)
    }, mc.cores = CPU, mc.allow.recursive = T) %>% Filter(Negate(is.null), .) %>% do.call(rbind, . ) %>%
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

#=============================================================================
# SINGLE-SITE GUARD
#=============================================================================
# IBD is a correlation between genetic and geographic distance ACROSS pairs of
# populations. With one site the pairwise loop below produces a single self
# comparison: either mantel() fails on an all-zero geographic distance matrix
# (tryCatch swallows it, leaving an empty result whose rename(V1, V2) then
# errors), or -- if coordinates happen to jitter within the site -- it returns a
# plausible-looking p-value for comparing a population with itself, which is
# worse. Note this is a count of SITES, distinct from the n() > 1 filters above,
# which count samples within a site.
#
# Write both declared outputs and exit 0: population stats are an opt-in extra
# and must not take the whole structure mode down with them.
n_sites <- length(unique(samples$site))
if (n_sites < 2) {
	message(paste0('WARNING: IBD skipped -- ', n_sites,
		       ' sampling site(s) with more than one sample, need >= 2. ',
		       'Isolation by distance is undefined without population pairs.'))
	data.table(Comparison = character(0),
		   IBD_pval = numeric(0),
		   IBD_stat = numeric(0)) %>%
		fwrite(paste0(TABLES_DIR, 'ibd_raw.tsv'), sep = '\t')
	data.table(pop1 = character(0), pop2 = character(0)) %>%
		fwrite(paste0(TABLES_DIR, 'ibd_pairs.tsv'), sep = '\t')
	quit(status = 0)
}

# Calculate IBD
ibd <- FUN_ibd(clusters, samples, CPU)

# Save
ibd$raw %>%
	fwrite(paste0(TABLES_DIR, 'ibd_raw.tsv'), sep = '\t')
ibd$not_ibd %>%
	do.call(rbind, . ) %>%
	as.data.table %>%
	dplyr::rename(pop1 = V1, pop2 = V2) %>%
	fwrite(paste0(TABLES_DIR, 'ibd_pairs.tsv'), sep = '\t')

#TODO plot heatmap
