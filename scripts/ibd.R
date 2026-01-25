library(dplyr)
library(data.table)
library(parallel)
library(vegan)
library(stringr)
library(geosphere)

args = commandArgs(trailingOnly=TRUE)
set.seed(42)
#######################
CLUSTERS = args[1]
SAMPLES = args[2]
CPU = args[3] %>% as.numeric
TABLES_DIR = args[4]   # output directory for tables
#######################

############################## Function
FUN_ibd <- function(clusters, samples, cpu) {

  message(paste0('INFO: RUN IBD analysis with: n = ',
                 samples$site %>% unique %>% length,
                 '; GeoMatrix = ', nrow(samples)))

  ibd <- mclapply(samples$site %>% unique, function(i) {
    mclapply(samples$site %>% unique, function(i2) {

      # Select populations
      clust_df = clusters[clusters$site %in% c(i, i2), ] %>% dplyr::select(starts_with('C'))
      geo_df = samples[samples$site %in% c(i, i2), ] %>% dplyr::select(latitude, longitude)

      clust_dist = vegdist(clust_df, method = 'jaccard')
      geo_dist = distm(geo_df, fun = distVincentyEllipsoid)

      # Perform Mantel test for checking IBD
      ibd = mantel(clust_dist, geo_dist,
                   method = 'pearson',
                   permutations = 999)

      res_dt <- data.table(Comparison = paste0(i, ' vs ', i2),
                           IBD_pval = ibd$signif,
                           IBD_stat = ibd$statistic)

      return(res_dt)
    }, mc.cores = cpu, mc.allow.recursive = TRUE) %>% do.call(rbind, .)
  }, mc.cores = cpu, mc.allow.recursive = TRUE) %>% do.call(rbind, .) %>%
    na.omit

  message("INFO: IBD calculation completed, process the result")
  message(ibd %>% str)

  # Look on not isolated:
  ibd.list <- ibd[ibd$IBD_pval > 0.05, ] %>%
    .$Comparison %>%
    str_split(' vs ') %>%
    lapply(sort) %>%
    unique

  return(list(raw = ibd,
              not_ibd = ibd.list))
}

############################ MAIN
# Load tables
# Filter out all populations which contain only 1 individual
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
  fwrite(paste0(TABLES_DIR, 'IBD_raw.tsv'), sep = '\t')

ibd$not_ibd %>%
  do.call(rbind, .) %>%
  as.data.table %>%
  dplyr::rename(pop1 = V1, pop2 = V2) %>%
  fwrite(paste0(TABLES_DIR, 'IBD_notIsolated.tsv'), sep = '\t')
