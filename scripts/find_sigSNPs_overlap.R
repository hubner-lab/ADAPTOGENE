library(data.table)
library(dplyr)
library(qvalue)
library(stringr)
library(magrittr)
library(GenomicRanges)
library(parallel)
library(purrr)
args = commandArgs(trailingOnly=TRUE)
options(scipen = 99999)
################
sigSNPs_vec = args[1] %>% str_split(',') %>% unlist # comma separated list of files
METHOD = args[2]
GAP = args[3] %>% as.numeric
PREDICTORS_SELECTED = args[4] %>% str_split(',') %>% unlist
################
message(sigSNPs_vec, METHOD, GAP, PREDICTORS_SELECTED)

################ Functions
FUN_make_gr <- function(dt){
  dt %>%
    dplyr::mutate(start = pos,
                  end = pos) %>%
    dplyr::select(chr, start, end, SNPID, trait, method) %>%
    GRanges()
}

######################################### Main
# Name by method
methods_vec = sigSNPs_vec %>% str_split('/') %>% sapply(function(x) x[2])
sigSNPs_vec %<>% setNames(methods_vec)

# Load tables with filter by only clim traits
sigSNPs_lst <- 
  lapply(sigSNPs_vec, function(x) 
    x %>% 
      fread %>% 
      dplyr::filter(trait %in% !!PREDICTORS_SELECTED) %>%
      dplyr::select(SNPID, chr, pos, trait, method, pvalue)
  ) %>%
  setNames(methods_vec)
                                                    # message(sigSNPs_lst %>% str)
if(METHOD == 'EMMAX'){
  snpIDs <-
	  sigSNPs_lst[['EMMAX']]
}

if(METHOD == 'LFMM'){
  snpIDs <-
	  sigSNPs_lst[['LFMM']]
}

if(METHOD == 'Sum'){
  snpIDs <-
    	rbind(sigSNPs_lst[['EMMAX']],
	      sigSNPs_lst[['LFMM']]
	)
}

if(METHOD == 'Overlap'){
  # convert to GR
  sigSNPs_lstGR <- 
    lapply(sigSNPs_lst, FUN_make_gr) %>% 
    setNames(methods_vec)
                                            # message(sigSNPs_lstGR[[1]] %>% as.data.table %>% str)
  snpIDs <-
    lapply(methods_vec, function(n1){
      lapply(methods_vec, function(n2){
        
        # Do not compare the same ranges
        if(n1 == n2){return(NULL)}
        
        overlaps_dt <-
          findOverlaps(sigSNPs_lstGR[[n1]],
                       sigSNPs_lstGR[[n2]], 
                       maxgap = GAP) %>% 
          as.data.table
        
        dt1 = sigSNPs_lst[[n1]]
	dt2 = sigSNPs_lst[[n2]]
	
	dt1_overlap = dt1[overlaps_dt$queryHits,]
	dt2_overlap = dt2[overlaps_dt$subjectHits,]	
       	
       	rbind(dt1_overlap,
	      dt2_overlap
      	)

      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) %>%
  unique
}

if(METHOD == 'PairOverlap'){
  # convert to GR
  sigSNPs_lstGR <- 
    lapply(sigSNPs_lst, FUN_make_gr) %>% 
    setNames(methods_vec)
  # message(sigSNPs_lstGR[[1]] %>% as.data.table %>% str)
  snpIDs <-
    lapply(methods_vec, function(n1){
      lapply(methods_vec, function(n2){
        lapply(PREDICTORS_SELECTED, function(bio){
          # Do not compare the same ranges
          if(n1 == n2){return(NULL)}
          
          gr1_bio = sigSNPs_lstGR[[n1]][sigSNPs_lstGR[[n1]]$trait == bio]
          gr2_bio = sigSNPs_lstGR[[n2]][sigSNPs_lstGR[[n2]]$trait == bio]
          
          if(length(gr1_bio) == 0 | length(gr2_bio) == 0){return(NULL)}
          
          overlaps_dt <-
            findOverlaps(gr1_bio,
                         gr2_bio, 
                         maxgap = GAP) %>% 
            as.data.table
          
    dt1 = sigSNPs_lst[[n1]] %>% dplyr::filter(trait == !!bio)
    dt2 = sigSNPs_lst[[n2]] %>% dplyr::filter(trait == !!bio)

    dt1_overlap = dt1[overlaps_dt$queryHits,]
    dt2_overlap = dt2[overlaps_dt$subjectHits,]	

    rbind(dt1_overlap,
	  dt2_overlap
    )
          
        }) %>% do.call(rbind, .) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) %>%
    unique
}


                                          # message(snpIDs)

# Save selected SNPIDs
snpIDs %>%
  as.data.table %>%
  dplyr::arrange(chr, pos) %>%
  fwrite(paste0('tables/Selected_SNPs_redundant_',
                METHOD,
                '_Gap', GAP,
                '.tsv'),
	sep = '\t' 
  	)
# Save uniq
snpIDs %>%
	as.data.table %>% 
	dplyr::arrange(chr, pos) %>%
	dplyr::filter(!duplicated(SNPID)) %>%
	dplyr::select(SNPID) %>%
  	fwrite(paste0('tables/Selected_SNPs_unique_',
                METHOD,
                '_Gap', GAP,
                '.list'),
	sep = '\t' 
  	)

