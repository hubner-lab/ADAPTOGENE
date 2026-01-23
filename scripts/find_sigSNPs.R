library(data.table)
library(dplyr)
library(qvalue)
library(stringr)
library(GenomicRanges)
library(parallel)
args = commandArgs(trailingOnly=TRUE)
################
ASSOC_TABLE = args[1] # pvalues	
ADJUST_LIST = args[2] %>% str_split(',') %>% unlist
#ADJUST = args[2] %>% str_split('_') %>% unlist %>% .[1] # bonf_PVAL or qval_PVAL or custom_PVAL
#PVAL_THRESHOLD = args[2] %>% str_split('_') %>% unlist %>% .[2] %>% as.numeric
DISTANCE = args[3] %>% as.numeric # Distance based on LD decay to count neighbor SNPs
METHOD = args[4]
CPU = args[5] %>% as.numeric
################

######################## Function
# Define the function to calculate FDR and return the max p-value with FDR < 0.05
FUN_max_pvalue_fdr <- function(pvalues, pval_threshold) {
  
  
  # Calculate the q-values (FDR adjusted p-values)
  qvalues_result <- qvalue(pvalues)
  
  # Filter the p-values where the q-value (FDR) is less than 0.05
  significant_pvalues <- pvalues[qvalues_result$qvalues < pval_threshold]
  
  # If there are any significant p-values, return the maximum one
  if (length(significant_pvalues) > 0) {
    return(max(significant_pvalues))
  } else {
    # If no significant p-values, return NA or a custom message
    return(NA)
  }
}

FUN_max_pvalue_top <- function(pvalues, topN){
         if (length(pvalues) < topN) {
                 stop("topN is larger than the number of available p-values")
                                }
         sorted_pvalues <- sort(pvalues, decreasing = FALSE) # Sort in ascending order
         return(max(sorted_pvalues[1:topN])) # Return the maximum in the topN lowest values
}

FUN_find_sigSNPs <- function(ASSOC_PVAL,  # pvalues only
			     adjustment = 'bonf',  # bonf or qval or topN
                             pval_threshold = 0.05, # for fdr\qvalue or for bonferonni base value
			     CPU = 1
                                      ){

  snps_assoc = fread(ASSOC_PVAL)
  										message(snps_assoc %>% str)
  traits = snps_assoc %>% dplyr::select(-SNPID, -chr, -pos) %>% colnames
										message(adjustment %>% str)  
  if(adjustment == 'bonf'){ pval_threshold = pval_threshold / nrow(snps_assoc) }
  if(adjustment == 'qval'){ pval_thresholds = sapply(snps_assoc %>% dplyr::select(-SNPID, -chr, -pos),
						    function(x) FUN_max_pvalue_fdr(x, pval_threshold) ) %>% setNames(traits) }
  if(adjustment == 'top'){ pval_thresholds =  sapply(snps_assoc %>% dplyr::select(-SNPID, -chr, -pos),
                                                    function(x) FUN_max_pvalue_top(x, pval_threshold) ) %>% setNames(traits)}

 # if custom don't change anything
  									
						
  snps_sig_dt <-
    lapply(traits, function(trait){
	
      if(adjustment %in% c('top', 'qval')) { pval_threshold = pval_thresholds[trait] } # choose pval_threshold for multi_thresholds options (qval, top)
      
      mask = snps_assoc[[trait]] < pval_threshold

      snps_assoc[mask, ] %>%
	      dplyr::select(SNPID, chr, pos, !!trait) %>%
	      setNames(c('SNPID', 'chr', 'pos', 'pvalue')) %>%
	      dplyr::mutate(pval_threshold = !!pval_threshold,
	      		    trait = !!trait)
    }) %>%
	      do.call(rbind, .)
										message(snps_sig_dt %>% str)
  # Add overlapping info
  overlaps_dt <-
	  mclapply(1:nrow(snps_sig_dt), function(i){
			 trait_name = snps_sig_dt[i,]$trait
			 current_snp_gr <-
				 snps_sig_dt[i,] %>%
				 dplyr::mutate(start = pos,
					       end = pos) %>%
				 dplyr::select(chr, start, end) %>%
				 GRanges()

			 other_traits_dt <-
				 snps_sig_dt %>%
				 	dplyr::filter(trait != !!trait_name)

			 other_traits_gr <-
				 other_traits_dt %>%
					dplyr::mutate(start = pos,
	                                               end = pos) %>%
	                                 dplyr::select(chr, start, end) %>%
	                                 GRanges()
			 overlaps  <-
				 findOverlaps(current_snp_gr, other_traits_gr, 
					      maxgap = DISTANCE, ignore.strand = T) %>%
				 as.data.frame

			 traits <-
				 other_traits_dt[overlaps$subjectHits,]$trait %>%
				 paste(collapse = ',')
			 snpids <-
				 paste0(other_traits_dt[overlaps$subjectHits,]$chr, ':',
					other_traits_dt[overlaps$subjectHits,]$pos) %>%
                                 paste(collapse = ',')

			 return(data.frame(overlap_traits = ifelse(traits == "", NA, traits),
					   overlap_snps = ifelse(snpids == ':', NA, snpids),
					   overlap_distance = DISTANCE))


    }, mc.cores = CPU) %>% do.call(rbind, .)

								message(snps_sig_dt %>% str)
  								message(overlaps_dt %>% str)
    snps_sig_overlap_dt <-
	    cbind(snps_sig_dt,
		  overlaps_dt)

      return(snps_sig_overlap_dt)
}

######################################### Main

lapply(ADJUST_LIST, function(adjust_pval){
	# Parse param
	ADJUST = adjust_pval %>% str_split('_') %>% unlist %>% .[1]
	PVAL_THRESHOLD = adjust_pval %>% str_split('_') %>% unlist %>% .[2]


	snps_sig <-
		FUN_find_sigSNPs(ASSOC_TABLE,
				ADJUST,
				PVAL_THRESHOLD %>% as.numeric,
				CPU)


	snps_sig %>%
		dplyr::mutate(method = METHOD) %>%
		dplyr::select(SNPID, chr, pos, pvalue, pval_threshold, method, 
			      trait, overlap_traits, overlap_snps, overlap_distance) %>%
		fwrite(paste0(ASSOC_TABLE %>% gsub('.tsv', '', .), # extract prefix
			      '_sigSNPs_',
			      ADJUST, '_', 
			      PVAL_THRESHOLD,
			      '.tsv'),
		       sep = '\t'
			)
		  })
