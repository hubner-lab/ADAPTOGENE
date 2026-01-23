library(dplyr)
library(parallel)
library(data.table)
library(CMplot)
library(stringr)
library(qvalue)
library(qs)
args = commandArgs(trailingOnly=TRUE)
################################
ASSOC_TABLE = args[1] # SNPID chr pos TRAITS...
ADJUST_LIST = args[2] %>% str_split(',') %>% unlist
Kbest = args[3] %>% as.numeric
METHOD = args[4]

									message(METHOD)
									message(Kbest)
									message(args[2])
									message(ADJUST_LIST %>% str)

#CUSTOM_SNP = fread(args[4], header = F)$V1 # should be list without header, NA if highlight be pval threshold
CUSTOM_SNP = NA
################################ FUNCTIONS
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
    return(NULL)
  }
}

FUN_max_pvalue_top <- function(pvalues, topN){
	 if (length(pvalues) < topN) {
   		 stop("topN is larger than the number of available p-values")
  				}
	 sorted_pvalues <- sort(pvalues, decreasing = FALSE) # Sort in ascending order
	 return(max(sorted_pvalues[1:topN])) # Return the maximum in the topN lowest values
}
############################### MAIN


# Prepare DFs

message('INFO: PREPARING GWAS1 DATA FRAME')
snps_assoc = fread(ASSOC_TABLE, sep = '\t', header = T) %>%
	dplyr::rename(SNP = SNPID, Chromosome = chr, Position = pos) %>% # rename for CMplot
	as.data.frame
# Define pval thresholds
                                                                                message(snps_assoc %>% str)
traits = snps_assoc %>% dplyr::select(-SNP, -Chromosome, -Position) %>% colnames


# Change wd for saving plots 
setwd(paste0('plots/', METHOD))

# Run each adjustament variant
lapply(ADJUST_LIST, function(ADJUST){ #TODO test mclapply here

	# Parse params
	adjustment = ADJUST %>% str_split('_') %>% unlist %>% .[1] # bonf_PVAL or qval_PVAL or top_N
	pval_threshold = ADJUST %>% str_split('_') %>% unlist %>% .[2] %>% as.numeric
													message(ADJUST)
													message(adjustment)
													message(pval_threshold)
	# Define thresholds       
	if(adjustment == 'bonf'){ pval_threshold = pval_threshold / nrow(snps_assoc) }
	if(adjustment == 'qval'){ pval_threshold = lapply(snps_assoc %>% dplyr::select(-SNP, -Chromosome, -Position),
                                                    function(x) FUN_max_pvalue_fdr(x, pval_threshold) ) %>% setNames(traits) }
	if(adjustment == 'top'){ pval_threshold =  lapply(snps_assoc %>% dplyr::select(-SNP, -Chromosome, -Position),
                                                    function(x) FUN_max_pvalue_top(x, pval_threshold) ) %>% setNames(traits)}
	# Add adjustment and pval to trait names
	colnames(snps_assoc)[4:ncol(snps_assoc)] = paste0(traits, '_K', Kbest, '_',  ADJUST)
										message(pval_threshold)
										message(snps_assoc %>% colnames)
	# Plot
	message('INFO: PLOTTING')

	#TODO how to save as object, probably remake with ggplot object!!!!

	CMplot(snps_assoc, # return some error to data.table 
       		plot.type="m", 
	       	LOG10=TRUE, 
		ylim=NULL,
	       	threshold= pval_threshold,
		threshold.lty=2,
	        threshold.lwd=1, 
		threshold.col="black", 
		amplify=T,
		bin.size=1e6,
		chr.den.col=c("darkgreen", "yellow", "red"),
		signal.col = 'red',
		signal.cex = 0.5,
		cex = 0.25,
       	
		#axis.cex = 0.85, # x/y tick labels

		file="pdf",
		dpi=600,
		multracks = F,
		file.output=TRUE,
		verbose=TRUE,
		chr.labels.angle=60
		,width=8,height=4
	
	)
})
