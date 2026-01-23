library(data.table)
library(dplyr)
library(ggplot2) 
library(CMplot)
library(stringr)
library(parallel)
args = commandArgs(trailingOnly=TRUE)

###################################
ASSOC_TABLE = args[1] # should be pval sorted with columns SNP, CHR, POS, pvalue
sigSNPs = args[2] # any table with SNP column with chr:POS IDs
TRAIT_PAIRs = args[3] %>% str_split(';') %>% unlist %>% str_split(',') # trait1,trait2;trait5,trait6
################################### MAIN

# LOAD GEA/GWAS results
#SNPID   chr     pos
assoc = fread(ASSOC_TABLE)

# Load sig SNPs table
# chr     pos     pvalue  pval_threshold  trait   overlap_traits  overlap_snps    overlap_distance
sigsnps = fread(sigSNPs) %>%
	dplyr::mutate(SNP = paste0(chr, ':', pos) ) %>%
	dplyr::select(SNP, everything())


lapply(TRAIT_PAIRs, function(trait_pair){
	
       	SNPs <- list(
		sigsnps$SNP[sigsnps$trait == trait_pair[1]],
		sigsnps$SNP[sigsnps$trait == trait_pair[2]]
	)

	thresholds = list(
		sigsnps$pval_threshold[sigsnps$trait == trait_pair[1]] %>% unique,
        	sigsnps$pval_threshold[sigsnps$trait == trait_pair[2]] %>% unique
        )

	assoc_subset <-
	       	assoc %>% 
		dplyr::select(SNPID, chr, pos, !!trait_pair)
										message(assoc_subset %>% str)

	# plot miami from example
	CMplot(assoc_subset, 
	       type="p",
	       plot.type="m", 
	       band=0.5, 
	       LOG10=TRUE, 
	       ylab="-log10(pvalue)",

	       threshold = thresholds,
	       threshold.lty=2, 
	       threshold.lwd=1, 
	       threshold.col="red", 
	       

	       amplify=TRUE,
	       cex=0.6,
	       chr.den.col=NULL, 
	       
	       file="pdf",
	       
	       dpi=300,
	       file.output=TRUE,
	       verbose=TRUE, 
	       multracks = TRUE,
       	       
	       chr.labels.angle=60,
	       highlight.type="l",
	       highlight=SNPs, 
#	       highlight.text=SNPs, 
	       highlight.col=NULL)
})
