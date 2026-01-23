library(LEA)
#library(parallel)
library(dplyr)
library(data.table)
library(WGCNA)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

#################################
LFMM_LD_imp = args[1]
LFMM_imp = args[2]
CLIMATE = args[3]
Kbest = args[4] %>% as.numeric # latent factor number
PREDICTORS_set = args[5] %>% str_split(',') %>% unlist # names of predictors to use separated by comma ex: bio_1,bio_2,bio_3 ... 
VCFSNP = args[6]

SEPARETELY=T # For future use, default performance run lfmm2 separetely on each predictor
OUTPREFIX = "tables/LFMM/LFMM_"
#CPU = args[6] # no needed
#################################

							message('INFO: READ LFMM_LD')
lfmm_ld_imp <- 
    fread(LFMM_LD_imp, sep = ' ', header = F)
							message(lfmm_ld_imp %>% str)

							message('INFO: READ PREDICTORS')

predictors <-
	      fread(CLIMATE, sep = '\t', header = T)	       %>%
		dplyr::select(!!PREDICTORS_set) # keep only predefined predictors
							message(predictors %>% str)
							message('INFO: SCALE PREDICTORS')
predictors[, names(predictors) := lapply(.SD, scale)] # inplace scaling

							message(predictors %>% str)
							message('INFO: READ LFMM')

lfmm_imp <-
        fread(LFMM_imp, sep = ' ', header = F)
							
							message('INFO: READ VCFSNP')
vcfsnp <-
	fread(VCFSNP, sep = ' ', header = F) %>%
		dplyr::select(V1, V2) %>%
		setNames(c('chr', 'pos')) %>%
		dplyr::mutate(SNPID = paste0(chr, ':', pos)) %>%
		dplyr::select(SNPID, chr, pos)

							 message('INFO: Run LFMM2 separately for each variable')
lapply(names(predictors), function(bio){
           					      message(paste0('INFO: Train LFMM model on LD dataset with ', bio, ' predictor'))
	# Build model
	lfmm.model <-
	    lfmm2(lfmm_ld_imp,
	          env = predictors[,..bio],
	          K = Kbest)
							message(lfmm.model %>% str)

# Save model
#saveRDS(lfmm.model, paste0(OUTPREFIX, '.Rds'))

							message('INFO: EXTRACT PVALUES FROM THE MODEL')

	lfmm.res <-
	      lfmm2.test(lfmm.model,
	                 input = lfmm_imp,
	                 env = predictors[,..bio])
							message(lfmm.res %>% str)

							message('INFO: SAVE PVALUES')
#	if(ncol(lfmm.res$pvalues) != ncol(predictors)){ # For more than one trait it needed to be transposed
#		lfmm.res$pvalues %>%
#		transposeBigData %>%
#		fwrite(paste0(OUTPREFIX, '_pvalues.tsv'), sep = '\t', col.names = T, row.names = F, quote = F)
#			 } else{
	return(lfmm.res$pvalues)
							message(paste0('INFO: FINISH with ', bio, ' predictor'))
}) %>% 
       do.call(cbind, . ) %>%
       as.data.table %>% 
      setNames(names(predictors)) %>%
      cbind(vcfsnp, .) %>% # add SNPID chr pos
      fwrite(paste0(OUTPREFIX, 'pvalues_K', Kbest, '.tsv'),
                                sep = '\t', col.names = T, row.names = F, quote = F)



