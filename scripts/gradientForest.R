library(dplyr)
library(data.table)
library(vegan)
library(gradientForest)
library(stringr)
library(qs)
args = commandArgs(trailingOnly=TRUE)
##############
LFMM = args[1]
SNPIDs = args[2] # list of SNPIDs
VCFSNP = args[3] # list of snps from filtred vcf
REMOVED = args[4] # removed snps after transform in lfmm
SAMPLES = args[5]
CLIM_PRESENT_SITE = args[6]
PREDICTORS_SELECTED = args[7] %>% str_split(',') %>% unlist
NTREE = args[8] %>% as.numeric # 1000
COR_THRESHOLD = args[9] %>% as.numeric # 0.5
PCNM = args[10] # with or without
GF_SUFFIX = args[11] # should be specified for not overwriting different runs
# For labelling only
#K_BEST = args[11]
#ADJUST = args[12]
##############
SUFFIX = paste0(GF_SUFFIX, '_', PCNM, 'PCNM')
##############

lfmm_dt <- fread(LFMM)
									message(lfmm_dt %>% ncol )
									message(SNPIDs)
sigsnps <- fread(SNPIDs, header = T)$SNPID
									message(sigsnps %>% str)

vcfsnp <- fread(VCFSNP, header = F) %>% dplyr::select(V1, V2) %>% setNames(c('chr', 'pos')) %>% dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>% .$chrpos

removed_dt = fread(REMOVED, header = F) 
if(is.null(removed_dt) != 0){ 
	removed <-
		removed_dt %>% dplyr::select(V1, V2) %>% setNames(c('chr', 'pos')) %>% dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>% .$chrpos

	mask_adaptive = (!vcfsnp %in% removed) & (vcfsnp %in% sigsnps)
} else {
	mask_adaptive = vcfsnp %in% sigsnps	
}

mask_random = sample(names(lfmm_dt), 
		     max(mask_adaptive %>% sum, 300)
		    ) #TODO can add additional parameter to increase random sampled SNPs for random model


lfmm_imp_adaptive <- lfmm_dt[,..mask_adaptive]
lfmm_imp_random <- lfmm_dt[,..mask_random]

						message(lfmm_imp_adaptive %>% str)
						message(lfmm_imp_random %>% str)
maxLevel <- log2(0.368 * nrow(lfmm_dt) / 2)

samples <- fread(SAMPLES,
                 colClasses = c("site" = "character", 
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

coords <- data.frame(long = samples$longitude,
                     lat = samples$latitude)

pcnm <- pcnm(dist(coords))  #this generates the PCNMs, you could stop here if you want all of them
keep <- round(length(which(pcnm$value > 0))/2)
pcnm <- vegan::scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors

							message(pcnm %>% str)
env.bio = fread(CLIM_PRESENT_SITE) %>% dplyr::select(!!PREDICTORS_SELECTED)

# Run GF adaptive

## Adaptive
if(PCNM == 'with'){
	input_matrix = cbind(env.bio,
	                     pcnm,
	                     lfmm_imp_adaptive)
								message(input_matrix %>% str)	
								message(lfmm_imp_adaptive %>% str)
								message('FUCK')

	gf <- gradientForest(input_matrix,
	                     predictor.vars=c(colnames(env.bio)
	                                      ,colnames(pcnm)
	                                      ),
	                     response.vars=colnames(lfmm_imp_adaptive),
	                     ntree=NTREE,
	                     maxLevel=maxLevel,
		                     trace=T,
                	     corr.threshold=COR_THRESHOLD)
	qsave(gf, paste0('intermediate/gradientForest_adaptive_', SUFFIX, '.qs'))

	#TODO mclapply to run in parallel?
	## Random
	input_matrix = cbind(env.bio,
	                     pcnm,
	                     lfmm_imp_random)
	gf_random <- gradientForest(input_matrix,
	                     predictor.vars=c(colnames(env.bio)
	                                      ,colnames(pcnm)
	                                      ),
	                     response.vars=colnames(lfmm_imp_random),
	                     ntree=NTREE,
	                     maxLevel=maxLevel,
	                     trace=T,
	                     corr.threshold=COR_THRESHOLD)
	
	qsave(gf_random, paste0('intermediate/gradientForest_random_', SUFFIX, '.qs'))

} else {
        input_matrix = cbind(env.bio,
                             lfmm_imp_adaptive)
								message(input_matrix %>% str)	
								message(lfmm_imp_adaptive %>% str)
								message('FUCK')
	gf <- gradientForest(input_matrix,
                             predictor.vars=colnames(env.bio),
                             response.vars=colnames(lfmm_imp_adaptive),
                             ntree=NTREE,
                             maxLevel=maxLevel,
                                     trace=T,
                             corr.threshold=COR_THRESHOLD)
        qsave(gf, paste0('intermediate/gradientForest_adaptive_', SUFFIX, '.qs'))

        #TODO mclapply to run in parallel?
        ## Random
        input_matrix = cbind(env.bio,
                             pcnm,
                             lfmm_imp_random)

        gf_random <- gradientForest(input_matrix,
                             predictor.vars=colnames(env.bio),
                             response.vars=colnames(lfmm_imp_random),
                             ntree=NTREE,
                             maxLevel=maxLevel,
                             trace=T,
                             corr.threshold=COR_THRESHOLD)

        qsave(gf_random, paste0('intermediate/gradientForest_random_', SUFFIX, '.qs'))

}
