#!/usr/bin/env Rscript
# LFMM association analysis (refactored)
# Trains LFMM model on LD-pruned data and tests on full dataset

library(LEA)
library(dplyr)
library(data.table)
library(WGCNA)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
#################################
LFMM_LD_IMP = args[1]   # LD-pruned imputed LFMM (for training)
LFMM_IMP = args[2]       # Full imputed LFMM (for testing)
CLIMATE = args[3]        # Climate data (scaled)
Kbest = args[4] %>% as.numeric
PREDICTORS_SET = args[5] %>% str_split(',') %>% unlist
VCFSNP = args[6]         # SNP positions
TABLES_DIR = args[7]     # Output directory
#################################

message('INFO: Starting LFMM analysis')
message(paste0('INFO: K = ', Kbest))
message(paste0('INFO: Predictors: ', paste(PREDICTORS_SET, collapse = ', ')))

# Load LD-pruned imputed LFMM (for training)
message('INFO: Read LFMM_LD (training data)')
lfmm_ld_imp <- fread(LFMM_LD_IMP, sep = ' ', header = F)
message(lfmm_ld_imp %>% str)

# Load predictors
message('INFO: Read predictors')
predictors <- fread(CLIMATE, sep = '\t', header = T) %>%
    dplyr::select(!!PREDICTORS_SET)
message(predictors %>% str)

# Scale predictors
message('INFO: Scale predictors')
predictors[, names(predictors) := lapply(.SD, scale)]
message(predictors %>% str)

# Load full imputed LFMM (for testing)
message('INFO: Read LFMM full (testing data)')
lfmm_imp <- fread(LFMM_IMP, sep = ' ', header = F)

# Load SNP positions
message('INFO: Read VCFSNP')
vcfsnp <- fread(VCFSNP, sep = ' ', header = F) %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chr = as.character(chr),
                  SNPID = paste0(chr, ':', pos)) %>%
    dplyr::select(SNPID, chr, pos)

# Run LFMM2 separately for each predictor
message('INFO: Run LFMM2 separately for each variable')
pval_list <- lapply(names(predictors), function(bio) {
    message(paste0('INFO: Train LFMM model on LD dataset with ', bio, ' predictor'))

    # Build model on LD data
    lfmm.model <- lfmm2(lfmm_ld_imp,
                        env = predictors[, ..bio],
                        K = Kbest)
    message(lfmm.model %>% str)

    # Test on full data
    message('INFO: Extract p-values from the model')
    lfmm.res <- lfmm2.test(lfmm.model,
                           input = lfmm_imp,
                           env = predictors[, ..bio])
    message(lfmm.res %>% str)

    message(paste0('INFO: Finished with ', bio, ' predictor'))
    return(lfmm.res$pvalues)
})

# Combine results
pval_dt <- pval_list %>%
    do.call(cbind, .) %>%
    as.data.table %>%
    setNames(names(predictors)) %>%
    cbind(vcfsnp, .)

# Save results
message('INFO: Save p-values')
pval_dt %>%
    fwrite(paste0(TABLES_DIR, 'LFMM_pvalues_K', Kbest, '.tsv'),
           sep = '\t', col.names = T, row.names = F, quote = F)

message('INFO: LFMM analysis complete')
