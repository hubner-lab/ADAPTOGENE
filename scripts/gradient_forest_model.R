library(dplyr)
library(data.table)
library(vegan)
library(gradientForest)
library(stringr)
library(qs)
args = commandArgs(trailingOnly=TRUE)
##############
LFMM = args[1]
SIGSNPS = args[2]       # Selected_SNPs.tsv with SNPID column
VCFSNP = args[3]        # .vcfsnp file
REMOVED = args[4]       # .removed file
SAMPLES = args[5]
CLIM_PRESENT_SITE = args[6]
PREDICTORS_SELECTED = args[7] %>% str_split(',') %>% unlist
NTREE = args[8] %>% as.numeric
COR_THRESHOLD = args[9] %>% as.numeric
PCNM = args[10]         # "with" or "without"
MODEL_TYPE = args[11]   # "adaptive" or "random"
OUTPUT = args[12]
##############

set.seed(42)
message(paste0('INFO: Building ', MODEL_TYPE, ' Gradient Forest model'))

# Load LFMM genotype matrix
lfmm_dt <- fread(LFMM)
message(paste0('INFO: LFMM matrix: ', nrow(lfmm_dt), ' samples x ', ncol(lfmm_dt), ' SNPs'))

# Load SNP IDs and removed SNPs
sigsnps <- fread(SIGSNPS, header = TRUE)$SNPID
vcfsnp <- fread(VCFSNP, header = FALSE) %>%
  dplyr::select(V1, V2) %>%
  setNames(c('chr', 'pos')) %>%
  dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>%
  .$chrpos

removed_dt <- fread(REMOVED, header = FALSE)
if (!is.null(removed_dt) && nrow(removed_dt) > 0) {
  removed <- removed_dt %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>%
    .$chrpos
  mask_adaptive <- (!vcfsnp %in% removed) & (vcfsnp %in% sigsnps)
} else {
  mask_adaptive <- vcfsnp %in% sigsnps
}

n_adaptive <- sum(mask_adaptive)
message(paste0('INFO: Adaptive SNPs: ', n_adaptive))

# Select SNPs based on model type
if (MODEL_TYPE == 'adaptive') {
  snp_subset <- lfmm_dt[, ..mask_adaptive]
  message(paste0('INFO: Using ', ncol(snp_subset), ' adaptive SNPs'))
} else if (MODEL_TYPE == 'random') {
  n_random <- min(max(n_adaptive, 300), ncol(lfmm_dt))
  random_cols <- sample(names(lfmm_dt), n_random)
  snp_subset <- lfmm_dt[, ..random_cols]
  message(paste0('INFO: Using ', ncol(snp_subset), ' random SNPs'))
} else {
  stop(paste0('ERROR: Unknown MODEL_TYPE: ', MODEL_TYPE))
}

# Compute maxLevel for gradient forest
maxLevel <- log2(0.368 * nrow(lfmm_dt) / 2)

# Load samples for PCNM
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# Load climate predictors
env.bio <- fread(CLIM_PRESENT_SITE) %>% dplyr::select(!!PREDICTORS_SELECTED)

# Build input matrix and predictor list
if (PCNM == 'with') {
  coords <- data.frame(long = samples$longitude, lat = samples$latitude)
  pcnm_result <- pcnm(dist(coords))
  keep <- round(length(which(pcnm_result$value > 0)) / 2)
  pcnm_scores <- vegan::scores(pcnm_result)[, 1:keep]
  message(paste0('INFO: Using ', ncol(pcnm_scores), ' PCNM axes'))

  input_matrix <- cbind(env.bio, pcnm_scores, snp_subset)
  predictor_vars <- c(colnames(env.bio), colnames(pcnm_scores))
} else {
  input_matrix <- cbind(env.bio, snp_subset)
  predictor_vars <- colnames(env.bio)
}

message(paste0('INFO: Input matrix: ', nrow(input_matrix), ' samples x ', ncol(input_matrix), ' variables'))
message(paste0('INFO: Predictor variables: ', paste(predictor_vars, collapse = ', ')))
message(paste0('INFO: Response variables (SNPs): ', ncol(snp_subset)))

# Run Gradient Forest
gf <- gradientForest(input_matrix,
                     predictor.vars = predictor_vars,
                     response.vars = colnames(snp_subset),
                     ntree = NTREE,
                     maxLevel = maxLevel,
                     trace = TRUE,
                     corr.threshold = COR_THRESHOLD)

# Save model
qsave(gf, OUTPUT)
message(paste0('INFO: Saved ', MODEL_TYPE, ' GF model to: ', OUTPUT))
