#!/usr/bin/env Rscript
# EMMAX association analysis (refactored)
# Runs EMMAX for all traits specified in PREDICTORS_SELECTED
# Uses pre-computed TPED/TFAM and kinship from separate Snakemake rules.

library(dplyr)
library(data.table)
library(magrittr)
library(purrr)
library(stringr)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
########################
VCF = args[1]
Kbest = args[2] %>% as.numeric
TRAIT_FILE = args[3]  # tsv file with climate/trait values
COVARIATES_FILE = args[4]  # eigenvec file from PCA
PREDICTORS_SELECTED = args[5] %>% str_split(',') %>% unlist
INTER_DIR = args[6]
SAMPLES_FILE = args[7]
TABLES_DIR = args[8]
TPED_PREFIX = args[9]     # Pre-computed TPED/TFAM prefix
KINSHIP_FILE = args[10]   # Pre-computed .aBN.kinf
########################

message("INFO: Starting EMMAX analysis")
message(paste0("INFO: K = ", Kbest))
message(paste0("INFO: Predictors: ", paste(PREDICTORS_SELECTED, collapse = ", ")))
message(paste0("INFO: TPED prefix: ", TPED_PREFIX))
message(paste0("INFO: Kinship file: ", KINSHIP_FILE))

# Load trait data
trait <- fread(TRAIT_FILE, sep = '\t', header = T)

# Load covariates (PCs from LEA projections or PLINK eigenvec)
# LEA format: space-separated, no header, cols = PCs (no FID/IID)
# PLINK format: space-separated, no header, cols = FID IID PC1 PC2 ...
cov_raw <- fread(COVARIATES_FILE, sep = ' ', header = F)
# Detect format: if first column is all same value (FID) and second is different (IID), it's PLINK format
if (ncol(cov_raw) > 2 && length(unique(cov_raw[[1]])) == 1 && length(unique(cov_raw[[2]])) == nrow(cov_raw)) {
    message("INFO: Detected PLINK eigenvec format")
    covariates <- cov_raw[, 3:(2 + Kbest)] %>% setNames(paste0('PC', 1:Kbest))
} else {
    message("INFO: Detected LEA projections format")
    covariates <- cov_raw[, 1:Kbest] %>% setNames(paste0('PC', 1:Kbest))
}

################################ Functions

FUN_emmax <- function(VCF, trait, covariates, OUT, TPED_PREFIX, KINSHIP_FILE) {

    tfam <- paste0(TPED_PREFIX, '.tfam')

    message("INFO: Read VCF and extract SNPID")
    snpid <- fread(cmd = paste('grep -v "##"', VCF, '| cut -f1-2')) %>%
        setNames(c('CHROM', 'POS')) %>%
        dplyr::mutate(SNPID = paste0(CHROM, ':', POS)) %>%
        .$SNPID

    message("INFO: Read VCF and extract SampleNames")
    SampleName <- fread(cmd = paste("head -1000", VCF, "| grep '#CHROM' | cut -f10- "), header = F) %>%
        as.character()

    # FILES produced
    traitname <- colnames(trait)[1]
    PHEN <- paste0(OUT, 'EMMAX_phenotype_', traitname, '.tsv')
    COVAR <- paste0(OUT, 'EMMAX_covariates_', traitname, '.tsv')
    EMMAXOUT <- paste0(OUT, 'EMMAX_OUT_', traitname)

    # Prepare phenotype file — use FID/IID from TFAM to ensure exact match
    tfam_ids <- fread(tfam, header = FALSE, select = list(character = 1:2),
                      col.names = c("FAMID", "INDID"))
    emmax.phen <- cbind(tfam_ids, trait) %T>%
        write.table(PHEN, sep = '\t', col.names = F, row.names = F, quote = F)

    if (!is.null(covariates)) {
        message('RUN emmax with PCA and kinship corrections')

        # Prepare CV file
        emmax.phen %>%
            dplyr::select(FAMID, INDID) %>%
            dplyr::mutate(smth = 1) %>%  # intercept
            cbind(covariates) %>%
            write.table(COVAR, sep = '\t', col.names = F, row.names = F, quote = F)

        # Run EMMAX
        system(paste('/pipeline/scripts/emmax-intel64 -v -d 10 -t', TPED_PREFIX,
                     '-p', PHEN,
                     '-k', KINSHIP_FILE,
                     '-c', COVAR,
                     '-o', EMMAXOUT))
    } else {
        message('RUN emmax without PCA correction')
        system(paste('/pipeline/scripts/emmax-intel64 -v -d 10 -t', TPED_PREFIX,
                     '-p', PHEN,
                     '-k', KINSHIP_FILE,
                     '-o', EMMAXOUT))
    }

    EMMAXOUT.ps <- paste0(EMMAXOUT, '.ps')

    # Load results
    df <- fread(EMMAXOUT.ps) %>%
        setNames(c('SNPID', 'beta', 'SE.beta', traitname)) %>%
        dplyr::mutate(
            SNPID = !!snpid,
            chr = SNPID %>% str_extract("(.*):", group = 1),
            pos = SNPID %>% str_extract(":(.*)", group = 1) %>% as.numeric
        ) %>%
        dplyr::select(SNPID, chr, pos, everything()) %>%
        dplyr::select(-beta, -SE.beta)

    return(df)
}

################################ Main

# Select only specified predictors
trait <- trait %>% dplyr::select(!!PREDICTORS_SELECTED)

# Add phenotype traits if present in metadata
samples <- fread(SAMPLES_FILE)
if (ncol(samples) > 4) {
    trait <- cbind(trait, samples %>% dplyr::select(-site, -sample, -latitude, -longitude))
}

# Run EMMAX for each trait
pval_dt <- lapply(1:ncol(trait), function(i) {
    FUN_emmax(VCF, trait[, ..i], covariates, INTER_DIR, TPED_PREFIX, KINSHIP_FILE)
}) %>%
    reduce(function(x, y) {
        left_join(x, y, by = c("SNPID", "chr", "pos"))
    })

message(pval_dt %>% str)

# Calculate q-values
qval_dt <- lapply(pval_dt %>% dplyr::select(-SNPID, -chr, -pos), function(biovec) {
    qvalue(biovec)$qvalues
}) %>%
    do.call(cbind, .) %>%
    cbind(pval_dt %>% dplyr::select(SNPID, chr, pos), .)

# Save results
pval_dt %>%
    fwrite(paste0(TABLES_DIR, "EMMAX_pvalues_K", Kbest, ".tsv"), sep = '\t')

qval_dt %>%
    fwrite(paste0(TABLES_DIR, "EMMAX_qvalues_K", Kbest, ".tsv"), sep = '\t')

message("INFO: EMMAX analysis complete")
