#!/usr/bin/env Rscript
# EMMAX association analysis (refactored)
# Runs EMMAX for all traits specified in PREDICTORS_SELECTED

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
########################

message("INFO: Starting EMMAX analysis")
message(paste0("INFO: K = ", Kbest))
message(paste0("INFO: Predictors: ", paste(PREDICTORS_SELECTED, collapse = ", ")))

# Load trait data
trait <- fread(TRAIT_FILE, sep = '\t', header = T)

# Load covariates (PCs from eigenvec)
covariates <- fread(COVARIATES_FILE, sep = ' ', header = F)[, 3:(2 + Kbest)] %>%
    setNames(paste0('PC', 1:Kbest))

################################ Functions

FUN_emmax <- function(VCF, trait, covariates, OUT, kinship = TRUE, force = FALSE) {

    f <- paste0(gsub('.vcf', '', basename(VCF)))
    tfam <- paste0(f, '.tfam')
    kinship.f <- paste0(f, '.aIBS.kinf')

    message("INFO: Read VCF and extract SNPID")
    snpid <- fread(cmd = paste('grep -v "##"', VCF, '| cut -f1-2')) %>%
        setNames(c('CHROM', 'POS')) %>%
        dplyr::mutate(SNPID = paste0(CHROM, ':', POS)) %>%
        .$SNPID

    message("INFO: Read VCF and extract SampleNames")
    SampleName <- fread(cmd = paste("head -1000", VCF, "| grep '#CHROM' | cut -f10- "), header = F) %>%
        as.character()

    if (!file.exists(kinship.f) | force) {
        message(paste0('EMMAX: Create ', f, ' file'))

        # Prepare TPED/TFAM files
        system(paste0('plink --vcf ', VCF,
                      ' --allow-extra-chr --recode12 transpose --output-missing-genotype 0 --out ', f))

        # Edit family and individual names in .tfam
        system(paste("cat", tfam, "| awk '{split($1, a, \"_\"); split($2, b, \"_\"); if (a[1] == b[1]) {$1 = a[1]; $2 = a[1]} print}' > tmp.tfam && mv tmp.tfam", tfam))

        message('EMMAX: Calculate kinship matrix')
        system(paste('../scripts/emmax-kin-intel64 -v -s -d 10 -x', f))
    }

    # FILES produced
    traitname <- colnames(trait)[1]
    PHEN <- paste0(OUT, 'EMMAX_phenotype_', traitname, '.tsv')
    COVAR <- paste0(OUT, 'EMMAX_covariates_', traitname, '.tsv')
    EMMAXOUT <- paste0(OUT, 'EMMAX_OUT_', traitname)

    # Prepare phenotype file
    emmax.phen <- trait %>%
        dplyr::mutate(FAMID = SampleName, INDID = SampleName) %>%
        dplyr::select(FAMID, INDID, everything()) %T>%
        write.table(PHEN, sep = '\t', col.names = F, row.names = F, quote = F)

    if (!is.null(covariates) & !is.null(kinship)) {
        message('RUN emmax with PCA and kinship corrections')

        # Prepare CV file
        emmax.phen %>%
            dplyr::select(FAMID, INDID) %>%
            dplyr::mutate(smth = 1) %>%  # intercept
            cbind(covariates) %>%
            write.table(COVAR, sep = '\t', col.names = F, row.names = F, quote = F)

        # Run EMMAX
        system(paste('../scripts/emmax-intel64 -v -d 10 -t', f,
                     '-p', PHEN,
                     '-k', kinship.f,
                     '-c', COVAR,
                     '-o', EMMAXOUT))
    }

    if (is.null(covariates) & !is.null(kinship)) {
        message('RUN emmax without PCA correction')
        system(paste('../scripts/emmax-intel64 -v -d 10 -t', f,
                     '-p', PHEN,
                     '-k', kinship.f,
                     '-o', EMMAXOUT))
    }

    if (is.null(covariates) & is.null(kinship)) {
        message('RUN emmax without PCA and kinship corrections')
        system(paste('../scripts/emmax-intel64 -v -d 10 -t', f,
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
    FUN_emmax(VCF, trait[, ..i], covariates, INTER_DIR)
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
