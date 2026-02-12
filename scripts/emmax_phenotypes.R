#!/usr/bin/env Rscript
# EMMAX association analysis for phenotype traits
# Receives pre-computed TPED, kinship, and PCA projections from separate rules.
# Formats inputs for EMMAX, runs per-trait, parses output, computes q-values.

library(data.table)
library(dplyr)
library(magrittr)
library(purrr)
library(stringr)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
########################
VCF              = args[1]                # Filtered VCF (possibly per-trait subset for DROP mode)
TPED_PREFIX      = args[2]                # Path prefix for pre-computed .tped/.tfam
KINSHIP_FILE     = args[3]                # Pre-computed .aIBS.kinf
PCA_PROJECTIONS  = args[4]                # LEA PCA projections (space-sep, no header)
K_BEST           = args[5] %>% as.numeric # Number of PCs as covariates
PHENOTYPE_FILE   = args[6]                # TSV: sample, trait1, [trait2, ...]
TABLES_DIR       = args[7]                # Output directory
SAMPLES_ORDER    = if (length(args) >= 8) args[8] else NULL  # Optional: full sample order file
OUTPUT_PVALUES   = if (length(args) >= 9) args[9] else NULL  # Optional: explicit output p-values path
########################

message("INFO: Starting EMMAX phenotype analysis")
message(paste0("INFO: VCF: ", VCF))
message(paste0("INFO: TPED prefix: ", TPED_PREFIX))
message(paste0("INFO: Kinship: ", KINSHIP_FILE))
message(paste0("INFO: PCA projections: ", PCA_PROJECTIONS))
message(paste0("INFO: K = ", K_BEST))
message(paste0("INFO: Phenotype file: ", PHENOTYPE_FILE))
message(paste0("INFO: Samples order: ", ifelse(is.null(SAMPLES_ORDER), "not provided", SAMPLES_ORDER)))

# Create output directory
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Extract SNP IDs from VCF ---
message("INFO: Extracting SNP IDs from VCF")
snpid <- fread(cmd = paste('grep -v "##"', VCF, '| cut -f1-2')) %>%
    setNames(c('CHROM', 'POS')) %>%
    dplyr::mutate(SNPID = paste0(CHROM, ':', POS)) %>%
    .$SNPID

# --- Extract sample names from VCF header ---
message("INFO: Extracting sample names from VCF")
SampleName <- fread(cmd = paste("head -1000", VCF, "| grep '#CHROM' | cut -f10-"), header = F) %>%
    as.character()
message(paste0("INFO: ", length(SampleName), " samples in VCF"))

# --- Load PCA projections and subset if needed ---
message("INFO: Loading PCA projections")
cov_raw <- fread(PCA_PROJECTIONS, sep = ' ', header = F)
message(paste0("INFO: PCA matrix: ", nrow(cov_raw), " rows x ", ncol(cov_raw), " cols"))

if (nrow(cov_raw) > length(SampleName) && !is.null(SAMPLES_ORDER)) {
    # DROP mode: PCA is from whole dataset, VCF is per-trait subset
    message("INFO: PCA has more rows than VCF samples — subsetting (DROP mode)")
    all_samples <- fread(SAMPLES_ORDER, header = FALSE)$V1
    message(paste0("INFO: Full sample order: ", length(all_samples), " samples"))
    keep_idx <- which(all_samples %in% SampleName)
    cov_raw <- cov_raw[keep_idx, ]
    message(paste0("INFO: Subsetted PCA to ", nrow(cov_raw), " rows"))
}

if (nrow(cov_raw) != length(SampleName)) {
    stop(paste0("ERROR: PCA row count (", nrow(cov_raw),
                ") does not match VCF sample count (", length(SampleName),
                "). Provide SAMPLES_ORDER for DROP mode."))
}

covariates <- cov_raw[, 1:K_BEST] %>% setNames(paste0('PC', 1:K_BEST))

# --- Read phenotype file ---
message("INFO: Reading phenotype file")
pheno <- fread(PHENOTYPE_FILE, sep = '\t', header = TRUE)
trait_cols <- colnames(pheno)[colnames(pheno) != 'sample']
message(paste0("INFO: Traits: ", paste(trait_cols, collapse = ', ')))
message(paste0("INFO: Phenotype samples: ", nrow(pheno)))

# Match phenotype sample order to VCF sample order
pheno_ordered <- data.table(sample = SampleName) %>%
    left_join(pheno, by = 'sample')

if (any(is.na(pheno_ordered$sample))) {
    stop("ERROR: Some VCF samples not found in phenotype file")
}

# --- Run EMMAX per trait ---
work_dir <- dirname(TPED_PREFIX)
results <- lapply(trait_cols, function(trait) {
    message(paste0("INFO: Processing trait: ", trait))

    trait_vals <- pheno_ordered[[trait]]
    if (all(is.na(trait_vals))) {
        message(paste0("WARNING: All values NA for trait ", trait, ", skipping"))
        return(NULL)
    }

    # EMMAX phenotype file: FID IID value (no header)
    phen_file <- file.path(work_dir, paste0('EMMAX_phenotype_', trait, '.tsv'))
    data.table(FID = SampleName, IID = SampleName, value = trait_vals) %>%
        fwrite(phen_file, sep = '\t', col.names = FALSE)

    # EMMAX covariate file: FID IID 1 PC1 PC2 ... (no header)
    covar_file <- file.path(work_dir, paste0('EMMAX_covariates_', trait, '.tsv'))
    data.table(FID = SampleName, IID = SampleName, intercept = 1) %>%
        cbind(covariates) %>%
        fwrite(covar_file, sep = '\t', col.names = FALSE)

    # Run EMMAX
    out_prefix <- file.path(work_dir, paste0('EMMAX_OUT_', trait))
    cmd <- paste('/pipeline/scripts/emmax-intel64 -v -d 10',
                 '-t', TPED_PREFIX,
                 '-p', phen_file,
                 '-k', KINSHIP_FILE,
                 '-c', covar_file,
                 '-o', out_prefix)
    message(paste0("INFO: Running: ", cmd))
    system(cmd)

    # Parse .ps output
    ps_file <- paste0(out_prefix, '.ps')
    if (!file.exists(ps_file)) {
        message(paste0("WARNING: EMMAX output not found for trait ", trait))
        return(NULL)
    }

    df <- fread(ps_file) %>%
        setNames(c('SNPID', 'beta', 'SE.beta', trait)) %>%
        dplyr::mutate(
            SNPID = !!snpid,
            chr = SNPID %>% str_extract("(.*):", group = 1),
            pos = SNPID %>% str_extract(":(.*)", group = 1) %>% as.numeric
        ) %>%
        dplyr::select(SNPID, chr, pos, all_of(trait))

    message(paste0("INFO: ", trait, " — ", nrow(df), " SNPs"))
    return(df)
})

# Remove NULLs and merge
results <- results[!sapply(results, is.null)]

if (length(results) == 0) {
    stop("ERROR: No traits produced results")
}

pval_dt <- results %>%
    reduce(function(x, y) left_join(x, y, by = c("SNPID", "chr", "pos")))

message(paste0("INFO: Combined p-value table: ", nrow(pval_dt), " SNPs x ", length(trait_cols), " traits"))

# --- Compute q-values ---
message("INFO: Computing q-values")
qval_dt <- lapply(pval_dt %>% dplyr::select(-SNPID, -chr, -pos), function(pvals) {
    tryCatch(
        qvalue(pvals)$qvalues,
        error = function(e) {
            message(paste0("WARNING: qvalue failed (", e$message, "), using BH adjustment"))
            p.adjust(pvals, method = 'BH')
        }
    )
}) %>%
    do.call(cbind, .) %>%
    as.data.table() %>%
    cbind(pval_dt %>% dplyr::select(SNPID, chr, pos), .)

# --- Write output ---
if (!is.null(OUTPUT_PVALUES)) {
    pval_path <- OUTPUT_PVALUES
    qval_path <- sub('pvalues', 'qvalues', OUTPUT_PVALUES)
} else {
    pval_path <- file.path(TABLES_DIR, paste0("EMMAX_phenotypes_pvalues_K", K_BEST, ".tsv"))
    qval_path <- file.path(TABLES_DIR, paste0("EMMAX_phenotypes_qvalues_K", K_BEST, ".tsv"))
}

fwrite(pval_dt, pval_path, sep = '\t')
fwrite(qval_dt, qval_path, sep = '\t')

message(paste0("INFO: Wrote p-values to ", pval_path))
message(paste0("INFO: Wrote q-values to ", qval_path))
message("INFO: EMMAX phenotype analysis complete")
