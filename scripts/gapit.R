#!/usr/bin/env Rscript
# GAPIT3 association analysis wrapper
# Runs GAPIT models (GLM, MLM, CMLM, ECMLM, SUPER, MLMM, FarmCPU, BLINK)
# Uses pre-computed GD/GM, kinship (BN from emmax-kin), and LEA PCA as covariates.
# Outputs standardized pvalue TSVs matching EMMAX/LFMM format.

library(data.table)
library(dplyr)
library(stringr)
library(GAPIT)

args <- commandArgs(trailingOnly = TRUE)
########################
GD_FILE         <- args[1]   # GAPIT numeric GD (Taxa + SNP columns)
GM_FILE         <- args[2]   # GAPIT numeric GM (Name, Chromosome, Position)
TRAIT_FILE      <- args[3]   # Climate site scaled or phenotype file
PCA_FILE        <- args[4]   # LEA PCA projections (space-sep, no header)
KINSHIP_FILE    <- args[5]   # .aBN.kinf from emmax-kin (n x n, no header/rownames)
K_BEST          <- as.integer(args[6])
MODELS          <- strsplit(args[7], ",")[[1]]  # Comma-separated model names
WORKDIR         <- args[8]   # Temp working directory for GAPIT native output
TABLES_DIR      <- args[9]   # Output dir for standardized pvalue TSVs
PREDICTORS      <- strsplit(args[10], ",")[[1]] # Trait column names
SAMPLES_FILE    <- args[11]  # Sample names file (for Taxa matching)
NATIVE_OUTDIR   <- args[12]  # Directory for organized GAPIT native output
SAMPLES_SUBSET  <- if (length(args) >= 13 && args[13] != "NULL") args[13] else NULL
OUTPUT_PREFIX   <- if (length(args) >= 14 && args[14] != "NULL") args[14] else NULL  # For DROP mode: trait name prefix
########################

message("INFO: Starting GAPIT analysis")
message("INFO: Models: ", paste(MODELS, collapse = ", "))
message("INFO: Traits: ", paste(PREDICTORS, collapse = ", "))
message("INFO: K = ", K_BEST)

# --- Load GD and GM ---
message("INFO: Loading GD and GM")
myGD <- fread(GD_FILE, header = TRUE)
myGM <- fread(GM_FILE, header = TRUE)
message("INFO: GD: ", nrow(myGD), " samples x ", ncol(myGD) - 1, " SNPs")
message("INFO: GM: ", nrow(myGM), " SNPs")

# --- Load traits ---
message("INFO: Loading trait data")
traits <- fread(TRAIT_FILE, sep = '\t', header = TRUE)

# --- Load PCA covariates ---
message("INFO: Loading PCA projections")
pca_raw <- fread(PCA_FILE, sep = ' ', header = FALSE)
pca_cols <- pca_raw[, 1:K_BEST]
setnames(pca_cols, paste0("PC", 1:K_BEST))

# --- Load kinship ---
message("INFO: Loading kinship matrix")
kin_raw <- fread(KINSHIP_FILE, header = FALSE)
message("INFO: Kinship dimensions: ", nrow(kin_raw), " x ", ncol(kin_raw))

# --- Load sample names for ordering ---
sample_names <- myGD$Taxa

# --- Subset if DROP mode ---
if (!is.null(SAMPLES_SUBSET)) {
    message("INFO: DROP mode — subsetting to sample list: ", SAMPLES_SUBSET)
    keep_samples <- fread(SAMPLES_SUBSET, header = FALSE)$V1

    # Subset GD
    keep_idx <- which(sample_names %in% keep_samples)
    myGD <- myGD[keep_idx, ]

    # Subset PCA
    pca_cols <- pca_cols[keep_idx, ]

    # Subset kinship (rows and columns)
    kin_raw <- kin_raw[keep_idx, ..keep_idx]

    sample_names <- myGD$Taxa
    message("INFO: Subsetted to ", length(sample_names), " samples")
}

# --- Build Y (phenotype matrix) ---
# GAPIT expects: Taxa, trait1, trait2, ...
message("INFO: Building phenotype matrix")
myY <- data.table(Taxa = sample_names)

# Match trait data to sample order
# Trait file might have a sample/Taxa column or be in the same order as metadata
if ("sample" %in% colnames(traits)) {
    traits_ordered <- data.table(sample = sample_names) %>%
        merge(traits, by = "sample", sort = FALSE)
    for (pred in PREDICTORS) {
        myY[[pred]] <- traits_ordered[[pred]]
    }
} else {
    # Traits are in same order as samples (from metadata)
    # Need to load samples file to match order
    all_samples <- fread(SAMPLES_FILE, header = TRUE)
    if ("sample" %in% colnames(all_samples)) {
        traits_with_sample <- cbind(data.table(sample = all_samples$sample), traits)
        traits_ordered <- data.table(sample = sample_names) %>%
            merge(traits_with_sample, by = "sample", sort = FALSE)
        for (pred in PREDICTORS) {
            myY[[pred]] <- traits_ordered[[pred]]
        }
    } else {
        # Assume same order
        for (pred in PREDICTORS) {
            myY[[pred]] <- traits[[pred]]
        }
    }
}

# --- Build CV (covariates) ---
# GAPIT expects: Taxa, PC1, PC2, ...
message("INFO: Building covariates matrix")
myCV <- data.table(Taxa = sample_names)
myCV <- cbind(myCV, pca_cols)

# --- Build KI (kinship) ---
# GAPIT expects: Taxa column + n x n matrix
message("INFO: Building kinship matrix for GAPIT")
myKI <- data.table(Taxa = sample_names)
myKI <- cbind(myKI, kin_raw)

# --- Create work directory ---
dir.create(WORKDIR, recursive = TRUE, showWarnings = FALSE)
old_wd <- getwd()
setwd(WORKDIR)

# --- Run GAPIT ---
message("INFO: Running GAPIT with models: ", paste(MODELS, collapse = ", "))
tryCatch({
    gapit_result <- GAPIT(
        Y = as.data.frame(myY),
        GD = as.data.frame(myGD),
        GM = as.data.frame(myGM),
        CV = as.data.frame(myCV),
        KI = as.data.frame(myKI),
        model = MODELS,
        PCA.total = 0,  # We provide our own PCA via CV
        file.output = TRUE
    )
}, error = function(e) {
    message("ERROR: GAPIT failed: ", e$message)
    setwd(old_wd)
    stop(e)
})

setwd(old_wd)

# --- Parse GAPIT output and create standardized pvalue tables ---
message("INFO: Parsing GAPIT results")

for (model in MODELS) {
    message("INFO: Processing model: ", model)

    # Collect p-values across traits
    pval_list <- list()

    for (trait in PREDICTORS) {
        # GAPIT output file pattern
        result_file <- file.path(WORKDIR, paste0("GAPIT.Association.GWAS_Results.", model, ".", trait, ".csv"))

        if (!file.exists(result_file)) {
            message("WARNING: Result file not found for ", model, "/", trait, ": ", result_file)
            next
        }

        gwas <- fread(result_file)

        # GAPIT result columns: SNP, Chr, Pos, P.value, MAF, nobs, ...
        pvals <- data.table(
            SNPID = gwas$SNP,
            chr = as.character(gwas$Chr),
            pos = as.integer(gwas$Pos)
        )
        pvals[[trait]] <- gwas$P.value

        pval_list[[trait]] <- pvals
    }

    if (length(pval_list) == 0) {
        message("WARNING: No results for model ", model, ", skipping")
        next
    }

    # Merge all traits into single table
    pval_dt <- pval_list[[1]]
    if (length(pval_list) > 1) {
        for (i in 2:length(pval_list)) {
            pval_dt <- merge(pval_dt, pval_list[[i]], by = c("SNPID", "chr", "pos"), all = TRUE)
        }
    }

    # Sort by chr, pos
    pval_dt <- pval_dt[order(chr, pos)]

    # Write standardized output
    model_dir <- file.path(TABLES_DIR, model)
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

    file_prefix <- if (!is.null(OUTPUT_PREFIX)) OUTPUT_PREFIX else model
    out_file <- file.path(model_dir, paste0(file_prefix, "_pvalues_K", K_BEST, ".tsv"))
    fwrite(pval_dt, out_file, sep = "\t")
    message("INFO: Wrote ", out_file, " (", nrow(pval_dt), " SNPs x ", length(PREDICTORS), " traits)")

    # Copy native GAPIT output files to organized directory
    native_model_dir <- file.path(NATIVE_OUTDIR, model)
    dir.create(native_model_dir, recursive = TRUE, showWarnings = FALSE)

    native_files <- list.files(WORKDIR, pattern = paste0("GAPIT.*", model), full.names = TRUE)
    if (length(native_files) > 0) {
        file.copy(native_files, native_model_dir, overwrite = TRUE)
        message("INFO: Copied ", length(native_files), " native GAPIT files to ", native_model_dir)
    }
}

message("INFO: GAPIT analysis complete")
