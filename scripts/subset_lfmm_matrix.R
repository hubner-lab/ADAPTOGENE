#!/usr/bin/env Rscript
# Row-subset a generic LFMM-format matrix (space-separated, no header, one row per
# sample in VCF order) to the coord-valid sample set (climate coordinate NA handling).
# Reused across GEA-LFMM (lfmm_imp, lfmm_imp_full) and Maladaptation
# (lfmm_full for gradient_forest_model, lfmm_imp_full for geometric_offset) since all
# of these positionally bind climate values (already coord-filtered by
# filter_coord_samples) to genotype-matrix rows.

library(data.table)

args = commandArgs(trailingOnly=TRUE)
################
LFMM_IN = args[1]              # Input LFMM matrix (space-separated, no header, VCF row order)
SAMPLES_ORDER = args[2]        # Full VCF sample order, one per line (matches LFMM_IN rows)
COORD_VALID_SAMPLES = args[3]  # plink --keep format (FID IID, no header) of coord-valid samples
LFMM_OUT = args[4]             # Output: row-subset LFMM matrix
################

message('INFO: Subsetting LFMM matrix to coord-valid samples')
message(paste0('INFO: Input: ', LFMM_IN))

all_samples <- fread(SAMPLES_ORDER, header = FALSE, colClasses = "character")$V1
keep_samples <- fread(COORD_VALID_SAMPLES, header = FALSE, colClasses = "character")$V2

lfmm <- fread(LFMM_IN, sep = ' ', header = FALSE)
message(paste0('INFO: Loaded matrix: ', nrow(lfmm), ' samples x ', ncol(lfmm), ' SNPs'))

if (nrow(lfmm) != length(all_samples)) {
    stop(paste0('ERROR: LFMM row count (', nrow(lfmm), ') does not match samples_order (',
                length(all_samples), '). Both must derive from the same VCF sample order.'))
}

keep_idx <- which(all_samples %in% keep_samples)
message(paste0('INFO: Keeping ', length(keep_idx), '/', nrow(lfmm), ' samples with valid coordinates'))

lfmm_sub <- lfmm[keep_idx, ]
fwrite(lfmm_sub, LFMM_OUT, sep = ' ', col.names = FALSE, quote = FALSE)

message(paste0('INFO: Wrote subset LFMM matrix (', nrow(lfmm_sub), ' samples) to ', LFMM_OUT))
message('INFO: Complete')
