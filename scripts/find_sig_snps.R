#!/usr/bin/env Rscript
# Find significant SNPs from association analysis

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(qvalue)
    library(stringr)
    library(parallel)
})

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/utils/io_selected_snps.R")
source("/pipeline/scripts/R/lib/sig_snps.R")

args = commandArgs(trailingOnly = TRUE)
################
ASSOC_TABLE       = args[1]  # p-values table
ADJUST            = args[2]  # e.g., "bonf_0.05"
CLUMPING_DISTANCE = suppressWarnings(as.numeric(args[3]))
if (is.na(CLUMPING_DISTANCE)) {
    message(paste0('INFO: snp_clumping_distance="', args[3], '" — using 1e6 bp fallback for overlap annotation'))
    CLUMPING_DISTANCE = 1e6
}
METHOD = args[4]  # EMMAX or LFMM
CPU    = as.numeric(args[5])
OUTPUT = args[6]
IS_WZA = length(args) >= 7 && args[7] == "--wza"
################

message(paste0('INFO: Finding significant SNPs for ', METHOD))
message(paste0('INFO: Adjustment: ', ADJUST))
message(paste0('INFO: SNP clumping distance: ', CLUMPING_DISTANCE))

adjustment <- str_split(ADJUST, '_')[[1]][1]
pval_value <- as.numeric(str_split(ADJUST, '_')[[1]][2])

wza_meta_cols <- if (IS_WZA) c("n_snps", "mean_maf") else character(0)

pvals <- fread(ASSOC_TABLE, colClasses = c(chr = "character"))
message(str(pvals))

result <- find_significant_snps_per_trait(pvals, adjustment, pval_value, CPU,
                                           exclude_cols = wza_meta_cols,
                                           is_wza = IS_WZA)

sig <- annotate_cross_trait_overlaps(result$sig_snps, as.integer(CLUMPING_DISTANCE))
sig[, method := METHOD]

final <- sig[, .(SNPID, chr, pos, pvalue, pval_threshold, method,
                  trait, overlap_traits, overlap_snps, overlap_distance)]
write_selected_snps(final, OUTPUT)

message(paste0('INFO: Saved ', nrow(final), ' significant SNPs to ', OUTPUT))
