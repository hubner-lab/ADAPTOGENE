#!/usr/bin/env Rscript
# Find significant LD blocks from WZA block-level p-values.
# Adjustment is applied over N_blocks (not N_SNPs), making it far less conservative.

library(data.table)
library(dplyr)
library(qvalue)
library(stringr)

source("/pipeline/scripts/R/utils/pval_threshold.R")

args <- commandArgs(trailingOnly = TRUE)
BLOCK_PVALUES <- args[1]  # {method}_block_pvalues_K{k}.tsv
ADJUST        <- args[2]  # e.g., "bonf_0.05"
METHOD        <- args[3]  # method label (for method column in output)
OUTPUT        <- args[4]

message(paste0("INFO: Finding significant blocks for ", METHOD))
message(paste0("INFO: Adjustment: ", ADJUST))

adjustment     <- str_split(ADJUST, "_")[[1]][1]
pval_threshold <- str_split(ADJUST, "_")[[1]][2] %>% as.numeric

blocks <- fread(BLOCK_PVALUES)
blocks$chr <- as.character(blocks$chr)

# Trait columns: those ending in _p (block p-values from WZA)
p_cols  <- names(blocks)[endsWith(names(blocks), "_p")]
traits  <- sub("_p$", "", p_cols)

message(paste0("INFO: Traits: ", paste(traits, collapse = ", ")))
message(paste0("INFO: Total blocks: ", nrow(blocks)))

sig_rows <- lapply(traits, function(trait) {
    p_vec <- blocks[[paste0(trait, "_p")]]
    valid <- is.finite(p_vec)
    if (sum(valid) == 0) {
        message(paste0("WARNING: No valid p-values for trait ", trait))
        return(NULL)
    }
    thresh <- compute_pval_threshold(p_vec[valid], adjustment, pval_threshold)
    if (is.na(thresh)) {
        message(paste0("WARNING: No significant blocks for trait ", trait))
        return(NULL)
    }
    sig <- blocks[valid][p_vec[valid] < thresh, ]
    if (nrow(sig) == 0) return(NULL)
    message(paste0("INFO: ", nrow(sig), " significant blocks for trait ", trait,
                   " (threshold=", signif(thresh, 4), ")"))
    sig %>%
        dplyr::select(block_id, chr, start, end, n_snps,
                      Z         = !!paste0(trait, "_Z"),
                      pvalue    = !!paste0(trait, "_p"),
                      nsnps_wza = !!paste0(trait, "_nsnps")) %>%
        dplyr::mutate(pval_threshold = thresh, trait = trait)
}) %>%
    Filter(Negate(is.null), .) %>%
    rbindlist(fill = TRUE)

if (is.null(sig_rows) || nrow(sig_rows) == 0) {
    message("WARNING: No significant blocks found")
    sig_rows <- data.table(
        block_id      = character(),
        chr           = character(),
        start         = integer(),
        end           = integer(),
        n_snps        = integer(),
        Z             = numeric(),
        pvalue        = numeric(),
        nsnps_wza     = integer(),
        pval_threshold = numeric(),
        trait         = character(),
        method        = character()
    )
} else {
    sig_rows[, method := METHOD]
}

fwrite(sig_rows, OUTPUT, sep = "\t")
message(paste0("INFO: Saved ", nrow(sig_rows), " significant block×trait rows to ", OUTPUT))
