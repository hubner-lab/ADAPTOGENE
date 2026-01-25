#!/usr/bin/env Rscript
# Manhattan plot for a single trait (refactored)
# Generates Manhattan plot using CMplot

library(dplyr)
library(data.table)
library(CMplot)
library(stringr)
library(qvalue)

args = commandArgs(trailingOnly=TRUE)
################################
ASSOC_TABLE = args[1]  # SNPID chr pos TRAITS...
ADJUST = args[2]       # e.g., "bonf_0.05" or "qval_0.05"
Kbest = args[3] %>% as.numeric
METHOD = args[4]       # EMMAX or LFMM
TRAIT = args[5]        # single trait name, e.g., "bio_1"
PLOT_DIR = args[6]     # output directory
################################

message(paste0('INFO: Manhattan plot for ', METHOD, ' - ', TRAIT))
message(paste0('INFO: K = ', Kbest))
message(paste0('INFO: Adjustment: ', ADJUST))

################################ Functions

FUN_max_pvalue_fdr <- function(pvalues, pval_threshold) {
    qvalues_result <- qvalue(pvalues)
    significant_pvalues <- pvalues[qvalues_result$qvalues < pval_threshold]
    if (length(significant_pvalues) > 0) {
        return(max(significant_pvalues))
    } else {
        return(NULL)
    }
}

FUN_max_pvalue_top <- function(pvalues, topN) {
    if (length(pvalues) < topN) {
        stop("topN is larger than the number of available p-values")
    }
    sorted_pvalues <- sort(pvalues, decreasing = FALSE)
    return(max(sorted_pvalues[1:topN]))
}

################################ Main

# Load association data
message('INFO: Loading association data')
snps_assoc <- fread(ASSOC_TABLE, sep = '\t', header = T)
message(snps_assoc %>% str)

# Check if trait exists
if (!(TRAIT %in% colnames(snps_assoc))) {
    stop(paste0("Trait '", TRAIT, "' not found in association table"))
}

# Parse adjustment parameters
adjustment <- ADJUST %>% str_split('_') %>% unlist %>% .[1]
pval_threshold <- ADJUST %>% str_split('_') %>% unlist %>% .[2] %>% as.numeric

message(paste0('INFO: Adjustment method: ', adjustment))
message(paste0('INFO: Threshold: ', pval_threshold))

# Calculate threshold based on adjustment method
if (adjustment == 'bonf') {
    pval_threshold <- pval_threshold / nrow(snps_assoc)
}
if (adjustment == 'qval') {
    pval_threshold <- FUN_max_pvalue_fdr(snps_assoc[[TRAIT]], pval_threshold)
}
if (adjustment == 'top') {
    pval_threshold <- FUN_max_pvalue_top(snps_assoc[[TRAIT]], pval_threshold)
}

message(paste0('INFO: Final p-value threshold: ', pval_threshold))

# Prepare data for CMplot
plot_data <- snps_assoc %>%
    dplyr::select(SNPID, chr, pos, !!TRAIT) %>%
    dplyr::rename(SNP = SNPID, Chromosome = chr, Position = pos)

# Rename trait column for plot
new_colname <- paste0(TRAIT, '_K', Kbest, '_', ADJUST)
colnames(plot_data)[4] <- new_colname

# Change to output directory
setwd(PLOT_DIR)

message('INFO: Generating Manhattan plot')

# Generate plot
CMplot(plot_data,
       plot.type = "m",
       LOG10 = TRUE,
       ylim = NULL,
       threshold = pval_threshold,
       threshold.lty = 2,
       threshold.lwd = 1,
       threshold.col = "black",
       amplify = TRUE,
       bin.size = 1e6,
       chr.den.col = c("darkgreen", "yellow", "red"),
       signal.col = 'red',
       signal.cex = 0.5,
       cex = 0.25,
       file = "pdf",
       dpi = 600,
       multracks = FALSE,
       file.output = TRUE,
       verbose = TRUE,
       chr.labels.angle = 60,
       width = 8,
       height = 4
)

# Rename output file to match expected pattern
# CMplot creates files with "Rect_Manhtn." prefix
old_name <- paste0("Rect_Manhtn.", new_colname, ".pdf")
new_name <- paste0("Manhattan_", TRAIT, "_K", Kbest, "_", ADJUST, ".pdf")

message(paste0('INFO: Looking for file: ', old_name))
message(paste0('INFO: Will rename to: ', new_name))

if (file.exists(old_name)) {
    file.rename(old_name, new_name)
    message('INFO: File renamed successfully')
} else {
    # Try alternative naming pattern
    alt_name <- paste0("Rectangular-Manhattan.", new_colname, ".pdf")
    if (file.exists(alt_name)) {
        file.rename(alt_name, new_name)
        message('INFO: File renamed successfully (alt pattern)')
    } else {
        message(paste0('WARNING: Could not find output file. Files in directory:'))
        message(paste(list.files(), collapse = ', '))
    }
}

message('INFO: Manhattan plot complete')
