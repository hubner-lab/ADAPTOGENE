#!/usr/bin/env Rscript
# Find significant SNPs from association analysis (refactored)
# Handles single adjustment parameter and outputs to specified path

library(data.table)
library(dplyr)
library(qvalue)
library(stringr)
library(GenomicRanges)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
################
ASSOC_TABLE = args[1]  # p-values table
ADJUST = args[2]       # e.g., "bonf_0.05"
DISTANCE = args[3] %>% as.numeric  # Distance for overlap detection
METHOD = args[4]       # EMMAX or LFMM
CPU = args[5] %>% as.numeric
OUTPUT = args[6]       # Output file path
################

message(paste0('INFO: Finding significant SNPs for ', METHOD))
message(paste0('INFO: Adjustment: ', ADJUST))
message(paste0('INFO: Distance: ', DISTANCE))

######################## Functions

FUN_max_pvalue_fdr <- function(pvalues, pval_threshold) {
    qvalues_result <- qvalue(pvalues)
    significant_pvalues <- pvalues[qvalues_result$qvalues < pval_threshold]
    if (length(significant_pvalues) > 0) {
        return(max(significant_pvalues))
    } else {
        return(NA)
    }
}

FUN_max_pvalue_top <- function(pvalues, topN) {
    if (length(pvalues) < topN) {
        stop("topN is larger than the number of available p-values")
    }
    sorted_pvalues <- sort(pvalues, decreasing = FALSE)
    return(max(sorted_pvalues[1:topN]))
}

FUN_find_sigSNPs <- function(ASSOC_PVAL, adjustment, pval_threshold, CPU) {

    snps_assoc <- fread(ASSOC_PVAL)
    snps_assoc$chr <- as.character(snps_assoc$chr)
    message(snps_assoc %>% str)

    traits <- snps_assoc %>% dplyr::select(-SNPID, -chr, -pos) %>% colnames
    message(paste0('INFO: Traits: ', paste(traits, collapse = ', ')))

    # Calculate thresholds based on adjustment method
    if (adjustment == 'bonf') {
        pval_threshold <- pval_threshold / nrow(snps_assoc)
    }
    if (adjustment == 'qval') {
        pval_thresholds <- sapply(snps_assoc %>% dplyr::select(-SNPID, -chr, -pos),
                                  function(x) FUN_max_pvalue_fdr(x, pval_threshold)) %>%
            setNames(traits)
    }
    if (adjustment == 'top') {
        pval_thresholds <- sapply(snps_assoc %>% dplyr::select(-SNPID, -chr, -pos),
                                  function(x) FUN_max_pvalue_top(x, pval_threshold)) %>%
            setNames(traits)
    }

    # Find significant SNPs for each trait
    snps_sig_dt <- lapply(traits, function(trait) {

        if (adjustment %in% c('top', 'qval')) {
            pval_threshold <- pval_thresholds[trait]
        }

        if (is.na(pval_threshold)) {
            message(paste0('WARNING: No significant SNPs for trait ', trait))
            return(NULL)
        }

        mask <- snps_assoc[[trait]] < pval_threshold

        snps_assoc[mask, ] %>%
            dplyr::select(SNPID, chr, pos, !!trait) %>%
            setNames(c('SNPID', 'chr', 'pos', 'pvalue')) %>%
            dplyr::mutate(pval_threshold = !!pval_threshold,
                          trait = !!trait)
    }) %>%
        do.call(rbind, .)

    if (is.null(snps_sig_dt) || nrow(snps_sig_dt) == 0) {
        message('WARNING: No significant SNPs found')
        return(data.table(
            SNPID = character(),
            chr = character(),
            pos = numeric(),
            pvalue = numeric(),
            pval_threshold = numeric(),
            trait = character(),
            overlap_traits = character(),
            overlap_snps = character(),
            overlap_distance = numeric()
        ))
    }

    message(snps_sig_dt %>% str)

    # Add overlapping info
    overlaps_dt <- mclapply(1:nrow(snps_sig_dt), function(i) {
        trait_name <- snps_sig_dt[i, ]$trait

        current_snp_gr <- snps_sig_dt[i, ] %>%
            dplyr::mutate(start = pos, end = pos) %>%
            dplyr::select(chr, start, end) %>%
            GRanges()

        other_traits_dt <- snps_sig_dt %>%
            dplyr::filter(trait != !!trait_name)

        if (nrow(other_traits_dt) == 0) {
            return(data.frame(
                overlap_traits = NA,
                overlap_snps = NA,
                overlap_distance = DISTANCE
            ))
        }

        other_traits_gr <- other_traits_dt %>%
            dplyr::mutate(start = pos, end = pos) %>%
            dplyr::select(chr, start, end) %>%
            GRanges()

        overlaps <- findOverlaps(current_snp_gr, other_traits_gr,
                                 maxgap = DISTANCE, ignore.strand = TRUE) %>%
            as.data.frame

        traits <- other_traits_dt[overlaps$subjectHits, ]$trait %>%
            paste(collapse = ',')

        snpids <- paste0(other_traits_dt[overlaps$subjectHits, ]$chr, ':',
                         other_traits_dt[overlaps$subjectHits, ]$pos) %>%
            paste(collapse = ',')

        return(data.frame(
            overlap_traits = ifelse(traits == "", NA, traits),
            overlap_snps = ifelse(snpids == ':', NA, snpids),
            overlap_distance = DISTANCE
        ))

    }, mc.cores = CPU) %>% do.call(rbind, .)

    message(snps_sig_dt %>% str)
    message(overlaps_dt %>% str)

    snps_sig_overlap_dt <- cbind(snps_sig_dt, overlaps_dt)

    return(snps_sig_overlap_dt)
}

######################################### Main

# Parse parameters
adjustment <- ADJUST %>% str_split('_') %>% unlist %>% .[1]
pval_threshold <- ADJUST %>% str_split('_') %>% unlist %>% .[2] %>% as.numeric

snps_sig <- FUN_find_sigSNPs(ASSOC_TABLE, adjustment, pval_threshold, CPU)

# Save results
snps_sig %>%
    dplyr::mutate(method = METHOD) %>%
    dplyr::select(SNPID, chr, pos, pvalue, pval_threshold, method,
                  trait, overlap_traits, overlap_snps, overlap_distance) %>%
    fwrite(OUTPUT, sep = '\t')

message(paste0('INFO: Saved significant SNPs to ', OUTPUT))
