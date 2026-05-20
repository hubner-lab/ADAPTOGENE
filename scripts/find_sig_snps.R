#!/usr/bin/env Rscript
# Find significant SNPs from association analysis (refactored)
# Handles single adjustment parameter and outputs to specified path

library(data.table)
library(dplyr)
library(qvalue)
library(stringr)
library(GenomicRanges)
library(parallel)

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/utils/io_selected_snps.R")

args = commandArgs(trailingOnly=TRUE)
################
ASSOC_TABLE       = args[1]  # p-values table
ADJUST            = args[2]  # e.g., "bonf_0.05"
CLUMPING_DISTANCE = suppressWarnings(as.numeric(args[3]))
if (is.na(CLUMPING_DISTANCE)) {
    message(paste0('INFO: snp_clumping_distance="', args[3], '" — using 1e6 bp fallback for overlap annotation'))
    CLUMPING_DISTANCE = 1e6
}
METHOD            = args[4]  # EMMAX or LFMM
CPU               = args[5] %>% as.numeric
OUTPUT            = args[6]  # Output file path
IS_WZA            = length(args) >= 7 && args[7] == "--wza"
################

message(paste0('INFO: Finding significant SNPs for ', METHOD))
message(paste0('INFO: Adjustment: ', ADJUST))
message(paste0('INFO: SNP clumping distance: ', CLUMPING_DISTANCE))
if (IS_WZA) message('INFO: WZA mode — threshold denominator = non-NA windows per trait')

######################## Functions

FUN_find_sigSNPs <- function(ASSOC_PVAL, adjustment, pval_value, CPU,
                             exclude_cols = character(0)) {

    snps_assoc <- fread(ASSOC_PVAL)
    snps_assoc$chr <- as.character(snps_assoc$chr)
    message(snps_assoc %>% str)

    traits <- setdiff(
        snps_assoc %>% dplyr::select(-SNPID, -chr, -pos) %>% colnames,
        exclude_cols
    )
    message(paste0('INFO: Traits: ', paste(traits, collapse = ', ')))

    # Find significant SNPs for each trait (threshold computed per-trait to
    # correctly handle NA windows in WZA output).
    snps_sig_dt <- lapply(traits, function(trait) {

        result <- compute_pval_threshold(snps_assoc[[trait]], adjustment, pval_value)

        message(paste0('INFO: Trait=', trait,
                       ' status=',        result$status,
                       ' threshold=',     result$threshold,
                       ' n_tested=',      result$n_tested,
                       ' n_na_dropped=',  result$n_na_dropped))

        if (result$status != "ok") {
            message(paste0('WARNING: No threshold for trait ', trait,
                           ' (status=', result$status, '). Writing empty result.'))
            return(NULL)
        }

        pval_threshold <- result$threshold

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
                overlap_distance = CLUMPING_DISTANCE
            ))
        }

        other_traits_gr <- other_traits_dt %>%
            dplyr::mutate(start = pos, end = pos) %>%
            dplyr::select(chr, start, end) %>%
            GRanges()

        overlaps <- findOverlaps(current_snp_gr, other_traits_gr,
                                 maxgap = CLUMPING_DISTANCE, ignore.strand = TRUE) %>%
            as.data.frame

        traits <- other_traits_dt[overlaps$subjectHits, ]$trait %>%
            paste(collapse = ',')

        snpids <- paste0(other_traits_dt[overlaps$subjectHits, ]$chr, ':',
                         other_traits_dt[overlaps$subjectHits, ]$pos) %>%
            paste(collapse = ',')

        return(data.frame(
            overlap_traits = ifelse(traits == "", NA, traits),
            overlap_snps = ifelse(snpids == ':', NA, snpids),
            overlap_distance = CLUMPING_DISTANCE
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
pval_value <- ADJUST %>% str_split('_') %>% unlist %>% .[2] %>% as.numeric

wza_meta_cols <- if (IS_WZA) c("n_snps", "mean_maf") else character(0)
snps_sig <- FUN_find_sigSNPs(ASSOC_TABLE, adjustment, pval_value, CPU,
                              exclude_cols = wza_meta_cols)

# Save results
snps_sig %>%
    dplyr::mutate(method = METHOD) %>%
    dplyr::select(SNPID, chr, pos, pvalue, pval_threshold, method,
                  trait, overlap_traits, overlap_snps, overlap_distance) %>%
    write_selected_snps(OUTPUT)

message(paste0('INFO: Saved significant SNPs to ', OUTPUT))
