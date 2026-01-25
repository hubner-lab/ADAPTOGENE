#!/usr/bin/env Rscript
# Find overlap between significant SNPs from different methods (refactored)
# Outputs to specified paths

library(data.table)
library(dplyr)
library(qvalue)
library(stringr)
library(magrittr)
library(GenomicRanges)
library(parallel)
library(purrr)

args = commandArgs(trailingOnly=TRUE)
options(scipen = 99999)
################
SIGSNPS_FILES = args[1]  # space-separated list of files
METHOD = args[2]         # EMMAX, LFMM, Sum, Overlap, or PairOverlap
GAP = args[3] %>% as.numeric
PREDICTORS_SELECTED = args[4] %>% str_split(',') %>% unlist
OUTPUT_REDUNDANT = args[5]
OUTPUT_UNIQUE = args[6]
################

message(paste0('INFO: Combining significant SNPs'))
message(paste0('INFO: Method: ', METHOD))
message(paste0('INFO: Gap: ', GAP))
message(paste0('INFO: Predictors: ', paste(PREDICTORS_SELECTED, collapse = ', ')))

################ Functions

FUN_make_gr <- function(dt) {
    dt %>%
        dplyr::mutate(start = pos, end = pos) %>%
        dplyr::select(chr, start, end, SNPID, trait, method) %>%
        GRanges()
}

######################################### Main

# Parse input files (space-separated)
sigSNPs_vec <- SIGSNPS_FILES %>% str_split(' ') %>% unlist
sigSNPs_vec <- sigSNPs_vec[sigSNPs_vec != ""]  # Remove empty strings

message(paste0('INFO: Input files: ', paste(sigSNPs_vec, collapse = ', ')))

# Extract method names from file paths
methods_vec <- sigSNPs_vec %>%
    str_split('/') %>%
    sapply(function(x) x[length(x) - 1])  # Get folder name (EMMAX or LFMM)

sigSNPs_vec %<>% setNames(methods_vec)

message(paste0('INFO: Methods detected: ', paste(methods_vec, collapse = ', ')))

# Load tables, filter by climate traits
sigSNPs_lst <- lapply(sigSNPs_vec, function(x) {
    dt <- fread(x)
    if (nrow(dt) == 0) return(dt)
    dt %>%
        dplyr::filter(trait %in% !!PREDICTORS_SELECTED) %>%
        dplyr::select(SNPID, chr, pos, trait, method, pvalue)
}) %>%
    setNames(methods_vec)

# Check if any SNPs found
total_snps <- sum(sapply(sigSNPs_lst, nrow))
if (total_snps == 0) {
    message('WARNING: No significant SNPs found in any method')

    # Create empty outputs
    empty_dt <- data.table(
        SNPID = character(),
        chr = character(),
        pos = numeric(),
        trait = character(),
        method = character(),
        pvalue = numeric()
    )
    empty_dt %>% fwrite(OUTPUT_REDUNDANT, sep = '\t')
    data.table(SNPID = character()) %>% fwrite(OUTPUT_UNIQUE, sep = '\t')

    quit(save = "no", status = 0)
}

# Combine based on METHOD
if (METHOD %in% names(sigSNPs_lst)) {
    # Single method (EMMAX or LFMM)
    snpIDs <- sigSNPs_lst[[METHOD]]

} else if (METHOD == 'Sum') {
    # Combine all methods
    snpIDs <- do.call(rbind, sigSNPs_lst)

} else if (METHOD == 'Overlap') {
    # Find overlapping SNPs between methods
    sigSNPs_lstGR <- lapply(sigSNPs_lst, FUN_make_gr) %>%
        setNames(methods_vec)

    snpIDs <- lapply(methods_vec, function(n1) {
        lapply(methods_vec, function(n2) {
            if (n1 == n2) return(NULL)
            if (length(sigSNPs_lstGR[[n1]]) == 0 || length(sigSNPs_lstGR[[n2]]) == 0) return(NULL)

            overlaps_dt <- findOverlaps(sigSNPs_lstGR[[n1]], sigSNPs_lstGR[[n2]],
                                        maxgap = GAP) %>%
                as.data.table

            if (nrow(overlaps_dt) == 0) return(NULL)

            dt1 <- sigSNPs_lst[[n1]]
            dt2 <- sigSNPs_lst[[n2]]

            dt1_overlap <- dt1[overlaps_dt$queryHits, ]
            dt2_overlap <- dt2[overlaps_dt$subjectHits, ]

            rbind(dt1_overlap, dt2_overlap)
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) %>%
        unique

} else if (METHOD == 'PairOverlap') {
    # Find pair-wise overlapping SNPs per predictor
    sigSNPs_lstGR <- lapply(sigSNPs_lst, FUN_make_gr) %>%
        setNames(methods_vec)

    snpIDs <- lapply(methods_vec, function(n1) {
        lapply(methods_vec, function(n2) {
            lapply(PREDICTORS_SELECTED, function(bio) {
                if (n1 == n2) return(NULL)
                if (length(sigSNPs_lstGR[[n1]]) == 0 || length(sigSNPs_lstGR[[n2]]) == 0) return(NULL)

                gr1_bio <- sigSNPs_lstGR[[n1]][sigSNPs_lstGR[[n1]]$trait == bio]
                gr2_bio <- sigSNPs_lstGR[[n2]][sigSNPs_lstGR[[n2]]$trait == bio]

                if (length(gr1_bio) == 0 | length(gr2_bio) == 0) return(NULL)

                overlaps_dt <- findOverlaps(gr1_bio, gr2_bio, maxgap = GAP) %>%
                    as.data.table

                if (nrow(overlaps_dt) == 0) return(NULL)

                dt1 <- sigSNPs_lst[[n1]] %>% dplyr::filter(trait == !!bio)
                dt2 <- sigSNPs_lst[[n2]] %>% dplyr::filter(trait == !!bio)

                dt1_overlap <- dt1[overlaps_dt$queryHits, ]
                dt2_overlap <- dt2[overlaps_dt$subjectHits, ]

                rbind(dt1_overlap, dt2_overlap)
            }) %>% do.call(rbind, .)
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) %>%
        unique

} else {
    stop(paste0("Unknown METHOD: ", METHOD))
}

# Handle empty results
if (is.null(snpIDs) || nrow(snpIDs) == 0) {
    message('WARNING: No SNPs found with specified method')

    empty_dt <- data.table(
        SNPID = character(),
        chr = character(),
        pos = numeric(),
        trait = character(),
        method = character(),
        pvalue = numeric()
    )
    empty_dt %>% fwrite(OUTPUT_REDUNDANT, sep = '\t')
    data.table(SNPID = character()) %>% fwrite(OUTPUT_UNIQUE, sep = '\t')

} else {
    # Save redundant (all SNPs with duplicates)
    snpIDs %>%
        as.data.table %>%
        dplyr::arrange(chr, pos) %>%
        fwrite(OUTPUT_REDUNDANT, sep = '\t')

    # Save unique (deduplicated by SNPID)
    snpIDs %>%
        as.data.table %>%
        dplyr::arrange(chr, pos) %>%
        dplyr::filter(!duplicated(SNPID)) %>%
        dplyr::select(SNPID) %>%
        fwrite(OUTPUT_UNIQUE, sep = '\t')
}

message(paste0('INFO: Saved redundant SNPs to ', OUTPUT_REDUNDANT))
message(paste0('INFO: Saved unique SNPs to ', OUTPUT_UNIQUE))
message('INFO: Complete')
