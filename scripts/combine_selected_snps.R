#!/usr/bin/env Rscript
# Combine significant SNPs from different methods into Selected_SNPs table
# Output: cleaned table with one row per unique SNPID (collapsed traits/methods)

library(data.table)
library(dplyr)
library(stringr)
library(magrittr)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
options(scipen = 99999)
################
SIGSNPS_FILES = args[1]  # space-separated list of files
METHOD = args[2]         # EMMAX, LFMM, Sum, Overlap, or PairOverlap
GAP = args[3] %>% as.numeric
PREDICTORS_SELECTED = args[4] %>% str_split(',') %>% unlist
OUTPUT = args[5]
################

message('INFO: Combining significant SNPs')
message(paste0('INFO: Method: ', METHOD))
message(paste0('INFO: Gap: ', GAP))

################ Functions

FUN_make_gr <- function(dt) {
    dt %>%
        dplyr::mutate(start = pos, end = pos) %>%
        dplyr::select(chr, start, end, SNPID, trait, method) %>%
        GRanges()
}

######################################### Main

# Parse input files
sigSNPs_vec <- SIGSNPS_FILES %>% str_split(' ') %>% unlist
sigSNPs_vec <- sigSNPs_vec[sigSNPs_vec != ""]

message(paste0('INFO: Input files: ', paste(sigSNPs_vec, collapse = ', ')))

# Extract method names from file paths
methods_vec <- sigSNPs_vec %>%
    str_split('/') %>%
    sapply(function(x) x[length(x) - 1])

sigSNPs_vec %<>% setNames(methods_vec)

# Load tables, filter by climate traits
sigSNPs_lst <- lapply(sigSNPs_vec, function(x) {
    dt <- fread(x)
    if (nrow(dt) == 0) return(dt)
    dt$chr <- as.character(dt$chr)
    dt %>%
        dplyr::filter(trait %in% !!PREDICTORS_SELECTED) %>%
        dplyr::select(SNPID, chr, pos, trait, method, pvalue)
}) %>%
    setNames(methods_vec)

# Check if any SNPs found
total_snps <- sum(sapply(sigSNPs_lst, nrow))
if (total_snps == 0) {
    message('WARNING: No significant SNPs found in any method')
    empty_dt <- data.table(
        SNPID = character(),
        chr = character(),
        pos = integer(),
        traits = character(),
        methods = character(),
        min_pvalue = numeric()
    )
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
    quit(save = "no", status = 0)
}

# Combine based on METHOD
if (METHOD %in% names(sigSNPs_lst)) {
    snpIDs <- sigSNPs_lst[[METHOD]]
} else if (METHOD == 'Sum') {
    snpIDs <- do.call(rbind, sigSNPs_lst)
} else if (METHOD == 'Overlap') {
    sigSNPs_lstGR <- lapply(sigSNPs_lst, FUN_make_gr) %>% setNames(methods_vec)
    snpIDs <- lapply(methods_vec, function(n1) {
        lapply(methods_vec, function(n2) {
            if (n1 == n2) return(NULL)
            if (length(sigSNPs_lstGR[[n1]]) == 0 || length(sigSNPs_lstGR[[n2]]) == 0) return(NULL)
            overlaps_dt <- findOverlaps(sigSNPs_lstGR[[n1]], sigSNPs_lstGR[[n2]], maxgap = GAP) %>% as.data.table
            if (nrow(overlaps_dt) == 0) return(NULL)
            dt1 <- sigSNPs_lst[[n1]]
            dt2 <- sigSNPs_lst[[n2]]
            rbind(dt1[overlaps_dt$queryHits, ], dt2[overlaps_dt$subjectHits, ])
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) %>% unique
} else if (METHOD == 'PairOverlap') {
    sigSNPs_lstGR <- lapply(sigSNPs_lst, FUN_make_gr) %>% setNames(methods_vec)
    snpIDs <- lapply(methods_vec, function(n1) {
        lapply(methods_vec, function(n2) {
            lapply(PREDICTORS_SELECTED, function(bio) {
                if (n1 == n2) return(NULL)
                if (length(sigSNPs_lstGR[[n1]]) == 0 || length(sigSNPs_lstGR[[n2]]) == 0) return(NULL)
                gr1_bio <- sigSNPs_lstGR[[n1]][sigSNPs_lstGR[[n1]]$trait == bio]
                gr2_bio <- sigSNPs_lstGR[[n2]][sigSNPs_lstGR[[n2]]$trait == bio]
                if (length(gr1_bio) == 0 | length(gr2_bio) == 0) return(NULL)
                overlaps_dt <- findOverlaps(gr1_bio, gr2_bio, maxgap = GAP) %>% as.data.table
                if (nrow(overlaps_dt) == 0) return(NULL)
                dt1 <- sigSNPs_lst[[n1]] %>% dplyr::filter(trait == !!bio)
                dt2 <- sigSNPs_lst[[n2]] %>% dplyr::filter(trait == !!bio)
                rbind(dt1[overlaps_dt$queryHits, ], dt2[overlaps_dt$subjectHits, ])
            }) %>% do.call(rbind, .)
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) %>% unique
} else {
    stop(paste0("Unknown METHOD: ", METHOD))
}

# Handle empty results
if (is.null(snpIDs) || nrow(snpIDs) == 0) {
    message('WARNING: No SNPs found with specified method')
    # Create empty table with method columns
    empty_dt <- data.table(
        SNPID = character(),
        chr = character(),
        pos = integer()
    )
    for (m in methods_vec) {
        empty_dt[[m]] <- character()
    }
    empty_dt$min_pvalue <- numeric()
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
} else {
    # Get unique SNPIDs from selected SNPs
    selected_snp_ids <- unique(snpIDs$SNPID)

    # Build per-method trait attribution for selected SNPs
    # For each SNP, list which traits it's significant for in each method
    method_traits <- lapply(methods_vec, function(m) {
        if (is.null(sigSNPs_lst[[m]]) || nrow(sigSNPs_lst[[m]]) == 0) {
            return(data.frame(SNPID = character(), traits = character(), stringsAsFactors = FALSE))
        }
        sigSNPs_lst[[m]] %>%
            dplyr::filter(SNPID %in% selected_snp_ids) %>%
            group_by(SNPID) %>%
            summarise(traits = paste(unique(trait), collapse = ','), .groups = 'drop') %>%
            as.data.frame()
    }) %>% setNames(methods_vec)

    # Build base result with unique SNPs
    result <- snpIDs %>%
        as.data.table %>%
        group_by(SNPID, chr, pos) %>%
        summarise(min_pvalue = min(pvalue), .groups = 'drop') %>%
        dplyr::arrange(chr, pos) %>%
        as.data.frame()

    # Add method-specific trait columns
    for (m in methods_vec) {
        mt <- method_traits[[m]]
        if (nrow(mt) > 0) {
            result <- result %>%
                left_join(mt %>% dplyr::rename(!!m := traits), by = "SNPID")
        } else {
            result[[m]] <- NA_character_
        }
        # Replace NA with empty string for cleaner output
        result[[m]][is.na(result[[m]])] <- ""
    }

    # Reorder columns: SNPID, chr, pos, <methods>, min_pvalue
    col_order <- c("SNPID", "chr", "pos", methods_vec, "min_pvalue")
    result <- result[, col_order]

    result %>% as.data.table %>% fwrite(OUTPUT, sep = '\t')
    message(paste0('INFO: Saved ', nrow(result), ' unique SNPs to ', OUTPUT))

    # Log summary per method
    for (m in methods_vec) {
        n_with_method <- sum(result[[m]] != "", na.rm = TRUE)
        message(paste0('INFO: ', m, ': ', n_with_method, ' SNPs with significant traits'))
    }
}

message('INFO: Complete')
