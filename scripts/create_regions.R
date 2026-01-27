#!/usr/bin/env Rscript
# Merge nearby significant SNPs into regions using single-linkage clustering
#
# REGION_DISTANCE defines the maximum gap between neighboring SNPs to be in the same region.
# SNPs are clustered using single-linkage: if SNP_A is within REGION_DISTANCE of SNP_B,
# and SNP_B is within REGION_DISTANCE of SNP_C, then all three belong to the same region,
# even if SNP_A and SNP_C are far apart.
#
# Output region_id format: CHR:START-END (e.g., "chr1:1000000-2500000")

library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
options(scipen = 99999)
################
SELECTED_SNPS = args[1]  # Selected_SNPs.tsv
REGION_DISTANCE = args[2] %>% as.numeric
OUTPUT = args[3]
################

message('INFO: Creating regions from selected SNPs using single-linkage clustering')
message(paste0('INFO: Maximum gap between neighboring SNPs: ', REGION_DISTANCE))

# Load selected SNPs
snps <- fread(SELECTED_SNPS)
if ('chr' %in% colnames(snps)) snps$chr <- as.character(snps$chr)

if (nrow(snps) == 0) {
    message('WARNING: No SNPs to process')
    empty_dt <- data.table(
        region_id = character(),
        chr = character(),
        start = integer(),
        end = integer(),
        length = integer(),
        snp_count = integer(),
        snp_ids = character(),
        traits = character(),
        methods = character(),
        min_pvalue = numeric()
    )
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
    quit(save = "no", status = 0)
}

message(paste0('INFO: Processing ', nrow(snps), ' SNPs'))

# Handle new Selected_SNPs format with method-specific columns (EMMAX, LFMM, etc.)
# Reconstruct 'traits' and 'methods' columns from method-specific columns
if (!'traits' %in% colnames(snps)) {
    message('INFO: Detected new format with method-specific columns, reconstructing traits/methods')

    # Identify method columns (everything except SNPID, chr, pos, min_pvalue)
    method_cols <- setdiff(colnames(snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
    message(paste0('INFO: Method columns: ', paste(method_cols, collapse = ', ')))

    # Build traits and methods for each SNP
    snps$traits <- sapply(1:nrow(snps), function(i) {
        all_traits <- c()
        for (m in method_cols) {
            val <- snps[[m]][i]
            if (!is.na(val) && val != '' && val != '""') {
                # Remove quotes and split by comma
                traits_m <- str_replace_all(val, '"', '') %>% str_split(',') %>% unlist
                all_traits <- c(all_traits, traits_m)
            }
        }
        paste(unique(all_traits[all_traits != '']), collapse = ',')
    })

    snps$methods <- sapply(1:nrow(snps), function(i) {
        active_methods <- c()
        for (m in method_cols) {
            val <- snps[[m]][i]
            if (!is.na(val) && val != '' && val != '""') {
                active_methods <- c(active_methods, m)
            }
        }
        paste(active_methods, collapse = ',')
    })
}

# Convert to GRanges for region merging
snps_gr <- snps %>%
    dplyr::mutate(start = pos, end = pos) %>%
    dplyr::select(chr, start, end, SNPID, traits, methods, min_pvalue) %>%
    GRanges()

# Single-linkage clustering: extend each SNP by REGION_DISTANCE/2 on both sides.
# Two SNPs merge if their extended ranges overlap, i.e., if they are within REGION_DISTANCE.
# The reduce() function then merges all overlapping ranges, achieving single-linkage clustering.
snps_extended <- snps_gr
start(snps_extended) <- pmax(1, start(snps_extended) - REGION_DISTANCE / 2)
end(snps_extended) <- end(snps_extended) + REGION_DISTANCE / 2

# Reduce to get merged regions per chromosome
regions_reduced <- reduce(snps_extended)

# Now assign SNPs to regions and compute actual region boundaries from SNP positions
regions_list <- lapply(seq_along(regions_reduced), function(i) {
    region <- regions_reduced[i]

    # Find overlapping SNPs
    overlaps <- findOverlaps(snps_extended, region)
    snp_indices <- queryHits(overlaps)

    if (length(snp_indices) == 0) return(NULL)

    # Get SNP data
    region_snps <- snps[snp_indices, ]

    # Region boundaries are defined by actual SNP positions, not extended ranges
    region_chr <- as.character(seqnames(region))
    region_start <- min(region_snps$pos)
    region_end <- max(region_snps$pos)

    data.table(
        region_id = paste0(region_chr, ':', region_start, '-', region_end),
        chr = region_chr,
        start = region_start,
        end = region_end,
        length = region_end - region_start,
        snp_count = nrow(region_snps),
        snp_ids = paste(region_snps$SNPID, collapse = ','),
        traits = paste(unique(unlist(str_split(region_snps$traits, ','))), collapse = ','),
        methods = paste(unique(unlist(str_split(region_snps$methods, ','))), collapse = ','),
        min_pvalue = min(region_snps$min_pvalue)
    )
})

# Combine all regions
regions <- do.call(rbind, regions_list)

if (is.null(regions) || nrow(regions) == 0) {
    message('WARNING: No regions created')
    empty_dt <- data.table(
        region_id = character(),
        chr = character(),
        start = integer(),
        end = integer(),
        length = integer(),
        snp_count = integer(),
        snp_ids = character(),
        traits = character(),
        methods = character(),
        min_pvalue = numeric()
    )
    empty_dt %>% fwrite(OUTPUT, sep = '\t')
} else {
    # Sort by chromosome and position
    regions <- regions %>%
        dplyr::arrange(chr, start)

    regions %>% fwrite(OUTPUT, sep = '\t')
    message(paste0('INFO: Created ', nrow(regions), ' regions from ', nrow(snps), ' SNPs'))
}

message('INFO: Complete')
