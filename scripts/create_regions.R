#!/usr/bin/env Rscript
# Merge nearby significant SNPs into regions using single-linkage clustering
#
# Creates TWO types of regions:
# 1. Per-trait regions: SNPs clustered separately for each trait (bio_1, bio_12, etc.)
# 2. Combined climate regions: All SNPs clustered together regardless of trait
#
# REGION_DISTANCE defines the maximum gap between neighboring SNPs to be in the same region.
# SNPs are clustered using single-linkage: if SNP_A is within REGION_DISTANCE of SNP_B,
# and SNP_B is within REGION_DISTANCE of SNP_C, then all three belong to the same region,
# even if SNP_A and SNP_C are far apart.
#
# Output region_id format: CHR:START-END_TRAIT (e.g., "1:1000000-2500000_bio_1")

library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
options(scipen = 99999)
################
SELECTED_SNPS = args[1]  # Selected_SNPs.tsv
REGION_DISTANCE = args[2] %>% as.numeric
OUTPUT_PER_TRAIT = args[3]  # Regions_per_trait.tsv
OUTPUT_COMBINED = args[4]   # Regions_climate_combined.tsv
################

message('INFO: Creating regions from selected SNPs using single-linkage clustering')
message(paste0('INFO: Maximum gap between neighboring SNPs: ', REGION_DISTANCE))
message('INFO: Producing two outputs: per-trait regions and combined climate regions')

######################## Functions

# Create regions from SNP data using single-linkage clustering
create_regions_from_snps <- function(snps_dt, trait_label = NULL) {
    if (nrow(snps_dt) == 0) {
        return(data.table(
            region_id = character(),
            trait = character(),
            chr = character(),
            start = integer(),
            end = integer(),
            length = integer(),
            snp_count = integer(),
            snp_ids = character(),
            methods = character(),
            min_pvalue = numeric()
        ))
    }

    # Convert to GRanges
    snps_gr <- snps_dt %>%
        dplyr::mutate(start = pos, end = pos) %>%
        dplyr::select(chr, start, end, SNPID, min_pvalue) %>%
        GRanges()

    # Extend SNPs by REGION_DISTANCE/2 for single-linkage clustering
    snps_extended <- snps_gr
    start(snps_extended) <- pmax(1, start(snps_extended) - REGION_DISTANCE / 2)
    end(snps_extended) <- end(snps_extended) + REGION_DISTANCE / 2

    # Reduce to get merged regions
    regions_reduced <- reduce(snps_extended)

    # Assign SNPs to regions
    regions_list <- lapply(seq_along(regions_reduced), function(i) {
        region <- regions_reduced[i]

        # Find overlapping SNPs
        overlaps <- findOverlaps(snps_extended, region)
        snp_indices <- queryHits(overlaps)

        if (length(snp_indices) == 0) return(NULL)

        # Get SNP data
        region_snps <- snps_dt[snp_indices, ]

        # Region boundaries: extend SNP range by REGION_DISTANCE on each side
        region_chr <- as.character(seqnames(region))
        region_start <- as.integer(pmax(1, min(region_snps$pos) - REGION_DISTANCE))
        region_end <- as.integer(max(region_snps$pos) + REGION_DISTANCE)

        # Extract methods from method-specific columns
        method_cols <- setdiff(colnames(region_snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
        active_methods <- c()
        for (m in method_cols) {
            vals <- region_snps[[m]]
            if (any(!is.na(vals) & vals != '' & vals != '""')) {
                active_methods <- c(active_methods, m)
            }
        }
        methods_str <- paste(unique(active_methods), collapse = ',')

        # Create region ID (uses extended boundaries)
        region_id <- if (!is.null(trait_label)) {
            paste0(region_chr, ':', region_start, '-', region_end, '_', trait_label)
        } else {
            paste0(region_chr, ':', region_start, '-', region_end)
        }

        data.table(
            region_id = region_id,
            trait = ifelse(!is.null(trait_label), trait_label, ''),
            chr = region_chr,
            start = region_start,
            end = region_end,
            length = region_end - region_start,
            snp_count = nrow(region_snps),
            snp_ids = paste(region_snps$SNPID, collapse = ','),
            methods = methods_str,
            min_pvalue = min(region_snps$min_pvalue)
        )
    })

    # Combine all regions
    regions <- do.call(rbind, regions_list)

    if (!is.null(regions)) {
        regions <- regions %>% dplyr::arrange(chr, start)
    }

    return(regions)
}

#################### Main

# Load selected SNPs
snps <- fread(SELECTED_SNPS)
if ('chr' %in% colnames(snps)) snps$chr <- as.character(snps$chr)

if (nrow(snps) == 0) {
    message('WARNING: No SNPs to process')
    empty_dt <- data.table(
        region_id = character(),
        trait = character(),
        chr = character(),
        start = integer(),
        end = integer(),
        length = integer(),
        snp_count = integer(),
        snp_ids = character(),
        methods = character(),
        min_pvalue = numeric()
    )
    empty_dt %>% fwrite(OUTPUT_PER_TRAIT, sep = '\t')
    empty_dt %>% fwrite(OUTPUT_COMBINED, sep = '\t')
    quit(save = "no", status = 0)
}

message(paste0('INFO: Processing ', nrow(snps), ' SNPs'))

# Extract unique traits from method-specific columns
method_cols <- setdiff(colnames(snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
all_traits <- c()
for (m in method_cols) {
    traits_m <- snps[[m]] %>%
        str_replace_all('"', '') %>%
        str_split(',') %>%
        unlist() %>%
        unique()
    traits_m <- traits_m[traits_m != '' & !is.na(traits_m)]
    all_traits <- c(all_traits, traits_m)
}
all_traits <- unique(all_traits)

message(paste0('INFO: Found ', length(all_traits), ' unique traits: ', paste(all_traits, collapse = ', ')))

#=============================================================================
# 1. CREATE PER-TRAIT REGIONS
#=============================================================================
message('INFO: Creating per-trait regions...')

per_trait_regions_list <- lapply(all_traits, function(trait) {
    message(paste0('INFO:   Processing trait: ', trait))

    # Filter SNPs for this trait
    trait_snps <- snps[0, ]  # Empty table with correct structure
    for (m in method_cols) {
        # Get SNPs where this method has this trait
        snps_with_trait <- snps[grepl(paste0('(^|,)', trait, '($|,)'), snps[[m]]), ]
        if (nrow(snps_with_trait) > 0) {
            trait_snps <- rbind(trait_snps, snps_with_trait)
        }
    }
    trait_snps <- unique(trait_snps)

    if (nrow(trait_snps) == 0) {
        message(paste0('INFO:     No SNPs found for ', trait))
        return(NULL)
    }

    message(paste0('INFO:     Found ', nrow(trait_snps), ' SNPs for ', trait))

    # Create regions for this trait
    regions <- create_regions_from_snps(trait_snps, trait_label = trait)

    if (!is.null(regions) && nrow(regions) > 0) {
        message(paste0('INFO:     Created ', nrow(regions), ' regions for ', trait))
    }

    return(regions)
})

# Combine all per-trait regions
per_trait_regions <- do.call(rbind, per_trait_regions_list)

if (is.null(per_trait_regions) || nrow(per_trait_regions) == 0) {
    message('WARNING: No per-trait regions created')
    per_trait_regions <- data.table(
        region_id = character(),
        trait = character(),
        chr = character(),
        start = integer(),
        end = integer(),
        length = integer(),
        snp_count = integer(),
        snp_ids = character(),
        methods = character(),
        min_pvalue = numeric()
    )
}

# Add cross-trait evidence: for each region, find SNPs from OTHER traits that fall within or near it
if (nrow(per_trait_regions) > 0) {
    message('INFO: Computing cross-trait evidence for per-trait regions...')
    other_traits_col <- character(nrow(per_trait_regions))
    other_snp_count_col <- integer(nrow(per_trait_regions))

    for (i in seq_len(nrow(per_trait_regions))) {
        region <- per_trait_regions[i, ]
        region_trait <- region$trait
        region_chr <- region$chr
        # Use pre-extended region boundaries as-is (already extended in create_regions_from_snps)
        snps_in_region <- snps[chr == region_chr & pos >= region$start & pos <= region$end, ]

        if (nrow(snps_in_region) == 0) {
            other_traits_col[i] <- ''
            other_snp_count_col[i] <- 0L
            next
        }

        # Collect traits from these SNPs (excluding the region's own trait)
        cross_traits <- c()
        cross_snp_ids <- c()
        for (m in method_cols) {
            for (j in seq_len(nrow(snps_in_region))) {
                val <- snps_in_region[[m]][j]
                if (!is.na(val) && val != '' && val != '""') {
                    traits_m <- str_replace_all(val, '"', '') %>% str_split(',') %>% unlist
                    traits_m <- traits_m[traits_m != '' & !is.na(traits_m)]
                    other_t <- traits_m[traits_m != region_trait]
                    if (length(other_t) > 0) {
                        cross_traits <- c(cross_traits, other_t)
                        cross_snp_ids <- c(cross_snp_ids, snps_in_region$SNPID[j])
                    }
                }
            }
        }

        other_traits_col[i] <- paste(unique(cross_traits), collapse = ',')
        other_snp_count_col[i] <- length(unique(cross_snp_ids))
    }

    per_trait_regions$other_traits <- other_traits_col
    per_trait_regions$other_snp_count <- other_snp_count_col
    message('INFO: Cross-trait evidence added')
}

# Save per-trait regions
per_trait_regions %>% fwrite(OUTPUT_PER_TRAIT, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(per_trait_regions), ' per-trait regions to ', OUTPUT_PER_TRAIT))

#=============================================================================
# 2. CREATE COMBINED CLIMATE REGIONS
#=============================================================================
message('INFO: Creating combined climate regions...')

combined_regions <- create_regions_from_snps(snps, trait_label = NULL)

# Add traits column showing all traits in each region
if (!is.null(combined_regions) && nrow(combined_regions) > 0) {
    # Compute traits, per-trait counts, and per-method counts for each combined region
    n_regions <- nrow(combined_regions)
    traits_col <- character(n_regions)
    row_data <- vector('list', n_regions)

    for (i in seq_len(n_regions)) {
        snp_ids <- str_split(combined_regions$snp_ids[i], ',')[[1]]
        region_traits <- c()
        trait_snp_map <- list()  # trait -> set of snp_ids
        method_snp_map <- list()  # method -> set of snp_ids

        for (snp_id in snp_ids) {
            snp_row <- snps[SNPID == snp_id, ]
            for (m in method_cols) {
                val <- snp_row[[m]]
                if (!is.na(val) && val != '' && val != '""') {
                    traits_m <- str_replace_all(val, '"', '') %>% str_split(',') %>% unlist
                    traits_m <- traits_m[traits_m != '' & !is.na(traits_m)]
                    region_traits <- c(region_traits, traits_m)
                    for (t in traits_m) {
                        trait_snp_map[[t]] <- c(trait_snp_map[[t]], snp_id)
                    }
                    method_snp_map[[m]] <- c(method_snp_map[[m]], snp_id)
                }
            }
        }

        traits_col[i] <- paste(unique(region_traits[region_traits != '']), collapse = ',')
        row_data[[i]] <- list(
            trait_counts = sapply(trait_snp_map, function(x) length(unique(x))),
            method_counts = sapply(method_snp_map, function(x) length(unique(x)))
        )
    }

    combined_regions$traits <- traits_col

    # Collect all trait and method names across all regions
    all_trait_names <- sort(unique(unlist(lapply(row_data, function(x) names(x$trait_counts)))))
    all_method_names <- sort(unique(unlist(lapply(row_data, function(x) names(x$method_counts)))))

    # Add per-trait count columns (e.g., bio_1_snps, bio_12_snps)
    for (t in all_trait_names) {
        combined_regions[[paste0(t, '_snps')]] <- sapply(row_data, function(x) {
            if (t %in% names(x$trait_counts)) x$trait_counts[[t]] else 0L
        })
    }

    # Add per-method count columns (e.g., EMMAX_snps, LFMM_snps)
    for (m in all_method_names) {
        combined_regions[[paste0(m, '_snps')]] <- sapply(row_data, function(x) {
            if (m %in% names(x$method_counts)) x$method_counts[[m]] else 0L
        })
    }

    message(paste0('INFO: Created ', nrow(combined_regions), ' combined climate regions'))
} else {
    message('WARNING: No combined regions created')
    combined_regions <- data.table(
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
}

# Remove empty trait column from combined (not needed)
combined_regions$trait <- NULL

# Save combined regions
combined_regions %>% fwrite(OUTPUT_COMBINED, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(combined_regions), ' combined regions to ', OUTPUT_COMBINED))

message('INFO: Complete')
