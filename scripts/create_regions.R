#!/usr/bin/env Rscript
# Merge nearby significant SNPs into regions using single-linkage clustering
#
# Creates TWO types of regions:
# 1. Per-trait regions: SNPs clustered separately for each trait (bio_1, bio_12, etc.)
# 2. Combined climate regions: All SNPs clustered together regardless of trait
#
# REGION_DISTANCE defines the maximum gap between neighboring SNPs to be in the
# same region. Can be a scalar (all chromosomes) or a named numeric vector
# (per-chromosome, keyed by chr name).
#
# Accepted region_distance_spec values:
#   auto_per_chromosome — invert Hill-Weir curve per chromosome at r2_threshold
#   auto_genome_wide    — invert Hill-Weir curve genome-wide at r2_threshold
#   auto                — alias for auto_genome_wide (backward compat)
#   <integer>           — fixed bp distance

library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
options(scipen = 99999)
################
SELECTED_SNPS        = args[1]
REGION_DISTANCE_SPEC = args[2]   # "auto_per_chromosome" | "auto_genome_wide" | integer
OUTPUT_PER_TRAIT     = args[3]
OUTPUT_COMBINED      = args[4]
LD_DECAY_PATH        = if (length(args) >= 5) args[5] else "NULL"
R2_THRESHOLD         = if (length(args) >= 6) as.numeric(args[6]) else 0.2
LD_DECAY_GROUP       = if (length(args) >= 7) args[7] else "All"
################

# Hill-Weir (1988) expected r^2 at distance d given parameter C and sample size n
hill_weir <- function(d, C, n) {
    x <- C * d
    ((10 + x) / ((2 + x) * (11 + x))) *
    (1 + ((3 + x) * (12 + 12 * x + x^2)) / (n * (2 + x) * (11 + x)))
}

# Invert Hill-Weir: find distance where r^2 = r2_target given C_hat, n, max_dist_bp
invert_hill_weir <- function(C_hat, n_samples, r2_target, max_dist_bp) {
    if (length(C_hat) == 0 || length(n_samples) == 0 || length(max_dist_bp) == 0 ||
        is.na(C_hat) || is.na(n_samples) || is.na(max_dist_bp)) return(NA_real_)
    tryCatch(
        uniroot(function(d) hill_weir(d, C_hat, n_samples) - r2_target,
                c(1, max_dist_bp))$root,
        error = function(e) {
            message(paste0("WARNING: Hill-Weir inversion failed (r2=", r2_target,
                           "): ", e$message))
            NA_real_
        }
    )
}

# For LOESS/failed rows: fall back to nearest precomputed column
loess_fallback <- function(row, r2_target) {
    half_val <- as.numeric(row[["half_decay_bp"]])
    r2_02_val <- as.numeric(row[["r2_02_bp"]])
    r2_intercept <- as.numeric(row[["r2_intercept"]])
    half_r2 <- r2_intercept / 2
    # nearest breakpoint: if target is closer to half-decay (r2_intercept/2) use that; else r2=0.2
    if (!is.na(half_val) && !is.na(r2_02_val)) {
        if (abs(r2_target - half_r2) < abs(r2_target - 0.2)) {
            message(paste0("WARNING: Using precomputed half_decay_bp (r2=", round(half_r2, 3),
                           ") as fallback for r2=", r2_target))
            return(half_val)
        } else {
            message(paste0("WARNING: Using precomputed r2_02_bp as fallback for r2=", r2_target))
            return(r2_02_val)
        }
    } else if (!is.na(r2_02_val)) {
        return(r2_02_val)
    } else if (!is.na(half_val)) {
        return(half_val)
    }
    return(NA_real_)
}

# Resolve distance for one LD-decay TSV row
resolve_row_distance <- function(row, r2_target) {
    method <- row[["method"]]
    if (method == "hill_weir") {
        C_hat      <- as.numeric(row[["C_hat"]])
        n_samples  <- as.numeric(row[["n_samples"]])
        max_dist   <- as.numeric(row[["max_dist_bp"]])
        # C_hat absent (old TSV without the column) — fall back to precomputed columns
        if (length(C_hat) == 0 || length(max_dist) == 0) {
            message("WARNING: C_hat/max_dist_bp missing from LD decay table — using precomputed fallback. Re-run mode=structure to regenerate.")
            return(loess_fallback(row, r2_target))
        }
        dist <- invert_hill_weir(C_hat, n_samples, r2_target, max_dist)
        if (!is.na(dist)) return(round(dist))
        # hill_weir inversion failed — try precomputed
        return(loess_fallback(row, r2_target))
    } else {
        return(loess_fallback(row, r2_target))
    }
}

#=============================================================================
# RESOLVE REGION_DISTANCE
#=============================================================================
spec_lower <- tolower(REGION_DISTANCE_SPEC)
auto_mode <- spec_lower %in% c("auto_per_chromosome", "auto_genome_wide", "auto")

if (auto_mode) {
    if (LD_DECAY_PATH == "NULL" || !file.exists(LD_DECAY_PATH)) {
        stop("region_distance='", REGION_DISTANCE_SPEC, "' but LD decay table not found: ", LD_DECAY_PATH)
    }
    ld_table <- fread(LD_DECAY_PATH)

    if (spec_lower %in% c("auto_genome_wide", "auto")) {
        if (spec_lower == "auto") {
            message("INFO: 'auto' region_distance treated as 'auto_genome_wide' (deprecated alias)")
        }
        row <- ld_table[group == LD_DECAY_GROUP & scope == "genome_wide"]
        if (nrow(row) == 0) stop("No genome-wide LD decay row found for group='", LD_DECAY_GROUP, "'")
        dist_val <- resolve_row_distance(as.list(row[1]), R2_THRESHOLD)
        if (is.na(dist_val)) stop("Could not resolve genome-wide region_distance from LD decay table")
        REGION_DISTANCE <- dist_val
        message(paste0("INFO: AUTO genome-wide region_distance (r2<", R2_THRESHOLD, ", group=",
                       LD_DECAY_GROUP, "): ", REGION_DISTANCE, " bp"))

    } else {
        # auto_per_chromosome: build named vector, one entry per chr
        chr_rows <- ld_table[group == LD_DECAY_GROUP & !scope %in% c("genome_wide")]
        gw_row   <- ld_table[group == LD_DECAY_GROUP & scope == "genome_wide"]
        if (nrow(gw_row) == 0 && nrow(chr_rows) == 0)
            stop("No LD decay rows found for group='", LD_DECAY_GROUP, "'")

        # genome-wide fallback distance
        gw_dist <- if (nrow(gw_row) > 0) {
            d <- resolve_row_distance(as.list(gw_row[1]), R2_THRESHOLD)
            if (!is.na(d)) d else NULL
        } else NULL

        dist_map <- setNames(numeric(0), character(0))
        for (i in seq_len(nrow(chr_rows))) {
            row <- as.list(chr_rows[i])
            chr_name <- row[["scope"]]
            d <- resolve_row_distance(row, R2_THRESHOLD)
            if (is.na(d)) {
                if (!is.null(gw_dist)) {
                    message(paste0("WARNING: chr ", chr_name, " r2 inversion failed; using genome-wide fallback: ", gw_dist, " bp"))
                    d <- gw_dist
                } else {
                    stop("r2 inversion failed for chr ", chr_name, " and no genome-wide fallback available")
                }
            }
            dist_map[chr_name] <- round(d)
        }
        REGION_DISTANCE <- dist_map
        message(paste0("INFO: AUTO per-chromosome region_distance (r2<", R2_THRESHOLD, ", group=",
                       LD_DECAY_GROUP, "):"))
        for (chr in names(REGION_DISTANCE)) {
            message(paste0("  chr ", chr, ": ", REGION_DISTANCE[chr], " bp"))
        }
        if (!is.null(gw_dist)) {
            message(paste0("  genome-wide fallback: ", round(gw_dist), " bp"))
        }
    }
} else {
    REGION_DISTANCE <- as.numeric(REGION_DISTANCE_SPEC)
    message(paste0("INFO: Fixed region_distance: ", REGION_DISTANCE, " bp"))
}

message('INFO: Creating regions from selected SNPs using single-linkage clustering')
message('INFO: Producing two outputs: per-trait regions and combined climate regions')

######################## Functions

# Get the region_distance for a given chromosome (scalar or named-vector)
get_dist <- function(chr_name) {
    if (is.numeric(REGION_DISTANCE) && length(REGION_DISTANCE) == 1 && is.null(names(REGION_DISTANCE))) {
        return(REGION_DISTANCE)
    }
    if (chr_name %in% names(REGION_DISTANCE)) {
        return(REGION_DISTANCE[[chr_name]])
    }
    # chromosome not in map — use genome-wide if we're in auto_per_chromosome mode
    # (we store gw_dist in the environment above only in auto_per_chr; fall back to max)
    fallback <- max(REGION_DISTANCE, na.rm = TRUE)
    message(paste0("WARNING: chr '", chr_name, "' not in per-chr distance map; using max = ", fallback, " bp"))
    return(fallback)
}

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

    # Process per chromosome so per-chr distances apply cleanly
    chrs <- unique(snps_dt$chr)
    all_regions <- lapply(chrs, function(this_chr) {
        dist <- get_dist(as.character(this_chr))
        chr_snps <- snps_dt[chr == this_chr]

        snps_gr <- chr_snps %>%
            dplyr::mutate(start = pos, end = pos) %>%
            dplyr::select(chr, start, end, SNPID, min_pvalue) %>%
            GRanges()

        # Extend by dist/2 on each side for single-linkage clustering
        snps_ext <- snps_gr
        start(snps_ext) <- pmax(1L, start(snps_ext) - as.integer(dist / 2))
        end(snps_ext)   <- end(snps_ext) + as.integer(dist / 2)

        regions_reduced <- reduce(snps_ext)

        lapply(seq_along(regions_reduced), function(i) {
            region <- regions_reduced[i]
            overlaps <- findOverlaps(snps_ext, region)
            snp_indices <- queryHits(overlaps)
            if (length(snp_indices) == 0) return(NULL)

            region_snps <- chr_snps[snp_indices, ]
            region_chr  <- as.character(seqnames(region))
            region_start <- as.integer(pmax(1L, min(region_snps$pos) - as.integer(dist)))
            region_end   <- as.integer(max(region_snps$pos) + as.integer(dist))

            method_cols <- setdiff(colnames(region_snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
            active_methods <- c()
            for (m in method_cols) {
                vals <- region_snps[[m]]
                if (any(!is.na(vals) & vals != '' & vals != '""')) {
                    active_methods <- c(active_methods, m)
                }
            }
            methods_str <- paste(unique(active_methods), collapse = ',')

            region_id <- if (!is.null(trait_label)) {
                paste0(region_chr, '_', region_start, '-', region_end, '_', trait_label)
            } else {
                paste0(region_chr, '_', region_start, '-', region_end)
            }

            data.table(
                region_id  = region_id,
                trait      = ifelse(!is.null(trait_label), trait_label, ''),
                chr        = region_chr,
                start      = region_start,
                end        = region_end,
                length     = region_end - region_start,
                snp_count  = nrow(region_snps),
                snp_ids    = paste(region_snps$SNPID, collapse = ','),
                methods    = methods_str,
                min_pvalue = min(region_snps$min_pvalue)
            )
        })
    })

    regions <- rbindlist(Filter(Negate(is.null), unlist(all_regions, recursive = FALSE)))
    if (nrow(regions) > 0) regions <- regions[order(chr, start)]
    return(regions)
}

#################### Main

# Load selected SNPs
snps <- fread(SELECTED_SNPS)
if ('chr' %in% colnames(snps)) snps$chr <- as.character(snps$chr)

if (nrow(snps) == 0) {
    message('WARNING: No SNPs to process')
    empty_dt <- data.table(
        region_id = character(), trait = character(), chr = character(),
        start = integer(), end = integer(), length = integer(),
        snp_count = integer(), snp_ids = character(),
        methods = character(), min_pvalue = numeric()
    )
    fwrite(empty_dt, OUTPUT_PER_TRAIT, sep = '\t')
    fwrite(empty_dt, OUTPUT_COMBINED,  sep = '\t')
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

    trait_snps <- snps[0, ]
    for (m in method_cols) {
        snps_with_trait <- snps[grepl(paste0('(^|,)', trait, '($|,)'), snps[[m]]), ]
        if (nrow(snps_with_trait) > 0) trait_snps <- rbind(trait_snps, snps_with_trait)
    }
    trait_snps <- unique(trait_snps)

    if (nrow(trait_snps) == 0) {
        message(paste0('INFO:     No SNPs found for ', trait))
        return(NULL)
    }

    message(paste0('INFO:     Found ', nrow(trait_snps), ' SNPs for ', trait))
    regions <- create_regions_from_snps(trait_snps, trait_label = trait)
    if (!is.null(regions) && nrow(regions) > 0)
        message(paste0('INFO:     Created ', nrow(regions), ' regions for ', trait))
    return(regions)
})

per_trait_regions <- rbindlist(Filter(Negate(is.null), per_trait_regions_list))

if (is.null(per_trait_regions) || nrow(per_trait_regions) == 0) {
    message('WARNING: No per-trait regions created')
    per_trait_regions <- data.table(
        region_id = character(), trait = character(), chr = character(),
        start = integer(), end = integer(), length = integer(),
        snp_count = integer(), snp_ids = character(),
        methods = character(), min_pvalue = numeric()
    )
}

# Cross-trait evidence: for each region, find SNPs from OTHER traits within its bounds
if (nrow(per_trait_regions) > 0) {
    message('INFO: Computing cross-trait evidence for per-trait regions...')
    other_traits_col    <- character(nrow(per_trait_regions))
    other_snp_count_col <- integer(nrow(per_trait_regions))

    for (i in seq_len(nrow(per_trait_regions))) {
        region       <- per_trait_regions[i, ]
        region_trait <- region$trait
        region_chr   <- region$chr
        snps_in_region <- snps[chr == region_chr & pos >= region$start & pos <= region$end, ]

        if (nrow(snps_in_region) == 0) {
            other_traits_col[i]    <- ''
            other_snp_count_col[i] <- 0L
            next
        }

        cross_traits  <- c()
        cross_snp_ids <- c()
        for (m in method_cols) {
            for (j in seq_len(nrow(snps_in_region))) {
                val <- snps_in_region[[m]][j]
                if (!is.na(val) && val != '' && val != '""') {
                    traits_m <- str_replace_all(val, '"', '') %>% str_split(',') %>% unlist
                    traits_m <- traits_m[traits_m != '' & !is.na(traits_m)]
                    other_t  <- traits_m[traits_m != region_trait]
                    if (length(other_t) > 0) {
                        cross_traits  <- c(cross_traits, other_t)
                        cross_snp_ids <- c(cross_snp_ids, snps_in_region$SNPID[j])
                    }
                }
            }
        }
        other_traits_col[i]    <- paste(unique(cross_traits), collapse = ',')
        other_snp_count_col[i] <- length(unique(cross_snp_ids))
    }
    per_trait_regions$other_traits    <- other_traits_col
    per_trait_regions$other_snp_count <- other_snp_count_col
    message('INFO: Cross-trait evidence added')
}

fwrite(per_trait_regions, OUTPUT_PER_TRAIT, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(per_trait_regions), ' per-trait regions to ', OUTPUT_PER_TRAIT))

#=============================================================================
# 2. CREATE COMBINED CLIMATE REGIONS
#=============================================================================
message('INFO: Creating combined climate regions...')

combined_regions <- create_regions_from_snps(snps, trait_label = NULL)

if (!is.null(combined_regions) && nrow(combined_regions) > 0) {
    n_regions  <- nrow(combined_regions)
    traits_col <- character(n_regions)
    row_data   <- vector('list', n_regions)

    for (i in seq_len(n_regions)) {
        snp_ids     <- str_split(combined_regions$snp_ids[i], ',')[[1]]
        region_traits   <- c()
        trait_snp_map   <- list()
        method_snp_map  <- list()

        for (snp_id in snp_ids) {
            snp_row <- snps[SNPID == snp_id, ]
            for (m in method_cols) {
                val <- snp_row[[m]]
                if (!is.na(val) && val != '' && val != '""') {
                    traits_m <- str_replace_all(val, '"', '') %>% str_split(',') %>% unlist
                    traits_m <- traits_m[traits_m != '' & !is.na(traits_m)]
                    region_traits <- c(region_traits, traits_m)
                    for (t in traits_m) trait_snp_map[[t]] <- c(trait_snp_map[[t]], snp_id)
                    method_snp_map[[m]] <- c(method_snp_map[[m]], snp_id)
                }
            }
        }

        traits_col[i] <- paste(unique(region_traits[region_traits != '']), collapse = ',')
        row_data[[i]] <- list(
            trait_counts  = sapply(trait_snp_map,  function(x) length(unique(x))),
            method_counts = sapply(method_snp_map, function(x) length(unique(x)))
        )
    }

    combined_regions$traits <- traits_col

    all_trait_names  <- sort(unique(unlist(lapply(row_data, function(x) names(x$trait_counts)))))
    all_method_names <- sort(unique(unlist(lapply(row_data, function(x) names(x$method_counts)))))

    for (t in all_trait_names) {
        combined_regions[[paste0(t, '_snps')]] <- sapply(row_data, function(x) {
            if (t %in% names(x$trait_counts)) x$trait_counts[[t]] else 0L
        })
    }
    for (m in all_method_names) {
        combined_regions[[paste0(m, '_snps')]] <- sapply(row_data, function(x) {
            if (m %in% names(x$method_counts)) x$method_counts[[m]] else 0L
        })
    }

    message(paste0('INFO: Created ', nrow(combined_regions), ' combined climate regions'))
} else {
    message('WARNING: No combined regions created')
    combined_regions <- data.table(
        region_id = character(), chr = character(),
        start = integer(), end = integer(), length = integer(),
        snp_count = integer(), snp_ids = character(),
        traits = character(), methods = character(), min_pvalue = numeric()
    )
}

combined_regions$trait <- NULL
fwrite(combined_regions, OUTPUT_COMBINED, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(combined_regions), ' combined regions to ', OUTPUT_COMBINED))

message('INFO: Complete')
