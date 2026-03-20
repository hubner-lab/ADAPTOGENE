#!/usr/bin/env Rscript
# Combine GEA (climate) and GWAS (phenotype) significant SNPs into unified analysis
#
# Inputs: Selected_SNPs from both association types
# Outputs:
#   1. Selected_SNPs_all.tsv - Union of all significant SNPs (climate + phenotype)
#   2. Regions_per_trait_all.tsv - Per-trait regions from combined SNPs
#   3. Regions_all_combined.tsv - All SNPs clustered together regardless of source/trait
#   4. Overlap_summary.tsv - Overlap statistics between extended climate and phenotype regions

library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
options(scipen = 99999)
################
GEA_SELECTED   = args[1]  # association/Selected_SNPs.tsv (climate)
GWAS_SELECTED  = args[2]  # association_phenotypes/Selected_SNPs.tsv (phenotype)
GEA_REGIONS    = args[3]  # association/Regions_climate_combined.tsv
GWAS_REGIONS   = args[4]  # association_phenotypes/Regions_phenotype_combined.tsv
OVERLAP_RDIST_RAW = args[5]  # OVERLAP_REGION_DISTANCE or "auto"
OUT_SELECTED   = args[6]  # Selected_SNPs_all.tsv
OUT_PER_TRAIT  = args[7]  # Regions_per_trait_all.tsv
OUT_COMBINED   = args[8]  # Regions_all_combined.tsv
OUT_OVERLAP    = args[9]  # Overlap_summary.tsv
LD_DECAY_PATH  = if (length(args) >= 10) args[10] else "NULL"

# Resolve "auto" region distance from LD decay table
if (tolower(OVERLAP_RDIST_RAW) == "auto") {
    if (LD_DECAY_PATH == "NULL" || !file.exists(LD_DECAY_PATH)) {
        stop("region_distance='auto' but LD decay table not found: ", LD_DECAY_PATH)
    }
    ld_table <- data.table::fread(LD_DECAY_PATH)
    auto_val <- ld_table[group == "All" & scope == "genome_wide", r2_02_bp]
    if (length(auto_val) == 0 || is.na(auto_val[1])) {
        stop("No genome-wide r2<0.2 distance found in LD decay table for 'All' group")
    }
    OVERLAP_RDIST <- as.numeric(auto_val[1])
    message(paste0('INFO: AUTO overlap_region_distance from LD decay: ', OVERLAP_RDIST, ' bp'))
} else {
    OVERLAP_RDIST <- as.numeric(OVERLAP_RDIST_RAW)
}
################

message('INFO: Combining GEA (climate) and GWAS (phenotype) significant SNPs')

# ─────────────────────────────────────────────────────────────────────────────
# Helper: Create regions from SNPs (same algorithm as create_regions.R)
# ─────────────────────────────────────────────────────────────────────────────
create_regions_from_snps <- function(snps_dt, rdist, trait_label = NULL) {
    if (nrow(snps_dt) == 0) {
        return(data.table(
            region_id = character(), trait = character(),
            chr = character(), start = integer(), end = integer(),
            length = integer(), snp_count = integer(), snp_ids = character(),
            methods = character(), min_pvalue = numeric()))
    }

    snps_gr <- snps_dt %>%
        dplyr::mutate(start = pos, end = pos) %>%
        dplyr::select(chr, start, end, SNPID, min_pvalue) %>%
        GRanges()

    snps_extended <- snps_gr
    start(snps_extended) <- pmax(1, start(snps_extended) - rdist / 2)
    end(snps_extended) <- end(snps_extended) + rdist / 2

    regions_reduced <- reduce(snps_extended)

    regions_list <- lapply(seq_along(regions_reduced), function(i) {
        region <- regions_reduced[i]
        overlaps <- findOverlaps(snps_extended, region)
        snp_indices <- queryHits(overlaps)
        if (length(snp_indices) == 0) return(NULL)

        region_snps <- snps_dt[snp_indices, ]
        region_chr <- as.character(seqnames(region))
        # Region boundaries: extend SNP range by rdist on each side
        region_start <- as.integer(pmax(1, min(region_snps$pos) - rdist))
        region_end <- as.integer(max(region_snps$pos) + rdist)

        # Extract methods from method-specific columns
        method_cols <- setdiff(colnames(region_snps), c('SNPID', 'chr', 'pos', 'min_pvalue', 'source'))
        active_methods <- c()
        for (m in method_cols) {
            vals <- region_snps[[m]]
            if (any(!is.na(vals) & vals != '' & vals != '""')) {
                active_methods <- c(active_methods, m)
            }
        }
        methods_str <- paste(unique(active_methods), collapse = ',')

        # Region ID uses extended boundaries
        region_id <- if (!is.null(trait_label)) {
            paste0(region_chr, '_', region_start, '-', region_end, '_', trait_label)
        } else {
            paste0(region_chr, '_', region_start, '-', region_end)
        }

        data.table(
            region_id = region_id, trait = ifelse(!is.null(trait_label), trait_label, ''),
            chr = region_chr, start = region_start, end = region_end,
            length = region_end - region_start, snp_count = nrow(region_snps),
            snp_ids = paste(region_snps$SNPID, collapse = ','),
            methods = methods_str, min_pvalue = min(region_snps$min_pvalue))
    })

    regions <- rbindlist(regions_list)
    if (nrow(regions) > 0) regions <- regions[order(chr, start)]
    return(regions)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load and merge Selected SNPs
# ─────────────────────────────────────────────────────────────────────────────
gea_snps <- fread(GEA_SELECTED)
gwas_snps <- fread(GWAS_SELECTED)

if ('chr' %in% names(gea_snps)) gea_snps$chr <- as.character(gea_snps$chr)
if ('chr' %in% names(gwas_snps)) gwas_snps$chr <- as.character(gwas_snps$chr)

message(paste0('INFO: GEA SNPs: ', nrow(gea_snps), ', GWAS SNPs: ', nrow(gwas_snps)))

# Get method columns from each source
gea_methods <- setdiff(names(gea_snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
gwas_methods <- setdiff(names(gwas_snps), c('SNPID', 'chr', 'pos', 'min_pvalue'))
all_method_cols <- unique(c(gea_methods, gwas_methods))

# Add source tag to track origin
gea_snps[, source := 'climate']
gwas_snps[, source := 'phenotype']

# Ensure both tables have all method columns
for (m in all_method_cols) {
    if (!m %in% names(gea_snps)) gea_snps[[m]] <- ""
    if (!m %in% names(gwas_snps)) gwas_snps[[m]] <- ""
}

# Select common columns for binding
common_cols <- c('SNPID', 'chr', 'pos', all_method_cols, 'min_pvalue', 'source')
merged <- rbind(gea_snps[, ..common_cols], gwas_snps[, ..common_cols])

# Collapse duplicates (same SNP significant in both GEA and GWAS)
if (any(duplicated(merged$SNPID))) {
    merged <- merged[, lapply(.SD, function(x) {
        if (is.numeric(x)) return(min(x, na.rm = TRUE))
        vals <- x[x != '' & x != '""' & !is.na(x)]
        if (length(vals) == 0) return('')
        paste(unique(unlist(str_split(vals, ','))), collapse = ',')
    }), by = .(SNPID, chr, pos), .SDcols = c(all_method_cols, 'min_pvalue', 'source')]
}

message(paste0('INFO: Merged unique SNPs: ', nrow(merged)))

# Save merged Selected_SNPs (without source column for compatibility)
out_snps <- copy(merged)
out_snps[, source := NULL]
fwrite(out_snps, OUT_SELECTED, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(out_snps), ' merged SNPs to ', OUT_SELECTED))

# ─────────────────────────────────────────────────────────────────────────────
# 2. Extract all traits and create per-trait + combined regions
# ─────────────────────────────────────────────────────────────────────────────

# Extract unique traits from method columns
all_traits <- c()
for (m in all_method_cols) {
    traits_m <- merged[[m]] %>%
        str_replace_all('"', '') %>%
        str_split(',') %>%
        unlist() %>%
        unique()
    traits_m <- traits_m[traits_m != '' & !is.na(traits_m)]
    all_traits <- c(all_traits, traits_m)
}
all_traits <- unique(all_traits)
message(paste0('INFO: All traits: ', paste(all_traits, collapse = ', ')))

# --- Per-trait regions ---
message('INFO: Creating per-trait regions from all SNPs...')
per_trait_list <- lapply(all_traits, function(trait) {
    trait_snps <- merged[0, ]
    for (m in all_method_cols) {
        snps_with_trait <- merged[grepl(paste0('(^|,)', trait, '($|,)'), merged[[m]]), ]
        if (nrow(snps_with_trait) > 0) trait_snps <- rbind(trait_snps, snps_with_trait)
    }
    trait_snps <- unique(trait_snps)
    if (nrow(trait_snps) == 0) return(NULL)
    message(paste0('INFO:   ', trait, ': ', nrow(trait_snps), ' SNPs'))
    create_regions_from_snps(trait_snps, OVERLAP_RDIST, trait_label = trait)
})
per_trait_regions <- rbindlist(per_trait_list)

# Add cross-trait evidence + source info
if (nrow(per_trait_regions) > 0) {
    other_traits_col <- character(nrow(per_trait_regions))
    other_snp_count_col <- integer(nrow(per_trait_regions))
    source_col <- character(nrow(per_trait_regions))

    for (i in seq_len(nrow(per_trait_regions))) {
        region <- per_trait_regions[i, ]
        # Use pre-extended region boundaries as-is
        snps_in_region <- merged[chr == region$chr &
                                 pos >= region$start &
                                 pos <= region$end, ]
        # Cross-trait evidence
        cross_traits <- c()
        for (m in all_method_cols) {
            for (j in seq_len(nrow(snps_in_region))) {
                val <- snps_in_region[[m]][j]
                if (!is.na(val) && val != '' && val != '""') {
                    traits_m <- str_replace_all(val, '"', '') %>% str_split(',') %>% unlist
                    cross_traits <- c(cross_traits, traits_m[traits_m != '' & traits_m != region$trait])
                }
            }
        }
        other_traits_col[i] <- paste(unique(cross_traits), collapse = ',')
        other_snp_count_col[i] <- length(unique(snps_in_region$SNPID)) - region$snp_count

        # Source: which analysis types contribute to this region
        region_snp_ids <- str_split(region$snp_ids, ',')[[1]]
        sources <- unique(merged[SNPID %in% region_snp_ids, source])
        source_col[i] <- paste(sort(sources), collapse = ',')
    }
    per_trait_regions$other_traits <- other_traits_col
    per_trait_regions$other_snp_count <- pmax(0L, other_snp_count_col)
    per_trait_regions$source <- source_col
}

fwrite(per_trait_regions, OUT_PER_TRAIT, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(per_trait_regions), ' per-trait regions to ', OUT_PER_TRAIT))

# --- Combined regions ---
message('INFO: Creating combined regions from all SNPs...')
combined_regions <- create_regions_from_snps(merged, OVERLAP_RDIST, trait_label = NULL)

if (nrow(combined_regions) > 0) {
    n_regions <- nrow(combined_regions)
    traits_col <- character(n_regions)
    source_col <- character(n_regions)
    row_data <- vector('list', n_regions)

    for (i in seq_len(n_regions)) {
        snp_ids <- str_split(combined_regions$snp_ids[i], ',')[[1]]
        region_traits <- c()
        trait_snp_map <- list()
        method_snp_map <- list()
        region_sources <- c()

        for (snp_id in snp_ids) {
            snp_row <- merged[SNPID == snp_id, ]
            if (nrow(snp_row) == 0) next
            region_sources <- c(region_sources, str_split(snp_row$source, ',')[[1]])
            for (m in all_method_cols) {
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
        source_col[i] <- paste(sort(unique(region_sources)), collapse = ',')
        row_data[[i]] <- list(
            trait_counts = sapply(trait_snp_map, function(x) length(unique(x))),
            method_counts = sapply(method_snp_map, function(x) length(unique(x))))
    }

    combined_regions$traits <- traits_col
    combined_regions$source <- source_col

    # Add per-trait and per-method count columns
    all_trait_names <- sort(unique(unlist(lapply(row_data, function(x) names(x$trait_counts)))))
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

    # Count climate vs phenotype SNPs per region
    # Climate traits start with "bio_", phenotype traits are the rest
    for (i in seq_len(nrow(combined_regions))) {
        snp_ids <- str_split(combined_regions$snp_ids[i], ',')[[1]]
        n_climate <- sum(merged[SNPID %in% snp_ids, grepl('climate', source)])
        n_phenotype <- sum(merged[SNPID %in% snp_ids, grepl('phenotype', source)])
        combined_regions[i, climate_snps := n_climate]
        combined_regions[i, phenotype_snps := n_phenotype]
    }
}

# Remove empty trait column
combined_regions$trait <- NULL

fwrite(combined_regions, OUT_COMBINED, sep = '\t', quote = FALSE)
message(paste0('INFO: Saved ', nrow(combined_regions), ' combined regions to ', OUT_COMBINED))

# ─────────────────────────────────────────────────────────────────────────────
# 3. Compute overlap between GEA and GWAS extended regions
# ─────────────────────────────────────────────────────────────────────────────
message('INFO: Computing overlap between GEA and GWAS regions...')

gea_regions <- fread(GEA_REGIONS)
gwas_regions <- fread(GWAS_REGIONS)
if ('chr' %in% names(gea_regions)) gea_regions$chr <- as.character(gea_regions$chr)
if ('chr' %in% names(gwas_regions)) gwas_regions$chr <- as.character(gwas_regions$chr)

overlap_results <- list()
gea_matched <- character(0)
gwas_matched <- character(0)

if (nrow(gea_regions) > 0 && nrow(gwas_regions) > 0) {
    # Regions already have extended boundaries (REGION_DISTANCE applied in create_regions.R)
    # Use start/end directly for overlap detection
    for (gi in seq_len(nrow(gea_regions))) {
        g <- gea_regions[gi]
        for (wi in seq_len(nrow(gwas_regions))) {
            w <- gwas_regions[wi]
            if (as.character(g$chr) != as.character(w$chr)) next
            if (g$start > w$end || g$end < w$start) next

            ov_start <- max(g$start, w$start)
            ov_end <- min(g$end, w$end)

            # Calculate overlap percentage relative to smaller region
            gea_len <- g$end - g$start
            gwas_len <- w$end - w$start
            ov_len <- ov_end - ov_start
            min_len <- min(gea_len, gwas_len)
            overlap_pct <- if (min_len > 0) round(100 * ov_len / min_len, 1) else 100

            # Count SNPs from each source in overlap area
            gea_snp_ids <- str_split(g$snp_ids, ',')[[1]]
            gwas_snp_ids <- str_split(w$snp_ids, ',')[[1]]
            shared <- intersect(gea_snp_ids, gwas_snp_ids)

            gea_traits_str <- if ('traits' %in% names(g)) as.character(g$traits) else ""
            gwas_traits_str <- if ('traits' %in% names(w)) as.character(w$traits) else ""

            overlap_results[[length(overlap_results) + 1]] <- data.table(
                gea_region_id = g$region_id,
                gwas_region_id = w$region_id,
                chr = as.character(g$chr),
                overlap_start = ov_start,
                overlap_end = ov_end,
                overlap_length = ov_len,
                overlap_pct = overlap_pct,
                gea_snp_count = g$snp_count,
                gwas_snp_count = w$snp_count,
                shared_snps = length(shared),
                gea_traits = gea_traits_str,
                gwas_traits = gwas_traits_str,
                gea_min_pvalue = g$min_pvalue,
                gwas_min_pvalue = w$min_pvalue
            )
            gea_matched <- c(gea_matched, g$region_id)
            gwas_matched <- c(gwas_matched, w$region_id)
        }
    }
}

if (length(overlap_results) > 0) {
    overlap_dt <- rbindlist(overlap_results)
} else {
    overlap_dt <- data.table(
        gea_region_id = character(), gwas_region_id = character(),
        chr = character(), overlap_start = integer(), overlap_end = integer(),
        overlap_length = integer(), overlap_pct = numeric(),
        gea_snp_count = integer(), gwas_snp_count = integer(),
        shared_snps = integer(), gea_traits = character(), gwas_traits = character(),
        gea_min_pvalue = numeric(), gwas_min_pvalue = numeric())
}

# Add non-overlapping region summaries
n_gea <- nrow(gea_regions)
n_gwas <- nrow(gwas_regions)
n_gea_overlap <- length(unique(gea_matched))
n_gwas_overlap <- length(unique(gwas_matched))

# Append summary rows
summary_rows <- data.table(
    gea_region_id = c('SUMMARY', 'SUMMARY', 'SUMMARY', 'SUMMARY', 'SUMMARY'),
    gwas_region_id = c('gea_total', 'gwas_total', 'gea_overlapping', 'gwas_overlapping', 'overlap_pairs'),
    chr = '', overlap_start = 0L, overlap_end = 0L,
    overlap_length = c(n_gea, n_gwas, n_gea_overlap, n_gwas_overlap, nrow(overlap_dt)),
    overlap_pct = c(
        if (n_gea > 0) round(100 * n_gea_overlap / n_gea, 1) else 0,
        if (n_gwas > 0) round(100 * n_gwas_overlap / n_gwas, 1) else 0,
        0, 0, 0),
    gea_snp_count = 0L, gwas_snp_count = 0L, shared_snps = 0L,
    gea_traits = '', gwas_traits = '',
    gea_min_pvalue = NA_real_, gwas_min_pvalue = NA_real_)

overlap_out <- rbind(overlap_dt, summary_rows)
fwrite(overlap_out, OUT_OVERLAP, sep = '\t', quote = FALSE)

message(paste0('INFO: GEA regions: ', n_gea, ', GWAS regions: ', n_gwas))
message(paste0('INFO: Overlapping pairs: ', nrow(overlap_dt)))
message(paste0('INFO: GEA overlapping: ', n_gea_overlap, '/', n_gea,
               ' (', if (n_gea > 0) round(100 * n_gea_overlap / n_gea, 1) else 0, '%)'))
message(paste0('INFO: GWAS overlapping: ', n_gwas_overlap, '/', n_gwas,
               ' (', if (n_gwas > 0) round(100 * n_gwas_overlap / n_gwas, 1) else 0, '%)'))
message(paste0('INFO: New combined regions: ', nrow(combined_regions)))
message('INFO: Complete')
