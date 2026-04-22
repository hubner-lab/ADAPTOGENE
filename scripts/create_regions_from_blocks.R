#!/usr/bin/env Rscript
# Reformat significant LD blocks as per-trait and combined region tables.
# Output schema is identical to create_regions.R, so all downstream scripts
# (find_genes_around_regions.R, enrichment, regionplot, haplotype) are unchanged.
#
# Region IDs follow the same convention: CHR_START-END_TRAIT (per-trait) and
# CHR_START-END (combined), using the block's own boundaries directly.

library(data.table)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
options(scipen = 99999)

SELECTED_BLOCKS <- args[1]  # selected_blocks.tsv
SELECTED_SNPS   <- args[2]  # selected_snps.tsv (within-block SNPs)
OUTPUT_PER_TRAIT <- args[3]
OUTPUT_COMBINED  <- args[4]

message("INFO: Creating regions from LD blocks")

blocks <- fread(SELECTED_BLOCKS)
blocks$chr <- as.character(blocks$chr)

snps <- fread(SELECTED_SNPS)
snps$chr <- as.character(snps$chr)

empty_per_trait <- data.table(
    region_id = character(), trait = character(),
    chr = character(), start = integer(), end = integer(),
    length = integer(), snp_count = integer(), snp_ids = character(),
    methods = character(), min_pvalue = numeric(),
    other_traits = character(), other_snp_count = integer()
)
empty_combined <- data.table(
    region_id = character(),
    chr = character(), start = integer(), end = integer(),
    length = integer(), snp_count = integer(), snp_ids = character(),
    traits = character(), methods = character(), min_pvalue = numeric()
)

if (nrow(blocks) == 0) {
    message("WARNING: No selected blocks")
    fwrite(empty_per_trait, OUTPUT_PER_TRAIT, sep = "\t", quote = FALSE)
    fwrite(empty_combined,  OUTPUT_COMBINED,  sep = "\t", quote = FALSE)
    quit(save = "no", status = 0)
}

# Method columns in selected_snps (SNPID, chr, pos, <method cols>, min_pvalue)
fixed_snp_cols <- c("SNPID", "chr", "pos", "min_pvalue")
method_cols <- setdiff(names(snps), fixed_snp_cols)

# SNP IDs within each block (by position overlap)
snps_in_block <- function(blk_chr, blk_start, blk_end) {
    snps[chr == blk_chr & pos >= blk_start & pos <= blk_end]
}

#=============================================================================
# 1. PER-TRAIT REGIONS
#=============================================================================
message("INFO: Building per-trait regions")

# Expand blocks long by trait
blocks_long <- str_split(blocks$traits, ",") %>%
    mapply(function(trts, i) {
        data.table(block_idx = i, trait = trimws(trts))
    }, ., seq_len(nrow(blocks)), SIMPLIFY = FALSE) %>%
    rbindlist()

if (nrow(blocks_long) == 0) {
    fwrite(empty_per_trait, OUTPUT_PER_TRAIT, sep = "\t", quote = FALSE)
} else {
    per_trait_rows <- lapply(seq_len(nrow(blocks_long)), function(i) {
        bi    <- blocks_long$block_idx[i]
        trait <- blocks_long$trait[i]
        blk   <- blocks[bi]

        snp_dt <- snps_in_block(blk$chr, blk$start, blk$end)
        # SNPs in selected_snps must carry this trait in at least one method column
        snp_has_trait <- rowSums(sapply(method_cols, function(m) {
            grepl(paste0("(^|,)", trait, "($|,)"), snp_dt[[m]], fixed = FALSE)
        })) > 0
        trait_snps <- snp_dt[snp_has_trait]

        # Other traits in this block (cross-trait evidence)
        other_trts <- setdiff(str_split(blk$traits, ",")[[1]], trait)
        other_snp_ids <- snps_in_block(blk$chr, blk$start, blk$end)$SNPID
        other_snp_ids <- setdiff(other_snp_ids, trait_snps$SNPID)

        region_id <- paste0(blk$chr, "_", blk$start, "-", blk$end, "_", trait)

        data.table(
            region_id       = region_id,
            trait           = trait,
            chr             = blk$chr,
            start           = blk$start,
            end             = blk$end,
            length          = blk$end - blk$start,
            snp_count       = nrow(trait_snps),
            snp_ids         = paste(trait_snps$SNPID, collapse = ","),
            methods         = blk$methods,
            min_pvalue      = blk$min_pvalue,
            other_traits    = paste(other_trts, collapse = ","),
            other_snp_count = length(other_snp_ids)
        )
    }) %>% rbindlist()

    per_trait_rows <- per_trait_rows %>% dplyr::arrange(chr, start)
    fwrite(per_trait_rows, OUTPUT_PER_TRAIT, sep = "\t", quote = FALSE)
    message(paste0("INFO: Saved ", nrow(per_trait_rows), " per-trait region rows"))
}

#=============================================================================
# 2. COMBINED REGIONS
#=============================================================================
message("INFO: Building combined regions")

all_trait_names   <- sort(unique(unlist(str_split(blocks$traits, ","))))
all_method_names  <- sort(unique(unlist(str_split(blocks$methods, ","))))

combined_rows <- lapply(seq_len(nrow(blocks)), function(i) {
    blk    <- blocks[i]
    snp_dt <- snps_in_block(blk$chr, blk$start, blk$end)
    region_id <- paste0(blk$chr, "_", blk$start, "-", blk$end)

    row <- data.table(
        region_id  = region_id,
        chr        = blk$chr,
        start      = blk$start,
        end        = blk$end,
        length     = blk$end - blk$start,
        snp_count  = nrow(snp_dt),
        snp_ids    = paste(snp_dt$SNPID, collapse = ","),
        traits     = blk$traits,
        methods    = blk$methods,
        min_pvalue = blk$min_pvalue
    )

    # Per-trait SNP counts
    for (t in all_trait_names) {
        snp_has_trait <- rowSums(sapply(method_cols, function(m) {
            grepl(paste0("(^|,)", t, "($|,)"), snp_dt[[m]], fixed = FALSE)
        })) > 0
        row[[paste0(t, "_snps")]] <- sum(snp_has_trait)
    }

    # Per-method SNP counts
    for (m in all_method_names) {
        if (m %in% method_cols) {
            row[[paste0(m, "_snps")]] <- sum(snp_dt[[m]] != "" & !is.na(snp_dt[[m]]))
        } else {
            row[[paste0(m, "_snps")]] <- 0L
        }
    }

    row
}) %>% rbindlist(fill = TRUE)

combined_rows <- combined_rows %>% dplyr::arrange(chr, start)
fwrite(combined_rows, OUTPUT_COMBINED, sep = "\t", quote = FALSE)
message(paste0("INFO: Saved ", nrow(combined_rows), " combined region rows"))
message("INFO: Complete")
