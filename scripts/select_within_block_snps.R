#!/usr/bin/env Rscript
# Select per-SNP significant entities within significant LD blocks.
# Produces selected_snps.tsv with the same schema as combine_selected_snps.R output,
# feeding Manhattan plots and Gradient Forest maladaptation.
#
# Within-block selection methods:
#   within_order : p < min_p_in_block * multiplier  (default, scientifically motivated)
#   threshold    : p < threshold
#   top_n        : top N SNPs per block per trait
#   lead         : only the lead (minimum p) SNP per block per trait
#   all          : all SNPs in the block (for debugging)

library(data.table)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
SELECTED_BLOCKS     <- args[1]  # selected_blocks.tsv
PVALUES_FILES       <- args[2]  # space-separated per-method pvalues TSVs
WITHIN_METHOD       <- args[3]  # within_order | threshold | top_n | lead | all
MULTIPLIER          <- args[4] %>% as.numeric
THRESHOLD           <- args[5] %>% as.numeric
TOP_N               <- args[6] %>% as.integer
PREDICTORS_SELECTED <- args[7] %>% str_split(",") %>% unlist
OUTPUT              <- args[8]

message("INFO: Selecting within-block SNPs")
message(paste0("INFO: Method: ", WITHIN_METHOD))

# Load selected blocks
blocks <- fread(SELECTED_BLOCKS)
blocks$chr <- as.character(blocks$chr)

if (nrow(blocks) == 0) {
    message("WARNING: No selected blocks — writing empty selected_snps.tsv")
    empty_dt <- data.table(SNPID = character(), chr = character(), pos = integer(),
                            traits = character(), methods = character(),
                            min_pvalue = numeric())
    fwrite(empty_dt, OUTPUT, sep = "\t")
    quit(save = "no", status = 0)
}

# Load per-method p-value tables
pvalues_vec <- str_split(PVALUES_FILES, " ")[[1]]
pvalues_vec <- pvalues_vec[pvalues_vec != ""]

methods_vec <- pvalues_vec %>%
    str_split("/") %>%
    sapply(function(x) x[length(x) - 1])

pvals_lst <- lapply(seq_along(pvalues_vec), function(i) {
    dt <- fread(pvalues_vec[i])
    dt$chr <- as.character(dt$chr)
    dt[, method_name := methods_vec[i]]
    dt
}) %>% setNames(methods_vec)

# Collect trait column names (all non-fixed columns across files)
fixed_cols <- c("SNPID", "chr", "pos")
trait_cols <- unique(unlist(lapply(pvals_lst, function(dt) {
    setdiff(names(dt), c(fixed_cols, "method_name"))
})))
trait_cols <- trait_cols[trait_cols %in% PREDICTORS_SELECTED]

if (length(trait_cols) == 0) {
    message("WARNING: No matching trait columns found")
    empty_dt <- data.table(SNPID = character(), chr = character(), pos = integer(),
                            traits = character(), methods = character(),
                            min_pvalue = numeric())
    fwrite(empty_dt, OUTPUT, sep = "\t")
    quit(save = "no", status = 0)
}

# For each block × method × trait, select SNPs
select_snps_in_block <- function(snp_pvals, method, trait, block_id) {
    p_vec <- snp_pvals[[trait]]
    if (all(is.na(p_vec))) return(NULL)
    min_p <- min(p_vec, na.rm = TRUE)

    keep <- switch(WITHIN_METHOD,
        within_order = p_vec <= min_p * MULTIPLIER,
        threshold    = p_vec < THRESHOLD,
        top_n        = p_vec <= sort(p_vec)[min(TOP_N, length(p_vec))],
        lead         = p_vec == min_p,
        all          = rep(TRUE, length(p_vec)),
        stop(paste0("Unknown within method: ", WITHIN_METHOD))
    )
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) return(NULL)
    snp_pvals[keep, .(SNPID, chr, pos)] %>%
        dplyr::mutate(pvalue = p_vec[keep], trait = trait, method = method,
                      block_id = block_id)
}

result_rows <- lapply(seq_len(nrow(blocks)), function(bi) {
    blk <- blocks[bi]
    blk_chr   <- blk$chr
    blk_start <- blk$start
    blk_end   <- blk$end
    blk_id    <- blk$block_id

    lapply(methods_vec, function(m) {
        dt <- pvals_lst[[m]]
        in_block <- dt$chr == blk_chr & dt$pos >= blk_start & dt$pos <= blk_end
        snps_in <- dt[in_block]
        if (nrow(snps_in) == 0) return(NULL)
        lapply(trait_cols, function(trait) {
            if (!trait %in% names(snps_in)) return(NULL)
            select_snps_in_block(snps_in, m, trait, blk_id)
        }) %>% Filter(Negate(is.null), .) %>% rbindlist(fill = TRUE)
    }) %>% Filter(Negate(is.null), .) %>% rbindlist(fill = TRUE)
}) %>%
    Filter(Negate(is.null), .) %>%
    rbindlist(fill = TRUE)

if (is.null(result_rows) || nrow(result_rows) == 0) {
    message("WARNING: No SNPs selected within blocks")
    empty_dt <- data.table(SNPID = character(), chr = character(), pos = integer(),
                            traits = character(), methods = character(),
                            min_pvalue = numeric())
    fwrite(empty_dt, OUTPUT, sep = "\t")
    quit(save = "no", status = 0)
}

# Collapse to one row per SNPID (same schema as combine_selected_snps.R)
all_snp_ids <- unique(result_rows$SNPID)

method_traits <- lapply(methods_vec, function(m) {
    result_rows %>%
        dplyr::filter(method == m, SNPID %in% all_snp_ids) %>%
        group_by(SNPID) %>%
        summarise(traits = paste(sort(unique(trait)), collapse = ","), .groups = "drop")
}) %>% setNames(methods_vec)

summary_dt <- result_rows %>%
    group_by(SNPID, chr, pos) %>%
    summarise(min_pvalue = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(chr, pos) %>%
    as.data.frame

for (m in methods_vec) {
    mt <- method_traits[[m]]
    if (nrow(mt) > 0) {
        summary_dt <- left_join(summary_dt, mt %>% dplyr::rename(!!m := traits), by = "SNPID")
    } else {
        summary_dt[[m]] <- NA_character_
    }
    summary_dt[[m]][is.na(summary_dt[[m]])] <- ""
}

col_order <- c("SNPID", "chr", "pos", methods_vec, "min_pvalue")
summary_dt <- summary_dt[, col_order]

fwrite(as.data.table(summary_dt), OUTPUT, sep = "\t")
message(paste0("INFO: Saved ", nrow(summary_dt), " within-block SNPs to ", OUTPUT))

for (m in methods_vec) {
    n <- sum(summary_dt[[m]] != "", na.rm = TRUE)
    message(paste0("INFO: ", m, ": ", n, " SNPs with significant traits"))
}
