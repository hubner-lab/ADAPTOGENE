#!/usr/bin/env Rscript
# Combine significant blocks from all methods into a single selected_blocks.tsv.
# Block IDs are unique keys (chr:start-end) shared across all methods,
# so combining is a set operation on block_id — no coordinate clustering needed.

library(data.table)
library(dplyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
SIGBLOCKS_FILES      <- args[1]  # space-separated list of sig_blocks files
METHOD               <- args[2]  # All, Overlap, MethodOverlap (or single method name)
GAP                  <- args[3] %>% as.numeric  # unused for block mode (kept for API parity)
PREDICTORS_SELECTED  <- args[4] %>% str_split(",") %>% unlist
OUTPUT               <- args[5]

message("INFO: Combining significant blocks")
message(paste0("INFO: Method: ", METHOD))

sigblocks_vec <- str_split(SIGBLOCKS_FILES, " ")[[1]]
sigblocks_vec <- sigblocks_vec[sigblocks_vec != ""]
message(paste0("INFO: Input files: ", paste(sigblocks_vec, collapse = ", ")))

methods_vec <- sigblocks_vec %>%
    str_split("/") %>%
    sapply(function(x) x[length(x) - 1])

# Load and filter to configured predictors
sig_lst <- lapply(seq_along(sigblocks_vec), function(i) {
    dt <- fread(sigblocks_vec[i])
    dt$chr <- as.character(dt$chr)
    dt %>% dplyr::filter(trait %in% PREDICTORS_SELECTED)
}) %>% setNames(methods_vec)

total_rows <- sum(sapply(sig_lst, nrow))
if (total_rows == 0) {
    message("WARNING: No significant blocks in any method")
    empty_dt <- data.table(
        block_id   = character(), chr = character(),
        start      = integer(),   end = integer(),
        traits     = character(), methods = character(),
        min_pvalue = numeric()
    )
    fwrite(empty_dt, OUTPUT, sep = "\t")
    quit(save = "no", status = 0)
}

# Select block_ids based on combine METHOD
if (METHOD %in% names(sig_lst)) {
    selected_ids <- unique(sig_lst[[METHOD]]$block_id)
} else if (METHOD %in% c("All", "Sum")) {
    selected_ids <- unique(rbindlist(sig_lst, fill = TRUE)$block_id)
} else if (METHOD == "Overlap") {
    # Blocks significant in at least two methods
    id_counts <- table(rbindlist(sig_lst, fill = TRUE)$block_id)
    selected_ids <- names(id_counts[id_counts >= 2])
} else if (METHOD %in% c("MethodOverlap", "PairOverlap")) {
    # Blocks significant for the same trait in at least two methods
    all_blks <- rbindlist(sig_lst, fill = TRUE)
    dup_key <- paste(all_blks$block_id, all_blks$trait)
    key_counts <- table(dup_key)
    valid_keys <- names(key_counts[key_counts >= 2])
    selected_ids <- unique(all_blks$block_id[paste(all_blks$block_id, all_blks$trait) %in% valid_keys])
} else {
    stop(paste0("Unknown METHOD: ", METHOD))
}

if (length(selected_ids) == 0) {
    message("WARNING: No blocks selected with specified method")
    empty_dt <- data.table(
        block_id   = character(), chr = character(),
        start      = integer(),   end = integer(),
        traits     = character(), methods = character(),
        min_pvalue = numeric()
    )
    fwrite(empty_dt, OUTPUT, sep = "\t")
    quit(save = "no", status = 0)
}

# Build combined table with one row per block_id
all_sig <- rbindlist(sig_lst, fill = TRUE) %>%
    dplyr::filter(block_id %in% selected_ids)

result <- all_sig %>%
    group_by(block_id, chr, start, end) %>%
    summarise(
        traits     = paste(sort(unique(trait)),  collapse = ","),
        methods    = paste(sort(unique(method)), collapse = ","),
        min_pvalue = min(pvalue, na.rm = TRUE),
        .groups    = "drop"
    ) %>%
    dplyr::arrange(chr, start) %>%
    as.data.table

fwrite(result, OUTPUT, sep = "\t")
message(paste0("INFO: Saved ", nrow(result), " selected blocks to ", OUTPUT))
