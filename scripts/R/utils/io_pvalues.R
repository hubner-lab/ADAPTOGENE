# io_pvalues.R — read/write helpers for p-value and q-value tables, and association data loading
# Usage: source("/pipeline/scripts/R/utils/io_pvalues.R")
# Requires: data.table, dplyr, stringr loaded by the sourcing script
# Requires: pval_threshold.R sourced BEFORE this file (for compute_pval_threshold)

# Read a p-values TSV produced by EMMAX/LFMM/GAPIT.
# chr column forced to character to avoid integer-inference on numeric chromosome names.
read_pvalues_tsv <- function(path, ...) {
    data.table::fread(path, sep = '\t', header = TRUE,
                      colClasses = c("chr" = "character"), ...)
}

# Write a p-values table (no quoting, tab-separated).
write_pvalues <- function(df, path) {
    data.table::fwrite(df, path, sep = '\t', quote = FALSE)
}

# Write a q-values table (no quoting, tab-separated).
write_qvalues <- function(df, path) {
    data.table::fwrite(df, path, sep = '\t', quote = FALSE)
}

# Load and threshold association data from an assoc_info list.
# assoc_info: named list produced by parse_assoc_files_str()
# traits: character vector of trait names
# panel_label: optional string added as a "panel" column (used by Miami plots)
# Returns list(data = data.table, thresholds = named list of -log10(threshold) per "method_trait")
# Requires compute_pval_threshold() from pval_threshold.R in the calling environment.
load_assoc_data <- function(assoc_info, traits, panel_label = NULL) {
    all_data        <- list()
    method_thresholds <- list()

    for (method_name in names(assoc_info)) {
        info <- assoc_info[[method_name]]
        if (!is.null(panel_label)) {
            message(paste0('INFO: Loading ', panel_label, ' ', method_name, ' from ', info$filepath))
        } else {
            message(paste0('INFO: Loading ', method_name, ' from ', info$filepath))
        }

        snps_assoc      <- read_pvalues_tsv(info$filepath)
        adjustment      <- stringr::str_split(info$adjust, '_')[[1]][1]
        pval_value      <- stringr::str_split(info$adjust, '_')[[1]][2] %>% as.numeric()

        for (trait in traits) {
            if (!(trait %in% colnames(snps_assoc))) {
                message(paste0('WARNING: Trait ', trait, ' not found in ', method_name))
                next
            }

            threshold_value <- compute_pval_threshold(snps_assoc[[trait]], adjustment, pval_value)
            threshold_log10 <- -log10(threshold_value)
            method_thresholds[[paste0(method_name, "_", trait)]] <- threshold_log10

            trait_data <- snps_assoc %>%
                dplyr::select(SNPID, chr, pos, !!dplyr::sym(trait)) %>%
                dplyr::rename(pvalue = !!dplyr::sym(trait)) %>%
                dplyr::mutate(
                    method         = method_name,
                    trait          = trait,
                    log10p         = -log10(pvalue),
                    is_significant = log10p >= threshold_log10
                )

            if (!is.null(panel_label)) {
                trait_data <- trait_data %>% dplyr::mutate(panel = panel_label)
            }

            n_sig <- sum(trait_data$is_significant)
            if (!is.null(panel_label)) {
                message(paste0('INFO: ', panel_label, ' ', method_name, ' - ', trait, ': ', n_sig, ' significant SNPs'))
            } else {
                message(paste0('INFO: ', method_name, ' - ', trait, ': ', n_sig, ' significant SNPs'))
            }

            all_data[[paste0(method_name, "_", trait)]] <- trait_data
        }
    }

    list(data = dplyr::bind_rows(all_data), thresholds = method_thresholds)
}
