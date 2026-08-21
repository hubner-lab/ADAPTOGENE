# pval_threshold.R — p-value threshold computation helpers
# Usage: source("/pipeline/scripts/R/utils/pval_threshold.R")
# Requires: qvalue package loaded by the sourcing script

# Return the maximum p-value still below the FDR threshold (NA if none pass).
# NA values in pvalues are silently dropped before qvalue computation.
max_pvalue_fdr <- function(pvalues, fdr) {
    pvalues <- pvalues[!is.na(pvalues)]
    if (length(pvalues) == 0) return(NA_real_)
    qvalues_result <- qvalue(pvalues)
    significant_pvalues <- pvalues[qvalues_result$qvalues < fdr]
    if (length(significant_pvalues) > 0) max(significant_pvalues) else NA_real_
}

# Return the maximum p-value among the top N most significant.
# NA values are dropped; topN is compared against the non-NA count.
max_pvalue_top <- function(pvalues, topN) {
    pvalues <- pvalues[!is.na(pvalues)]
    if (length(pvalues) < topN) {
        stop(paste0("topN (", topN, ") is larger than the number of non-NA p-values (",
                    length(pvalues), ")"))
    }
    sorted_pvalues <- sort(pvalues, decreasing = FALSE)
    max(sorted_pvalues[1:topN])
}

# Compute the raw p-value cutoff for a given adjustment method.
#
# NA values in pvalues are silently dropped before threshold computation — this
# is essential for WZA where some windows have NA p for a trait when no SNPs
# passed the MAF filter within that window.
#
# Returns a named list:
#   threshold    numeric (NA when status != "ok")
#   n_tested     integer — number of non-NA p-values used
#   n_na_dropped integer — number of NA entries removed
#   status       character: "ok" | "too_few_tests" | "no_tests"
#
# "too_few_tests" is returned (not a silent BH fallback) so the caller can
# surface an actionable warning. No automatic fallback is applied.
#
# THE RETURNED CUTOFF IS INCLUSIVE. Callers must select with `p <= threshold`.
# For 'qval' and 'top' the returned value is itself a member of the intended
# call set (max_pvalue_fdr() returns the largest passing p; max_pvalue_top()
# returns the N-th smallest p), so a strict `<` drops the boundary observation
# and every observation tied with it. Under 'top' the inclusive rule can
# therefore return MORE than N rows when the N-th smallest p has ties — that is
# correct, not a defect.
compute_pval_threshold <- function(pvalues, adjustment, value) {
    n_total      <- length(pvalues)
    pvalues      <- pvalues[!is.na(pvalues)]
    n_tested     <- length(pvalues)
    n_na_dropped <- n_total - n_tested

    if (n_tested == 0L) {
        return(list(threshold = NA_real_, n_tested = 0L,
                    n_na_dropped = n_na_dropped, status = "no_tests"))
    }

    if (adjustment == 'bonf') {
        return(list(threshold    = value / n_tested,
                    n_tested     = n_tested,
                    n_na_dropped = n_na_dropped,
                    status       = "ok"))
    }

    if (adjustment == 'qval') {
        # qvalue needs enough tests to estimate pi0; empirically fails below ~10.
        if (n_tested < 10L) {
            message(paste0("WARNING: qvalue requires >=10 tests; only ", n_tested,
                           " non-NA p-values available. Use bonf or top instead."))
            return(list(threshold = NA_real_, n_tested = n_tested,
                        n_na_dropped = n_na_dropped, status = "too_few_tests"))
        }
        t <- tryCatch(
            max_pvalue_fdr(pvalues, value),
            error = function(e) {
                message(paste0("WARNING: qvalue() failed: ", e$message,
                               ". Use bonf or top instead."))
                NA_real_
            }
        )
        if (is.na(t)) {
            return(list(threshold = NA_real_, n_tested = n_tested,
                        n_na_dropped = n_na_dropped, status = "too_few_tests"))
        }
        return(list(threshold = t, n_tested = n_tested,
                    n_na_dropped = n_na_dropped, status = "ok"))
    }

    if (adjustment == 'top') {
        topN <- as.integer(value)
        if (n_tested < topN) {
            message(paste0("WARNING: top N=", topN, " > n_tested=", n_tested,
                           ". Use a smaller N or switch to bonf."))
            return(list(threshold = NA_real_, n_tested = n_tested,
                        n_na_dropped = n_na_dropped, status = "too_few_tests"))
        }
        return(list(threshold    = max_pvalue_top(pvalues, topN),
                    n_tested     = n_tested,
                    n_na_dropped = n_na_dropped,
                    status       = "ok"))
    }

    if (adjustment == 'custom') {
        custom_threshold <- suppressWarnings(as.numeric(value))
        if (is.na(custom_threshold) || custom_threshold <= 0) {
            message(paste0("WARNING: custom threshold value '", value,
                           "' invalid (must be a positive number). Returning no threshold."))
            return(list(threshold = NA_real_, n_tested = n_tested,
                        n_na_dropped = n_na_dropped, status = "too_few_tests"))
        }
        return(list(threshold    = custom_threshold,
                    n_tested     = n_tested,
                    n_na_dropped = n_na_dropped,
                    status       = "ok"))
    }

    stop(paste0("Unknown adjustment method: ", adjustment))
}

# Compute q-values with tryCatch fallback to BH adjustment.
# The safer version for scripts that already have a per-SNP pvalue vector
# and want soft q-values (not a significance threshold). For hard threshold
# computation use compute_pval_threshold() instead.
compute_qvalues_safe <- function(pvalues) {
    tryCatch(
        qvalue(pvalues)$qvalues,
        error = function(e) {
            message(paste0("WARNING: qvalue failed (", e$message, "), using BH adjustment"))
            p.adjust(pvalues, method = 'BH')
        }
    )
}
