# pval_threshold.R — p-value threshold computation helpers
# Usage: source("/pipeline/scripts/R/utils/pval_threshold.R")
# Requires: qvalue package loaded by the sourcing script

# Return the maximum p-value still below the FDR threshold (NA if none pass).
# Unified version: always returns NA (not NULL) when no SNPs pass.
max_pvalue_fdr <- function(pvalues, fdr) {
    qvalues_result <- qvalue(pvalues)
    significant_pvalues <- pvalues[qvalues_result$qvalues < fdr]
    if (length(significant_pvalues) > 0) {
        return(max(significant_pvalues))
    } else {
        return(NA_real_)
    }
}

# Return the maximum p-value among the top N most significant.
max_pvalue_top <- function(pvalues, topN) {
    if (length(pvalues) < topN) {
        stop("topN is larger than the number of available p-values")
    }
    sorted_pvalues <- sort(pvalues, decreasing = FALSE)
    return(max(sorted_pvalues[1:topN]))
}

# Compute the raw p-value cutoff for a given adjustment method.
# adjustment: one of "bonf", "qval", "top"
# value: the numeric parameter (alpha for bonf/qval, N for top)
# Returns the raw p-value threshold (not -log10).
# For "qval" with no passing SNPs, falls back silently to Bonferroni.
compute_pval_threshold <- function(pvalues, adjustment, value) {
    if (adjustment == 'bonf') {
        return(value / length(pvalues))
    } else if (adjustment == 'qval') {
        t <- max_pvalue_fdr(pvalues, value)
        if (is.na(t)) {
            message("WARNING: No significant SNPs at q-value threshold, falling back to Bonferroni")
            return(value / length(pvalues))
        }
        return(t)
    } else if (adjustment == 'top') {
        return(max_pvalue_top(pvalues, as.numeric(value)))
    } else {
        stop(paste0("Unknown adjustment method: ", adjustment))
    }
}

# Compute q-values with tryCatch fallback to BH adjustment.
# The safer version from emmax_phenotypes.R — use everywhere instead of bare qvalue().
compute_qvalues_safe <- function(pvalues) {
    tryCatch(
        qvalue(pvalues)$qvalues,
        error = function(e) {
            message(paste0("WARNING: qvalue failed (", e$message, "), using BH adjustment"))
            p.adjust(pvalues, method = 'BH')
        }
    )
}
