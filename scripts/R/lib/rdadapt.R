# rdadapt.R — robust-Mahalanobis RDA outlier test (Capblancq et al. 2018).
# Usage: source("/pipeline/scripts/R/lib/rdadapt.R")
# Requires: robust, qvalue loaded by the sourcing script.
#
# Extracted from scripts/rda.R (the GEA-mode RDA scan) so it can be shared
# with preGEA's RDA-setup block (scripts/pregea_rda_setup.R) without
# duplicating the statistic — both need the SAME candidate-calling math, just
# applied to different Condition()-PC ladders. Verbatim port from
# Capblancq/RDA-genome-scan/Script_RDA.R; one deviation: drop=FALSE on the
# loadings subset (behavior-preserving given the K>=2 floor callers must
# enforce; guards against silent vector coercion if that floor is ever
# relaxed).

#' Robust-Mahalanobis RDA candidate test.
#'
#' @param rda_obj A fitted vegan::rda() (or partial RDA) object.
#' @param K Number of retained constrained axes to test over (must be >= 2 —
#'   covRob's pairwiseGK estimator needs >= 2 loading columns; callers must
#'   enforce this floor themselves, verified empirically: K=1 crashes).
#' @return list(p.values, q.values, gif_lambda, qvalue_method, distance) —
#'   distance is returned so candidate reporting never recomputes covRob
#'   (CLAUDE.md Rule 1: compute once, reuse downstream).
#' @noRd
rdadapt <- function(rda_obj, K) {
    loadings    <- rda_obj$CCA$v[, 1:as.numeric(K), drop = FALSE]
    resscale    <- apply(loadings, 2, scale)
    resmaha     <- covRob(resscale, distance = TRUE, na.action = na.omit,
                          estim = "pairwiseGK")$dist
    lambda      <- median(resmaha) / qchisq(0.5, df = K)
    reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
    # qvalue fallback chain — storey -> qvalue(lambda=0) -> BH.
    qval_result <- tryCatch(
        list(qvalues = qvalue(reschi2test)$qvalues, method = "storey"),
        error = function(e) tryCatch(
            list(qvalues = qvalue(reschi2test, lambda = 0)$qvalues, method = "storey_lambda0"),
            error = function(e2) list(qvalues = p.adjust(reschi2test, "BH"), method = "BH")
        )
    )
    list(p.values = reschi2test, q.values = qval_result$qvalues,
         gif_lambda = lambda, qvalue_method = qval_result$method,
         distance = resmaha)
}
