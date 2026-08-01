#!/usr/bin/env Rscript
# =============================================================================
# eval_detection.R — score ONE GEA p-value table against the causal-locus truth.
#
# Implements Axis 2 of the two-axis evaluation harness specified in
# docs/laruson_automation_roadmap.md Phase 5, with three corrections the spec
# pre-dates (see docs/laruson_dataset_notes.md and the LD-decay measurement in
# {PROJECT}_results/pipeline_summary.tsv):
#
#   1. Exact-position matching is PRIMARY. The spec's window rule is a
#      sensitivity view, requested via --tp-window-kb (comma list, 0 = exact).
#   2. Recall is reported against TWO denominators: causal loci actually present
#      in the tested SNP set (`recall_testable`) and all simulated causal loci
#      (`recall_simulated`). The second exposes the MAF filter ceiling — at the
#      full Láruson cohort only 28 of 102 causal loci survive MAF 0.05.
#   3. Hits on `linked_neutral` positions are reported separately as
#      `expected_linked` — they are physically linked to true QTNs and are NOT
#      pipeline errors. Only `background_neutral` hits count as false positives.
#
# Significance calling reuses compute_pval_threshold() from the pipeline's own
# scripts/R/utils/pval_threshold.R, so --mode=threshold reproduces
# find_sig_snps.R exactly. Scoring logic lives in benchmarks/lib_detection.R and
# is shared with sweep_thresholds.R — do not duplicate it here.
#
# Usage:
#   Rscript eval_detection.R --mode=threshold --pvalues=FILE --truth=FILE \
#       --label=STR --adjust=bonf|qval|top|custom --value=NUM \
#       [--tp-window-kb=0,0.6,5] [--traits=all|a,b] [--out=FILE]
#
#   Rscript eval_detection.R --mode=rank --pvalues=FILE --truth=FILE \
#       --label=STR [--tp-window-kb=0] [--traits=all|a,b] [--out=FILE]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(qvalue)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/pval_threshold.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args <- parse_kv_args(commandArgs(trailingOnly = TRUE))
req  <- function(k) {
    if (is.null(args[[k]]) || !nzchar(args[[k]])) stop("Missing required --", gsub("_", "-", k))
    args[[k]]
}
opt  <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]

MODE      <- req("mode")
PVAL_FILE <- req("pvalues")
LABEL     <- req("label")
TRUTH_F   <- req("truth")
OUT_TSV   <- opt("out", file.path(PIPELINE_ROOT, "benchmarks/laruson_eval/benchmark_eval.tsv"))
WINDOWS   <- as.numeric(strsplit(opt("tp_window_kb", "0"), ",", fixed = TRUE)[[1]])

if (!MODE %in% c("threshold", "rank")) stop("--mode must be 'threshold' or 'rank'")
if (any(is.na(WINDOWS)) || any(WINDOWS < 0)) stop("--tp-window-kb must be non-negative numbers")

truth <- load_truth(TRUTH_F)
lp    <- load_pvalues(PVAL_FILE, opt("traits", "all"))
pv    <- lp$pv; trait_cols <- lp$trait_cols

# --genomic-control=yes recalibrates before scoring. PreGEA writes RAW p-values
# (genomic_control = FALSE) while production mode=gea applies GC, so this is how
# a ladder rung gets scored under production calling behaviour.
if (tolower(opt("genomic_control", "no")) %in% c("yes", "true", "1")) {
    pv <- apply_genomic_control(pv, trait_cols)
    LABEL <- paste0(LABEL, "+gc")
}

testable   <- testable_causal(pv, trait_cols, truth)
n_testable <- length(testable)
if (n_testable == 0L) stop("No causal loci present in the tested SNP set — check chr naming in ", PVAL_FILE)

E <- new_emitter()
E$emit(LABEL, "truth", "n_causal_simulated", truth$n_simulated)
E$emit(LABEL, "truth", "n_causal_testable",  n_testable)
E$emit(LABEL, "truth", "n_snps_tested",      nrow(pv))

# ------------------------------------------------------------- mode=threshold
if (MODE == "threshold") {
    ADJUST <- req("adjust"); VALUE <- req("value")
    res    <- call_by_threshold(pv, trait_cols, ADJUST, VALUE, quiet = FALSE)
    for (tc in names(res$info)) {
        E$emit(LABEL, paste0("trait_", tc), "threshold", res$info[[tc]]$threshold)
        E$emit(LABEL, paste0("trait_", tc), "n_tested",  res$info[[tc]]$n_tested)
        if (!identical(res$info[[tc]]$status, "ok")) {
            message("NOTE: trait ", tc, " threshold status = ", res$info[[tc]]$status,
                    " — contributed 0 hits")
        }
    }
    config <- paste0(LABEL, "|", ADJUST, "_", VALUE)
    for (w in WINDOWS) {
        s <- score_calls(res$called, w, truth, testable)
        E$emit_scored(config, s, w)
        message(sprintf(
            "%-46s w=%-4skb called=%-5d TP=%-3d (any=%-3d) linked=%-4d FP_bg=%-4d prec=%.3f rec_test=%.3f rec_sim=%.3f F1=%.3f",
            config, w, s$n_called, s$tp, s$tp_any_causal, s$expected_linked, s$fp_background,
            ifelse(is.na(s$precision_strict), 0, s$precision_strict),
            s$recall_testable, s$recall_simulated, ifelse(is.na(s$f1), 0, s$f1)))
    }
}

# ------------------------------------------------------------------ mode=rank
if (MODE == "rank") {
    rk     <- rank_keys(pv, trait_cols)
    ranked <- rk$keys
    n_rank <- length(ranked)
    w      <- WINDOWS[1]
    E$emit(LABEL, "rank", "n_ranked", n_rank)

    is_causal <- ranked %in% truth$causal$key
    is_bg     <- ranked %in% truth$bg_keys
    tp_cum <- cumsum(is_causal); fp_cum <- cumsum(is_bg)
    prec   <- tp_cum / (tp_cum + fp_cum); prec[is.nan(prec)] <- 1
    rec    <- tp_cum / n_testable
    f1v    <- ifelse((prec + rec) > 0, 2 * prec * rec / (prec + rec), 0)

    # AUC-PR by the interpolation-free step rule: mean precision at each TP,
    # with un-retrieved causal loci contributing zero.
    auc_pr <- sum(prec[is_causal]) / n_testable
    best   <- which.max(f1v)
    E$emit(LABEL, "rank", "auc_pr",          auc_pr)
    E$emit(LABEL, "rank", "best_f1",         f1v[best])
    E$emit(LABEL, "rank", "best_f1_n",       best)
    E$emit(LABEL, "rank", "best_f1_tp",      tp_cum[best])
    E$emit(LABEL, "rank", "best_f1_prec",    prec[best])
    E$emit(LABEL, "rank", "best_f1_recall",  rec[best])

    grid <- c(10, 20, 30, 50, 100, 200, 500)
    grid <- grid[grid <= n_rank]
    for (n in grid) {
        s <- score_calls(ranked[seq_len(n)], w, truth, testable)
        cfg <- paste0(LABEL, "|top", n)
        E$emit(cfg, paste0("window_", w, "kb"), "tp",              s$tp)
        E$emit(cfg, paste0("window_", w, "kb"), "fp_background",   s$fp_background)
        E$emit(cfg, paste0("window_", w, "kb"), "expected_linked", s$expected_linked)
        E$emit(cfg, paste0("window_", w, "kb"), "recall_testable", s$recall_testable)
    }
    message(sprintf("%-46s n_ranked=%-6d testable=%-3d AUC-PR=%.4f  bestF1=%.3f @n=%d (TP=%d)",
                    LABEL, n_rank, n_testable, auc_pr, f1v[best], best, tp_cum[best]))
    message(sprintf("%-46s TP@top: %s", LABEL,
                    paste(sapply(grid, function(n) paste0(n, "=", sum(is_causal[seq_len(n)]))),
                          collapse = " ")))
}

append_eval(E$collect(), OUT_TSV)
