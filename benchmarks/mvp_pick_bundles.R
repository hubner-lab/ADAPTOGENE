#!/usr/bin/env Rscript
# mvp_pick_bundles.R -- choose the two mixed-rung bundles for one seed and print their
# sweep_thresholds.R specs to stdout, one per line:  kind <TAB> axis <TAB> truth <TAB> spec
#
#   Rscript benchmarks/mvp_pick_bundles.R --seed=1232548 [--cells=5]
#
# THE TWO BUNDLES, and why there are two
# --------------------------------------
# Per-cell bundles hold every method at the SAME rung index, which is an artifact of how the
# runs were amortized, not a configuration anyone would choose. The interesting bundles let
# each method sit at its own rung:
#
#   oracle  each method at the rung with the highest AUC-PR against the truth table.
#           This is an UPPER BOUND and nothing else. Selecting a hyperparameter by scoring it
#           against the answer key is not something a user can do -- reporting it as "what the
#           pipeline achieves" would reduce the whole claim to "knowing the answer helps you
#           find the answer". It is here to bound how much the achievable arm leaves on the
#           table.
#
#   pregea  each method at the rung mode=pregea recommended, read out of
#           PreGEA/tables/pregea_recommendations.tsv. No truth table is consulted. THIS is the
#           configuration a user following the pipeline's own workflow actually reaches, and
#           it is the number the paper claim rests on.
#
# BLINK has no PreGEA recommendation of its own (pregea_recommend.R emits LFMM.K, EMMAX.n_pcs,
# RDA.condition_pcs, RDA.axes, RDA.predictor_set only). It inherits EMMAX's n_pcs, since both
# parameters are literally the same quantity -- the number of PCA covariates entered as fixed
# effects. Flagged in the output so the report can state it rather than bury it.
#
# A recommendation that falls between two rungs snaps to the nearest one and is recorded in
# pregea_offladder.tsv; a bundle built on a snapped rung is not the recommended configuration
# and must not be reported as one without saying so.

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || is.na(x) || identical(x, "")) y else x
A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
    a <- sub("^--", "", a)
    A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}
SEED  <- A$seed %||% stop("Missing --seed")
CELLS <- as.integer(if (is.null(A$cells)) 5L else A$cells)
ROOT  <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
PROJ  <- paste0("MVP", SEED)

UNIVARIATE  <- c("EMMAX", "LFMM", "BLINK")
ALL_METHODS <- c("EMMAX", "LFMM", "RDA", "BLINK")

sweep_dir  <- file.path(ROOT, "benchmarks/mvp_eval/sweep", PROJ)
params_dir <- file.path("/pipeline/benchmarks/mvp_eval/params", PROJ)   # container-side
truth_dir  <- file.path("/pipeline/data/mvp", PROJ)

n_causal <- function(f) {
    tr <- fread(file.path(ROOT, "data/mvp", PROJ, f), colClasses = c(chr = "character"))
    sum(tr$category == "causal")
}
axes <- list(
    temp = list(truth = "truth_temp.tsv", methods = UNIVARIATE,  suffix = "_pvalues_bio_1.tsv"),
    sal  = list(truth = "truth_sal.tsv",  methods = UNIVARIATE,  suffix = "_pvalues_bio_2.tsv"),
    any  = list(truth = "truth_any.tsv",  methods = ALL_METHODS, suffix = "_pvalues.tsv"))
axes <- axes[vapply(axes, function(a) n_causal(a$truth) > 0L, logical(1))]

# One line per bundle, read by mvp_score_sweep.sh with `IFS=$'\t' read`. NOT cat(fill=TRUE):
# that wraps at getOption("width") and would split the method spec across two lines.
emit <- function(kind, axis, spec)
    cat(paste(kind, axis, file.path(truth_dir, axes[[axis]]$truth), spec, sep = "\t"), "\n", sep = "")
spec_of <- function(axis, cell_of) paste(vapply(axes[[axis]]$methods, function(m)
    sprintf("%s=%s/c%d/%s%s", m, params_dir, cell_of[[m]], m, axes[[axis]]$suffix),
    character(1)), collapse = ",")

# ------------------------------------------------------------------- oracle bundle
rank_f <- file.path(sweep_dir, "rank.tsv")
if (file.exists(rank_f)) {
    rk <- fread(rank_f)
    rk <- rk[metric == "auc_pr"]
    # config is the --label: "MVP{seed}|{method}|c{i}|{axis}". The parsed axis goes into
    # ax_parsed, NOT axis: eval_detection.R already emits an `axis` column (always
    # "detection"), and assigning a mix of new and existing names on a `:=` LHS is a good
    # way to filter on something other than what you think you are filtering on.
    rk[, c("proj", "method", "cell", "ax_parsed") := tstrsplit(config, "|", fixed = TRUE)]
    rk[, cell_i := as.integer(sub("^c", "", cell))]
    for (ax in names(axes)) {
        sub_rk <- rk[ax_parsed == ax & method %in% axes[[ax]]$methods]
        if (!nrow(sub_rk)) next
        best <- sub_rk[order(-as.numeric(value))][, .SD[1], by = method]
        cell_of <- setNames(as.list(best$cell_i), best$method)
        if (!all(axes[[ax]]$methods %in% names(cell_of))) next
        emit("oracle", ax, spec_of(ax, cell_of))
    }
} else {
    message("WARN: no rank.tsv for ", PROJ, " — oracle bundle skipped")
}

# ------------------------------------------------------------------- pregea bundle
rec_f   <- file.path(ROOT, sprintf("%s_results/PreGEA/tables/pregea_recommendations.tsv", PROJ))
cells_f <- file.path(ROOT, "benchmarks/mvp_sweep_cells.tsv")
if (file.exists(rec_f) && file.exists(cells_f)) {
    rec   <- fread(rec_f)
    ladder <- fread(cells_f, colClasses = c(seed = "character"))[seed == SEED]
    ladder[, cell_i := as.integer(sub("^c", "", cell))]

    want <- c(EMMAX = "EMMAX.n_pcs", LFMM = "LFMM.K", RDA = "RDA.condition_pcs", BLINK = "EMMAX.n_pcs")
    off  <- list()
    cell_of <- list()
    for (m in ALL_METHODS) {
        key <- want[[m]]
        row <- rec[paste0(method, ".", param) == key]
        rung <- ladder[method == m]
        if (!nrow(row) || !nrow(rung)) { cell_of[[m]] <- NA_integer_; next }
        v  <- suppressWarnings(as.numeric(row$recommended_value[1]))
        d  <- abs(rung$value - v)
        pick <- rung$cell_i[which.min(d)]
        cell_of[[m]] <- pick
        if (min(d) > 0) off[[length(off) + 1L]] <- data.table(
            seed = SEED, method = m, param = key, recommended = v,
            snapped_to = rung$value[which.min(d)], cell = pick,
            inherited = (m == "BLINK"))
        else if (m == "BLINK") off[[length(off) + 1L]] <- data.table(
            seed = SEED, method = m, param = key, recommended = v,
            snapped_to = v, cell = pick, inherited = TRUE)
    }
    if (length(off)) fwrite(rbindlist(off), file.path(sweep_dir, "pregea_offladder.tsv"), sep = "\t")

    # The full recommended-rung mapping, not just the snapped rows. Downstream needs it to
    # report the ACHIEVABLE arm per method -- comparing the published baseline against our
    # default rung alone would understate the pipeline, since the default is measurably the
    # wrong rung for EMMAX and RDA on this landscape.
    pc <- rbindlist(lapply(ALL_METHODS, function(m) data.table(
        seed = SEED, method = m,
        cell = if (is.na(cell_of[[m]])) NA_character_ else paste0("c", cell_of[[m]]),
        param = want[[m]], inherited = (m == "BLINK"))))
    fwrite(pc, file.path(sweep_dir, "pregea_cells.tsv"), sep = "\t")
    for (ax in names(axes)) {
        need <- axes[[ax]]$methods
        if (any(is.na(unlist(cell_of[need])))) {
            message("WARN: ", PROJ, " axis ", ax, " — missing PreGEA recommendation, bundle skipped")
            next
        }
        emit("pregea", ax, spec_of(ax, cell_of))
    }
} else {
    message("WARN: no pregea_recommendations.tsv for ", PROJ, " — achievable bundle unavailable")
}
