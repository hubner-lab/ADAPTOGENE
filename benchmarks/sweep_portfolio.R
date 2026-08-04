#!/usr/bin/env Rscript
# =============================================================================
# sweep_portfolio.R -- STAGE 2 of the GEA calling search: heterogeneous per-method
# operating points, over restricted method subsets, at several agreement windows.
#
# WHY THIS IS A SEPARATE SCRIPT FROM sweep_thresholds.R
# ----------------------------------------------------
# sweep_thresholds.R answers "what does each method do across the whole threshold
# grid, and what happens if every method is forced onto the SAME rule". That shared
# rule is the thing being questioned here: methods differ enormously in p-value
# calibration (RDA's pmax gives correct ORDERING but not a correct p-value,
# MVP_README.md:283-306), so forcing one rule on all of them handicaps exactly the
# combine rules it is meant to evaluate. The pipeline never required it either --
# GEA.configs is a list of {method, adjust, threshold} and has been per-method since
# common.smk:395-409.
#
# Keeping the two stages apart also keeps journal 05's script auditable: its only
# journal-06 edits are the capped grid, the combine window and calls_per_tp.
#
# THE COMBINATORICS, AND HOW THEY ARE BOUNDED
# -------------------------------------------
# 11 methods x ~17 grid points is not enumerable as a cross-product. Instead:
#
#   Stage 1 (sweep_thresholds.R, already run) gives each method's full curve.
#   Here:  reduce each method to <= 3 CANDIDATE operating points --
#            best-F1 | max-recall at precision_strict >= HP_PRECISION | bonf_0.05
#          then enumerate all pairs and triples drawn from the top-N methods by solo
#          F1, plus the all-methods at_least_k rules, each method sitting at each of
#          its own candidates.
#
# --rank-source is what keeps the three configuration groups (default / pregea /
# oracle) comparable: the method pool is ranked ONCE, on the default cell, and that
# same pool is reused for the others. Re-ranking per group would silently give each
# group a different candidate pool -- the same failure journal 05 sec.10 called out for
# reading combine evidence off cell c3.
#
# Usage:
#   Rscript sweep_portfolio.R --methods="EMMAX=/p/a.tsv,LFMM=/p/b.tsv,..." \
#       --truth=FILE --tag=STR --sweep=/p/{TAG}_threshold_sweep.tsv \
#       [--rank-source=FILE] [--outdir=DIR] [--combine-window-kb=0,1,2.5,5] \
#       [--top-max=100] [--pool=6] [--hp-precision=0.9] [--tp-window-kb=0]
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

TAG          <- req("tag")
TRUTH_F      <- req("truth")
SWEEP_F      <- req("sweep")
RANK_F       <- opt("rank_source", SWEEP_F)
OUTDIR       <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/portfolio"))
WINDOW       <- as.numeric(opt("tp_window_kb", "0"))
POOL_N       <- as.integer(opt("pool", "6"))
HP_PRECISION <- as.numeric(opt("hp_precision", "0.9"))
TOP_MAX      <- as.numeric(opt("top_max", "100"))
# Optional. When given, a SECOND table {TAG}_transfer.tsv is written holding fixed schemes
# imported from another replicate. Without the flag nothing changes and {TAG}_portfolio.tsv
# is byte-identical -- the same regression discipline journal 06 used for the config writer.
SCHEMES_F    <- opt("schemes", "")
TOP_GRID     <- c(10, 20, 30, 50, 100, 200, 500, 1000)
TOP_GRID     <- TOP_GRID[TOP_GRID <= TOP_MAX]
COMBINE_WINDOWS <- as.numeric(strsplit(opt("combine_window_kb", "0,1,2.5,5"), ",", fixed = TRUE)[[1]])
if (any(is.na(COMBINE_WINDOWS))) stop("--combine-window-kb must be a comma-separated numeric list")

spec <- strsplit(req("methods"), ",", fixed = TRUE)[[1]]
methods <- list()
for (s in spec) {
    kv <- strsplit(s, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Bad --methods entry (expected NAME=PATH): ", s)
    if (!file.exists(kv[2])) stop("Method p-value table not found: ", kv[2])
    methods[[kv[1]]] <- kv[2]
}
if (length(methods) < 2) stop("--methods needs at least two methods to combine")

truth <- load_truth(TRUTH_F)
message("INFO: truth = ", truth$n_simulated, " simulated causal loci")

# ------------------------------------------------------------------- load once
loaded <- list()
for (m in names(methods)) {
    lp <- load_pvalues(methods[[m]], "all")
    tc <- testable_causal(lp$pv, lp$trait_cols, truth)
    loaded[[m]] <- list(pv = lp$pv, trait_cols = lp$trait_cols,
                        testable = tc, rank = rank_keys(lp$pv, lp$trait_cols))
    message(sprintf("INFO: %-8s %d SNPs, traits=[%s], testable causal=%d",
                    m, nrow(lp$pv), paste(lp$trait_cols, collapse = ","), length(tc)))
}

# --------------------------------------------------- candidate operating points
# Three per method, deduplicated. The bonf_0.05 baseline is kept even when it is
# dominated: it is what the pipeline does today, so every table needs a row that says
# what the user gets without touching anything.
sw <- fread(SWEEP_F)
if (!"kind" %in% names(sw)) stop("--sweep table has no `kind` column: ", SWEEP_F)
solo <- sw[kind %in% names(loaded) & !grepl("^combine_", kind)]
if (!nrow(solo)) stop("No per-method rows found in ", SWEEP_F)

pick_candidates <- function(m) {
    d <- solo[kind == m]
    out <- list()
    bf <- d[!is.na(f1)][order(-f1, n_called)]
    if (nrow(bf)) out[["best_f1"]] <- bf[1, .(adjust, value)]
    hp <- d[!is.na(precision) & precision >= HP_PRECISION & n_called > 0][order(-recall_testable, n_called)]
    if (nrow(hp)) out[["hp_recall"]] <- hp[1, .(adjust, value)]
    out[["default"]] <- data.table(adjust = "bonf", value = "0.05")
    cand <- rbindlist(out, idcol = "role")
    unique(cand, by = c("adjust", "value"))
}
CAND <- rbindlist(lapply(names(loaded), function(m)
    cbind(method = m, pick_candidates(m))))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
fwrite(CAND, file.path(OUTDIR, paste0(TAG, "_candidates.tsv")), sep = "\t")

# ------------------------------------------------------------- called-set cache
# key: "METHOD|adjust_value". Every downstream combination reuses these, so each
# threshold is applied exactly once per method regardless of how many subsets use it.
CALLED <- list()
for (i in seq_len(nrow(CAND))) {
    m <- CAND$method[i]; a <- CAND$adjust[i]; v <- CAND$value[i]
    k <- paste0(m, "|", a, "_", v)
    if (!is.null(CALLED[[k]])) next
    CALLED[[k]] <- call_by_threshold(loaded[[m]]$pv, loaded[[m]]$trait_cols, a, v, quiet = TRUE)$called
}

# ----------------------------------------------------------------- method pool
# Ranked on --rank-source, which defaults to this tag's own sweep but is pointed at the
# DEFAULT cell's sweep when scoring pregea/oracle, so all three groups search the same
# subsets. See the header.
rk <- fread(RANK_F)
rk <- rk[kind %in% names(loaded) & !is.na(f1), .(best_f1 = max(f1)), by = kind][order(-best_f1)]
POOL <- head(rk$kind, POOL_N)
message("INFO: subset pool (top ", length(POOL), " by solo F1 on ", basename(RANK_F), "): ",
        paste(POOL, collapse = ", "))

subsets <- c(
    combn(POOL, 2, simplify = FALSE),
    if (length(POOL) >= 3) combn(POOL, 3, simplify = FALSE) else list(),
    list(names(loaded))                    # the all-methods bundle, for at_least_k
)

# --------------------------------------------------------------------- scoring
rows <- list()
add <- function(subset, assign_lbl, rule, k, cw, s, n_testable) {
    rows[[length(rows) + 1L]] <<- data.table(
        tag = TAG, subset = paste(subset, collapse = "+"), n_methods = length(subset),
        assignment = assign_lbl, rule = rule, k = k, combine_window_kb = cw,
        n_called = s$n_called, tp = s$tp, tp_any_causal = s$tp_any_causal,
        expected_linked = s$expected_linked, fp_background = s$fp_background,
        fn_testable = s$fn_testable,
        precision = s$precision_strict, precision_all = s$precision_all,
        recall_testable = s$recall_testable, recall_simulated = s$recall_simulated,
        f1 = s$f1, calls_per_tp = s$calls_per_tp, n_testable = n_testable)
}

for (ss in subsets) {
    # Union of what ANY method in the subset could have found. Using one method's
    # testable set would misreport the subset's recall.
    testable_ss <- unique(unlist(lapply(ss, function(m) loaded[[m]]$testable)))
    cands <- lapply(ss, function(m) CAND[method == m])
    names(cands) <- ss

    # Assignment enumeration is bounded by subset size, not by taste. For pairs and
    # triples the full cross-product is <= 3^3 = 27 and is enumerated. For the
    # all-methods bundle it would be 3^11 = 177 147, so only the three COHERENT
    # assignments are used: every method at its best-F1 point, every method at its
    # high-precision point, every method at the pipeline default. Those are the three
    # a user could actually state as a policy; the 177 145 mixtures in between are not
    # a recommendation anyone could follow.
    grid <- if (length(ss) <= 3) {
        do.call(expand.grid, c(lapply(cands, function(d) seq_len(nrow(d))),
                               list(KEEP.OUT.ATTRS = FALSE)))
    } else {
        roles <- unique(unlist(lapply(cands, `[[`, "role")))
        as.data.table(do.call(rbind, lapply(roles, function(rl)
            vapply(ss, function(m) {
                i <- which(cands[[m]]$role == rl)
                if (length(i)) i[1] else which(cands[[m]]$role == "default")[1]
            }, integer(1)))))
    }

    for (g in seq_len(nrow(grid))) {
        parts <- list(); lbl <- character(0)
        for (j in seq_along(ss)) {
            m <- ss[j]; r <- cands[[m]][as.integer(grid[[j]][g])]
            parts[[m]] <- CALLED[[paste0(m, "|", r$adjust, "_", r$value)]]
            lbl <- c(lbl, sprintf("%s:%s_%s", m, r$adjust, r$value))
        }
        assign_lbl <- paste(lbl, collapse = ";")

        for (cw in COMBINE_WINDOWS) {
            support <- combine_support(parts, cw * 1000)
            ks <- if (length(ss) > 3) c(1L, 2L, 3L, 4L, length(ss)) else seq_len(length(ss))
            for (k in unique(ks)) {
                sel  <- support[n_methods >= k, snp]
                rule <- if (k == 1L) "union" else if (k == length(ss)) "intersection"
                        else paste0("at_least_", k)
                add(ss, assign_lbl, rule, k, cw,
                    score_calls(sel, WINDOW, truth, testable_ss), length(testable_ss))
            }
        }
    }

    # ---- rank_sum: threshold-free, so it has no assignment and no window. Kept because
    # it is the only combine rule that does not inherit RDA's uncalibrated p-values --
    # journal 05 sec.10 showed any intersection containing RDA degenerates to zero.
    rt <- rbindlist(lapply(ss, function(m) {
        r <- loaded[[m]]$rank$keys
        data.table(snp = r, rank = seq_along(r), method = m)
    }))
    worst <- rt[, .(worst = max(rank)), by = method]
    wide  <- dcast(rt, snp ~ method, value.var = "rank")
    for (m in ss) set(wide, which(is.na(wide[[m]])), m, worst[method == m, worst] + 1L)
    wide[, rank_sum := rowSums(.SD), .SDcols = ss]
    setorder(wide, rank_sum)
    for (n in TOP_GRID) {
        if (n > nrow(wide)) next
        add(ss, "rank_sum", "rank_sum", NA_integer_, NA_real_,
            score_calls(wide$snp[seq_len(n)], WINDOW, truth, testable_ss), length(testable_ss))
        rows[[length(rows)]][, n_called_target := n]
    }
}

# ---- solo rows at each candidate point, so "best single method" and "best combination"
# are compared inside ONE table rather than across two files with different denominators.
for (m in names(loaded)) {
    d <- CAND[method == m]
    for (i in seq_len(nrow(d))) {
        add(m, sprintf("%s:%s_%s", m, d$adjust[i], d$value[i]), "solo", 1L, NA_real_,
            score_calls(CALLED[[paste0(m, "|", d$adjust[i], "_", d$value[i])]],
                        WINDOW, truth, loaded[[m]]$testable),
            length(loaded[[m]]$testable))
    }
}

out <- rbindlist(rows, fill = TRUE)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(OUTDIR, paste0(TAG, "_portfolio.tsv"))
fwrite(out, out_f, sep = "\t")
message("INFO: wrote ", out_f, " (", nrow(out), " operating points)")

message("\n=== best F1 overall (ties broken by fewer calls) ===")
print(head(out[!is.na(f1)][order(-f1, n_called),
    .(subset, rule, cw = combine_window_kb, assignment, n_called, tp,
      precision = round(precision, 3), recall = round(recall_testable, 3),
      f1 = round(f1, 3), per_tp = round(calls_per_tp, 1))], 15))

message("\n=== best F1 per rule ===")
print(out[!is.na(f1)][order(-f1, n_called), .SD[1], by = rule][
    , .(rule, subset, cw = combine_window_kb, n_called, tp,
        precision = round(precision, 3), recall = round(recall_testable, 3),
        f1 = round(f1, 3))])

# =============================================================================
# TRANSFER: fixed schemes imported from another replicate
# =============================================================================
# Everything above is an IN-SAMPLE search -- subset, assignment and window are all chosen on
# the same causal loci the resulting F1 is then reported on. That number cannot answer "does
# combining generalise"; it can only answer "what is the best combination here".
#
# This block answers the first question by scoring schemes that were fixed on a DIFFERENT
# replicate. Two modes per scheme, and the difference between them is the point:
#
#   verbatim  the anchor's subset, rule, window AND per-method operating points, unchanged.
#             The strict test. Note that an operating point is not a portable object:
#             `top_100` selects 100 SNPs whether the replicate has 43 causal loci or 395, so
#             a verbatim loss on a polygenic seed is expected and is not evidence about the
#             combine structure.
#   recal     the anchor's subset, rule and window, but each method's operating point re-
#             picked from THIS replicate's own candidate set. This isolates the transferable
#             claim -- that the STRUCTURE of the combination generalises -- from the
#             calibration, which demonstrably does not.
#
# Neither mode chooses the subset on this replicate, so neither inherits the search bias that
# {TAG}_portfolio.tsv carries by construction.
if (nzchar(SCHEMES_F)) {
    if (!file.exists(SCHEMES_F)) stop("--schemes given but not found: ", SCHEMES_F)
    sch <- fread(SCHEMES_F)
    for (cn in c("scheme_id", "subset", "rule", "combine_window_kb", "assignment"))
        if (!cn %in% names(sch)) stop("--schemes table has no `", cn, "` column")

    called_of <- function(m, a, v) {
        k <- paste0(m, "|", a, "_", v)
        if (is.null(CALLED[[k]]))
            CALLED[[k]] <<- call_by_threshold(loaded[[m]]$pv, loaded[[m]]$trait_cols,
                                              a, v, quiet = TRUE)$called
        CALLED[[k]]
    }
    # An operating point is written "<adjust>_<value>"; split at the FIRST underscore only,
    # or `custom_1e-04` loses its exponent and `at_least_2` would not round-trip.
    split_point <- function(s) {
        i <- regexpr("_", s, fixed = TRUE)
        if (i < 1L) stop("Bad operating point (expected <adjust>_<value>): ", s)
        c(substr(s, 1L, i - 1L), substr(s, i + 1L, nchar(s)))
    }
    parse_assignment <- function(a, ss) {
        out <- list()
        for (t in strsplit(a, ";", fixed = TRUE)[[1]]) {
            kv <- strsplit(t, ":", fixed = TRUE)[[1]]
            if (length(kv) != 2L) stop("Bad assignment token: ", t)
            pt <- split_point(kv[2])
            if (kv[1] == "*") for (m in ss) out[[m]] <- pt else out[[kv[1]]] <- pt
        }
        if (length(setdiff(ss, names(out))))
            stop("Assignment does not cover ", paste(setdiff(ss, names(out)), collapse = ","))
        out[ss]
    }
    rule_k <- function(rule, n) {
        if (rule %in% c("solo", "union")) 1L
        else if (rule == "intersection") as.integer(n)
        else if (grepl("^at_least_[0-9]+$", rule)) as.integer(sub("at_least_", "", rule))
        else stop("Unknown rule: ", rule)
    }
    eval_scheme <- function(ss, pts, rule, cw) {
        # Union of what ANY member could have found -- using one method's testable set would
        # misreport the subset's recall, same convention as the main enumeration above.
        testable_ss <- unique(unlist(lapply(ss, function(m) loaded[[m]]$testable)))
        parts <- lapply(ss, function(m) called_of(m, pts[[m]][1], pts[[m]][2]))
        names(parts) <- ss
        sel <- if (length(ss) == 1L) parts[[1]] else {
            combine_support(parts, cw * 1000)[n_methods >= rule_k(rule, length(ss)), snp]
        }
        list(s = score_calls(sel, WINDOW, truth, testable_ss), n_testable = length(testable_ss))
    }
    # This replicate's own candidate points, enumerated exactly as the main search does:
    # full cross-product for <= 3 methods, the three COHERENT role assignments beyond that.
    recal_grid <- function(ss) {
        cands <- lapply(ss, function(m) CAND[method == m]); names(cands) <- ss
        if (length(ss) <= 3L) {
            g <- do.call(expand.grid, c(lapply(cands, function(d) seq_len(nrow(d))),
                                        list(KEEP.OUT.ATTRS = FALSE)))
            lapply(seq_len(nrow(g)), function(i) {
                p <- lapply(ss, function(m) {
                    r <- cands[[m]][as.integer(g[[m]][i])]; c(r$adjust, r$value) })
                names(p) <- ss; p })
        } else {
            roles <- unique(unlist(lapply(cands, `[[`, "role")))
            lapply(roles, function(rl) {
                p <- lapply(ss, function(m) {
                    d <- cands[[m]]
                    j <- which(d$role == rl)
                    # Falling back to "default" is NOT safe: pick_candidates() dedups by
                    # (adjust, value), so a method whose high-precision point happens to BE
                    # bonf_0.05 loses its `default` row to the dedup and keeps only
                    # `hp_recall`. The old fallback then indexed with integer(0) -> NA, and
                    # an NA operating point reaches compute_pval_threshold(), which halts the
                    # whole script -- taking that cell's entire transfer table with it. Seen
                    # on MVP1231858 c1/c5, where MLM/CMLM/ECMLM/MLMM all lack a default row.
                    if (!length(j)) j <- which(d$role == "default")
                    if (!length(j)) j <- 1L
                    c(d$adjust[j[1]], d$value[j[1]]) })
                names(p) <- ss; p })
        }
    }
    lbl_of <- function(ss, pts) paste(vapply(ss, function(m)
        sprintf("%s:%s_%s", m, pts[[m]][1], pts[[m]][2]), character(1)), collapse = ";")

    trows <- list()
    tadd <- function(row, ss, rule, cw, mode, pts, r) {
        trows[[length(trows) + 1L]] <<- data.table(
            tag = TAG, scheme_id = row$scheme_id,
            label = if ("label" %in% names(row)) row$label else NA_character_,
            subset = paste(ss, collapse = "+"), n_methods = length(ss),
            rule = rule, combine_window_kb = cw, mode = mode,
            assignment = lbl_of(ss, pts),
            n_called = r$s$n_called, tp = r$s$tp, fp_background = r$s$fp_background,
            precision = r$s$precision_strict, precision_all = r$s$precision_all,
            recall_testable = r$s$recall_testable, f1 = r$s$f1,
            calls_per_tp = r$s$calls_per_tp, n_testable = r$n_testable)
    }

    for (i in seq_len(nrow(sch))) {
        row <- sch[i]
        ss  <- if (row$subset == "ALL") names(loaded)
               else strsplit(row$subset, "+", fixed = TRUE)[[1]]
        miss <- setdiff(ss, names(loaded))
        if (length(miss)) {
            message("SKIP ", row$scheme_id, ": method(s) not available here: ",
                    paste(miss, collapse = ","))
            next
        }
        cw <- as.numeric(row$combine_window_kb)

        pts <- tryCatch(parse_assignment(row$assignment, ss), error = function(e) e)
        if (inherits(pts, "error")) {
            message("SKIP ", row$scheme_id, ": ", conditionMessage(pts)); next
        }
        tadd(row, ss, row$rule, cw, "verbatim", pts, eval_scheme(ss, pts, row$rule, cw))

        # recalibrated: same structure, this replicate's own points, best F1 (ties -> fewer calls)
        best <- NULL; best_pts <- NULL
        for (p in recal_grid(ss)) {
            r <- eval_scheme(ss, p, row$rule, cw)
            if (is.na(r$s$f1)) next
            if (is.null(best) || r$s$f1 > best$s$f1 ||
                (r$s$f1 == best$s$f1 && r$s$n_called < best$s$n_called)) {
                best <- r; best_pts <- p
            }
        }
        if (!is.null(best)) tadd(row, ss, row$rule, cw, "recal", best_pts, best)
    }

    if (length(trows)) {
        tout <- rbindlist(trows, fill = TRUE)
        tf <- file.path(OUTDIR, paste0(TAG, "_transfer.tsv"))
        fwrite(tout, tf, sep = "\t")
        message("INFO: wrote ", tf, " (", nrow(tout), " transferred scheme evaluations)")
        message("\n=== transferred schemes ===")
        print(tout[, .(scheme_id, mode, n_called, tp, precision = round(precision, 3),
                       recall = round(recall_testable, 3), f1 = round(f1, 3))])
    } else message("WARNING: --schemes given but no scheme was evaluable here")
}
