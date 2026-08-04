#!/usr/bin/env Rscript
# =============================================================================
# mvp_crossseed_report.R -- fold 14 replicates x 2 MAF arms x 5 structure-correction rungs
# into the tables journal 07 reports.
#
#   Rscript mvp_crossseed_report.R --seeds=CSV [--cells=5] [--outdir=DIR]
#
# WHAT THIS IS FOR, AND THE ONE DISTINCTION EVERYTHING RESTS ON
# ------------------------------------------------------------
# Journal 06 measured, on ONE replicate, that combining methods beats the best single method.
# That measurement is IN-SAMPLE on both sides: the winning subset, the per-method operating
# points and the agreement window were all chosen on the same causal loci the resulting F1 was
# then reported on, and so was the "best solo" it was compared against. With ~30 subsets x up
# to 27 assignments x 4 windows searched against 43 causal loci, some of that gain is
# selection, and nothing inside one replicate can say how much.
#
# So this script keeps two families of numbers strictly apart and never averages them:
#
#   TRANSFER   (portfolio07*/..._transfer.tsv) -- schemes fixed on the ANCHOR replicate
#              (1232548, journal 06) and applied unchanged to replicates that had no part in
#              choosing them. This is the only family that can support a general claim.
#   RE-SEARCH  (portfolio07*/..._portfolio.tsv) -- the full search repeated per replicate.
#              This is the per-dataset recommendation, and it is what a user would actually
#              do; it is NOT evidence that combining generalises.
#
# Their difference, per replicate, IS the selection bias -- reported as its own table rather
# than buried, because it is the honest error bar on journal 06's headline.
#
# WHY ABSOLUTE F1 IS NOT COMPARED ACROSS SEEDS
# --------------------------------------------
# TOP_MAX is per seed (benchmarks/mvp_seed_topmax.tsv): 100 where the replicate has tens of
# causal loci, 1000-2000 where it has hundreds. Without that, a fixed cap of 100 would hold
# recall below 0.19 on the polygenic replicates by construction and the comparison would
# measure the cap. The cost is that F1 levels are not commensurable between seeds. Every
# cross-seed statement here is therefore a WITHIN-seed contrast -- a delta against that same
# seed's own baseline, or a rank -- never a pooled mean of F1.
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt  <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]

SEEDS  <- strsplit(opt("seeds", ""), ",", fixed = TRUE)[[1]]
if (!length(SEEDS)) stop("Missing required --seeds")
CELLS  <- as.integer(opt("cells", "5"))
OUTDIR <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/report07"))
# --eval-dir exists so this script can be exercised against a fixture tree before the real
# 28-hour run finishes. Defaults to the real location; the driver never passes it.
EVAL   <- opt("eval_dir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ARMS <- data.table(arm = c("base", "m05"), tag = c("", "_m05"),
                   sfx = c("", "M05"), maf = c("0.01", "0.05"))

MAN <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
             colClasses = c(seed = "character"))
TM  <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seed_topmax.tsv"),
             colClasses = c(seed = "character"))
STRATA <- MAN[, .(seed, arch = arch_level, demog = demog_level, n_traits,
                  k_best, meanFst, r2_pc1_temp, n_causal_maf01)]
STRATA <- merge(STRATA, TM[, .(seed, top_max)], by = "seed", all.x = TRUE)
# The two Est-Clines replicates are the framework's own falsification control: R2(PC1~temp)
# 0.92-0.94, i.e. structure and the predictor are nearly the same variable. MVP_README
# predicts they collapse once structure correction is switched on. They are labelled here and
# excluded from every cross-seed vote, because they are not samples from the same population
# of scenarios as the SS-Mtn replicates.
STRATA[, control := demog == "Est-Clines"]
# Replicates with under ~15 causal loci cannot carry a "scheme A beats scheme B" statement:
# one locus moves F1 by more than the differences being compared. Flagged, reported, and kept
# out of the votes -- not silently dropped.
STRATA[, underpowered := n_causal_maf01 < 15]

rd <- function(f) if (file.exists(f)) fread(f) else NULL

# ------------------------------------------------------------------ collect
transfer <- list(); research <- list(); solo <- list(); redun <- list()
for (a in seq_len(nrow(ARMS))) {
    tag <- ARMS$tag[a]; sfx <- ARMS$sfx[a]; arm <- ARMS$arm[a]
    for (s in SEEDS) {
        proj <- paste0("MVP", s, sfx)
        r <- rd(file.path(EVAL, paste0("redundancy07", tag),
                          paste0(proj, "_c3_method_clusters.tsv")))
        if (!is.null(r)) redun[[length(redun) + 1L]] <- cbind(seed = s, arm = arm, r)
        for (i in seq_len(CELLS)) {
            base <- file.path(EVAL, paste0("portfolio07", tag))
            tf <- rd(file.path(base, sprintf("%s_c%d_any_transfer.tsv", proj, i)))
            pf <- rd(file.path(base, sprintf("%s_c%d_any_portfolio.tsv", proj, i)))
            if (!is.null(tf))
                transfer[[length(transfer) + 1L]] <- cbind(seed = s, arm = arm, cell = i, tf)
            if (is.null(pf)) next
            # In-sample winner and in-sample best solo, from the SAME table, so the two share
            # a denominator. Ties broken toward fewer calls, as in journal 06.
            cmb <- pf[rule != "solo" & !is.na(f1)][order(-f1, n_called)]
            slo <- pf[rule == "solo"  & !is.na(f1)][order(-f1, n_called)]
            if (nrow(cmb)) research[[length(research) + 1L]] <-
                cbind(seed = s, arm = arm, cell = i, cmb[1])
            if (nrow(slo)) solo[[length(solo) + 1L]] <-
                cbind(seed = s, arm = arm, cell = i, slo[1])
        }
    }
}
TR <- if (length(transfer)) rbindlist(transfer, fill = TRUE) else NULL
RS <- if (length(research)) rbindlist(research, fill = TRUE) else NULL
SO <- if (length(solo))     rbindlist(solo,     fill = TRUE) else NULL
RD <- if (length(redun))    rbindlist(redun,    fill = TRUE) else NULL
if (is.null(TR) || is.null(RS)) stop("Nothing to aggregate -- run the score stage first")

message(sprintf("INFO: %d transfer rows, %d re-search rows, %d/%d seeds present",
                nrow(TR), nrow(RS), uniqueN(RS$seed), length(SEEDS)))

# ------------------------------------------- 1. transfer, with its two baselines
# Two baselines per (seed, arm, cell), and they answer different questions:
#   best_solo_f1     the best single method AFTER searching this seed's own grid. The
#                    conservative comparison -- an in-sample baseline against an out-of-sample
#                    scheme, so any positive delta is a lower bound.
#   default_solo_f1  LFMM at bonf 0.05, scheme B1. What the pipeline gives a user today with
#                    one method and no tuning. The user-facing comparison.
BASE <- merge(
    SO[, .(seed, arm, cell, best_solo_f1 = f1, best_solo_subset = subset,
           best_solo_called = n_called)],
    TR[scheme_id == "B1" & mode == "verbatim",
       .(seed, arm, cell, default_solo_f1 = f1, default_solo_called = n_called)],
    by = c("seed", "arm", "cell"), all = TRUE)

TR <- merge(TR, BASE, by = c("seed", "arm", "cell"), all.x = TRUE)
TR <- merge(TR, STRATA, by = "seed", all.x = TRUE)
# A baseline that calls nothing (or calls without a single true positive) has recall 0, hence
# F1 0 -- score_calls returns NaN because precision is undefined at n_called = 0. Scoring it as
# NA instead would DROP those cells, and they are not random: 26 of them, concentrated on the
# highly-polygenic replicates where LFMM at bonf 0.05 calls 0-1 SNPs. Dropping them would
# delete exactly the cases where the pipeline default performs worst and bias every "gain over
# the default" figure in the default's favour.
TR[is.na(default_solo_f1), default_solo_f1 := 0]
TR[is.na(best_solo_f1),    best_solo_f1    := 0]
TR[, `:=`(gain_vs_best_solo    = f1 - best_solo_f1,
          gain_vs_default_solo = f1 - default_solo_f1)]
setorder(TR, arm, seed, cell, scheme_id, mode)
fwrite(TR, file.path(OUTDIR, "transfer.tsv"), sep = "\t")

# --------------------------------------- 2. does a fixed scheme transfer, and where
# Aggregated over cells WITHIN a seed first (median across the 5 rungs), so a seed whose
# ladder happens to have more cells cannot outvote one that has fewer.
TRS <- TR[!(scheme_id %in% c("B1")) & !is.na(f1),
          .(f1 = median(f1, na.rm = TRUE),
            gain_best = median(gain_vs_best_solo, na.rm = TRUE),
            gain_default = median(gain_vs_default_solo, na.rm = TRUE),
            n_called = median(n_called, na.rm = TRUE),
            precision = median(precision, na.rm = TRUE),
            recall = median(recall_testable, na.rm = TRUE)),
          by = .(scheme_id, label, mode, arm, seed, arch, control, underpowered)]
fwrite(TRS, file.path(OUTDIR, "transfer_by_seed.tsv"), sep = "\t")

vote_pool <- TRS[control == FALSE & underpowered == FALSE]
# n_seeds vs n_seeds_with_baseline: `gain_default` is NA whenever scheme B1 (LFMM @ bonf 0.05)
# did not land for that seed -- which happens if LFMM is a declared exclusion in
# mvp_method_exclusions.tsv. The `n_beat_*` counters use na.rm, so without the second column a
# missing baseline would deflate the numerator while uniqueN(seed) kept the denominator full.
# (That table is currently empty, so the two columns should be equal; if they ever differ, the
# beat-counts are being read against the wrong denominator.)
SUM <- vote_pool[, .(
    n_seeds      = uniqueN(seed),
    n_seeds_with_baseline = uniqueN(seed[!is.na(gain_default)]),
    n_beat_best  = sum(gain_best > 0, na.rm = TRUE),
    n_beat_deflt = sum(gain_default > 0, na.rm = TRUE),
    med_gain_best    = round(median(gain_best, na.rm = TRUE), 4),
    med_gain_default = round(median(gain_default, na.rm = TRUE), 4),
    med_f1       = round(median(f1, na.rm = TRUE), 4),
    med_called   = round(median(n_called, na.rm = TRUE), 1),
    med_precision = round(median(precision, na.rm = TRUE), 3)),
    by = .(scheme_id, label, mode, arm)]
setorder(SUM, arm, mode, -med_gain_default)
fwrite(SUM, file.path(OUTDIR, "transfer_summary.tsv"), sep = "\t")

# per architecture stratum -- the deliverable the user asked for: which scheme for which data
by_arch <- function(pool) pool[, .(
    n_seeds = uniqueN(seed),
    n_beat_default = sum(gain_default > 0, na.rm = TRUE),
    med_gain_default = round(median(gain_default, na.rm = TRUE), 4),
    med_f1 = round(median(f1, na.rm = TRUE), 4)),
    by = .(arch, arm, mode, scheme_id, label)]

BYARCH <- by_arch(vote_pool)
setorder(BYARCH, arch, arm, mode, -med_gain_default)
fwrite(BYARCH, file.path(OUTDIR, "transfer_by_arch.tsv"), sep = "\t")

# Every one of the four oligogenic replicates has under 15 causal loci, so `vote_pool` -- which
# excludes underpowered seeds -- contains NO oligogenic stratum at all. Reporting nothing for
# it would read as "no scheme works there" rather than "this benchmark cannot tell". So the
# stratum table is emitted a SECOND time with the underpowered filter lifted, carrying the flag
# so the two are never confused. Nothing from this table feeds a recommendation.
BYARCH_ALL <- by_arch(TRS[control == FALSE])
BYARCH_ALL <- merge(BYARCH_ALL,
                    unique(TRS[control == FALSE, .(arch, underpowered)])[
                        , .(any_underpowered = any(underpowered)), by = arch],
                    by = "arch", all.x = TRUE)
setorder(BYARCH_ALL, arch, arm, mode, -med_gain_default)
fwrite(BYARCH_ALL, file.path(OUTDIR, "transfer_by_arch_incl_underpowered.tsv"), sep = "\t")

# the recommendation table: best transferable scheme per (stratum, arm), recal mode.
# seq_len(min(3, .N)) not 1:3 -- `.SD[1:3]` pads a short group with NA rows, which would put
# NA-filled lines into a table whose entire purpose is to be acted on.
BEST_BY_ARCH <- BYARCH[mode == "recal"][order(arch, arm, -med_gain_default),
                                        .SD[seq_len(min(3L, .N))], by = .(arch, arm)]
fwrite(BEST_BY_ARCH, file.path(OUTDIR, "recommended_schemes.tsv"), sep = "\t")

# ------------------------------------------------ 3. re-search, and selection bias
RS2 <- merge(RS, BASE, by = c("seed", "arm", "cell"), all.x = TRUE)
RS2 <- merge(RS2, STRATA, by = "seed", all.x = TRUE)
RS2[, `:=`(gain_vs_best_solo = f1 - best_solo_f1,
           gain_vs_default_solo = f1 - default_solo_f1)]
setorder(RS2, arm, seed, cell)
fwrite(RS2, file.path(OUTDIR, "research.tsv"), sep = "\t")

# Selection bias = (best combination found BY searching this seed) - (best ANCHOR scheme,
# recalibrated on this seed). Both are combinations, both use this seed's own calibration;
# the only difference is whether the SUBSET/RULE/WINDOW was chosen here or imported. So the
# gap is what the search bought, and equivalently what journal 06's headline over-states.
TR_BEST <- TR[mode == "recal" & scheme_id != "B1" & !is.na(f1),
              .(transfer_best_f1 = max(f1),
                transfer_best_scheme = scheme_id[which.max(f1)]),
              by = .(seed, arm, cell)]
BIAS <- merge(RS2[, .(seed, arm, cell, arch, control, underpowered,
                      research_f1 = f1, research_subset = subset, research_rule = rule,
                      best_solo_f1, default_solo_f1)],
              TR_BEST, by = c("seed", "arm", "cell"), all.x = TRUE)
BIAS[, selection_bias := research_f1 - transfer_best_f1]
setorder(BIAS, arm, seed, cell)
fwrite(BIAS, file.path(OUTDIR, "selection_bias.tsv"), sep = "\t")

# -------------------------------------------------------- 4. the MAF arm contrast
# Paired WITHIN seed and rung. Recall levels are not comparable between arms (different
# denominators by construction: fewer causal loci survive MAF 0.05), so the columns that
# carry the argument are precision, calls_per_tp and the threshold-free AUC-PR below.
ARMCMP <- dcast(RS2[, .(seed, arm, cell, arch, control, f1, precision, n_called,
                        calls_per_tp, recall_testable)],
                seed + cell + arch + control ~ arm,
                value.var = c("f1", "precision", "n_called", "calls_per_tp", "recall_testable"))
if (all(c("f1_base", "f1_m05") %in% names(ARMCMP))) {
    ARMCMP[, `:=`(d_f1 = f1_m05 - f1_base,
                  d_precision = precision_m05 - precision_base,
                  d_calls_per_tp = calls_per_tp_m05 - calls_per_tp_base)]
}
setorder(ARMCMP, seed, cell)
fwrite(ARMCMP, file.path(OUTDIR, "arm_compare.tsv"), sep = "\t")

# AUC-PR is threshold-free, so it is the one quantity the softer Bonferroni denominator under
# MAF 0.05 (n_tested falls, so the same nominal level is a laxer absolute threshold) cannot
# explain away. Journal 06 found it rose for all 11 methods on the anchor.
if (!is.null(RD) && "auc_pr" %in% names(RD)) {
    AUC <- dcast(RD[, .(seed, arm, method, auc_pr)], seed + method ~ arm, value.var = "auc_pr")
    if (all(c("base", "m05") %in% names(AUC))) {
        AUC[, ratio := m05 / base]
        setnames(AUC, c("base", "m05"), c("aucpr_maf01", "aucpr_maf05"))
    }
    AUC <- merge(AUC, STRATA[, .(seed, arch, control)], by = "seed", all.x = TRUE)
    setorder(AUC, method, seed)
    fwrite(AUC, file.path(OUTDIR, "aucpr_by_arm.tsv"), sep = "\t")

    RANK <- RD[, .(med_aucpr = median(auc_pr, na.rm = TRUE),
                   n_seeds = uniqueN(seed)), by = .(arm, method)]
    RANK[, rank := frank(-med_aucpr, ties.method = "min"), by = arm]
    setorder(RANK, arm, rank)
    fwrite(RANK, file.path(OUTDIR, "method_rank.tsv"), sep = "\t")
}

# --------------------------------------------------------------------- 5. plots
theme_j07 <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))

p1 <- ggplot(vote_pool[mode == "recal"],
             aes(x = reorder(scheme_id, gain_default, median), y = gain_default,
                 fill = arm)) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.85) +
    coord_flip() + facet_wrap(~ arch) +
    labs(x = NULL, y = "F1(scheme) - F1(LFMM @ bonf 0.05), per seed",
         fill = "Filter.maf") + theme_j07
ggsave(file.path(OUTDIR, "transfer_gain.png"), p1, width = 10, height = 5, dpi = 150)
ggsave(file.path(OUTDIR, "transfer_gain.svg"), p1, width = 10, height = 5)

if (nrow(BIAS[!is.na(selection_bias)])) {
    p2 <- ggplot(BIAS[control == FALSE & underpowered == FALSE & !is.na(selection_bias)],
                 aes(x = arm, y = selection_bias, fill = arm)) +
        geom_hline(yintercept = 0, linewidth = 0.4) +
        geom_boxplot(outlier.size = 0.7, alpha = 0.85, show.legend = FALSE) +
        facet_wrap(~ arch) +
        labs(x = NULL, y = "F1(searched here) - F1(best imported scheme)") + theme_j07
    ggsave(file.path(OUTDIR, "selection_bias.png"), p2, width = 8, height = 4, dpi = 150)
    ggsave(file.path(OUTDIR, "selection_bias.svg"), p2, width = 8, height = 4)
}

if (exists("AUC") && "ratio" %in% names(AUC)) {
    p3 <- ggplot(AUC[control == FALSE & is.finite(ratio)],
                 aes(x = reorder(method, ratio, median), y = ratio)) +
        geom_hline(yintercept = 1, linewidth = 0.4) +
        geom_boxplot(outlier.size = 0.7, alpha = 0.85) +
        coord_flip() + scale_y_log10() +
        labs(x = NULL, y = "AUC-PR(maf 0.05) / AUC-PR(maf 0.01), per seed") + theme_j07
    ggsave(file.path(OUTDIR, "aucpr_arm_ratio.png"), p3, width = 7, height = 5, dpi = 150)
    ggsave(file.path(OUTDIR, "aucpr_arm_ratio.svg"), p3, width = 7, height = 5)
}

# ------------------------------------------------------------------- 6. console
message("\n=== transferred schemes, pooled over non-control, adequately powered seeds ===")
print(SUM[mode == "recal"][order(arm, -med_gain_default)])

message("\n=== recommended scheme per architecture stratum (recal mode) ===")
print(BEST_BY_ARCH[, .(arch, arm, scheme_id, n_seeds, n_beat_default,
                       med_gain_default, med_f1)])

message("\n=== selection bias: what searching this seed bought over importing a scheme ===")
print(BIAS[control == FALSE & underpowered == FALSE,
           .(median_bias = round(median(selection_bias, na.rm = TRUE), 4),
             q25 = round(quantile(selection_bias, .25, na.rm = TRUE), 4),
             q75 = round(quantile(selection_bias, .75, na.rm = TRUE), 4),
             n = .N), by = .(arm, arch)])

message("\nINFO: wrote ", OUTDIR)
