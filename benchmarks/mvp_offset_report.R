#!/usr/bin/env Rscript
# =============================================================================
# mvp_offset_report.R -- fold benchmarks/eval_offset.R's per-(seed, method, set) scores into
# the tables and figures journal 08 reports.
#
# THE ONE DISTINCTION EVERYTHING RESTS ON: a SNP set here is one of three things, and they
# are never pooled --
#   candidate  best / union / solo   -- the things a user could actually build
#   ceiling    truth                 -- the causal loci; what the candidates are chasing
#   floor      rand_*                -- size-matched background-neutral draws; the null
#
# Every candidate claim is therefore reported as a WITHIN-SEED contrast against that seed's
# own ceiling and floor, never as a pooled absolute number: rho levels differ across seeds
# with architecture (6 to 395 causal loci) and with the strength of the climate-distance
# null (measured 0.016-0.738 across the 14 scenarios), so pooled means would mostly measure
# which seeds are in the pool.
#
# Est-Clines control seeds (R^2(PC1~temp) 0.92-0.94) are excluded from every aggregate and
# reported in their own table, exactly as journals 06-07 do.
#
# Usage:  Rscript mvp_offset_report.R [--outdir=DIR]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))   # parse_kv_args()

args   <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt    <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
OUTDIR <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset08"))

SC <- fread(file.path(OUTDIR, "scores.tsv"), colClasses = c("seed" = "character"))
TS <- fread(file.path(OUTDIR, "target_summary.tsv"), colClasses = c("seed" = "character"))
ID <- fread(file.path(OUTDIR, "fitness_identity.tsv"), colClasses = c("seed" = "character"))
SC <- merge(SC, TS[, .(seed, null_r2_scenario = null_r2_clim_dist,
                       oracle_r2 = oracle_r2_axis_weighted, headroom,
                       delta_logf_med, n_demes_gaining)], by = "seed", all.x = TRUE)

# Selection regime, read off the fitness identity's own coefficients rather than from the
# architecture label: it is what decides whether the unweighted Euclidean null is already the
# right metric. On equal-S replicates (b1 = b2) it IS right, so no SNP set can beat it and the
# replicate carries no information about set choice -- pooling those in would dilute the very
# contrast the journal is measuring.
SC <- merge(SC, ID[, .(seed, b1, b2)], by = "seed", all.x = TRUE)
SC[, regime := fifelse(abs(b2) < 1e-5, "1-trait",
               fifelse(abs(abs(b1) - abs(b2)) < 1e-6, "equal-S", "unequal-S"))]
SC[, discriminating := headroom > 0.05]

SET_ORDER <- c("best", "union", "solo", "truth",
               "rand_best1", "rand_best2", "rand_best3", "rand_union1", "rand_solo1")
SC[, set := factor(set, levels = intersect(SET_ORDER, unique(set)))]
SC[, method := factor(method, levels = c("gradient_forest", "geometric_offset", "rda_offset"))]

MAIN <- SC[control == FALSE]
CTRL <- SC[control == TRUE]

# ---------------------------------------------------------------- 1. by set x method
# Primary aggregate is over the DISCRIMINATING replicates only (headroom > 0.05). The
# equal-S replicates are reported separately in by_regime.tsv, not silently dropped.
BY_SET_ALL <- MAIN[, .(n_seeds = uniqueN(seed), rho_med = median(rho),
                       delta_r2_med = median(delta_r2)),
                   by = .(method, set, role)]
fwrite(BY_SET_ALL, file.path(OUTDIR, "by_set_all_seeds.tsv"), sep = "\t")

BY_REGIME <- MAIN[, .(n_seeds = uniqueN(seed), headroom_med = median(headroom),
                      rho_med = median(rho), delta_r2_med = median(delta_r2),
                      partial_rho_med = median(partial_rho)),
                  by = .(regime, method, set, role)]
setorder(BY_REGIME, regime, method, -rho_med)
fwrite(BY_REGIME, file.path(OUTDIR, "by_regime.tsv"), sep = "\t")

MAIN <- MAIN[discriminating == TRUE]
BY_SET <- MAIN[, .(n_seeds       = uniqueN(seed),
                   rho_med       = median(rho),
                   rho_q25       = quantile(rho, 0.25),
                   rho_q75       = quantile(rho, 0.75),
                   delta_r2_med  = median(delta_r2),
                   partial_rho_med = median(partial_rho),
                   n_snps_med    = median(n_snps, na.rm = TRUE)),
               by = .(method, set, role)]
setorder(BY_SET, method, -rho_med)
fwrite(BY_SET, file.path(OUTDIR, "by_set.tsv"), sep = "\t")

# ------------------------------------------------- 2. candidate vs its own floor/ceiling
# Floor for a candidate = the median rho of the random draws matched to ITS size:
#   best -> rand_best1..3, union -> rand_union1, solo -> rand_solo1
FLOOR_MAP <- list(best = c("rand_best1", "rand_best2", "rand_best3"),
                  union = "rand_union1", solo = "rand_solo1")

cand <- MAIN[role == "candidate"]
gap_rows <- list()
for (i in seq_len(nrow(cand))) {
    r  <- cand[i]
    fl <- MAIN[seed == r$seed & method == r$method & set %in% FLOOR_MAP[[as.character(r$set)]], rho]
    ce <- MAIN[seed == r$seed & method == r$method & set == "truth", rho]
    gap_rows[[i]] <- data.table(
        seed = r$seed, arch = r$arch, method = r$method, set = r$set,
        n_snps = r$n_snps, n_causal_in_set = r$n_causal, n_causal_total = r$n_causal_maf01,
        rho = r$rho, rho_floor = if (length(fl)) median(fl) else NA_real_,
        rho_ceiling = if (length(ce)) ce[1] else NA_real_,
        delta_r2 = r$delta_r2,
        null_r2_scenario = r$null_r2_scenario)
}
GAP <- rbindlist(gap_rows)
GAP[, beats_floor   := rho > rho_floor]
GAP[, gap_to_ceiling := rho_ceiling - rho]
fwrite(GAP, file.path(OUTDIR, "candidate_vs_refs.tsv"), sep = "\t")

WIN <- GAP[, .(n_seeds = .N,
               n_beat_floor = sum(beats_floor, na.rm = TRUE),
               rho_med = median(rho),
               floor_med = median(rho_floor, na.rm = TRUE),
               ceiling_med = median(rho_ceiling, na.rm = TRUE),
               gap_to_ceiling_med = median(gap_to_ceiling, na.rm = TRUE)),
           by = .(method, set)]
setorder(WIN, method, -rho_med)
fwrite(WIN, file.path(OUTDIR, "candidate_win_table.tsv"), sep = "\t")

# ------------------------------------------- 3. by architecture (discriminating seeds only)
# MAIN was narrowed above, so this table -- like by_set.tsv and candidate_win_table.tsv --
# covers only replicates where an axis weighting can beat the null. by_regime.tsv is the one
# table that spans all of them.
BY_ARCH <- MAIN[, .(n_seeds = uniqueN(seed), rho_med = median(rho),
                    delta_r2_med = median(delta_r2)),
                by = .(arch, method, set, role)]
setorder(BY_ARCH, arch, method, -rho_med)
fwrite(BY_ARCH, file.path(OUTDIR, "by_arch.tsv"), sep = "\t")

# ---------------------------------------------------------------- 4. per-seed wide views
for (m in levels(MAIN$method)) {
    sub <- SC[method == m]
    if (!nrow(sub)) next
    fwrite(dcast(sub, seed + arch + control ~ set, value.var = "rho"),
           file.path(OUTDIR, paste0("rho_by_seed_", m, ".tsv")), sep = "\t")
}
if (nrow(CTRL)) fwrite(CTRL, file.path(OUTDIR, "controls.tsv"), sep = "\t")

# ---------------------------------------------------------------- 5. figures
theme_j08 <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
ROLE_COL <- c(candidate = "#0B775E", ceiling = "#35274A", floor = "#F2300F")

p1 <- ggplot(MAIN, aes(set, rho, fill = role)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.6) +
    facet_wrap(~ method) +
    scale_fill_manual(values = ROLE_COL) +
    labs(x = NULL, y = "Spearman rho (offset vs realized fitness decline)", fill = "role",
         title = "Genetic offset vs fitness decline",
         subtitle = paste0("discriminating replicates only (headroom > 0.05): n = ",
                           uniqueN(MAIN$seed))) +
    theme_j08
ggsave(file.path(OUTDIR, "rho_by_set.png"), p1, width = 11, height = 5, dpi = 150)
ggsave(file.path(OUTDIR, "rho_by_set.svg"), p1, width = 11, height = 5)

p2 <- ggplot(MAIN, aes(set, delta_r2, fill = role)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.6) +
    facet_wrap(~ method) +
    scale_fill_manual(values = ROLE_COL) +
    labs(x = NULL, y = expression(Delta*R^2~"over the climate-distance null"),
         fill = "role", title = "Increment over the Euclidean climate-distance null") +
    theme_j08
ggsave(file.path(OUTDIR, "delta_r2_by_set.png"), p2, width = 11, height = 5, dpi = 150)
ggsave(file.path(OUTDIR, "delta_r2_by_set.svg"), p2, width = 11, height = 5)

p3 <- ggplot(GAP[!is.na(rho_floor)], aes(rho_floor, rho, colour = method, shape = set)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(size = 2, alpha = 0.85) +
    labs(x = "rho of the size-matched random floor", y = "rho of the candidate set",
         title = "Every point above the diagonal is a candidate beating its own null") +
    theme_j08 + theme(axis.text.x = element_text(angle = 0))
ggsave(file.path(OUTDIR, "candidate_vs_floor.png"), p3, width = 8, height = 5.5, dpi = 150)
ggsave(file.path(OUTDIR, "candidate_vs_floor.svg"), p3, width = 8, height = 5.5)

message("=== by set x method (non-control seeds) ===")
print(BY_SET)
message("\n=== candidates vs their own floor and ceiling ===")
print(WIN)
message("\nWrote tables and figures to ", OUTDIR)
