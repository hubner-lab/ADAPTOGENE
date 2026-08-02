#!/usr/bin/env Rscript
# mvp_select_seeds.R -- build the MVP benchmark seed manifest.
#
# Inputs (both persisted under data/mvp/selection/):
#   summary_20220428_20220726.csv  -- BCO-DMO doi:10.26008/1912/bco-dmo.889769.1, one row per seed
#   0b-final_params-20220428.txt   -- ModelValidationProgram/MVP-NonClinalAF src/, seed -> design cell
#
# Output:
#   benchmarks/mvp_seeds.tsv  -- 14 rows: 12 primary + 2 degenerate controls
#
# Selection rule (docs/gea_simulation_benchmarks.md Sec 8.2):
#   primary  : meanFst >= 0.05 AND cor_PC1_temp^2 in [0.30, 0.60] AND cor_PC1_sal^2 <= 0.20
#              -> 55 seeds, all SS-Mtn. The 12 used here are the stratified subset in Sec 8.3.
#   controls : the INVERSE on the alignment axis (cor_PC1_temp^2 > 0.9, i.e. degenerate = the
#              Laruson failure mode) while holding demography and architecture constant against
#              two of the primaries, so LANDSCAPE is the only thing that differs.
#
# The 12 primary seeds are transcribed from the doc, then re-verified here against the raw
# summary table. A mismatch means the doc and the deposit disagree -- fail loudly, do not proceed.

suppressPackageStartupMessages({
  library(data.table)
})

ROOT      <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
SEL_DIR   <- file.path(ROOT, "data/mvp/selection")
OUT_TSV   <- file.path(ROOT, "benchmarks/mvp_seeds.tsv")

# Seeds proposed in docs/gea_simulation_benchmarks.md Sec 8.3 (stratified over architecture).
PRIMARY_SEEDS <- c(1231288, 1232218, 1231993, 1232398,   # oligogenic
                   1232578, 1232548, 1231858, 1233208,   # moderately polygenic
                   1232923, 1231798, 1232908, 1233133)   # highly polygenic

message("== reading selection inputs ==")
summ <- fread(file.path(SEL_DIR, "summary_20220428_20220726.csv"))
# The params file was written by R with row.names -- header carries one fewer field than the
# data rows, so read.table's row-name inference is what we want here (fread would misalign).
par  <- as.data.table(read.table(file.path(SEL_DIR, "0b-final_params-20220428.txt"),
                                 header = TRUE, stringsAsFactors = FALSE))

message(sprintf("  summary: %d seeds x %d cols", nrow(summ), ncol(summ)))
message(sprintf("  params : %d seeds x %d cols", nrow(par),  ncol(par)))

keep_par <- c("seed", "arch", "arch_level", "arch_level_sub", "demog_level", "demog_level_sub",
              "N_traits", "ispleiotropy", "MIG_x", "MIG_y")
stopifnot("K" %in% names(summ))
d <- merge(summ, par[, ..keep_par], by = "seed")
stopifnot(nrow(d) == nrow(summ))
message(sprintf("  joined : %d / %d seeds", nrow(d), nrow(summ)))

# The deposit ships Kendall correlations between structure and environment; the screening
# criteria are stated on R^2, so square them here rather than anywhere downstream.
d[, r2_pc1_temp := cor_PC1_temp^2]
d[, r2_pc1_sal  := cor_PC1_sal^2]

# ---------------------------------------------------------------- primary set ----
message("== primary set ==")
missing_seeds <- setdiff(PRIMARY_SEEDS, d$seed)
if (length(missing_seeds)) {
  stop("Seeds from docs Sec 8.3 not present in the deposit summary: ",
       paste(missing_seeds, collapse = ", "))
}
primary <- d[seed %in% PRIMARY_SEEDS]

# Re-verify the doc's own selection rule holds for every transcribed seed. This is the check
# that the doc was not mis-transcribed, and it costs nothing.
bad <- primary[!(meanFst >= 0.05 & r2_pc1_temp >= 0.30 & r2_pc1_temp <= 0.60 & r2_pc1_sal <= 0.20)]
if (nrow(bad)) {
  print(bad[, .(seed, meanFst, r2_pc1_temp, r2_pc1_sal)])
  stop("Primary seeds violate the Sec 8.2 selection rule -- doc and deposit disagree.")
}
if (length(unique(primary$demog_level)) != 1L || primary$demog_level[1] != "SS-Mtn") {
  stop("Expected all primary seeds in the SS-Mtn landscape; got: ",
       paste(unique(primary$demog_level), collapse = ", "))
}
message(sprintf("  %d seeds verified against the Sec 8.2 rule, all %s / %s",
                nrow(primary), primary$demog_level[1],
                paste(unique(primary$demog_level_sub), collapse = ",")))

# ---------------------------------------------------------------- controls ------
# Degenerate arm: same demography sub-level and architecture as a primary seed, but the
# Est-Clines landscape, where R^2(PC1~temp) ~ 0.96. Expectation when these are run: structure
# correction destroys the signal, reproducing the Laruson collapse on purpose.
message("== degenerate controls ==")
fst_lo <- min(primary$meanFst); fst_hi <- max(primary$meanFst)
cand <- d[demog_level == "Est-Clines" &
          demog_level_sub == unique(primary$demog_level_sub)[1] &
          r2_pc1_temp > 0.90 &
          meanFst >= fst_lo & meanFst <= fst_hi]
message(sprintf("  %d Est-Clines candidates in the primary set's Fst range [%.3f, %.3f]",
                nrow(cand), fst_lo, fst_hi))

# One control per architecture level that carries the most primaries, highest Fst within each.
target_arch <- c("moderately-polygenic", "highly-polygenic")
target_arch <- intersect(target_arch, unique(cand$arch_level))
if (length(target_arch) < 2L) {
  # Fall back to whatever architecture levels are available, still 2 controls.
  target_arch <- head(unique(cand$arch_level), 2L)
}
controls <- rbindlist(lapply(target_arch, function(a) {
  cand[arch_level == a][order(-meanFst)][1L]
}))
if (nrow(controls) != 2L || anyNA(controls$seed)) {
  print(cand[, .(seed, arch_level, meanFst, r2_pc1_temp, r2_pc1_sal)])
  stop("Could not pick 2 degenerate controls -- inspect the candidate table above.")
}
message(sprintf("  picked %s", paste(controls$seed, collapse = ", ")))

# ---------------------------------------------------------------- manifest ------
mk <- function(x, arm) {
  data.table(
    seed            = as.integer(x$seed),
    arm             = arm,
    project         = paste0("MVP", x$seed),
    architecture    = x$arch,
    arch_level      = x$arch_level,
    demog_level     = x$demog_level,
    demog_level_sub = x$demog_level_sub,
    n_traits        = x$N_traits,
    ispleiotropy    = x$ispleiotropy,
    # "number of populations used in analyses" -- the authors' own per-seed K, shipped as a
    # column. It is the defensible source for sNMF.k_best here because sNMF cross-entropy has
    # no interior minimum on this landscape (100 demes on a continuous stepping-stone grid;
    # verified on seed 1231288, where it decreases monotonically from K=2 to K=10). Floored at
    # 2 for the pipeline, which needs >= 2 clusters for structure plots and pop_diff_test.
    K_authors       = as.integer(x$K),
    k_best          = pmax(2L, as.integer(x$K)),
    meanFst         = round(x$meanFst, 4),
    r2_pc1_temp     = round(x$r2_pc1_temp, 4),
    r2_pc1_sal      = round(x$r2_pc1_sal, 4),
    n_snps          = x$nSNPs,
    n_causal_pre    = x$num_causal_prefilter,
    n_causal_maf01  = x$num_causal_postfilter,
    n_causal_temp   = x$num_causal_temp,
    n_causal_sal    = x$num_causal_sal,
    final_LA        = round(x$final_LA, 4)
  )
}
manifest <- rbind(mk(primary, "primary"), mk(controls, "control_degenerate"))
setorder(manifest, arm, -meanFst)

dir.create(dirname(OUT_TSV), showWarnings = FALSE, recursive = TRUE)
fwrite(manifest, OUT_TSV, sep = "\t")
message(sprintf("== wrote %s (%d seeds) ==", OUT_TSV, nrow(manifest)))
print(manifest[, .(seed, arm, arch_level, meanFst, r2_pc1_temp, r2_pc1_sal,
                   n_snps, n_causal_maf01, K_authors, k_best)])
