#!/usr/bin/env Rscript
# probe_rda_cost.R -- where does the RDA wall-clock actually go?
#
#   Rscript benchmarks/probe_rda_cost.R [--snps=4000] [--perms=99]
#
# Diagnostic only. Writes nothing into any {PROJECT}_results/ tree and touches no pipeline
# state. It reconstructs the same call shape scripts/rda.R uses -- rda(Y ~ preds +
# Condition(PCs), scale = TRUE) on the imputed genotype matrix -- and times the ordination
# separately from each of the three anova.cca calls, because the pipeline emits no log line
# between them and a stalled run is therefore indistinguishable between "still in rda()" and
# "still in the permutation tests".
#
# Scaling matters as much as the split: cost is reported per-SNP so the full-set (10-13k SNP)
# figure can be extrapolated from a cheap pruned-set (~4k SNP) run.

suppressPackageStartupMessages({ library(data.table); library(vegan) })

A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
    a <- sub("^--", "", a)
    A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}
NSNP  <- as.integer(if (is.null(A$snps))  4000L else A$snps)
PERMS <- as.integer(if (is.null(A$perms)) 99L   else A$perms)
ROOT  <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
RES   <- file.path(ROOT, "MVP1232548_results")

cat("OPENBLAS_NUM_THREADS=", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n", sep = "")

G <- as.matrix(fread(file.path(RES, "_intermediate/climate_subset/lfmm_imp_climate.lfmm"),
                     header = FALSE))
if (ncol(G) > NSNP) G <- G[, seq_len(NSNP)]
env <- fread(file.path(RES, "climate/tables/present/climate_present_site_scaled.tsv"))
pcs <- as.matrix(fread(file.path(RES,
    "_work/maf0.01_miss0.1_smiss0.5/ld0.2_win10_step2/MVP1232548.pca/MVP1232548.projections"),
    header = FALSE))[, 1:5]
colnames(pcs) <- paste0("PC", 1:5)

# Env is per site; the genotype matrix is per individual. Any row-aligned surrogate is fine
# for a COST probe -- the timing depends on matrix dimensions and rank, not on which numbers
# sit in the predictor columns. Nothing here is scored.
set.seed(42)
n <- nrow(G)
V <- data.frame(bio_1 = rnorm(n), bio_2 = rnorm(n), pcs[seq_len(n), , drop = FALSE])
cat(sprintf("matrix: %d individuals x %d SNPs\n", nrow(G), ncol(G)))

t_rda <- system.time(
    m <- rda(G ~ bio_1 + bio_2 + Condition(PC1 + PC2 + PC3 + PC4 + PC5), data = V, scale = TRUE)
)[["elapsed"]]
cat(sprintf("rda()                       %8.1f s\n", t_rda))

t_full <- system.time(anova.cca(m, permutations = PERMS, parallel = 8))[["elapsed"]]
cat(sprintf("anova.cca full (parallel=8) %8.1f s   (perms=%d)\n", t_full, PERMS))

t_axis <- system.time(anova.cca(m, by = "axis", permutations = PERMS))[["elapsed"]]
cat(sprintf("anova.cca by=axis (serial)  %8.1f s   (perms=%d)\n", t_axis, PERMS))

t_marg <- system.time(anova.cca(m, by = "margin", permutations = PERMS))[["elapsed"]]
cat(sprintf("anova.cca by=margin (serial)%8.1f s   (perms=%d)\n", t_marg, PERMS))

anova_total <- t_full + t_axis + t_marg
cat(sprintf("\nanova total                 %8.1f s  (%.0f%% of one fit's cost)\n",
            anova_total, 100 * anova_total / (t_rda + anova_total)))
cat(sprintf("at perms=999 the anova block would be ~%.1f s\n", anova_total * 999 / PERMS))
cat(sprintf("per-SNP: rda %.4f s, anova %.4f s\n", t_rda / ncol(G), anova_total / ncol(G)))
