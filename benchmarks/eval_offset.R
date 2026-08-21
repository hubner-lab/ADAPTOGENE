#!/usr/bin/env Rscript
# =============================================================================
# eval_offset.R -- score predicted genetic offset against realized fitness decline.
#
# Companion to benchmarks/eval_detection.R, which scores p-values against causal-locus
# truth. This one scores the OFFSET arm: for each (seed, method, SNP set), how well does
# the predicted per-deme offset track the fitness the demes actually lose when the climate
# moves, and does it add anything over simply measuring how far the climate moved.
#
# INPUTS
#   {PROJECT}_results/Maladaptation/tables/{method}/{set}_{spatial}/genetic_offset_site.tsv
#       columns: site, sample, genetic_offset  (one row per SAMPLE; all samples in a deme
#       share the deme's value, so it is averaged by site here)
#   {outdir}/target_{seed}.tsv   from benchmarks/mvp_make_future_env.R
#       columns: site, delta_logf, clim_dist, ...
#
# METRICS, per (seed, method, set), over the demes
#   decline    = -delta_logf   (positive = the deme loses fitness; sign flipped so that a
#                               better offset predictor gives a POSITIVE correlation)
#   rho        Spearman(offset, decline) -- rank agreement, the primary metric. Rank-based
#              because the three methods' offsets are on incomparable scales (GF: turnover
#              units; geometric: LEA gap units; RDA: weighted axis distance) and nothing in
#              the question depends on the scale.
#   r2         Pearson R^2 of decline ~ offset
#   null_r2    R^2 of decline ~ clim_dist -- the Euclidean climate-transfer distance null.
#              This is the bar: a genomic offset that cannot beat "how far did the climate
#              move" adds nothing. Reported for the Capblancq/Gain comparison (0.33-0.77 and
#              ~0.45 respectively).
#   full_r2    R^2 of decline ~ clim_dist + offset
#   delta_r2   full_r2 - null_r2 -- the increment the genomics buys over the null
#   partial_rho  Spearman of the two residual series after regressing both decline and
#              offset on clim_dist -- the rank-based version of the same increment
#
# Usage:
#   Rscript eval_offset.R [--seeds=all|CSV] [--outdir=DIR] [--spatial=nospatial]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))   # parse_kv_args()

args    <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt     <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
SEEDS_A <- opt("seeds", "all")
OUTDIR  <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset08"))
SPATIAL <- opt("spatial", "nospatial")
METHODS <- strsplit(opt("methods", "gradient_forest,geometric_offset,rda_offset"), ",")[[1]]

MAN   <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
               colClasses = c("seed" = "character"))
SEEDS <- if (SEEDS_A == "all") MAN$seed else strsplit(SEEDS_A, ",", fixed = TRUE)[[1]]

SETSUM <- fread(file.path(OUTDIR, "snp_sets_summary.tsv"), colClasses = c("seed" = "character"))

# Which sets are candidates, which are the ceiling, which are the floor. Written into every
# row so no downstream reading can mistake `truth` for a recommendation.
set_role <- function(x) fifelse(x %in% c("best", "union", "solo"), "candidate",
                        fifelse(x == "truth", "ceiling",
                        fifelse(grepl("^rand", x), "floor", "other")))

rows <- list()
inv_rows <- list()   # per (seed, method): does the SNP set change the deme RANKING at all?

for (s in SEEDS) {
    proj <- paste0("MVP", s)
    tgt_f <- file.path(OUTDIR, paste0("target_", s, ".tsv"))
    if (!file.exists(tgt_f)) { message(proj, ": no target table -- skipped"); next }
    TG <- fread(tgt_f, colClasses = c("site" = "character"))
    TG[, decline := -delta_logf]

    res_dir <- file.path(PIPELINE_ROOT, paste0(proj, "_results"), "Maladaptation", "tables")
    for (m in METHODS) {
        mdir <- file.path(res_dir, m)
        if (!dir.exists(mdir)) next
        offs_by_set <- list()   # site-level offset vectors, for the invariance check below
        for (dn in list.dirs(mdir, full.names = FALSE, recursive = FALSE)) {
            if (!grepl(paste0("_", SPATIAL, "$"), dn)) next
            set_name <- sub(paste0("_", SPATIAL, "$"), "", dn)
            f <- file.path(mdir, dn, "genetic_offset_site.tsv")
            if (!file.exists(f)) next

            GO <- fread(f, colClasses = c("site" = "character", "sample" = "character"))
            # One value per deme. The three scripts all write one row per sample carrying
            # its deme's value; mean() is exact when they agree and would expose it loudly
            # in n_distinct if they ever did not.
            GS <- GO[, .(offset = mean(genetic_offset, na.rm = TRUE),
                         offset_sd_within = sd(genetic_offset, na.rm = TRUE)), by = site]

            D <- merge(TG, GS, by = "site")
            D <- D[is.finite(offset) & is.finite(decline)]
            if (nrow(D) < 10) {
                message(sprintf("%s %s %s: only %d demes -- skipped", proj, m, set_name, nrow(D)))
                next
            }

            ct   <- suppressWarnings(cor.test(D$offset, D$decline, method = "spearman"))
            r2   <- suppressWarnings(cor(D$offset, D$decline)^2)
            nullm <- lm(decline ~ clim_dist, data = D)
            fullm <- lm(decline ~ clim_dist + offset, data = D)
            null_r2 <- summary(nullm)$r.squared
            full_r2 <- summary(fullm)$r.squared
            # rank-based increment: both series stripped of the climate-distance signal first
            r_dec <- residuals(lm(decline ~ clim_dist, data = D))
            r_off <- residuals(lm(offset  ~ clim_dist, data = D))
            prho  <- suppressWarnings(cor(r_dec, r_off, method = "spearman"))

            ns <- SETSUM[seed == s & set == set_name, n_snps]
            nc <- SETSUM[seed == s & set == set_name, n_causal]

            rows[[length(rows) + 1L]] <- data.table(
                seed = s, project = proj, method = m, set = set_name,
                role = set_role(set_name),
                n_snps  = if (length(ns)) ns[1] else NA_integer_,
                n_causal = if (length(nc)) nc[1] else NA_integer_,
                n_demes = nrow(D),
                rho = unname(ct$estimate), rho_p = ct$p.value,
                r2 = r2, null_r2 = null_r2, full_r2 = full_r2,
                delta_r2 = full_r2 - null_r2, partial_rho = prho,
                offset_min = min(D$offset), offset_max = max(D$offset),
                offset_within_site_sd_max = max(D$offset_sd_within, na.rm = TRUE)
            )
            offs_by_set[[set_name]] <- setNames(D$offset, D$site)
        }

        # ---- does the SNP set change anything? -------------------------------------
        # If two sets' offset vectors are exact scalar multiples of each other, the SNP set
        # cannot change the deme ranking and the whole question is moot FOR THAT METHOD --
        # measured, not assumed. min_pearson near 1 with a ratio spread at machine precision
        # is the signature. (Observed on rda_offset: ratio sd ~ 1e-15 across all nine sets.)
        if (length(offs_by_set) >= 2L) {
            sites_common <- Reduce(intersect, lapply(offs_by_set, names))
            M <- sapply(offs_by_set, function(v) v[sites_common])
            pear <- cor(M);  spear <- cor(M, method = "spearman")
            ratio_spread <- max(vapply(seq_len(ncol(M))[-1], function(j) {
                r <- M[, j] / M[, 1]; diff(range(r)) / mean(r)
            }, numeric(1)))
            inv_rows[[length(inv_rows) + 1L]] <- data.table(
                seed = s, method = m, n_sets = ncol(M),
                min_pearson = min(pear[lower.tri(pear)]),
                min_spearman = min(spear[lower.tri(spear)]),
                max_rel_ratio_spread = ratio_spread)
        }
    }
}

if (!length(rows)) stop("No offset outputs found -- has mode=maladaptation run?")

SC <- rbindlist(rows)
SC <- merge(SC, MAN[, .(seed, arch = arch_level, demog = demog_level, n_traits,
                        n_causal_maf01, r2_pc1_temp, meanFst)],
            by = "seed", all.x = TRUE)
SC[, control := demog == "Est-Clines"]
setorder(SC, seed, method, -rho)

fwrite(SC, file.path(OUTDIR, "scores.tsv"), sep = "\t")
message("\nWrote ", file.path(OUTDIR, "scores.tsv"), "  (", nrow(SC), " rows)")

if (length(inv_rows)) {
    INV <- rbindlist(inv_rows)
    fwrite(INV, file.path(OUTDIR, "set_invariance.tsv"), sep = "\t")
    message("\n=== does the SNP set change the deme ranking at all? (per method, over seeds) ===")
    print(INV[, .(n_seeds = .N,
                  min_pearson_med  = median(min_pearson),
                  min_spearman_med = median(min_spearman),
                  ratio_spread_med = median(max_rel_ratio_spread)), by = method])
}

message("\n=== Spearman rho (offset vs fitness decline), median over non-control seeds ===")
print(dcast(SC[control == FALSE, .(rho = median(rho)), by = .(method, set)],
            set ~ method, value.var = "rho"))
message("\n=== delta_R2 over the climate-distance null, median over non-control seeds ===")
print(dcast(SC[control == FALSE, .(d = median(delta_r2)), by = .(method, set)],
            set ~ method, value.var = "d"))
