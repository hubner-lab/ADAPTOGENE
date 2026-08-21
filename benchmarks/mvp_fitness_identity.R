#!/usr/bin/env Rscript
# =============================================================================
# mvp_fitness_identity.R -- recover the MVP fitness function per seed, as an identity.
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# The MVP deposit ships fitness "of the individual in the deme where it was sampled"
# (home-deme fitness only), so an offset benchmark needs TRANSPLANT fitness -- fitness of
# the same genotypes under a different environmental optimum. docs/gea_simulation_benchmarks.md
# and benchmarks/capblancq_eval/FINDINGS.md both recorded this as "re-derive Lotterhos 2023
# Eq. 3, then validate the re-derivation on a deposit that ships fitness".
#
# That re-derivation is unnecessary. Fitness in these sims is a deterministic function of
# the phenotype-optimum mismatch, and the deposit ships phenotypes (phen_temp, phen_sal),
# deme optima (bio_1, bio_2) and realized fitness for the same 1000 individuals. So the
# function can be READ OUT of the shipped data by regressing
#
#     log(fitness_i) ~ (phen_temp_i - opt1_j)^2 + (phen_sal_i - opt2_j)^2
#
# If the model is right, this is not a fit but an identity: R^2 = 1 to machine precision and
# the coefficients are the negative inverse-doubled selection variances. It is then applied
# to ANY optimum map to get transplant fitness, with no assumption imported from the paper.
#
# THE COEFFICIENTS ARE PER SEED -- this is the load-bearing result. Three regimes appear:
#   equal-S 2-trait    b1 = -1.000, b2 = -1.000
#   unequal-S 2-trait  b1 = -0.125, b2 = -1.000   (sigma differs between traits)
#   1-trait            b1 = -2.000, b2 =  0.000   (phen_sal is constant, the axis is inert)
# A single hardcoded fitness function would fabricate a selection penalty on the inert axis
# for the 1-trait seeds and use the wrong selection strength on the unequal-S seeds, in both
# cases silently.
#
# Usage:
#   Rscript mvp_fitness_identity.R [--seeds=all|CSV] [--outdir=DIR]
#
# Output:
#   {outdir}/fitness_identity.tsv  -- one row per seed: n, b0, b1, b2, r2, max_abs_resid,
#                                     sd_phen_sal, n_traits, arch, gate
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))   # parse_kv_args()

args    <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt     <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
SEEDS_A <- opt("seeds", "all")
OUTDIR  <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset08"))
DATADIR <- file.path(PIPELINE_ROOT, "data/mvp")

# R^2 gate. The claim is an identity, so anything below this is a modelling failure, not a
# noisy fit -- fail loudly rather than propagate a wrong fitness target downstream.
R2_GATE <- 1 - 1e-6

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MAN <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
             colClasses = c("seed" = "character"))
SEEDS <- if (SEEDS_A == "all") MAN$seed else strsplit(SEEDS_A, ",", fixed = TRUE)[[1]]

out <- vector("list", length(SEEDS))

for (i in seq_along(SEEDS)) {
    s    <- SEEDS[i]
    proj <- paste0("MVP", s)
    d    <- file.path(DATADIR, proj)

    md_f  <- file.path(d, paste0(proj, "_metadata.tsv"))
    env_f <- file.path(d, paste0(proj, "_env_present.tsv"))
    for (f in c(md_f, env_f)) if (!file.exists(f)) stop("Missing input: ", f)

    md  <- fread(md_f,  colClasses = c("site" = "character", "sample" = "character"))
    env <- fread(env_f, colClasses = c("site" = "character"))

    dt <- merge(md, env, by = "site", all.x = TRUE)
    if (anyNA(dt$bio_1) || anyNA(dt$bio_2))
        stop(proj, ": individuals whose deme has no optimum row -- metadata/env key mismatch")

    dt[, d1 := (phen_temp - bio_1)^2]
    dt[, d2 := (phen_sal  - bio_2)^2]
    dt[, logf := log(fitness)]

    fit   <- lm(logf ~ d1 + d2, data = dt)
    cf    <- coef(fit)
    resid <- residuals(fit)
    r2    <- summary(fit)$r.squared

    # 1-trait seeds: phen_sal is constant, so d2 has (near) zero variance and its
    # coefficient is not identified -- lm returns NA. That is the correct answer (the axis
    # carries no selection), recorded as 0 rather than propagated as NA.
    b2 <- unname(cf["d2"])
    if (is.na(b2)) b2 <- 0

    sd_sal   <- sd(dt$phen_sal)
    n_traits <- MAN[seed == s, n_traits]

    gate <- r2 >= R2_GATE
    # 1-trait seeds must come back with no penalty on the inert axis; a nonzero b2 there
    # means the target column is wrong, not that selection exists on that axis.
    # Tolerance is 1e-5, the scale of the OLS residuals themselves (max|resid| ~ 8e-6 on every
    # seed): on the 1-trait seeds phen_sal is constant but the OPTIMUM bio_2 still varies across
    # demes, so d2 is identified and returns b2 ~ 1e-8 -- numerically zero, not exactly zero.
    if (length(n_traits) == 1L && !is.na(n_traits) && n_traits == 1L && abs(b2) > 1e-5)
        gate <- FALSE

    out[[i]] <- data.table(
        seed          = s,
        project       = proj,
        n             = nrow(dt),
        b0            = unname(cf["(Intercept)"]),
        b1            = unname(cf["d1"]),
        b2            = b2,
        r2            = r2,
        max_abs_resid = max(abs(resid)),
        sd_phen_sal   = sd_sal,
        n_traits      = if (length(n_traits) == 1L) n_traits else NA_integer_,
        arch          = MAN[seed == s, arch_level][1],
        gate          = gate
    )

    message(sprintf("%s  n=%d  b1=%+.6f  b2=%+.6f  R2=%.12f  max|resid|=%.3e  %s",
                    proj, nrow(dt), unname(cf["d1"]), b2, r2, max(abs(resid)),
                    if (gate) "OK" else "GATE FAIL"))
}

OUT <- rbindlist(out)
fwrite(OUT, file.path(OUTDIR, "fitness_identity.tsv"), sep = "\t")
message("\nWrote ", file.path(OUTDIR, "fitness_identity.tsv"))

if (any(!OUT$gate)) {
    print(OUT[gate == FALSE])
    stop("Fitness identity gate failed on ", sum(!OUT$gate), " seed(s) -- ",
         "the transplant-fitness target cannot be trusted on those. Stop here.")
}
message("All ", nrow(OUT), " seeds pass the identity gate (R2 >= ", format(R2_GATE), ").")
