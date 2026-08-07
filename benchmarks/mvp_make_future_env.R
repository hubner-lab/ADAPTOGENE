#!/usr/bin/env Rscript
# =============================================================================
# mvp_make_future_env.R -- build the future-climate table and the offset scoring target.
#
# WHAT IT WRITES, PER SEED
#   1. data/mvp/MVP{seed}/MVP{seed}_env_future.tsv   -- future deme optima, in the same
#      schema as the shipped *_env_present.tsv (site, bio_1, bio_2). This is what
#      Climate.custom.future_table points at, so it MUST live under Input.dir and carry the
#      same key/columns as the present table (stage_custom_climate.R reads both the same way).
#   2. {outdir}/target_{seed}.tsv  -- the scoring target, one row per deme.
#
# THE SCENARIO: BOTH axes move, with spatially independent amplitudes
# --------------------------------------------------------------------
#   delta1_j = D1MIN + (D1MAX - D1MIN) * (lon_j - lon_min)/(lon_max - lon_min)   [bio_1]
#   delta2_j = a fixed permutation of seq(-D2AMP, +D2AMP) over demes, seed D2SEED [bio_2]
#   bio_1_future = bio_1 + delta1 ;  bio_2_future = bio_2 + delta2
#
# WHY BOTH AXES -- this is the correction of the first version of this script, and it is
# the strongest methodological finding of journal 08. The first version shifted bio_1 only.
# Measured consequence on the anchor replicate: EVERY quadratic-form offset became a scalar
# multiple of the transfer distance, because with movement confined to one axis the
# quadratic form collapses to sqrt(w11)*|delta1| no matter what the weights are:
#
#     rda_offset        corr(offset, clim_dist) = 1.000000   -- identical for all 9 SNP sets
#     geometric_offset  corr(offset, clim_dist) = 0.990267   -- identical for all 9 SNP sets
#
# So the SNP set could not matter BY CONSTRUCTION for two of the three methods, and the
# benchmark could not answer its own question. A single-axis climate shift is a trap for
# any genomic-offset benchmark, not a property of these particular methods.
#
# WHY delta1 IS ALL-POSITIVE AND MONOTONE IN LONGITUDE. In these sims bio_1 is a perfect
# linear function of latitude (corr = 1.000) and bio_2 is a perfect V of longitude
# (corr(bio_2, -|lon - mid|) = -1.000, the SS-Mtn "mountain"). A latitude-graded delta1
# would therefore make bio_1_future an affine function of bio_1 -- transfer distance would
# be a relabelled present temperature. A longitude-graded delta1 that is CENTRED on zero is
# no better: the distance uses |delta1|, which is then a V of longitude and so collinear
# with bio_2 (measured R^2 = 0.478). An all-positive monotone delta1 keeps |delta1| monotone
# and lands orthogonal to both axes.
#
# WHY delta2 IS A DRAW, NOT A GRADIENT. After delta1 takes the one longitude-monotone
# direction, every remaining smooth spatial pattern on this 10x10 grid is collinear with
# bio_1 (latitude) or bio_2 (folded longitude). A fixed pseudo-random assignment over demes
# is independent of both by construction -- the same device as the "novel climate" gardens
# of the MVP-offsets pipeline. The permutation is keyed on the SORTED DEME NAMES with a
# fixed seed, so all 14 replicates get the SAME scenario (they share one 10x10 deme grid)
# and cross-replicate contrasts are not contaminated by scenario differences.
#
# Measured on the built scenario (anchor replicate, and identical by construction on the
# rest): corr(delta1, delta2) = -0.015, R^2(clim_dist ~ bio_1) = 0.002,
# R^2(clim_dist ~ bio_2) = 0.002, and the share of movement on axis 1 varies 0.13-0.98
# across demes -- which is what gives an axis-weighting model something to get right.
#
#   A spatially UNIFORM delta was rejected: it makes transfer distance constant, so the
#   Euclidean-climate-distance null (the bar every offset method must clear -- 0.33-0.77 on
#   the Capblancq deposit, ~0.45 in Gain et al. 2023) has zero variance and cannot be
#   computed at all.
#
# NOTE ON THE NULL. clim_dist = sqrt(delta1^2 + delta2^2) is the UNWEIGHTED Euclidean
# distance -- the standard null, and deliberately the wrong metric for the 1-trait seeds
# (b2 = 0, movement on bio_2 costs no fitness) and the unequal-S seeds (b1 = -0.125,
# b2 = -1.000, movement on bio_2 costs 8x more). That mismatch is the headroom: an offset
# model that learns the right axis weighting from the genotypes can beat it, and one that
# just measures how far the climate moved cannot.
#
# THE TARGET: realized per-deme log-fitness decline
# -------------------------------------------------
# Transplant fitness is exact, not modelled -- see benchmarks/mvp_fitness_identity.R:
#     log f_i(opt) = b1_s * (phen_temp_i - opt1)^2 + b2_s * (phen_sal_i - opt2)^2
# with the SEED'S OWN (b1, b2) read out of the shipped data at R^2 = 1. Applying it to the
# future optimum map gives instantaneous post-shift fitness (no re-adaptation, no generations
# elapsed -- the same construction as the fitness_pred column the Capblancq deposit ships).
#
#     delta_logf_j = mean_i in deme j [ log f_i(future) - log f_i(present) ]
#
# Individuals are taken from the PIPELINE's own climate-valid sample list
# ({PROJECT}_results/_intermediate/samples/metadata_climate.tsv), not from the raw deposit
# metadata: filtering can drop samples, and scoring a deme mean over individuals the model
# never saw would silently mismatch.
#
# Usage:
#   Rscript mvp_make_future_env.R [--seeds=all|CSV] [--outdir=DIR]
#                                 [--d1min=0.1] [--d1max=0.6] [--d2amp=0.5] [--d2seed=7]
#                                 [--force=true]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))   # parse_kv_args()

args    <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt     <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
SEEDS_A <- opt("seeds", "all")
OUTDIR  <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset08"))
D1MIN   <- as.numeric(opt("d1min", "0.1"))    # bio_1 shift at the western edge
D1MAX   <- as.numeric(opt("d1max", "0.6"))    # bio_1 shift at the eastern edge
D2AMP   <- as.numeric(opt("d2amp", "0.5"))    # bio_2 draw spans [-D2AMP, +D2AMP]
D2SEED  <- as.integer(opt("d2seed", "7"))     # fixes ONE scenario shared by all replicates
FORCE   <- tolower(opt("force", "false")) %in% c("true", "1", "yes")
DATADIR <- file.path(PIPELINE_ROOT, "data/mvp")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

IDENT_F <- file.path(OUTDIR, "fitness_identity.tsv")
if (!file.exists(IDENT_F))
    stop("Run benchmarks/mvp_fitness_identity.R first -- missing ", IDENT_F)
IDENT <- fread(IDENT_F, colClasses = c("seed" = "character"))
if (any(!IDENT$gate)) stop("fitness_identity.tsv contains failed gates -- refusing to build a target")

MAN   <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
               colClasses = c("seed" = "character"))
SEEDS <- if (SEEDS_A == "all") MAN$seed else strsplit(SEEDS_A, ",", fixed = TRUE)[[1]]

summary_rows <- vector("list", length(SEEDS))

for (i in seq_along(SEEDS)) {
    s    <- SEEDS[i]
    proj <- paste0("MVP", s)
    d    <- file.path(DATADIR, proj)

    env_f  <- file.path(d, paste0(proj, "_env_present.tsv"))
    fut_f  <- file.path(d, paste0(proj, "_env_future.tsv"))
    meta_f <- file.path(d, paste0(proj, "_metadata.tsv"))
    pipe_f <- file.path(PIPELINE_ROOT, paste0(proj, "_results"),
                        "_intermediate/samples/metadata_climate.tsv")
    for (f in c(env_f, meta_f, pipe_f)) if (!file.exists(f)) stop(proj, ": missing ", f)

    env  <- fread(env_f,  colClasses = c("site" = "character"))
    meta <- fread(meta_f, colClasses = c("site" = "character", "sample" = "character"))
    pipe <- fread(pipe_f, colClasses = c("site" = "character", "sample" = "character"))

    b1 <- IDENT[seed == s, b1]; b2 <- IDENT[seed == s, b2]
    if (length(b1) != 1L) stop(proj, ": no fitness-identity row")

    # ---- deme coordinates (one row per deme, taken from the raw metadata: coordinates are a
    # property of the deme, not of which individuals survived filtering)
    coords <- unique(meta[, .(site, latitude, longitude)])
    if (nrow(coords) != nrow(env))
        stop(proj, ": ", nrow(coords), " demes with coordinates vs ", nrow(env), " in env table")

    E <- merge(env, coords, by = "site")
    lon_rng <- range(E$longitude)
    if (diff(lon_rng) <= 0) stop(proj, ": zero longitude span -- amplitude gradient undefined")

    setorder(E, site)   # the delta2 permutation is keyed on sorted deme names
    E[, delta1 := D1MIN + (D1MAX - D1MIN) * (longitude - lon_rng[1]) / diff(lon_rng)]
    set.seed(D2SEED)
    E[, delta2 := sample(seq(-D2AMP, D2AMP, length.out = .N))]
    E[, bio_1_f := bio_1 + delta1]
    E[, bio_2_f := bio_2 + delta2]
    E[, clim_dist := sqrt(delta1^2 + delta2^2)]

    # Scenario contract. These are the properties the benchmark rests on, so they are
    # checked on every build rather than asserted in a comment.
    cor_dd  <- cor(E$delta1, E$delta2)
    r2_db1  <- cor(E$clim_dist, E$bio_1)^2
    r2_db2  <- cor(E$clim_dist, E$bio_2)^2
    share1  <- abs(E$delta1) / E$clim_dist
    # Hard gate on the two properties that are geometry-independent: the two shift axes must
    # be mutually independent, and the direction of movement must vary enough across demes
    # for an axis weighting to be discriminable at all.
    if (abs(cor_dd) > 0.1)
        stop(proj, ": the two shift axes are not independent -- corr(d1,d2)=", round(cor_dd, 3))
    # Correlation of the transfer distance with the PRESENT environment is a property of the
    # replicate's own landscape, not of the scenario: the SS-Mtn seeds put Env2 in a folded
    # longitude, the Est-Clines controls do not, so one shared scenario cannot be orthogonal
    # to both. Reported per replicate (r2_clim_dist_vs_opt1/2 below) and warned about rather
    # than enforced -- the alternative, re-optimising the scenario per replicate, would make
    # every cross-replicate contrast a comparison of scenarios as well.
    if (max(r2_db1, r2_db2) > 0.1)
        message(sprintf("  WARNING %s: transfer distance is correlated with the present env ",
                        proj), sprintf("(R2 vs bio_1 %.3f, vs bio_2 %.3f)", r2_db1, r2_db2))
    if (sd(share1) < 0.05)
        stop(proj, ": movement direction barely varies across demes (sd of axis-1 share = ",
             round(sd(share1), 3), ") -- axis weighting cannot be discriminated")

    # ---- future climate table, same schema as the present one
    FUT <- E[, .(site, bio_1 = bio_1_f, bio_2 = bio_2_f)]
    setorder(FUT, site)
    if (file.exists(fut_f) && !FORCE) {
        old <- fread(fut_f, colClasses = c("site" = "character"))
        setorder(old, site)
        if (!isTRUE(all.equal(as.data.frame(old), as.data.frame(FUT), tolerance = 1e-12)))
            stop(proj, ": ", fut_f, " already exists with DIFFERENT content. Refusing to ",
                 "overwrite an input file -- pass --force if the scenario really changed.")
        message(proj, ": future env table already present and identical -- kept")
    } else {
        fwrite(FUT, fut_f, sep = "\t")
    }

    # ---- target: per-deme realized log-fitness decline over the pipeline's own samples
    P <- merge(pipe[, .(site, sample, phen_temp, phen_sal)],
               E[, .(site, opt1 = bio_1, opt2 = bio_2, opt1_f = bio_1_f, opt2_f = bio_2_f,
                     delta1, delta2, clim_dist, latitude, longitude)],
               by = "site", all.x = TRUE)
    if (anyNA(P$opt1)) stop(proj, ": pipeline samples in demes absent from the env table")

    P[, logf_pres := b1 * (phen_temp - opt1)^2   + b2 * (phen_sal - opt2)^2]
    P[, logf_fut  := b1 * (phen_temp - opt1_f)^2 + b2 * (phen_sal - opt2_f)^2]

    TG <- P[, .(n_ind      = .N,
                latitude   = latitude[1],
                longitude  = longitude[1],
                opt1       = opt1[1],
                opt2       = opt2[1],
                opt1_f     = opt1_f[1],
                opt2_f     = opt2_f[1],
                delta1     = delta1[1],
                delta2     = delta2[1],
                clim_dist  = clim_dist[1],
                # adaptive lag: how far the deme's mean phenotype already sits from its own
                # optimum. It is the ONLY genomic quantity that enters delta_logf besides the
                # deterministic delta^2 term, so it is carried through for diagnostics.
                lag_temp   = mean(phen_temp) - opt1[1],
                logf_pres  = mean(logf_pres),
                logf_fut   = mean(logf_fut),
                delta_logf = mean(logf_fut - logf_pres),
                fit_pres   = mean(exp(logf_pres)),
                fit_fut    = mean(exp(logf_fut))),
            by = site]
    setorder(TG, site)
    TG[, seed := s]

    fwrite(TG, file.path(OUTDIR, paste0("target_", s, ".tsv")), sep = "\t")

    n_gain <- TG[delta_logf > 0, .N]
    null_r2 <- summary(lm(delta_logf ~ clim_dist, data = TG))$r.squared
    # How much better an ORACLE axis weighting would do than the unweighted null. This is the
    # headroom a genomically-informed offset is competing for; if it is ~0 on a replicate,
    # no SNP set can distinguish itself there and that must be visible, not inferred.
    orc_r2 <- summary(lm(delta_logf ~ I(abs(delta1)) + I(abs(delta2)), data = TG))$r.squared
    message(sprintf(paste0("%s  demes=%d  ind=%d  d1 %.2f-%.2f d2 %.2f-%.2f  ",
                           "dlogf med=%.3f [%.3f,%.3f]  gaining=%d  null R2=%.3f  oracle R2=%.3f"),
                    proj, nrow(TG), sum(TG$n_ind), min(TG$delta1), max(TG$delta1),
                    min(TG$delta2), max(TG$delta2),
                    median(TG$delta_logf), min(TG$delta_logf), max(TG$delta_logf),
                    n_gain, null_r2, orc_r2))

    summary_rows[[i]] <- data.table(
        seed = s, project = proj, n_demes = nrow(TG), n_ind = sum(TG$n_ind),
        b1 = b1, b2 = b2, d1min = D1MIN, d1max = D1MAX, d2amp = D2AMP, d2seed = D2SEED,
        delta_logf_med = median(TG$delta_logf),
        delta_logf_min = min(TG$delta_logf), delta_logf_max = max(TG$delta_logf),
        n_demes_gaining = n_gain,
        lag_sd = sd(TG$lag_temp),
        null_r2_clim_dist = null_r2,
        oracle_r2_axis_weighted = orc_r2,
        headroom = orc_r2 - null_r2,
        # Scenario contract, recomputed per replicate so the journal can show it was checked
        # everywhere rather than on one replicate.
        cor_d1_d2 = cor(TG$delta1, TG$delta2),
        r2_clim_dist_vs_opt1 = summary(lm(clim_dist ~ opt1, data = TG))$r.squared,
        r2_clim_dist_vs_opt2 = summary(lm(clim_dist ~ opt2, data = TG))$r.squared,
        share_axis1_sd = sd(abs(TG$delta1) / TG$clim_dist)
    )
}

SUM <- rbindlist(summary_rows)
fwrite(SUM, file.path(OUTDIR, "target_summary.tsv"), sep = "\t")
message("\nWrote ", file.path(OUTDIR, "target_summary.tsv"))
print(SUM[, .(seed, n_demes, delta_logf_med, n_demes_gaining,
              null_r2 = round(null_r2_clim_dist, 3),
              oracle_r2 = round(oracle_r2_axis_weighted, 3),
              headroom = round(headroom, 3),
              cor_d1d2 = round(cor_d1_d2, 3),
              r2_dist_b1 = round(r2_clim_dist_vs_opt1, 4),
              r2_dist_b2 = round(r2_clim_dist_vs_opt2, 4))])
