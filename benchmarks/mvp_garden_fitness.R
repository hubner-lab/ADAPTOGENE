#!/usr/bin/env Rscript
# =============================================================================
# mvp_garden_fitness.R -- build the in-silico common gardens of Lind & Lotterhos 2025
# and the population-mean-fitness matrix they validate genomic offset against.
#
# NOTHING HERE IS OUR OWN DESIGN. Every choice below was read out of the authors' own
# code (ModelValidationProgram/MVP-offsets) or their paper, and the source is named inline.
#
# GARDENS
#   within-landscape : all 100 demes, each deme's own (bio_1, bio_2) is a garden.
#       "full factorial in silico reciprocal transplant design", 100 gardens per simulation.
#   climate-novelty  : 12 gardens at deviations from the climate centre, taken verbatim from
#       02_analysis/04_outlier_climate_runs/06_multivariate_calculate_climate_outlier_fitness.ipynb:
#           deviations = [0, 1.72, 2.35, 2.74, 3.13, 3.53, 3.92, 4.31, 4.70, 5.09, 5.48, 5.88]
#           outlier_vals[dev][env] = envdata[env].mean() + (dev * np.std(envdata[env]))
#       -- i.e. BOTH environmental axes are moved to mean + dev*SD simultaneously (the loop runs
#       `for env in envs`), and dev = 0 (the climate centre itself) is included as a garden.
#       np.std is the population SD (ddof = 0), so sd() here uses ddof = 0 to match.
#
# FITNESS
#   MVP_climate_outlier_fitness_calculator.R computes
#       fits = dmvnorm(phenos, opts, sigma = fitness_varcov) / fitness_norm
#       fitness[1, transplant_ID] = mean(fits)
#   i.e. a normalised multivariate-normal density at each individual's phenotype relative to the
#   garden optimum, then the ARITHMETIC MEAN over the individuals of a source population. With a
#   diagonal covariance that density ratio is exactly exp(b1*d1^2 + b2*d2^2) with b_i = -1/(2*sigma_i^2),
#   which is the identity benchmarks/mvp_fitness_identity.R recovered from the shipped data at
#   R^2 = 1.000000. So this script needs no distributional assumption of its own -- it applies the
#   seed's own (b1, b2) to a new optimum.
#
#   NOTE the mean is over FITNESS, not over log-fitness (journal 08 averaged log-fitness). Kept as
#   they do it, because their Kendall's tau is computed on that quantity.
#
# Usage:
#   Rscript mvp_garden_fitness.R [--seeds=all|CSV] [--outdir=DIR]
#
# Outputs, per seed:
#   {outdir}/gardens_{seed}.tsv        -- garden_id, type, dev, opt1, opt2
#   {outdir}/garden_fitness_{seed}.tsv -- garden_id x source_site population mean fitness (long)
#   {outdir}/garden_qc.tsv             -- local-adaptation checks, all seeds
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))   # parse_kv_args()

args    <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt     <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
SEEDS_A <- opt("seeds", "all")
OUTDIR  <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset09"))
IDENT_F <- opt("identity", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset08/fitness_identity.tsv"))

# Verbatim from the authors' notebook.
DEVIATIONS <- c(0, 1.72, 2.35, 2.74, 3.13, 3.53, 3.92, 4.31, 4.70, 5.09, 5.48, 5.88)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(IDENT_F)) stop("Missing fitness identity table: ", IDENT_F)
IDENT <- fread(IDENT_F, colClasses = c("seed" = "character"))
if (any(!IDENT$gate)) stop("fitness_identity.tsv has failed gates -- refusing to build gardens")

MAN   <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
               colClasses = c("seed" = "character"))
SEEDS <- if (SEEDS_A == "all") MAN$seed else strsplit(SEEDS_A, ",", fixed = TRUE)[[1]]

# population SD, ddof = 0, to match numpy's np.std default
sd_pop <- function(x) sqrt(mean((x - mean(x))^2))

qc <- list()

for (s in SEEDS) {
    proj <- paste0("MVP", s)
    env_f  <- file.path(PIPELINE_ROOT, "data/mvp", proj, paste0(proj, "_env_present.tsv"))
    pipe_f <- file.path(PIPELINE_ROOT, paste0(proj, "_results"),
                        "_intermediate/samples/metadata_climate.tsv")
    for (f in c(env_f, pipe_f)) if (!file.exists(f)) stop(proj, ": missing ", f)

    env  <- fread(env_f,  colClasses = c("site" = "character"))
    pipe <- fread(pipe_f, colClasses = c("site" = "character", "sample" = "character"))
    setorder(env, site)

    b1 <- IDENT[seed == s, b1]; b2 <- IDENT[seed == s, b2]
    if (length(b1) != 1L) stop(proj, ": no fitness-identity row")

    # ---- gardens ------------------------------------------------------------
    G_land <- env[, .(garden_id = site, type = "landscape", dev = NA_real_,
                      opt1 = bio_1, opt2 = bio_2)]

    c1 <- mean(env$bio_1); s1 <- sd_pop(env$bio_1)
    c2 <- mean(env$bio_2); s2 <- sd_pop(env$bio_2)
    G_nov <- data.table(garden_id = sprintf("novel_%s", format(DEVIATIONS, trim = TRUE)),
                        type = "novelty", dev = DEVIATIONS,
                        opt1 = c1 + DEVIATIONS * s1,
                        opt2 = c2 + DEVIATIONS * s2)

    G <- rbind(G_land, G_nov)
    fwrite(G, file.path(OUTDIR, paste0("gardens_", s, ".tsv")), sep = "\t")

    # ---- population mean fitness in every garden ----------------------------
    # Vectorised over gardens: for each garden, mean over individuals of the source deme.
    P <- pipe[, .(site, sample, phen_temp, phen_sal)]
    out <- vector("list", nrow(G))
    for (g in seq_len(nrow(G))) {
        o1 <- G$opt1[g]; o2 <- G$opt2[g]
        # The per-individual fitness must be a COLUMN, not an external vector: inside
        # DT[, .(...), by=] an outside vector is not subset by group, so `mean(f)` would
        # return the same grand mean for every population (silently, and it did).
        P[, fit := exp(b1 * (phen_temp - o1)^2 + b2 * (phen_sal - o2)^2)]
        out[[g]] <- P[, .(garden_id = G$garden_id[g], type = G$type[g],
                          n_ind = .N, fitness = mean(fit)), by = .(source_site = site)]
    }
    FM <- rbindlist(out)
    FM[, seed := s]
    setcolorder(FM, c("seed", "garden_id", "type", "source_site", "n_ind", "fitness"))
    fwrite(FM, file.path(OUTDIR, paste0("garden_fitness_", s, ".tsv")), sep = "\t")

    # ---- QC: local adaptation ----------------------------------------------
    # Sympatric = source population in its OWN deme's garden. Under local adaptation it should be
    # at or near the top of that population's row, and mean(sympatric) - mean(allopatric) should
    # reproduce the deposit's own final_LA for this seed.
    L <- FM[type == "landscape"]
    symp  <- L[garden_id == source_site, fitness]
    allo  <- L[garden_id != source_site, fitness]
    home_rank <- L[, {
        ord <- order(-fitness)
        .(home_rank = which(garden_id[ord] == source_site[1]))
    }, by = source_site]

    la_ours <- mean(symp) - mean(allo)
    la_dep  <- MAN[seed == s, final_LA]

    qc[[length(qc) + 1L]] <- data.table(
        seed = s, n_gardens = nrow(G), n_landscape = nrow(G_land), n_novelty = nrow(G_nov),
        n_sources = uniqueN(FM$source_site),
        mean_sympatric = mean(symp), mean_allopatric = mean(allo),
        LA_ours = la_ours, LA_deposit = la_dep, LA_diff = la_ours - la_dep,
        median_home_rank = median(home_rank$home_rank),
        pct_home_in_top10 = 100 * mean(home_rank$home_rank <= 10),
        env_centre_bio1 = c1, env_sd_bio1 = s1, env_centre_bio2 = c2, env_sd_bio2 = s2,
        novelty_opt1_max = max(G_nov$opt1), novelty_opt2_max = max(G_nov$opt2)
    )

    message(sprintf(paste0("%s  gardens=%d (%d landscape + %d novelty)  sources=%d  ",
                           "LA ours=%.4f deposit=%.4f diff=%+.4f  median home rank=%.0f"),
                    proj, nrow(G), nrow(G_land), nrow(G_nov), uniqueN(FM$source_site),
                    la_ours, la_dep, la_ours - la_dep, median(home_rank$home_rank)))
}

QC <- rbindlist(qc)
fwrite(QC, file.path(OUTDIR, "garden_qc.tsv"), sep = "\t")
message("\nWrote ", file.path(OUTDIR, "garden_qc.tsv"))
print(QC[, .(seed, n_gardens, LA_ours = round(LA_ours, 4), LA_deposit = round(LA_deposit, 4),
             LA_diff = round(LA_diff, 4), median_home_rank, pct_home_in_top10)])

# Hard gates.
#  1. The home deme must be a near-best garden for its own population. Not rank 1 exactly: demes sit
#     on a 10x10 grid at 0.2 deg spacing, so a population's immediate neighbours have almost the same
#     optimum and routinely tie with home. Observed median home rank is 2-7 out of 100 gardens.
#  2. Our sympatric-minus-allopatric difference must reproduce the deposit's own final_LA. This is the
#     real check that the garden fitness matrix is the authors' quantity and not something of ours.
bad <- QC[median_home_rank > 15]
if (nrow(bad)) { print(bad); stop("Sympatric fitness is not near-top-ranked on some seeds") }
bad2 <- QC[abs(LA_diff) > 0.05]
if (nrow(bad2)) { print(bad2); stop("Local adaptation does not reproduce the deposit's final_LA") }
message("Local-adaptation QC passed on all ", nrow(QC), " seeds ",
        "(max |LA - final_LA| = ", round(max(abs(QC$LA_diff)), 4),
        ", max median home rank = ", max(QC$median_home_rank), ").")
