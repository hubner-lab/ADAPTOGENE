#!/usr/bin/env Rscript
# =============================================================================
# eval_offset_lind.R -- score the garden sweep with the metric of Lind & Lotterhos 2025.
#
# NO METRIC OF OUR OWN. Both performance flavours are taken from their validation code,
# 01_src/MVP_03_validate_gradient_forests.py::calculate_performance:
#
#   garden_performance   "correlation of fitness and offset within gardens across transplants"
#                        offset.corrwith(fitness, axis=1, method='kendall')
#                        -> ONE Kendall's tau per garden, across the 100 source populations
#   source_performance   "correlation of fitness and offset for transplants across gardens"
#                        -> ONE tau per source population, across gardens
#   slopes               linregress(offset, fitness) at the same two levels
#
# Sign convention is theirs: fitness falls as offset rises, so a GOOD model gives a NEGATIVE
# tau. They invert the axis in figures; nothing is inverted in this table.
#
# Marker-set vocabulary is mapped to theirs so the comparison table can be read directly:
#   all -> all | best -> adaptive | rand_best1 -> neutral | truth -> (ours, not in their design)
#
# Inputs:
#   {outdir}/gardens_{seed}.tsv, {outdir}/garden_fitness_{seed}.tsv   (mvp_garden_fitness.R)
#   {outdir}/gardens/{seed}/{garden}/{method}__{set}.tsv              (mvp_garden_run.sh)
#
# Usage:  Rscript eval_offset_lind.R [--seeds=all|CSV] [--outdir=DIR]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args    <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt     <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
SEEDS_A <- opt("seeds", "all")
OUTDIR  <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/offset09"))

# Panel -> vocabulary. THEIR three terms are mapped onto the panels that match their
# definitions, which is not the naive name match:
#   their "adaptive" = SNPs at QTNs (Methods: "mean N = 188")      -> our `truth`
#   their "neutral"  = SNPs on linkage groups without any QTNs      -> our `neutral_all`
#   their "all"      = every marker                                 -> our `all`
# Our GEA-derived panels have NO counterpart in their design -- they never tested a
# GEA-identified candidate set -- so they carry their own labels and are reported separately.
MARKER_SET <- c(all         = "all",
                neutral_all = "neutral",
                truth       = "adaptive",
                best        = "gea_best",
                union       = "gea_union",
                intersect3  = "gea_strict",
                rand_best1  = "random_matched",
                # Single-method panels: same scan, same verbatim operating points as
                # best/union/intersect3, agreement step removed. They are what turns
                # "combining methods helps" into a measured contrast rather than an
                # assumption. A panel absent from this map is SILENTLY SKIPPED by the
                # loop below, so anything new must be registered here.
                solo_lfmm   = "gea_lfmm_only",
                solo_rda    = "gea_rda_only",
                solo_emmax  = "gea_emmax_only")
METHOD_LABEL <- c(gradient_forest = "GFoffset", geometric_offset = "LFMM2offset",
                  rda_offset = "RDA-uncorrected", rda_corrected = "RDA-corrected")

MAN   <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
               colClasses = c("seed" = "character"))
SEEDS <- if (SEEDS_A == "all") MAN$seed else strsplit(SEEDS_A, ",", fixed = TRUE)[[1]]

kend <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok], method = "kendall"))
}
slope <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    # A source population is scored across gardens, and a garden across sources; either series
    # can be constant (e.g. one garden harvested so far), which makes the slope undefined.
    if (sum(ok) < 3L || sd(x[ok]) == 0) return(NA_real_)
    unname(coef(lm(y[ok] ~ x[ok]))[2])
}

garden_rows <- list(); source_rows <- list(); miss <- list()

for (s in SEEDS) {
    gdir <- file.path(OUTDIR, "gardens", s)
    if (!dir.exists(gdir)) { message("MVP", s, ": no harvested gardens -- skipped"); next }

    GARD <- fread(file.path(OUTDIR, paste0("gardens_", s, ".tsv")))
    FIT  <- fread(file.path(OUTDIR, paste0("garden_fitness_", s, ".tsv")),
                  colClasses = c("seed" = "character"))
    setkey(FIT, garden_id, source_site)

    for (gd in list.dirs(gdir, full.names = FALSE, recursive = FALSE)) {
        files <- list.files(file.path(gdir, gd), pattern = "\\.tsv$", full.names = TRUE)
        for (f in files) {
            base <- sub("\\.tsv$", "", basename(f))
            parts <- strsplit(base, "__", fixed = TRUE)[[1]]
            if (length(parts) != 2L) next
            meth <- parts[1]; set_name <- parts[2]
            if (!set_name %in% names(MARKER_SET)) next

            GO <- fread(f, colClasses = c("site" = "character", "sample" = "character"))
            # one value per source population (all individuals of a deme share its value)
            OS <- GO[, .(offset = mean(genetic_offset, na.rm = TRUE)), by = .(source_site = site)]
            D  <- merge(OS, FIT[garden_id == gd, .(source_site, fitness)], by = "source_site")
            if (nrow(D) < 10 || sd(D$offset) == 0) {
                miss[[length(miss) + 1L]] <- data.table(seed = s, garden = gd, method = meth,
                                                        set = set_name, n = nrow(D),
                                                        reason = "too few sources or constant offset")
                next
            }
            garden_rows[[length(garden_rows) + 1L]] <- data.table(
                seed = s, garden_id = gd,
                garden_type = GARD[garden_id == gd, type][1],
                dev = GARD[garden_id == gd, dev][1],
                method = meth, method_label = METHOD_LABEL[[meth]],
                set = set_name, marker_set = MARKER_SET[[set_name]],
                n_sources = nrow(D),
                tau = kend(D$offset, D$fitness),
                slope = slope(D$offset, D$fitness),
                # Dispersion of the predicted offset ACROSS source populations within this
                # garden. A garden where every deme gets nearly the same offset cannot rank
                # them, so tau there measures noise rather than method quality -- this is the
                # quantity that decides whether a weak tau is a floor effect or a real failure.
                offset_mean = mean(D$offset), offset_sd = sd(D$offset),
                offset_cv = sd(D$offset) / abs(mean(D$offset)),
                offset_iqr = IQR(D$offset),
                fitness_cv = sd(D$fitness) / abs(mean(D$fitness)))
        }
    }

    # ---- source performance: per source population, across gardens -----------
    # Guard: a seed can have a gardens/ directory whose files are still being written by a
    # running lane, so there may be no scored rows for it yet.
    if (!length(garden_rows)) next
    GR <- rbindlist(garden_rows)[seed == s]
    if (!nrow(GR)) next
    for (meth in unique(GR$method)) for (st in unique(GR$set)) {
        long <- rbindlist(lapply(unique(GR[method == meth & set == st, garden_id]), function(gd) {
            f <- file.path(gdir, gd, paste0(meth, "__", st, ".tsv"))
            if (!file.exists(f)) return(NULL)
            GO <- fread(f, colClasses = c("site" = "character", "sample" = "character"))
            GO[, .(garden_id = gd, offset = mean(genetic_offset, na.rm = TRUE)),
               by = .(source_site = site)]
        }))
        if (!nrow(long)) next
        D <- merge(long, FIT[, .(garden_id, source_site, fitness)],
                   by = c("garden_id", "source_site"))
        src <- D[, .(n_gardens = .N, tau = kend(offset, fitness),
                     slope = slope(offset, fitness)), by = source_site]
        src[, `:=`(seed = s, method = meth, method_label = METHOD_LABEL[[meth]],
                   set = st, marker_set = MARKER_SET[[st]])]
        source_rows[[length(source_rows) + 1L]] <- src
    }
}

if (!length(garden_rows)) stop("No harvested offsets found under ", file.path(OUTDIR, "gardens"))

GP <- rbindlist(garden_rows)
GP <- merge(GP, MAN[, .(seed, arch = arch_level, demog = demog_level, final_LA,
                        n_causal_maf01, r2_pc1_temp)], by = "seed", all.x = TRUE)
GP[, control := demog == "Est-Clines"]
fwrite(GP, file.path(OUTDIR, "garden_performance.tsv"), sep = "\t")

SP <- rbindlist(source_rows)
if (nrow(SP)) {
    SP <- merge(SP, MAN[, .(seed, arch = arch_level, demog = demog_level, final_LA)],
                by = "seed", all.x = TRUE)
    SP[, control := demog == "Est-Clines"]
    fwrite(SP, file.path(OUTDIR, "source_performance.tsv"), sep = "\t")
}
if (length(miss)) fwrite(rbindlist(miss), file.path(OUTDIR, "scoring_skipped.tsv"), sep = "\t")

message("\nWrote garden_performance.tsv (", nrow(GP), " rows) and source_performance.tsv (",
        nrow(SP), " rows)")
message("\n=== median Kendall's tau per method x marker set (landscape gardens, non-control) ===")
print(dcast(GP[garden_type == "landscape" & control == FALSE,
                .(tau = median(tau, na.rm = TRUE)), by = .(method_label, marker_set)],
            marker_set ~ method_label, value.var = "tau"))
message("\n=== novelty gardens: median tau by deviation (non-control, adaptive set) ===")
print(dcast(GP[garden_type == "novelty" & control == FALSE & marker_set == "adaptive",
                .(tau = median(tau, na.rm = TRUE)), by = .(dev, method_label)],
            dev ~ method_label, value.var = "tau"))
