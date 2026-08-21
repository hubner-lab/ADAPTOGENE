#!/usr/bin/env Rscript
# =============================================================================
# mvp_write_offset_configs.R -- generate config_MVP{seed}_off.yaml for the offset arm.
#
# The generated file is the seed's BASE config plus exactly two additions:
#   1. Climate.custom.future_table  -> the scenario table written by mvp_make_future_env.R
#   2. a Maladaptation: block       -> the three offset methods, nospatial, snp_sets: all
#
# Everything else -- project_name, Input.*, Filter.*, LD.*, sNMF.*, Climate.*, PreGEA.*,
# GEA.*, GWAS.* -- is copied through byte-for-byte. That is the point of patching the base
# config's own text instead of re-templating: `project_name` is unchanged, so the run lands
# in the SAME {PROJECT}_results/ tree and reuses every upstream product (imputed matrices,
# PCA, kinship, staged present climate) rather than recomputing it. CLAUDE.md's
# "two configs, one project_name" hazard is being used deliberately here, and it is safe
# only because the two configs are identical upstream by construction -- which is exactly
# what generating rather than hand-authoring guarantees.
#
# Method settings and why:
#   spatial_correction: without  -- geometric_offset and rda_offset are nospatial-only by
#     wildcard constraint (workflow/rules/maladaptation.smk), and no MVP project has run
#     mode=climate, so climate/tables/{spatial,varpart}/ are empty and a spatial GF would
#     have no dbMEM inputs. Uniform nospatial also keeps the three methods comparable.
#   random_model: false -- GF's random model feeds only the cumulative-importance and
#     importance plots; it never reaches the offset. The random comparison in this benchmark
#     is a SNP SET (rand_*), run through all three methods, not this flag.
#   snp_sets: all -- resolve_active_snp_sets() globs _intermediate/snp_sets/, so the roster
#     adapts per seed automatically (a set that came out empty simply has no directory).
#
# Usage:
#   Rscript mvp_write_offset_configs.R [--seeds=all|CSV] [--ntree=500] [--suffix=_off]
# =============================================================================

suppressPackageStartupMessages(library(data.table))

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))   # parse_kv_args()

args    <- parse_kv_args(commandArgs(trailingOnly = TRUE))
opt     <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]
SEEDS_A <- opt("seeds", "all")
NTREE   <- opt("ntree", "500")
SUFFIX  <- opt("suffix", "_off")

MAN   <- fread(file.path(PIPELINE_ROOT, "benchmarks/mvp_seeds.tsv"),
               colClasses = c("seed" = "character"))
SEEDS <- if (SEEDS_A == "all") MAN$seed else strsplit(SEEDS_A, ",", fixed = TRUE)[[1]]

mala_block <- function(seed) c(
    "",
    "#-----------------------------------------------------------------------------",
    "# MALADAPTATION -- genetic-offset arm (journal 08). Added by",
    "# benchmarks/mvp_write_offset_configs.R; do not hand-edit, regenerate.",
    "#",
    "# All three methods run nospatial: geometric_offset and rda_offset are nospatial-only",
    "# by wildcard constraint, and mode=climate (dbMEM/varpart) has never been run for the",
    "# MVP projects, so a spatial GF would have no inputs.",
    "#",
    "# random_model: false -- GF's neutral model only feeds the importance plots and never",
    "# reaches the offset. The random control here is a SNP SET (rand_best*/rand_union1/",
    "# rand_solo1, built by benchmarks/mvp_build_snp_sets.R from background-neutral loci",
    "# only), scored through all three methods like any other set.",
    "#-----------------------------------------------------------------------------",
    "Maladaptation:",
    "  methods:",
    "    gradient_forest:",
    sprintf("      ntree: '%s'", NTREE),
    "      cor_threshold: '0.5'",
    "      spatial_correction: without",
    "      random_model: false",
    "    geometric_offset:",
    "      scale: true",
    "      k: ''                                  # empty -> uses sNMF.k_best for this seed",
    "    rda_offset:",
    "      axes: auto",
    "      axis_alpha: 0.05",
    "      permutations: 999",
    "      condition_pcs: 0                       # canonical (docs/rda_research.md B7)",
    "      seed: 42",
    "  snp_sets: all                              # glob _intermediate/snp_sets/"
)

written <- character(0)
for (s in SEEDS) {
    base_f <- file.path(PIPELINE_ROOT, sprintf("config_MVP%s.yaml", s))
    out_f  <- file.path(PIPELINE_ROOT, sprintf("config_MVP%s%s.yaml", s, SUFFIX))
    if (!file.exists(base_f)) stop("Missing base config: ", base_f)

    L <- readLines(base_f, warn = FALSE)

    if (any(grepl("^Maladaptation:", L)))
        stop(base_f, " already has a Maladaptation block -- patching would duplicate it")
    if (any(grepl("future_table", L) & !grepl("^\\s*#", L)))
        stop(base_f, " already sets future_table -- refusing to patch")

    # ---- 1. future_table, inserted inside Climate.custom, right after present_table
    i <- grep("^\\s*present_table:", L)
    if (length(i) != 1L) stop(base_f, ": expected exactly one present_table: line, got ", length(i))
    indent <- sub("^(\\s*).*$", "\\1", L[i])
    fut_line <- sprintf("%sfuture_table: \"MVP%s_env_future.tsv\"   # journal 08 scenario",
                        indent, s)
    L <- append(L, fut_line, after = i)

    # ---- 2. Maladaptation block at the end
    L <- c(L, mala_block(s))

    writeLines(L, out_f)
    written <- c(written, out_f)
    message("wrote ", basename(out_f))
}

message("\n", length(written), " config(s) written. Diff check against the base config:")
for (f in written) {
    s <- sub(".*config_MVP([0-9]+).*", "\\1", basename(f))
    base_f <- file.path(PIPELINE_ROOT, sprintf("config_MVP%s.yaml", s))
    d <- system2("diff", c(shQuote(base_f), shQuote(f)), stdout = TRUE, stderr = TRUE)
    added   <- sum(grepl("^> ", d))
    removed <- sum(grepl("^< ", d))
    if (removed != 0L) stop(basename(f), ": diff removes ", removed,
                            " line(s) from the base config -- upstream is NOT identical")
    message(sprintf("  %-32s +%d lines, -0", basename(f), added))
}
