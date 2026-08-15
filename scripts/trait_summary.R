#!/usr/bin/env Rscript
# Per-trait summary statistics for the Phenotypic Factors module (mode=traits).
#
# Reads the project metadata (site, sample, latitude, longitude, then one
# column per trait) and writes one row per trait:
#   trait, n, n_missing, pct_missing, mean, sd, min, median, max
#
# Trait columns are every numeric column that is not sample/site/coordinate
# metadata — the same exclusion rule as plot_density.R:27-31 and
# check_climate_variance.R's `traits` mode.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
META_FILE <- args[1]
OUT_FILE  <- args[2]

meta <- fread(META_FILE, sep = "\t", header = TRUE,
              colClasses = c("site" = "character", "sample" = "character"))

meta_cols <- c("sample", "site", "latitude", "longitude", "lat", "lon", "ID")
cand <- colnames(meta)[!tolower(colnames(meta)) %in% tolower(meta_cols)]
trait_cols <- cand[vapply(cand, function(x) is.numeric(meta[[x]]), logical(1))]

if (length(trait_cols) == 0) {
    message("WARNING: no numeric trait columns found in ", META_FILE,
            " — writing header-only summary")
    fwrite(data.table(trait = character(), n = integer(), n_missing = integer(),
                      pct_missing = numeric(), mean = numeric(), sd = numeric(),
                      min = numeric(), median = numeric(), max = numeric()),
           OUT_FILE, sep = "\t", quote = FALSE)
    quit(status = 0)
}

summary_dt <- rbindlist(lapply(trait_cols, function(tr) {
    vals <- suppressWarnings(as.numeric(meta[[tr]]))
    ok   <- vals[!is.na(vals)]
    data.table(
        trait       = tr,
        n           = length(ok),
        n_missing   = sum(is.na(vals)),
        pct_missing = round(100 * sum(is.na(vals)) / length(vals), 2),
        mean        = if (length(ok) > 0) round(mean(ok), 4) else NA_real_,
        sd          = if (length(ok) > 1) round(stats::sd(ok), 4) else NA_real_,
        min         = if (length(ok) > 0) round(min(ok), 4) else NA_real_,
        median      = if (length(ok) > 0) round(stats::median(ok), 4) else NA_real_,
        max         = if (length(ok) > 0) round(max(ok), 4) else NA_real_
    )
}))

fwrite(summary_dt, OUT_FILE, sep = "\t", quote = FALSE, na = "NA")
message("INFO: Wrote summary for ", nrow(summary_dt), " trait(s): ",
        paste(summary_dt$trait, collapse = ", "))
