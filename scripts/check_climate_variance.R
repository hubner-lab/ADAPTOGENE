#!/usr/bin/env Rscript
# Detect invariant (zero-variance or all-NA) factors at sample sites.
# Writes a TSV with columns: predictor, reason
# Header-only file when all factors are valid.
#
# Serves two modules (Phase 3):
#   COLUMNS=bio     climate predictors, "^bio_" prefixed (mode=climate)
#   COLUMNS=traits  phenotypic traits — every numeric column that is not
#                   sample/site/coordinate metadata (mode=traits), the same
#                   exclusion rule as plot_density.R:27-31

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
SITE_FILE  <- args[1]
OUT_FILE   <- args[2]
COLUMNS    <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "bio"

site <- fread(SITE_FILE, sep = '\t', header = TRUE)

meta_cols <- c("sample", "site", "latitude", "longitude", "lat", "lon", "ID")

if (identical(COLUMNS, "bio")) {
    target_cols <- setdiff(grep("^bio_", colnames(site), value = TRUE), meta_cols)
    label <- "climate predictors"
} else if (identical(COLUMNS, "traits")) {
    cand <- colnames(site)[!tolower(colnames(site)) %in% tolower(meta_cols)]
    target_cols <- cand[vapply(cand, function(x) is.numeric(site[[x]]), logical(1))]
    label <- "traits"
} else {
    stop("check_climate_variance.R: unknown column selection '", COLUMNS,
         "' (expected 'bio' or 'traits')")
}

invariant <- data.table(predictor = character(), reason = character())

for (col in target_cols) {
    vals <- site[[col]]
    # Coerce to numeric — empty strings become NA
    vals <- suppressWarnings(as.numeric(as.character(vals)))

    if (all(is.na(vals))) {
        invariant <- rbind(invariant, data.table(predictor = col, reason = "all_na"))
        message("WARNING: ", col, " — all values are NA/empty")
    } else if (sd(vals, na.rm = TRUE) == 0) {
        invariant <- rbind(invariant, data.table(predictor = col, reason = "zero_variance"))
        message("WARNING: ", col, " — zero variance across all samples")
    }
}

fwrite(invariant, OUT_FILE, sep = '\t', quote = FALSE)

if (nrow(invariant) == 0) {
    message("INFO: All ", length(target_cols), " ", label, " have valid variance")
} else {
    message("INFO: Found ", nrow(invariant), " invariant ", label, ": ",
            paste(invariant$predictor, collapse = ", "))
}
