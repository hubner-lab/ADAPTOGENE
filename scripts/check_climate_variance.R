#!/usr/bin/env Rscript
# Detect invariant (zero-variance or all-NA) bioclimatic predictors at sample sites.
# Writes a TSV with columns: predictor, reason
# Header-only file when all predictors are valid.

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
SITE_FILE  <- args[1]
OUT_FILE   <- args[2]

site <- fread(SITE_FILE, sep = '\t', header = TRUE)

# Bio columns only (exclude sample/site/lat/lon metadata)
meta_cols <- c("sample", "site", "latitude", "longitude", "ID")
bio_cols  <- setdiff(grep("^bio_", colnames(site), value = TRUE), meta_cols)

invariant <- data.table(predictor = character(), reason = character())

for (col in bio_cols) {
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
    message("INFO: All ", length(bio_cols), " climate predictors have valid variance")
} else {
    message("INFO: Found ", nrow(invariant), " invariant predictor(s): ",
            paste(invariant$predictor, collapse = ", "))
}
