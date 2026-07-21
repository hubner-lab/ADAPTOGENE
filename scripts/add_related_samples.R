#!/usr/bin/env Rscript
# add_related_samples.R — one-off dev utility, NOT a Snakemake rule.
#
# Injects near-duplicate ("clone") samples into data/SIMDATA.vcf and
# data/SIMDATA_metadata.tsv so the Processing relatedness/IBS filter has
# something real to catch during the alpha test. For each donor sample below,
# appends an exact genotype copy under a new "<donor>_DUP" sample name:
#   - Same GT at every variant -> IBS ~= 1.0 with its donor (survives LD
#     pruning regardless of pruning parameters; a near-duplicate with allele
#     flips risks dropping below the 0.99 threshold on this small dataset).
#   - A small extra fraction of genotypes are set to missing (on top of
#     whatever was already missing) so the duplicate has a distinctly higher
#     F_MISS than its donor -- this exercises the relatedness_filter tie-break
#     rule ("drop the higher-missingness member of each over-threshold pair").
#     The extra missingness rate is chosen well below Filter.sample_miss
#     (0.5) so the duplicate survives sample-missingness filtering and
#     actually reaches the relatedness step.
#
# One donor per site (Negev/TelAviv/Galilee) so all three populations get a
# duplicate cluster in the relatedness MDS. Idempotent: aborts if any
# "*_DUP" sample already exists in either file, rather than double-injecting.
#
# Usage: Rscript /pipeline/scripts/add_related_samples.R
# (run from the repo root / mounted /pipeline, so the relative data/ paths resolve)

suppressPackageStartupMessages(library(data.table))

set.seed(42)

VCF_PATH      <- "data/SIMDATA.vcf"
METADATA_PATH <- "data/SIMDATA_metadata.tsv"

# One donor per site/population.
DONORS <- c("NEG01", "TAV05", "GAL03")
EXTRA_MISSING_RATE <- 0.10   # additional fraction of non-missing calls flipped to './.'

#=============================================================================
# 1. VCF — read, locate donor columns, append duplicate columns
#=============================================================================
lines <- readLines(VCF_PATH)
header_idx <- which(startsWith(lines, "#CHROM"))
stopifnot(length(header_idx) == 1)

chrom_cols <- strsplit(lines[header_idx], "\t")[[1]]
sample_cols <- chrom_cols[-(1:9)]  # columns 1-9 are CHROM..FORMAT

dup_names <- paste0(DONORS, "_DUP")
if (any(dup_names %in% sample_cols)) {
    stop("Duplicate samples already present (", paste(intersect(dup_names, sample_cols), collapse = ", "),
         ") -- add_related_samples.R already ran on this VCF. Aborting to avoid double injection.")
}

donor_idx <- match(DONORS, sample_cols)
if (any(is.na(donor_idx))) {
    stop("Donor sample(s) not found in VCF: ", paste(DONORS[is.na(donor_idx)], collapse = ", "))
}
# Column index in the full tab-split line (offset by the 9 fixed VCF columns)
donor_line_idx <- donor_idx + 9

# New header: original + duplicate sample names
new_header <- paste(c(chrom_cols, dup_names), collapse = "\t")

data_idx <- (header_idx + 1):length(lines)
n_variants <- length(data_idx)

new_data_lines <- character(n_variants)
for (i in seq_along(data_idx)) {
    fields <- strsplit(lines[data_idx[i]], "\t")[[1]]
    donor_gt <- fields[donor_line_idx]

    dup_gt <- donor_gt
    # Flip a subset of currently-non-missing calls to missing (per donor, independently)
    non_missing <- which(dup_gt != "./.")
    flip_mask <- non_missing[runif(length(non_missing)) < EXTRA_MISSING_RATE]
    dup_gt[flip_mask] <- "./."

    new_data_lines[i] <- paste(c(fields, dup_gt), collapse = "\t")
}

writeLines(c(lines[1:(header_idx - 1)], new_header, new_data_lines), VCF_PATH)

message("INFO: Added ", length(dup_names), " duplicate sample(s) to ", VCF_PATH,
        " (", n_variants, " variants updated): ", paste(dup_names, collapse = ", "))

#=============================================================================
# 2. Metadata — append one row per duplicate, copying the donor's site/coords/traits
#=============================================================================
meta <- fread(METADATA_PATH, colClasses = c(site = "character", sample = "character"))

if (any(dup_names %in% meta$sample)) {
    stop("Duplicate sample rows already present in metadata -- aborting to avoid double injection.")
}

donor_rows <- meta[match(DONORS, sample)]
if (any(is.na(donor_rows$sample))) {
    stop("Donor sample(s) not found in metadata: ", paste(DONORS[is.na(donor_rows$sample)], collapse = ", "))
}

dup_rows <- copy(donor_rows)
dup_rows[, sample := dup_names]

fwrite(rbind(meta, dup_rows), METADATA_PATH, sep = "\t", quote = FALSE, na = "NA")

message("INFO: Added ", nrow(dup_rows), " duplicate metadata row(s) to ", METADATA_PATH, ": ",
        paste(dup_names, collapse = ", "))
message("INFO: SIMDATA now has ", nrow(meta) + nrow(dup_rows), " samples (was ", nrow(meta), ").")
