#!/usr/bin/env Rscript
# select_haplotype_regions.R
# Select top regions for haplotype analysis from association/phenotype/overlap results
#
# Args:
#   1  regions_file     - Combined regions TSV (from any association mode)
#   2  output_file      - Output TSV with selected regions + suggested MGmin
#   3  top_n            - Number of top regions to select
#   4  source           - Source label (association, association_phenotypes, overlapping)

args <- commandArgs(trailingOnly = TRUE)
regions_file   <- args[1]
output_file    <- args[2]
top_n          <- as.integer(args[3])
source_label   <- args[4]

library(data.table)

message("=== Select Haplotype Regions ===")
message("Regions file: ", regions_file)
message("Top N: ", top_n)
message("Source: ", source_label)

regions <- fread(regions_file)
message("Total regions: ", nrow(regions))

# Standardize column names (different sources may use different names)
# Expected columns: region_id, chr, start, end, snp_count (at minimum)
if (!"region_id" %in% names(regions) && "combined_region_id" %in% names(regions)) {
  setnames(regions, "combined_region_id", "region_id")
}

# Sort by snp_count descending to pick the most informative regions
if ("snp_count" %in% names(regions)) {
  setorder(regions, -snp_count)
} else if ("n_snps" %in% names(regions)) {
  setnames(regions, "n_snps", "snp_count")
  setorder(regions, -snp_count)
}

# Select top N regions
selected <- head(regions, top_n)

# Add suggested MGmin: max(10, floor(snp_count / 50))
if ("snp_count" %in% names(selected)) {
  selected[, suggested_mgmin := pmax(10L, as.integer(floor(snp_count / 50)))]
} else {
  selected[, suggested_mgmin := 10L]
}

# Add source label
selected[, source := source_label]

# Ensure required columns exist
required <- c("region_id", "chr", "start", "end")
for (col in required) {
  if (!col %in% names(selected)) {
    message("WARNING: Column '", col, "' not found in regions file")
  }
}

message("Selected ", nrow(selected), " regions")
fwrite(selected, output_file, sep = "\t", quote = FALSE)
message("Done.")
