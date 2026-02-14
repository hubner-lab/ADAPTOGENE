#!/usr/bin/env Rscript
# run_haplotype_scan.R
# Run crosshap haplotyping across multiple epsilon values for selected regions
# Produces HapObject .qs files and clustree visualization plots
#
# Args:
#   1  selected_regions  - TSV with region_id, chr, start, end, snp_count
#   2  filtered_vcf      - Full filtered VCF (real genotypes)
#   3  ld_vcf            - Imputed VCF for LD computation (or same as filtered)
#   4  metadata_file     - Aligned metadata TSV (site, sample, lat, lon, traits...)
#   5  metadata_type     - 'site' or 'cluster_K{N}'
#   6  clusters_file     - Cluster assignments TSV (or 'NULL')
#   7  mgmin             - DBSCAN MinPts
#   8  minhap            - Minimum haplotype group size
#   9  epsilon_range     - Comma-separated epsilon values
#  10  intermediate_dir  - Output dir for .qs HapObject files
#  11  plots_dir         - Output dir for clustree plots

args <- commandArgs(trailingOnly = TRUE)
selected_regions <- args[1]
filtered_vcf     <- args[2]
ld_vcf           <- args[3]
metadata_file    <- args[4]
metadata_type    <- args[5]
clusters_file    <- args[6]
mgmin            <- as.integer(args[7])
minhap           <- as.integer(args[8])
epsilon_range    <- as.numeric(strsplit(args[9], ",")[[1]])
intermediate_dir <- args[10]
plots_dir        <- args[11]

library(data.table)
library(crosshap)
library(qs)
library(ggplot2)

message("=== Haplotype Scan ===")
message("Regions: ", selected_regions)
message("Filtered VCF: ", filtered_vcf)
message("LD VCF: ", ld_vcf)
message("Metadata type: ", metadata_type)
message("MGmin: ", mgmin, "  minHap: ", minhap)
message("Epsilon range: ", paste(epsilon_range, collapse = ", "))

# Read regions
regions <- fread(selected_regions)
message("Processing ", nrow(regions), " regions")

# Read metadata for phenotype and grouping
meta <- fread(metadata_file)
sample_col <- names(meta)[2]  # sample column
site_col   <- names(meta)[1]  # site column

# Phenotype: use first numeric trait (column 5+)
trait_cols <- names(meta)[5:ncol(meta)]
first_trait <- trait_cols[1]
message("Using trait '", first_trait, "' for crosshap phenotype input")

# Prepare metadata grouping
if (metadata_type == "site") {
  meta_groups <- meta[, .(ind = get(sample_col), group = get(site_col))]
} else if (grepl("^cluster", metadata_type) && clusters_file != "NULL") {
  clusters <- fread(clusters_file)
  # clusters table has sample and assigned_cluster columns
  clust_col <- grep("assigned|cluster", names(clusters), value = TRUE, ignore.case = TRUE)
  if (length(clust_col) == 0) clust_col <- names(clusters)[ncol(clusters)]
  meta_groups <- data.table(
    ind = clusters[[names(clusters)[grep("sample", names(clusters), ignore.case = TRUE)[1]]]],
    group = as.character(clusters[[clust_col[1]]])
  )
} else {
  meta_groups <- meta[, .(ind = get(sample_col), group = get(site_col))]
}

# Process each region
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  rid <- reg$region_id
  chr <- reg$chr
  start_pos <- reg$start
  end_pos <- reg$end

  message("\n--- Region ", rid, ": ", chr, ":", start_pos, "-", end_pos, " ---")

  # Create temp directory for this region
  tmp_dir <- file.path(intermediate_dir, paste0("tmp_region_", rid))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  region_vcf <- file.path(tmp_dir, "region.vcf")
  region_ld_vcf <- file.path(tmp_dir, "region_ld.vcf")

  # Extract region from VCF using awk (works with uncompressed VCFs without indexing)
  # Keeps header lines + data lines matching chr and within position range
  extract_region_vcf <- function(input_vcf, output_vcf, chr, start, end) {
    cmd <- paste0(
      "awk -v chr='", chr, "' -v start=", start, " -v end=", end,
      " 'BEGIN{OFS=\"\\t\"} /^#/{print; next} $1==chr && $2>=start && $2<=end{print}' ",
      input_vcf, " > ", output_vcf
    )
    system(cmd)
  }

  message("Extracting region VCF: ", chr, ":", start_pos, "-", end_pos)
  extract_region_vcf(filtered_vcf, region_vcf, chr, start_pos, end_pos)

  # Count data lines in extracted VCF
  vcf_lines <- readLines(region_vcf)
  data_lines <- vcf_lines[!grepl("^#", vcf_lines)]
  n_snps <- length(data_lines)
  if (n_snps < 3) {
    message("WARNING: Only ", n_snps, " SNPs in region ", rid, ", skipping")
    unlink(tmp_dir, recursive = TRUE)
    next
  }
  message("Extracted ", n_snps, " SNPs")

  # Extract region from LD VCF (imputed, for LD computation)
  extract_region_vcf(ld_vcf, region_ld_vcf, chr, start_pos, end_pos)

  # Compute pairwise LD with plink
  ld_prefix <- file.path(tmp_dir, "region_ld")
  cmd_plink <- paste0(
    "plink --vcf ", region_ld_vcf,
    " --const-fid 0",
    " --r2 square --out ", ld_prefix,
    " --allow-extra-chr"
  )
  message("Computing LD matrix: ", cmd_plink)
  system(cmd_plink)

  ld_file <- paste0(ld_prefix, ".ld")
  if (!file.exists(ld_file)) {
    message("WARNING: LD computation failed for region ", rid, ", skipping")
    unlink(tmp_dir, recursive = TRUE)
    next
  }

  # Prepare phenotype file (two columns, no header: ind | value)
  pheno_file <- file.path(tmp_dir, "pheno.txt")
  pheno_data <- meta[, .(ind = get(sample_col), value = get(first_trait))]
  pheno_data <- pheno_data[!is.na(value)]
  fwrite(pheno_data, pheno_file, sep = "\t", col.names = FALSE, quote = FALSE)

  # Prepare metadata file (two columns, no header: ind | group)
  meta_grp_file <- file.path(tmp_dir, "metadata.txt")
  fwrite(meta_groups, meta_grp_file, sep = "\t", col.names = FALSE, quote = FALSE)

  # Run crosshap
  tryCatch({
    message("Reading VCF for crosshap...")
    vcf_obj <- read_vcf(region_vcf)

    # Ensure unique SNP IDs (crosshap requires this)
    if (any(duplicated(vcf_obj$ID)) || all(vcf_obj$ID == ".")) {
      vcf_obj <- dplyr::mutate(vcf_obj, ID = paste0("SNP_", POS))
      message("Renamed SNP IDs to SNP_<POS> for uniqueness")
    }

    message("Reading LD matrix...")
    ld_obj <- read_LD(ld_file, vcf_obj)

    message("Reading phenotype...")
    pheno_obj <- read_pheno(pheno_file)

    message("Reading metadata...")
    meta_obj <- read_metadata(meta_grp_file)

    # Adjust MGmin if needed (don't exceed n_snps / 3)
    local_mgmin <- min(mgmin, max(3L, as.integer(n_snps / 3)))
    if (local_mgmin != mgmin) {
      message("Adjusted MGmin from ", mgmin, " to ", local_mgmin, " (region has ", n_snps, " SNPs)")
    }

    message("Running haplotyping with ", length(epsilon_range), " epsilon values...")
    hap_obj <- run_haplotyping(
      vcf = vcf_obj,
      LD = ld_obj,
      pheno = pheno_obj,
      metadata = meta_obj,
      epsilon = epsilon_range,
      MGmin = local_mgmin,
      minHap = minhap
    )

    # crosshap HapObject is a named list of per-epsilon results:
    #   e.g., list(Haplotypes_MGmin8_E0.4 = ..., Haplotypes_MGmin8_E0.6 = ...)
    # Check if any epsilon produced valid results
    if (is.null(hap_obj) || length(hap_obj) == 0) {
      message("WARNING: No marker groups found for region ", rid,
              " (SNPs may be too sparse or not in LD). Skipping.")
      unlink(tmp_dir, recursive = TRUE)
      next
    }
    n_valid_eps <- sum(sapply(hap_obj, function(x) !is.null(x) && length(x) > 0))
    message("Found results for ", n_valid_eps, "/", length(epsilon_range), " epsilon values")

    # Save HapObject
    hap_file <- file.path(intermediate_dir, paste0("Region_", rid, "_HapObject.qs"))
    qsave(hap_obj, hap_file)
    message("Saved HapObject: ", hap_file)

    # Generate clustree visualization (requires >= 2 epsilon results)
    if (n_valid_eps >= 2) {
      tryCatch({
        message("Generating clustree visualizations...")

        # Marker Group clustree
        p_mg <- clustree_viz(hap_obj, type = "MG")
        clustree_mg_file <- file.path(plots_dir, paste0("Region_", rid, "_clustree_MG"))
        qsave(p_mg, paste0(clustree_mg_file, ".qs"))
        ggsave(paste0(clustree_mg_file, ".png"), p_mg, width = 10, height = 8, dpi = 150)
        ggsave(paste0(clustree_mg_file, ".svg"), p_mg, width = 10, height = 8)
        message("Saved MG clustree: ", clustree_mg_file)

        # Haplotype clustree
        p_hap <- clustree_viz(hap_obj, type = "hap")
        clustree_hap_file <- file.path(plots_dir, paste0("Region_", rid, "_clustree_hap"))
        qsave(p_hap, paste0(clustree_hap_file, ".qs"))
        ggsave(paste0(clustree_hap_file, ".png"), p_hap, width = 10, height = 8, dpi = 150)
        ggsave(paste0(clustree_hap_file, ".svg"), p_hap, width = 10, height = 8)
        message("Saved hap clustree: ", clustree_hap_file)
      }, error = function(e) {
        message("WARNING: Clustree visualization failed for region ", rid, ": ", e$message)
      })
    } else {
      message("NOTE: Only ", n_valid_eps, " epsilon(s) produced results; ",
              "clustree requires >= 2. Skipping clustree visualization.")
    }

  }, error = function(e) {
    # All crosshap internal errors are non-fatal — log as WARNING and continue
    message("WARNING: Haplotyping failed for region ", rid, ": ", e$message)
  })

  # Clean up temp files
  unlink(tmp_dir, recursive = TRUE)
}

message("\n=== Haplotype scan complete ===")
