#!/usr/bin/env Rscript
# add_pregea_sites.R — one-off dev utility, NOT a Snakemake rule.
#
# SIMDATA originally has 3 sites (Negev/TelAviv/Galilee) — enough for GEA/GWAS
# testing but degenerate for the preGEA varpart/dbMEM block, whose whole
# purpose is spatial-eigenvector analysis (adespatial::dbmem() + ordiR2step +
# vegan::varpart). 3 points give at most 1-2 dbMEM vectors and no real
# geographic gradient to detect; dbMEM's hard floor of 3 sites (hardcoded in
# scripts/pregea_dbmem.R, not a config key) would push a bare-3-site SIMDATA
# run right up against the "too_few_sites" skip branch, leaving that code
# effectively untested on the pipeline's primary dataset.
#
# Adds 6 new sites (3 samples each -> 18 new samples, 9 sites total) INSIDE
# the existing coordinate bounding box (lat 30.854-33.0128, lon
# 34.7817-35.4985 -- see Map.zoom_extent / the 3 existing sites) so
# Climate.climate_extent: auto (data-derived) recomputes to the SAME bounding
# box and download_climate_present's WorldClim raster crop is reused, not
# re-downloaded.
#
# Genotypes are NOT random per new site -- that would give zero spatial
# autocorrelation, defeating the point of the exercise. Each new site's
# per-SNP allele frequency is inverse-distance-weighted (IDW) from the 3
# existing sites' empirical frequencies, producing a smooth geographic
# gradient (exactly what dbMEM/Moran's-eigenvector analysis is designed to
# detect). Phenotypes (height/flowering_time/disease_score) are interpolated
# the same way, with per-sample noise matched to the existing data's spread.
#
# Idempotent: aborts if any new site name already exists in metadata, rather
# than double-injecting.
#
# Usage: Rscript /pipeline/scripts/add_pregea_sites.R
# (run from the repo root / mounted /pipeline, so the relative data/ paths resolve)

suppressPackageStartupMessages(library(data.table))

set.seed(43)   # distinct from generate_simdata.R's 42 and add_related_samples.R's 42

VCF_PATH      <- "data/SIMDATA.vcf"
METADATA_PATH <- "data/SIMDATA_metadata.tsv"

# 6 new sites, all inside the existing bounding box
# (lat 30.854-33.0128, lon 34.7817-35.4985), well separated from each other
# and from the 3 existing sites (Negev 30.854/34.7826, TelAviv 32.0837/34.7817,
# Galilee 33.0128/35.4985).
NEW_SITES <- list(
    list(site = "Carmel",       lat = 32.70,   lon = 35.05,   n = 3, prefix = "CRM"),
    list(site = "Jerusalem",    lat = 31.7683, lon = 35.2137, n = 3, prefix = "JRS"),
    list(site = "DeadSea",      lat = 31.50,   lon = 35.40,   n = 3, prefix = "DSE"),
    list(site = "Samaria",      lat = 32.20,   lon = 35.25,   n = 3, prefix = "SAM"),
    list(site = "JordanValley", lat = 32.40,   lon = 35.45,   n = 3, prefix = "JOR"),
    list(site = "Sharon",       lat = 32.30,   lon = 34.85,   n = 3, prefix = "SHR")
)

#=============================================================================
# 0. Guards — read existing files, abort on re-run
#=============================================================================
meta <- fread(METADATA_PATH, colClasses = c(site = "character", sample = "character"))
new_site_names <- vapply(NEW_SITES, `[[`, character(1), "site")
if (any(new_site_names %in% meta$site)) {
    stop("Some new site name(s) already present in metadata (",
         paste(intersect(new_site_names, meta$site), collapse = ", "),
         ") -- add_pregea_sites.R already ran on this dataset. Aborting to avoid double injection.")
}

lines <- readLines(VCF_PATH)
header_idx <- which(startsWith(lines, "#CHROM"))
stopifnot(length(header_idx) == 1)
chrom_cols  <- strsplit(lines[header_idx], "\t")[[1]]
sample_cols <- chrom_cols[-(1:9)]

# Existing samples with valid, non-duplicate coordinates -- these define the
# 3 known site centroids the IDW interpolation anchors to. *_DUP samples
# share their donor's coordinates exactly, so including them would just
# double-weight those 3 points; excluded for a clean anchor set.
meta_coord <- meta[!is.na(latitude) & !is.na(longitude) & !grepl("_DUP$", sample)]
existing_sites <- unique(meta_coord[, .(site, latitude, longitude)])
stopifnot(nrow(existing_sites) == 3)   # Negev, TelAviv, Galilee -- fails loudly if the base fixture changed

message("INFO: Anchor sites for IDW interpolation:")
for (i in seq_len(nrow(existing_sites))) {
    message(sprintf("  %-10s lat=%.4f lon=%.4f", existing_sites$site[i],
                    existing_sites$latitude[i], existing_sites$longitude[i]))
}

# ── Inverse-distance weights: new site -> named vector over the 3 anchor sites
idw_weights <- function(lat, lon, anchors) {
    d <- sqrt((anchors$latitude - lat)^2 + (anchors$longitude - lon)^2)
    d[d == 0] <- 1e-6   # guard exact coincidence (shouldn't happen given the coords chosen)
    w <- 1 / d^2
    stats::setNames(w / sum(w), anchors$site)
}

#=============================================================================
# 1. VCF — per-existing-site allele frequency per SNP, IDW-interpolate for
#    each new site, generate new sample genotype columns
#=============================================================================
data_idx   <- (header_idx + 1):length(lines)
n_variants <- length(data_idx)
split_lines <- strsplit(lines[data_idx], "\t")

site_of_sample <- stats::setNames(meta_coord$site, meta_coord$sample)
site_sample_idx <- lapply(existing_sites$site, function(s) {
    which(sample_cols %in% names(site_of_sample)[site_of_sample == s]) + 9   # +9 fixed VCF cols
})
names(site_sample_idx) <- existing_sites$site

allele_count <- function(gt) {
    # "0/0"->0, "0/1"|"1/0"->1, "1/1"->2, "./."->NA. Biallelic SNPs only (SIMDATA convention).
    if (gt == "./.") return(NA_real_)
    a <- as.integer(strsplit(gt, "/", fixed = TRUE)[[1]])
    sum(a)
}

new_sample_names <- unlist(lapply(NEW_SITES, function(s) paste0(s$prefix, sprintf("%02d", seq_len(s$n)))))
new_header <- paste(c(chrom_cols, new_sample_names), collapse = "\t")

new_data_lines <- character(n_variants)
for (i in seq_len(n_variants)) {
    fields <- split_lines[[i]]

    site_freq <- vapply(existing_sites$site, function(s) {
        gts <- fields[site_sample_idx[[s]]]
        ac  <- vapply(gts, allele_count, numeric(1))
        ac  <- ac[!is.na(ac)]
        if (length(ac) == 0) 0.5 else mean(ac) / 2   # fallback: uninformative if a site is all-missing here
    }, numeric(1))

    new_gts <- character(0)
    for (s in NEW_SITES) {
        w <- idw_weights(s$lat, s$lon, existing_sites)
        p_alt <- sum(w[existing_sites$site] * site_freq[existing_sites$site])
        p_alt <- min(max(p_alt, 0.01), 0.99)   # keep away from fixation so covRob/rda never see zero-variance here
        for (k in seq_len(s$n)) {
            a1 <- rbinom(1, 1, p_alt)
            a2 <- rbinom(1, 1, p_alt)
            gt <- paste0(a1, "/", a2)
            if (runif(1) < 0.05) gt <- "./."   # background missingness, matching generate_simdata.R's normal rate
            new_gts <- c(new_gts, gt)
        }
    }
    new_data_lines[i] <- paste(c(fields, new_gts), collapse = "\t")
}

writeLines(c(lines[1:(header_idx - 1)], new_header, new_data_lines), VCF_PATH)
message("INFO: Added ", length(new_sample_names), " sample(s) across ", length(NEW_SITES),
       " new site(s) to ", VCF_PATH, " (", n_variants, " variants updated).")

#=============================================================================
# 2. Metadata — IDW-interpolated phenotypes + per-sample noise
#=============================================================================
site_pheno <- meta_coord[, .(height = mean(height, na.rm = TRUE),
                             flowering_time = mean(flowering_time, na.rm = TRUE),
                             disease_score  = mean(disease_score,  na.rm = TRUE)),
                         by = site]
# Per-trait SD across existing samples -> the noise scale for new samples,
# so the new sites' phenotype spread looks like the rest of the dataset
# rather than being suspiciously tight around the interpolated mean.
pheno_sd <- meta_coord[, .(height = stats::sd(height, na.rm = TRUE),
                           flowering_time = stats::sd(flowering_time, na.rm = TRUE),
                           disease_score  = stats::sd(disease_score,  na.rm = TRUE))]

new_rows <- list()
for (s in NEW_SITES) {
    w <- idw_weights(s$lat, s$lon, existing_sites)
    interp <- function(trait) sum(w[site_pheno$site] * site_pheno[[trait]][match(existing_sites$site, site_pheno$site)])
    mu <- c(height = interp("height"), flowering_time = interp("flowering_time"),
           disease_score = interp("disease_score"))
    for (k in seq_len(s$n)) {
        new_rows[[length(new_rows) + 1]] <- data.table(
            site = s$site, sample = paste0(s$prefix, sprintf("%02d", k)),
            latitude = s$lat, longitude = s$lon,
            height         = round(mu["height"]         + stats::rnorm(1, 0, pheno_sd$height),         1),
            flowering_time = round(mu["flowering_time"] + stats::rnorm(1, 0, pheno_sd$flowering_time), 1),
            disease_score  = round(mu["disease_score"]  + stats::rnorm(1, 0, pheno_sd$disease_score),  1)
        )
    }
}
new_meta <- rbindlist(new_rows)
# Match the existing dataset's occasional NA phenotype cells (not coords —
# these are the coord-VALID subset dbMEM/varpart needs) so the new sites
# exercise the same missing-value handling as the rest of the pipeline.
if (nrow(new_meta) >= 6) {
    na_idx <- sample(seq_len(nrow(new_meta)), 2)
    new_meta$disease_score[na_idx[1]] <- NA
    new_meta$flowering_time[na_idx[2]] <- NA
}

fwrite(rbind(meta, new_meta), METADATA_PATH, sep = "\t", quote = FALSE, na = "NA")
message("INFO: Added ", nrow(new_meta), " new metadata row(s) across ", length(NEW_SITES),
       " new site(s) to ", METADATA_PATH, ".")
message("INFO: SIMDATA now has ", nrow(meta) + nrow(new_meta), " samples (was ", nrow(meta), "), ",
       nrow(existing_sites) + length(NEW_SITES), " distinct sites (was ", nrow(existing_sites), ").")
