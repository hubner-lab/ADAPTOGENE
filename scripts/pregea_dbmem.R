#!/usr/bin/env Rscript
# pregea_dbmem.R — SHARED spatial-vector artifact (dbMEM eigenvectors),
# reused by preGEA's varpart block AND (follow-up, see docs/rda_research.md
# A20/C4) Gradient Forest, replacing gradient_forest_model.R's ad hoc
# raw-Euclidean-dist + "keep half the positive eigenvalues" heuristic.
#
# Ports rda_varpart_lasky.Rmd's PCNM block (geodesic distance via terra,
# MST-truncation threshold, keep-all-positive-eigenvalues), upgraded from
# vegan::pcnm() to adespatial::dbmem() per the spec.
#
# Inputs are COORDINATES ONLY (deliberate) — cheap enough to be pulled into
# maladaptation mode on demand once the Gradient Forest swap lands, without
# dragging genotype work into that DAG.
#
# SPATIAL_LEVEL='auto': site level (avoids the pseudoreplication warned about
# in rda_varpart_lasky.Rmd — climate/spatial predictors are constant within a
# site, so sample-level analysis replicates each value ~N times per site and
# artificially inflates the shared variance fraction) when the dataset has
# enough distinct sites; else sample level with a loud diagnostic note.
# Site-level MEMs are BROADCAST back to samples by coordinate/site match, so
# the output artifact is ALWAYS sample-indexed regardless of level chosen —
# consumers (varpart, future GF) must match() on `sample`, never assume row
# order (a wider coord-valid sample set may be joining against this
# climate-valid one).

suppressPackageStartupMessages({
    library(data.table)
    library(vegan)
    library(adespatial)
    library(terra)
    library(ggplot2)
})

source("/pipeline/scripts/R/utils/theme_adaptogene.R")

args <- commandArgs(trailingOnly = TRUE)
################################################################################
METADATA_CLIMATE_VALID <- args[1]
CLIMATE_VALID_SAMPLES  <- args[2]
SPATIAL_LEVEL          <- args[3]            # "auto" | "site" | "sample"
MIN_SITES              <- as.integer(args[4])
OUT_VECTORS_TSV        <- args[5]
OUT_DIAG_TSV           <- args[6]
OUT_PNG                <- args[7]
INTER_DIR              <- args[8]
################################################################################

for (d in c(dirname(OUT_VECTORS_TSV), dirname(OUT_DIAG_TSV), dirname(OUT_PNG), INTER_DIR))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
OUT_SVG <- sub("\\.png$", ".svg", OUT_PNG)

diag <- list()
add_diag <- function(k, v) diag[[k]] <<- v

meta <- fread(METADATA_CLIMATE_VALID, colClasses = c(site = "character", sample = "character"))
sample_order <- fread(CLIMATE_VALID_SAMPLES, header = FALSE, colClasses = "character")$V2
meta <- meta[match(sample_order, sample)]   # enforce CLIMATE_VALID_SAMPLES order
n_rows <- nrow(meta)

coords_key <- paste(meta$latitude, meta$longitude)
n_unique_coords <- length(unique(coords_key))
add_diag("n_rows", n_rows)
add_diag("n_unique_coords", n_unique_coords)

level <- SPATIAL_LEVEL
if (level == "auto") level <- if (n_unique_coords >= MIN_SITES) "site" else "sample"
add_diag("spatial_level", level)
message("INFO: preGEA dbMEM — level=", level, " n_rows=", n_rows, " n_unique_coords=", n_unique_coords,
       " (min_sites=", MIN_SITES, ")")

write_skip <- function(status) {
    add_diag("status", status)
    add_diag("mst_threshold_m", NA_real_); add_diag("mst_threshold_km", NA_real_)
    add_diag("n_mem_total", 0L); add_diag("n_mem_positive", 0L)
    out <- meta[, .(sample, site)]
    fwrite(out, OUT_VECTORS_TSV, sep = "\t", quote = FALSE)
    diag_dt <- data.table(names(diag), sapply(diag, as.character)); setnames(diag_dt, c("key", "value"))
    fwrite(diag_dt, OUT_DIAG_TSV, sep = "\t", quote = FALSE)
    # Empty-state placeholder — this IS the entire plot content (CLAUDE.md Rule 8 exception).
    g <- ggplot() + annotate("text", x = 0, y = 0, label = paste0(
            "dbMEM skipped: ", status, "\n(", n_unique_coords, " distinct coordinate(s), need >= ", MIN_SITES, ")")) +
        theme_void()
    ggsave(OUT_PNG, g, width = 6, height = 4, dpi = 150)
    ggsave(OUT_SVG, g, width = 6, height = 4, device = svglite::svglite, bg = "transparent")
    message("INFO: dbMEM skipped (", status, ") — wrote 0-MEM vectors file")
}

if (n_unique_coords < MIN_SITES) {
    write_skip("skipped_too_few_sites")
    quit(status = 0)
}

################################################################################
# Build the point set at the chosen level
################################################################################
if (level == "site") {
    site_coords <- unique(meta[, .(site, latitude, longitude)])
    pt_ids <- site_coords$site
    lon <- site_coords$longitude; lat <- site_coords$latitude
} else {
    pt_ids <- meta$sample
    lon <- meta$longitude; lat <- meta$latitude
}
n_pts <- length(pt_ids)

################################################################################
# Geodesic distance (terra::distance on EPSG:4326 — NOT raw Euclidean, unlike
# gradient_forest_model.R's current heuristic) + MST truncation
################################################################################
pts_sv <- terra::vect(cbind(lon, lat), crs = "EPSG:4326")
d_mat  <- as.matrix(terra::distance(pts_sv))
diag(d_mat) <- 0
rownames(d_mat) <- colnames(d_mat) <- pt_ids

sp_tree   <- vegan::spantree(as.dist(d_mat))
threshold <- max(sp_tree$dist)
add_diag("mst_threshold_m", threshold)
add_diag("mst_threshold_km", threshold / 1000)
message("INFO: MST truncation threshold = ", round(threshold / 1000, 2), " km")

################################################################################
# adespatial::dbmem() — keep ALL positive-eigenvalue MEMs
################################################################################
mem_obj <- tryCatch(
    adespatial::dbmem(as.dist(d_mat), thresh = threshold, MEM.autocor = "positive", silent = TRUE),
    error = function(e) { message("WARNING: dbmem() failed: ", e$message); NULL }
)

if (is.null(mem_obj) || ncol(as.data.frame(mem_obj)) == 0) {
    write_skip("degenerate_rank")
    quit(status = 0)
}

mem_df <- as.data.frame(mem_obj)
n_mem <- ncol(mem_df)
mem_eigs <- attr(mem_obj, "values")
n_mem_positive <- if (!is.null(mem_eigs)) sum(mem_eigs > 0) else n_mem
colnames(mem_df) <- paste0("MEM", seq_len(n_mem))
mem_df$.pt_id <- pt_ids

add_diag("n_mem_total", n_mem)
add_diag("n_mem_positive", n_mem_positive)
add_diag("status", "ok")
message("INFO: dbMEM produced ", n_mem, " positive-autocorrelation MEM(s)")

################################################################################
# Broadcast site-level MEMs back to samples (always sample-indexed output)
################################################################################
if (level == "site") {
    out <- merge(meta[, .(sample, site)], mem_df, by.x = "site", by.y = ".pt_id", all.x = TRUE)
    setcolorder(out, c("sample", "site", paste0("MEM", seq_len(n_mem))))
} else {
    out <- merge(meta[, .(sample, site)], mem_df, by.x = "sample", by.y = ".pt_id", all.x = TRUE)
}
out <- out[match(sample_order, sample)]   # restore CLIMATE_VALID_SAMPLES order
fwrite(out, OUT_VECTORS_TSV, sep = "\t", quote = FALSE)

diag_dt <- data.table(names(diag), sapply(diag, as.character)); setnames(diag_dt, c("key", "value"))
fwrite(diag_dt, OUT_DIAG_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote ", nrow(out), " sample rows x ", n_mem, " MEM(s) to ", OUT_VECTORS_TSV)

################################################################################
# Screeplot
################################################################################
eig_dt <- data.table(mem = seq_len(n_mem), eigenvalue = if (!is.null(mem_eigs)) mem_eigs[seq_len(n_mem)] else rep(NA_real_, n_mem))
g <- ggplot(eig_dt, aes(x = factor(mem), y = eigenvalue)) +
    geom_col(fill = ADAPT_RETAINED) +
    labs(x = "MEM", y = "Eigenvalue (Moran's I)",
        title = "dbMEM spatial eigenvector screeplot",
        subtitle = sprintf("level=%s | %d positive-autocorrelation MEM(s) | MST threshold=%.1f km",
                           level, n_mem_positive, threshold / 1000)) +
    theme_adaptogene()
ggsave(OUT_PNG, g, width = 7, height = 5, dpi = 300)
ggsave(OUT_SVG, g, width = 7, height = 5, device = svglite::svglite, bg = "transparent")

message("INFO: preGEA dbMEM complete.")
