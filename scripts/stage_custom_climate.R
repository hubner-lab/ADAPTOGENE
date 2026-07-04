#!/usr/bin/env Rscript
# stage_custom_climate.R — Stage user-supplied environmental tables into pipeline-standard climate outputs.
# Enables Climate.source: custom (bypass WorldClim/CMIP6 download).
#
# Usage — present mode (9 positional args):
#   Rscript stage_custom_climate.R present \
#       METADATA ENV_TABLE COLUMNS KEY GRID_RES \
#       OUT_RASTER OUT_ALL OUT_SITE OUT_SITE_SCALED
#
# Usage — future mode (10 positional args):
#   Rscript stage_custom_climate.R future \
#       METADATA ENV_TABLE COLUMNS KEY \
#       PRESENT_RASTER PRESENT_ALL \
#       OUT_RASTER OUT_ALL OUT_SITE
#
# Args:
#   1  MODE            "present" | "future"
#   2  METADATA        Pipeline metadata TSV (O['metadata']): site sample latitude longitude ...
#   3  ENV_TABLE       User env table: site (or sample) + one column per predictor (named bio_<N>)
#   4  COLUMNS         Comma-separated predictor column names (must match env table headers)
#   5  KEY             Join column linking metadata to env table ("site" | "sample")
#   6  GRID_RES        Raster cell size in degrees (present mode only; future inherits present grid)
#   --- present mode continues ---
#   7  OUT_RASTER      Output .tif path (W['climate_raster'])
#   8  OUT_ALL         Output landscape table path (O['climate_all'])
#   9  OUT_SITE        Output per-sample table path (O['climate_site'])
#   10 OUT_SITE_SCALED Output per-sample scaled table path (O['climate_site_scaled'])
#   --- future mode: arg 6 = PRESENT_RASTER, arg 7 = PRESENT_ALL, then ---
#   8  OUT_RASTER      W['climate_future_raster']
#   9  OUT_ALL         O['climate_future_all']
#   10 OUT_SITE        O['climate_future_site']
#
# Output schemas (byte-compatible with download_climate_present.R / merge_climate_future.R):
#   _site.tsv         sample <tab> bio_* columns (one row per sample, metadata row order)
#   _site_scaled.tsv  sample <tab> bio_* columns, column z-scored via base::scale()
#   _all.tsv          ID <tab> bio_* columns (one row per non-NA raster cell, ID = cell index)
#   raster.tif        SpatRaster, one named layer per predictor, BAND/LZW compressed

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(terra)
})

args <- commandArgs(trailingOnly = TRUE)

MODE <- args[1]
if (!MODE %in% c("present", "future")) {
    stop("ERROR: MODE must be 'present' or 'future', got: ", MODE)
}

METADATA <- args[2]
ENV_TABLE <- args[3]
COLUMNS_STR <- args[4]
KEY <- args[5]

if (MODE == "present") {
    GRID_RES    <- as.numeric(args[6])
    OUT_RASTER  <- args[7]
    OUT_ALL     <- args[8]
    OUT_SITE    <- args[9]
    OUT_SCALED  <- args[10]
} else {
    PRESENT_RASTER <- args[6]
    PRESENT_ALL    <- args[7]
    OUT_RASTER     <- args[8]
    OUT_ALL        <- args[9]
    OUT_SITE       <- args[10]
}

# Parse predictor column list
COLS <- trimws(strsplit(COLUMNS_STR, ",")[[1]])
COLS <- COLS[nchar(COLS) > 0]
if (length(COLS) == 0) stop("ERROR: COLUMNS is empty — no predictor columns specified")
message("INFO: Predictors: ", paste(COLS, collapse = ", "))

#=============================================================================
# Read metadata (canonical sample order, lat/lon coords)
#=============================================================================
meta <- fread(METADATA,
              colClasses = c(site      = "character",
                             sample    = "character",
                             latitude  = "numeric",
                             longitude = "numeric"))
message("INFO: Metadata loaded — ", nrow(meta), " samples")

#=============================================================================
# Read and validate user env table
#=============================================================================
env <- fread(ENV_TABLE, colClasses = c(site = "character", sample = "character"))
message("INFO: Env table loaded — ", nrow(env), " rows, columns: ", paste(names(env), collapse = ", "))

# Validate that all predictor columns exist in env table
missing_cols <- setdiff(COLS, names(env))
if (length(missing_cols) > 0) {
    stop("ERROR: Predictor column(s) not found in env table: ",
         paste(missing_cols, collapse = ", "),
         "\nAvailable columns: ", paste(names(env), collapse = ", "))
}

# Validate KEY column exists in both tables
if (!KEY %in% names(meta)) stop("ERROR: KEY column '", KEY, "' not found in metadata")
if (!KEY %in% names(env))  stop("ERROR: KEY column '", KEY, "' not found in env table")

#=============================================================================
# Build per-sample site table (LEFT JOIN from metadata → env, metadata order)
# This preserves the row order that emmax.R / lfmm.R expect against the VCF.
#=============================================================================
# Unique per-key env rows (in case env has one row per site but metadata one per sample)
env_unique <- unique(env[, c(KEY, COLS), with = FALSE])

# Index: metadata key → env row
key_col_meta <- meta[[KEY]]
key_col_env  <- env_unique[[KEY]]
idx <- match(key_col_meta, key_col_env)

unmatched <- which(is.na(idx))
if (length(unmatched) > 0) {
    stop("ERROR: Metadata ", KEY, " values not found in env table: ",
         paste(unique(key_col_meta[unmatched]), collapse = ", "))
}

site_vals <- env_unique[idx, COLS, with = FALSE]

message("INFO: Site table built — ", nrow(site_vals), " rows")
if (nrow(site_vals) != nrow(meta)) {
    stop("ERROR: Site table row count (", nrow(site_vals),
         ") != metadata rows (", nrow(meta), ")")
}

# Write site table (unscaled): sample + bio_* columns, matching download_climate_present.R:196
site_dt <- cbind(data.table(sample = meta$sample), site_vals)

#=============================================================================
# Build site-coordinate lookup for raster construction (unique sites)
#=============================================================================
coords_dt <- unique(meta[, .(key = get(KEY), longitude, latitude)])
message("INFO: Distinct sites for raster: ", nrow(coords_dt))

#=============================================================================
# Raster construction
#=============================================================================

if (MODE == "present") {
    # ── Present raster: build geographic grid from site bounding box ──────────
    lon_range <- range(coords_dt$longitude)
    lat_range <- range(coords_dt$latitude)
    buf <- max(GRID_RES, 0.5)  # buffer >= 0.5° so sites don't fall on boundary
    ext_v <- terra::ext(lon_range[1] - buf, lon_range[2] + buf,
                        lat_range[1] - buf, lat_range[2] + buf)
    r0 <- terra::rast(ext_v, resolution = GRID_RES, crs = "EPSG:4326")
    message("INFO: Raster grid: ", nrow(r0), " rows x ", ncol(r0), " cols = ",
            terra::ncell(r0), " cells")

    # Map distinct sites to cell indices
    cells <- terra::cellFromXY(r0, cbind(coords_dt$longitude, coords_dt$latitude))
    if (anyNA(cells)) {
        stop("ERROR: Some site coordinates fall outside the raster extent — adjust grid_resolution or check coordinates")
    }
    if (length(unique(cells)) != nrow(coords_dt)) {
        stop("ERROR: Two distinct sites map to the same raster cell. Decrease Climate.custom.grid_resolution. ",
             "Site count: ", nrow(coords_dt), ", unique cells: ", length(unique(cells)))
    }

    # Build one layer per predictor with site values; all other cells = NA
    # Site value = env value keyed by site; env_unique already unique per key
    env_site_keyed <- setNames(as.data.frame(env_unique[match(coords_dt$key, env_unique[[KEY]]),
                                                         COLS, with = FALSE]),
                               COLS)

    layers <- lapply(seq_along(COLS), function(i) {
        lyr <- r0
        lyr[] <- NA_real_
        lyr[cells] <- env_site_keyed[[COLS[i]]]
        names(lyr) <- COLS[i]
        lyr
    })
    clim <- terra::rast(layers)
    names(clim) <- COLS

    message("INFO: Raster layers built: ", paste(names(clim), collapse = ", "))

    # Write raster (byte-identical compression to download_climate_present.R:186-189)
    dir.create(dirname(OUT_RASTER), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(clim, OUT_RASTER, overwrite = TRUE,
                       gdal = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
    message("INFO: Raster written → ", OUT_RASTER)

    # Build _all table: ID = cell index, predictor columns, non-NA cells only
    # Mirrors download_climate_present.R:124-127
    all_dt <- terra::as.data.frame(clim, cells = TRUE, na.rm = TRUE)
    names(all_dt)[names(all_dt) == "cell"] <- "ID"
    setcolorder(all_dt, c("ID", COLS))
    all_dt <- as.data.table(all_dt)

    message("INFO: _all table: ", nrow(all_dt), " non-NA cells")

} else {
    # ── Future raster: inherit geometry from present raster ───────────────────
    clim_pres <- terra::rast(PRESENT_RASTER)
    message("INFO: Present raster geometry: ", nrow(clim_pres), "r x ", ncol(clim_pres), "c = ",
            terra::ncell(clim_pres), " cells")

    # Remap site future values onto same cell positions as present raster
    # (cell positions come from present _all table ID column)
    pres_all <- fread(PRESENT_ALL, colClasses = c(ID = "integer"))
    message("INFO: Present _all table: ", nrow(pres_all), " cells")

    # Site values for future keyed by KEY → match to coords_dt
    env_site_keyed_fut <- as.data.frame(
        env_unique[match(coords_dt$key, env_unique[[KEY]]), COLS, with = FALSE])

    cells_pres <- terra::cellFromXY(clim_pres,
                                    cbind(coords_dt$longitude, coords_dt$latitude))
    if (anyNA(cells_pres)) {
        stop("ERROR: Site coordinates outside present raster extent in future staging")
    }

    # Build future raster with same geometry
    layers_fut <- lapply(seq_along(COLS), function(i) {
        lyr <- clim_pres[[1]]  # borrow geometry; values will be replaced
        lyr[] <- NA_real_
        lyr[cells_pres] <- env_site_keyed_fut[[COLS[i]]]
        names(lyr) <- COLS[i]
        lyr
    })
    clim <- terra::rast(layers_fut)
    names(clim) <- COLS

    # Explicit cell-count check (mirrors merge_climate_future.R:69-71)
    if (terra::ncell(clim) != terra::ncell(clim_pres)) {
        stop("ERROR: Future raster ncell (", terra::ncell(clim),
             ") != present raster ncell (", terra::ncell(clim_pres),
             ") — geometry must be identical")
    }

    message("INFO: Future raster ncell == present ncell: ", terra::ncell(clim), " ✓")

    # Write future raster
    dir.create(dirname(OUT_RASTER), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(clim, OUT_RASTER, overwrite = TRUE,
                       gdal = c("INTERLEAVE=BAND", "COMPRESS=LZW"))
    message("INFO: Future raster written → ", OUT_RASTER)

    # Build future _all: same row order as present _all (positional match for offset scripts)
    # geometric_offset.R:125-132 and gradient_forest_offset.R:37-41 use positional subtraction
    fut_vals <- terra::extract(clim, pres_all$ID)
    fut_vals$ID <- pres_all$ID
    setcolorder(fut_vals, c("ID", COLS))
    all_dt <- as.data.table(fut_vals)

    message("INFO: Future _all table: ", nrow(all_dt), " cells (row-matched to present _all)")
}

#=============================================================================
# Write outputs
#=============================================================================
dir.create(dirname(OUT_ALL),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_SITE), recursive = TRUE, showWarnings = FALSE)

# _all.tsv — landscape table (same schema for present/future)
fwrite(all_dt, OUT_ALL, sep = '\t', quote = FALSE)
message("INFO: _all written → ", OUT_ALL)

# _site.tsv — per-sample table (unscaled)
fwrite(site_dt, OUT_SITE, sep = '\t', quote = FALSE)
message("INFO: _site written → ", OUT_SITE)

if (MODE == "present") {
    # _site_scaled.tsv — column z-score via dplyr::mutate_all(scale)
    # Replicates download_climate_present.R:200 exactly: cbind(sample, mutate_all(scale))
    dir.create(dirname(OUT_SCALED), recursive = TRUE, showWarnings = FALSE)
    site_scaled_dt <- cbind(
        data.table(sample = meta$sample),
        as.data.table(site_vals %>% dplyr::mutate_all(scale))
    )
    fwrite(site_scaled_dt, OUT_SCALED, sep = '\t', quote = FALSE)
    message("INFO: _site_scaled written → ", OUT_SCALED)

    # Validate scaling
    bio_cols <- COLS
    col_means <- sapply(bio_cols, function(b) mean(site_scaled_dt[[b]], na.rm = TRUE))
    col_sds   <- sapply(bio_cols, function(b) sd(site_scaled_dt[[b]],   na.rm = TRUE))
    if (any(abs(col_means) > 1e-10) || any(abs(col_sds - 1) > 1e-10)) {
        message("WARNING: Scaled table validation — max |mean| = ",
                max(abs(col_means)), ", max |sd-1| = ", max(abs(col_sds - 1)))
    } else {
        message("INFO: Scaling validated — column means ≈ 0, sds ≈ 1 ✓")
    }
}

message("INFO: stage_custom_climate.R (", MODE, ") complete")
