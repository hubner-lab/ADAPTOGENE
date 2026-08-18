library(LEA)
library(data.table)
library(dplyr)
library(terra)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
##############################################################################
# Arg order (17 positional):
#  1  LFMM_IMP_FULL   — full imputed genotype matrix (.lfmm, rows=samples cols=SNPs)
#  2  VCFSNP          — SNP chr/pos table (.vcfsnp, no header, V1=chr V2=pos)
#  3  REMOVED         — removed SNPs table (.removed, no header)
#  4  CANDIDATE_SNPS  — selected_snps.tsv with SNPID column (chr:pos format)
#  5  ENV_SITE_PRES   — present climate at sampling sites (climate_present_site.tsv)
#  6  ENV_SITE_FUT    — future  climate at sampling sites (climate_future_site.tsv)
#  7  ENV_ALL_PRES    — present climate for all raster cells (climate_present_all.tsv)
#  8  ENV_ALL_FUT     — future  climate for all raster cells (climate_future_all.tsv)
#  9  PRES_RASTER     — present climate raster (.tif, used as spatial template)
# 10  SAMPLES         — metadata TSV (site, sample, latitude, longitude, ...)
# 11  PREDICTORS      — comma-separated predictor names
# 12  K               — number of latent factors (integer, or '' to error-check)
# 13  SCALE           — 'TRUE' or 'FALSE' — standardise environment before fitting
# 14  OUT_SITE_TSV    — output: genetic_offset_site.tsv
# 15  OUT_MAP_TSV     — output: genetic_offset_map.tsv (matrix, no header)
# 16  OUT_RASTER_TIF  — output: offset raster (.tif)
# 17  OUT_IMPORTANCE  — output: predictor importance plot (.png)
##############################################################################
LFMM_IMP_FULL  <- args[1]
VCFSNP         <- args[2]
REMOVED        <- args[3]
CANDIDATE_SNPS <- args[4]
ENV_SITE_PRES  <- args[5]
ENV_SITE_FUT   <- args[6]
ENV_ALL_PRES   <- args[7]
ENV_ALL_FUT    <- args[8]
PRES_RASTER    <- args[9]
SAMPLES        <- args[10]
PREDICTORS     <- args[11] %>% str_split(',') %>% unlist
K              <- args[12] %>% as.integer
SCALE          <- as.logical(args[13])
OUT_SITE_TSV   <- args[14]
OUT_MAP_TSV    <- args[15]
OUT_RASTER_TIF <- args[16]
OUT_IMPORTANCE <- args[17]
##############################################################################

set.seed(42)
message('INFO: Geometric Genetic Offset (LEA::genetic.gap)')
message(paste0('INFO: K=', K, '  scale=', SCALE))

# ── 1. Load full imputed genotype matrix ─────────────────────────────────────
message('INFO: Loading full imputed LFMM matrix')
lfmm_imp <- as.matrix(fread(LFMM_IMP_FULL))
message(paste0('INFO: LFMM matrix: ', nrow(lfmm_imp), ' samples x ', ncol(lfmm_imp), ' SNPs'))

# ── 2. Build chr:pos vector for all columns ───────────────────────────────────
vcfsnp_dt <- fread(VCFSNP, header = FALSE) %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chrpos = paste0(chr, ':', pos))
vcfsnp_full <- vcfsnp_dt$chrpos

# ── 3. CRITICAL ALIGNMENT ASSERTION ──────────────────────────────────────────
# The full imputed matrix columns MUST correspond 1-to-1 with vcfsnp_full rows.
# A silent column reorder would corrupt every candidate.loci index.
if (ncol(lfmm_imp) != nrow(vcfsnp_dt)) {
    stop(sprintf(
        paste0('FATAL: Column count mismatch - lfmm_imp has %d columns but vcfsnp has %d rows. ',
               'Imputation must preserve SNP order.'),
        ncol(lfmm_imp), nrow(vcfsnp_dt)
    ))
}
# Spot-check: first SNP in vcfsnp_full must match first column name (if columns named)
if (!is.null(colnames(lfmm_imp)) && colnames(lfmm_imp)[1] != '' &&
        !startsWith(colnames(lfmm_imp)[1], 'V')) {
    expected <- vcfsnp_full[1]
    actual   <- colnames(lfmm_imp)[1]
    if (!grepl(sub('.*:', '', expected), actual, fixed = TRUE)) {
        stop(sprintf(
            'FATAL: First column name "%s" does not match first SNP "%s". Column order mismatch.',
            actual, expected
        ))
    }
}
message('INFO: Alignment assertion passed')

# ── 4. Build removed set ──────────────────────────────────────────────────────
removed_dt <- fread(REMOVED, header = FALSE)
if (!is.null(removed_dt) && nrow(removed_dt) > 0) {
    removed <- removed_dt %>%
        dplyr::select(V1, V2) %>%
        setNames(c('chr', 'pos')) %>%
        dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>%
        .$chrpos
} else {
    removed <- character(0)
}

# ── 5. Load candidate SNPs and derive INTEGER index vector ────────────────────
# CRITICAL: pass the FULL matrix to genetic.gap() — restrict via candidate.loci
# integer indices, never by subsetting the matrix (would break latent factor estimation).
sig_dt   <- fread(CANDIDATE_SNPS, header = TRUE,
                  colClasses = c('SNPID' = 'character'))
sig_snps <- sig_dt$SNPID   # expected format: chr:pos

candidate <- which((!vcfsnp_full %in% removed) & (vcfsnp_full %in% sig_snps))
message(paste0('INFO: Candidate loci: ', length(candidate),
               ' of ', ncol(lfmm_imp), ' total SNPs'))
if (length(candidate) == 0) {
    stop('FATAL: No candidate loci found. Check SNPID format in ', CANDIDATE_SNPS,
         ' (expected chr:pos) vs vcfsnp file.')
}

# ── 6. Load climate data ──────────────────────────────────────────────────────
# Both _site tables carry a leading `sample` column (download_climate_present.R /
# merge_climate_future.R / stage_custom_climate.R all write it). Keep the raw tables:
# the sample IDs are the only key linking these rows to a sample identity, and section 8
# writes the output against them instead of trusting positional alignment.
env_site_pres_raw <- fread(ENV_SITE_PRES, colClasses = c('sample' = 'character'))
env_site_fut_raw  <- fread(ENV_SITE_FUT,  colClasses = c('sample' = 'character'))
for (nm in c('ENV_SITE_PRES', 'ENV_SITE_FUT')) {
    tbl <- if (nm == 'ENV_SITE_PRES') env_site_pres_raw else env_site_fut_raw
    if (!'sample' %in% names(tbl)) {
        stop('FATAL: ', get(nm), ' has no `sample` column (found: ',
             paste(utils::head(names(tbl), 5), collapse = ', '),
             '). Expected schema: sample <tab> one column per predictor, one row per sample. ',
             'This is likely a stale table from an older pipeline version -- re-run mode=structure.')
    }
}

env_site_pres <- env_site_pres_raw %>% dplyr::select(!!PREDICTORS)
env_site_fut  <- env_site_fut_raw  %>% dplyr::select(!!PREDICTORS)
env_all_pres  <- fread(ENV_ALL_PRES)  %>% dplyr::select(ID, !!PREDICTORS)
env_all_fut   <- fread(ENV_ALL_FUT)   %>% dplyr::select(ID, !!PREDICTORS)

# NOTE: climate_{present,future}_site.tsv are one row per SAMPLE (each site's climate values
# repeated for every sample at that site), so this is a sample count, not a site count -- it is
# the number of leading rows of gap$offset that belong to sampled individuals rather than to
# landscape raster cells. The sample IDs in column 1 are what section 8 keys the output on.
n_sample <- nrow(env_site_pres)
message(paste0('INFO: Sample rows: ', n_sample,
               '  Raster rows: ', nrow(env_all_pres)))

# Defense-in-depth: env_site_pres/env_site_fut and lfmm_imp are all supposed to already be
# row-aligned to the same climate-valid sample set (subset_lfmm_climate + the climate_site
# tables narrowed by download_climate_present/merge_climate_future). This should never fire --
# if it does, something upstream desynced the sample sets and feeding mismatched rows into
# LEA::genetic.gap() would silently bind the wrong sample's genotypes to the wrong climate
# values. Fail loud instead.
if (nrow(env_site_fut) != n_sample || nrow(lfmm_imp) != n_sample) {
    stop(paste0('FATAL: Sample-count mismatch between climate and genotype inputs -- ',
                'env_site_pres has ', n_sample, ' rows, env_site_fut has ', nrow(env_site_fut),
                ' rows, lfmm_imp has ', nrow(lfmm_imp), ' rows. These must match (all three ',
                'should already be subset to the same climate-valid sample set). Re-run ',
                'filter_climate_valid_samples/subset_lfmm_climate and check for a stale ',
                'intermediate file.'))
}

# The count guard above cannot see a REORDERING. offset_i is computed from the pair
# (env_site_pres[i], env_site_fut[i]) -- if the two tables list the same samples in a
# different order, every offset is a present-of-one-sample vs future-of-another delta:
# a wrong NUMBER, not merely a wrong label, and silent. Both tables carry the sample key,
# so check it. (lfmm_imp has no IDs -- its alignment remains count-checked only, guaranteed
# upstream by subset_lfmm_climate operating on the same climate-valid sample list.)
if (!identical(env_site_pres_raw$sample, env_site_fut_raw$sample)) {
    mismatched <- which(env_site_pres_raw$sample != env_site_fut_raw$sample)
    stop('FATAL: ', ENV_SITE_PRES, ' and ', ENV_SITE_FUT, ' list the same number of rows but ',
         'not in the same sample order (', length(mismatched), ' position(s) differ, first at ',
         'row ', mismatched[1], ': "', env_site_pres_raw$sample[mismatched[1]], '" vs "',
         env_site_fut_raw$sample[mismatched[1]], '"). genetic.gap() pairs these two tables ',
         'positionally, so proceeding would compute each offset from one sample\'s present ',
         'climate and another sample\'s future climate. Both tables must be written in the ',
         'row order of metadata_climate_valid.tsv -- re-run mode=structure.')
}

# Drop raster cells with any NA across present+future predictors
ok <- complete.cases(env_all_pres[, -1]) & complete.cases(env_all_fut[, -1])
env_all_pres_ok <- env_all_pres[ok, ]
env_all_fut_ok  <- env_all_fut[ok, ]
message(paste0('INFO: Non-NA raster cells: ', sum(ok)))

# Stack: sampled individuals first (top n_sample rows), then landscape cells
new.env  <- as.matrix(rbind(env_site_pres,      env_all_pres_ok[, -1]))
pred.env <- as.matrix(rbind(env_site_fut,       env_all_fut_ok[,  -1]))
env      <- as.matrix(env_site_pres)             # training environment = present at sites

# ── 7. Compute genetic.gap() ──────────────────────────────────────────────────
message('INFO: Running LEA::genetic.gap() — this may take several minutes on large datasets')
gap <- LEA::genetic.gap(
    input          = lfmm_imp,
    env            = env,
    new.env        = new.env,
    pred.env       = pred.env,
    K              = K,
    scale          = SCALE,
    candidate.loci = candidate
)

message(paste0('INFO: genetic.gap complete. Offset length: ', length(gap$offset)))

# Slice: first n_sample rows = sampled individuals, remainder = landscape
offset_site      <- gap$offset[seq_len(n_sample)]
offset_landscape <- gap$offset[(n_sample + 1):length(gap$offset)]

message(paste0('INFO: Site offset range: ',
               round(min(offset_site, na.rm = TRUE), 4), ' – ',
               round(max(offset_site, na.rm = TRUE), 4)))

# ── 8. Write site TSV (byte-compatible with GF output) ───────────────────────
samples <- fread(SAMPLES,
                 colClasses = c('site' = 'character', 'sample' = 'character',
                                'latitude' = 'numeric',  'longitude' = 'numeric'))

# Key the offsets on the sample IDs that came with the climate table they were computed from
# (column 1 of climate_present_site.tsv), then JOIN site from the metadata. The previous
# version cbind()-ed offset_site onto unique(samples[, c('site','sample')]) -- a purely
# positional bind whose only protection was the row-count check above, so any reordering of
# the metadata relative to the climate table would have silently attached each offset to the
# wrong sample while every count still matched.
go_site <- data.table(sample = env_site_pres_raw$sample,
                      genetic_offset = offset_site)
go_site <- merge(go_site, unique(samples[, c('site', 'sample')]),
                 by = 'sample', all.x = TRUE, sort = FALSE)
if (anyNA(go_site$site)) {
    stop('FATAL: ', sum(is.na(go_site$site)), ' sample ID(s) from ', ENV_SITE_PRES,
         ' have no row in ', SAMPLES, ' (e.g. ',
         paste(utils::head(go_site$sample[is.na(go_site$site)], 3), collapse = ', '),
         '). The climate table and the metadata must describe the same sample set.')
}
go_site <- go_site %>%
    dplyr::select(site, sample, genetic_offset) %>%
    dplyr::arrange(desc(genetic_offset))

fwrite(go_site, OUT_SITE_TSV, sep = '\t')
message(paste0('INFO: Saved site GO values: ', OUT_SITE_TSV))

# ── 9. Write raster (mirror gradient_forest_offset.R:47-56) ──────────────────
clim_pres   <- rast(PRES_RASTER)
rast_offset <- clim_pres[[PREDICTORS[1]]]

# Assign offset to the valid cell IDs
rast_offset[] <- NA
rast_offset[env_all_pres_ok$ID] <- offset_landscape

writeRaster(rast_offset,
            filename  = OUT_RASTER_TIF,
            overwrite = TRUE,
            gdal      = c('INTERLEAVE=BAND', 'COMPRESS=LZW'))
message(paste0('INFO: Saved offset raster: ', OUT_RASTER_TIF))

# ── 10. Write map TSV (mirror gradient_forest_offset.R:58-62) ─────────────────
GO_matrix <- matrix(values(rast_offset), ncol = ncol(rast_offset), byrow = TRUE)
GO_matrix %>%
    as.data.table() %>%
    fwrite(OUT_MAP_TSV, sep = '\t', col.names = FALSE)
message(paste0('INFO: Saved map GO values: ', OUT_MAP_TSV))

# ── 11. Write importance plot (per-predictor contribution from eigenvalues) ────
# Produces a bar plot of predictor contributions; writes a placeholder on failure.
tryCatch({
    # gap$eigenvalues: per-K eigenvalues; gap$vectors: loadings matrix (predictors x K)
    # Contribution = sum of squared loadings weighted by eigenvalue
    evals <- gap$eigenvalues
    vecs  <- gap$vectors
    contrib <- colSums(t(vecs^2) * evals)
    names(contrib) <- PREDICTORS

    png(OUT_IMPORTANCE, width = 800, height = 500)
    par(mar = c(8, 4, 3, 1))
    barplot(sort(contrib, decreasing = TRUE),
            las    = 2,
            col    = '#1B7A6E',
            main   = 'Geometric Offset: Predictor Contribution',
            ylab   = 'Weighted loading (sum λ·v²)',
            cex.names = 0.85)
    dev.off()
    message(paste0('INFO: Saved importance plot: ', OUT_IMPORTANCE))
}, error = function(e) {
    # Write placeholder so downstream piemap/summary targets always resolve
    message(paste0('WARNING: Importance plot failed (', conditionMessage(e),
                   '). Writing placeholder.'))
    png(OUT_IMPORTANCE, width = 400, height = 200)
    par(mar = c(1, 1, 2, 1))
    plot.new()
    text(0.5, 0.5, 'Importance plot unavailable\n(see log for details)',
         cex = 1.2, col = 'grey40')
    dev.off()
})

message('INFO: Geometric genetic offset complete')
