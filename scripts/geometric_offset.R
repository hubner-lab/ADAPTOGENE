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
# 18  OUT_DIAGNOSTICS — output: geometric_offset_diagnostics.tsv (conditioning of
#                       the environment covariance + per-eigen-direction share of
#                       the realised offset; see section 10b)
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
OUT_DIAGNOSTICS <- args[18]
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

# ── 10b. Conditioning diagnostic for the offset and its contribution plot ─────
# WHY THIS EXISTS.  LEA computes the offset as a quadratic form in the predictor
# space (genetic_gap.R):
#
#     gg  = rowSums(((X.new - X.pred) %*% t(B))^2) / nrow(B)
#     eig = eigen(cov(B))
#
# and B, the LFMM2 effect-size matrix, is fitted with a RIDGE-REGULARISED inverse
# of the environment cross-product (lfmm2.R, lambda defaults to 1e-5):
#
#     B = (t(Ys - W) %*% Xs) %*% solve(t(Xs) %*% Xs + diag(lambda))
#
# So an environment direction whose singular value^2 is far below lambda is
# amplified by up to ~1/lambda = 1e5 in B, and therefore in the offset. With more
# predictors than distinct site environments the environment covariance is
# rank-deficient and this fires hard. Measured on Trifolium86 (18 predictors,
# 11 sites, condition number ~1e6): 93-99% of the realised offset sat on ONE
# eigenvector with 99.85% of its mass on environment PC10, a direction carrying
# 7.9e-7 of the environmental variance. Both the site RANKING and the
# contribution bar plot below are then statements about the conditioning of the
# predictor block, not about the markers -- and nothing in any output said so.
#
# WHAT THIS DECOMPOSES, EXACTLY.  gap$eigenvalues/gap$vectors are the eigen-pair
# of cov(B) -- COLUMN-CENTRED and divided by L-1 -- while the offset's own matrix
# is the uncentred t(B) %*% B / L. They differ by the rank-1 column-mean term of
# B. genetic.gap() does not return B, and recomputing it would mean a second
# LFMM2 fit on what is already the critical path, so we decompose cov(B) and
# REPORT how well that reproduces gap$offset. If it does not reproduce it, the
# per-direction rows are withheld rather than published as if they described this
# offset. Expect a small non-zero error (~1e-3), not exact agreement.
ENV_VAR_TOL <- 1e-3   # share of total environmental variance below which a direction is "null"
RECON_TOL   <- 1e-2   # max relative reconstruction error at which the decomposition is publishable
diag_dt <- data.table(quantity = character(), value = character())
add_gdiag <- function(k, v) diag_dt <<- rbind(diag_dt, data.table(quantity = k, value = as.character(v)))
dir_dt <- data.table()

tryCatch({
    evals <- as.numeric(gap$eigenvalues)         # length = n predictors (NOT K)
    vecs  <- as.matrix(gap$vectors)              # predictors x predictors, orthonormal
    D_raw <- as.matrix(pred.env) - as.matrix(new.env)
    sdv   <- apply(as.matrix(env), 2, sd)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    # LEA centres and scales X/X.new/X.pred by env's OWN mean and sd when
    # scale = TRUE (genetic_gap.R). A centre cancels in a difference, so only the
    # sd matters. Rather than trust SCALE, pick whichever space reproduces
    # gap$offset and report the error.
    cands <- list(raw = D_raw, scaled = sweep(D_raw, 2, sdv, '/'))
    errs  <- vapply(cands, function(D) {
        recon <- as.numeric(((D %*% vecs)^2) %*% evals)
        max(abs(recon - gap$offset)) / max(abs(gap$offset))
    }, numeric(1))
    space <- names(which.min(errs))
    D     <- cands[[space]]
    E     <- if (space == 'raw') as.matrix(env) else sweep(as.matrix(env), 2, sdv, '/')
    add_gdiag('decomposition_matrix', 'cov(B) (column-centred); offset uses t(B)%*%B/L')
    add_gdiag('decomposition_space', space)
    add_gdiag('reconstruction_rel_error', signif(min(errs), 3))
    add_gdiag('reconstruction_tolerance', RECON_TOL)

    S      <- stats::cov(E)
    s_eig  <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    cond   <- max(s_eig) / max(min(s_eig), .Machine$double.eps)
    n_env_rows <- nrow(unique(as.data.table(E)))
    add_gdiag('n_predictors', ncol(E))
    add_gdiag('n_rows_training_env', nrow(E))
    add_gdiag('n_distinct_env_rows', n_env_rows)
    add_gdiag('env_cov_rank_numeric', sum(s_eig > max(s_eig) * 1e-8))
    add_gdiag('env_cov_condition_number', signif(cond, 4))

    if (min(errs) > RECON_TOL) {
        # The eigen-pair does not describe THIS offset well enough to attribute
        # shares to directions. Say so instead of publishing a decomposition of
        # something else.
        add_gdiag('status', paste0('reconstruction failed (rel. error ', signif(min(errs), 3),
                                   ' > ', RECON_TOL, ') - per-direction decomposition withheld'))
        message('WARNING: could not reconstruct gap$offset from (eigenvalues, vectors) to better ',
                'than rel. error ', signif(min(errs), 3), ' (tolerance ', RECON_TOL,
                '). Per-direction offset shares are NOT reported for this run; the environment ',
                'conditioning numbers above are still valid.')
    } else {
        proj      <- D %*% vecs                                  # rows x predictors
        realised  <- evals * colSums(proj^2)                     # whole map (sites + raster cells)
        realised_site <- evals * colSums(proj[seq_len(n_site), , drop = FALSE]^2)
        env_var   <- vapply(seq_len(ncol(vecs)),
                            function(j) as.numeric(t(vecs[, j]) %*% S %*% vecs[, j]), numeric(1))
        env_share <- env_var / sum(diag(S))
        dir_dt <- data.table(
            direction         = seq_along(evals),
            eigenvalue        = evals,
            env_var_share     = env_share,
            offset_share      = realised / sum(realised),
            offset_share_site = realised_site / sum(realised_site),
            top_predictor     = PREDICTORS[apply(abs(vecs), 2, which.max)],
            top_loading_frac  = apply(vecs^2, 2, max)
        )[order(-offset_share)]

        null_share      <- sum(dir_dt$offset_share[dir_dt$env_var_share < ENV_VAR_TOL])
        null_share_site <- sum(dir_dt$offset_share_site[dir_dt$env_var_share < ENV_VAR_TOL])
        add_gdiag('offset_share_from_null_directions', signif(null_share, 4))
        add_gdiag('offset_share_from_null_directions_site', signif(null_share_site, 4))
        add_gdiag('env_var_tolerance', ENV_VAR_TOL)
        add_gdiag('status', 'ok')
        if (max(null_share, null_share_site) > 0.5) {
            message('WARNING: ', round(100 * null_share, 1), '% of the realised geometric offset ',
                    '(', round(100 * null_share_site, 1), '% at the sampled sites) comes from ',
                    'eigen-directions carrying < ', 100 * ENV_VAR_TOL, '% of the environmental ',
                    'variance each. The environment covariance is ill-conditioned (', ncol(E),
                    ' predictors, ', n_env_rows, ' distinct environment rows, condition number ',
                    signif(cond, 3), '), so the ridge-regularised inverse of t(X)X inside ',
                    'LEA::lfmm2() (lambda = 1e-5) is amplifying numerically null directions. ',
                    'Both the offset RANKING and the predictor-contribution plot are then set by ',
                    'the environment, not by the markers. Reduce/decorrelate Climate.predictors ',
                    'before interpreting either.')
        }
    }
}, error = function(e) {
    message('WARNING: offset conditioning diagnostic failed (', conditionMessage(e), ')')
    add_gdiag('status', paste0('failed: ', conditionMessage(e)))
})

if (nrow(dir_dt) > 0) {
    diag_dt <- rbind(diag_dt, data.table(
        quantity = paste0('direction', dir_dt$direction,
                          '_[eigenvalue|env_var_share|offset_share|offset_share_site|top_predictor]'),
        value = paste(signif(dir_dt$eigenvalue, 4), signif(dir_dt$env_var_share, 4),
                      signif(dir_dt$offset_share, 4), signif(dir_dt$offset_share_site, 4),
                      dir_dt$top_predictor, sep = ' | ')))
}
fwrite(diag_dt, OUT_DIAGNOSTICS, sep = '\t')
message(paste0('INFO: Saved offset conditioning diagnostics: ', OUT_DIAGNOSTICS))

# ── 11. Write importance plot (per-predictor contribution from eigenvalues) ────
# Produces a bar plot of predictor contributions; writes a placeholder on failure.
# NOTE: this is diag(cov(B)) — the per-predictor spread of LFMM2 effect sizes
# across candidate loci, i.e. the offset produced by a unit shift in one
# predictor, summed over EVERY eigen-direction including the numerically null
# ones. Read it together with geometric_offset_diagnostics.tsv (section 10b) —
# when offset_share_from_null_directions is large these bars are an artefact of
# an ill-conditioned predictor block, not a marker result.
tryCatch({
    # gap$eigenvalues / gap$vectors: eigen-pair of cov(B), so both are sized by the
    # number of PREDICTORS (not by K). vecs is predictors x predictors, orthonormal.
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
