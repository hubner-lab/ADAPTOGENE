library(LEA)
library(data.table)
library(dplyr)
library(terra)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
##############################################################################
# Arg order (17 positional, + optional 18th):
#  1  LFMM_IMP_FULL   — full imputed genotype matrix (.lfmm, rows=samples cols=SNPs)
#  2  VCFSNP          — SNP chr/pos table (.vcfsnp, no header, V1=chr V2=pos)
#  3  REMOVED         — removed SNPs table (.removed, no header)
#  4  CANDIDATE_SNPS  — selected_snps.tsv with SNPID column (chr:pos format)
#  5  ENV_SITE_PRES   — present climate at sampling sites (climate_present_site.tsv)
#  6  ENV_SITE_FUT    — future  climate at sites      — COMMA-SEPARATED, one per scenario
#  7  ENV_ALL_PRES    — present climate for all raster cells (climate_present_all.tsv)
#  8  ENV_ALL_FUT     — future  climate for all cells — COMMA-SEPARATED, one per scenario
#  9  PRES_RASTER     — present climate raster (.tif, used as spatial template)
# 10  SAMPLES         — metadata TSV (site, sample, latitude, longitude, ...)
# 11  PREDICTORS      — comma-separated predictor names
# 12  K               — number of latent factors (integer, or '' to error-check)
# 13  SCALE           — 'TRUE' or 'FALSE' — standardise environment before fitting
# 14  OUT_SITE_TSV    — output: genetic_offset_site.tsv  — COMMA-SEPARATED, one per scenario
# 15  OUT_MAP_TSV     — output: genetic_offset_map.tsv   — COMMA-SEPARATED, one per scenario
# 16  OUT_RASTER_TIF  — output: offset raster (.tif)     — COMMA-SEPARATED, one per scenario
# 17  OUT_IMPORTANCE  — output: predictor importance plot (.png) — SINGLE, scenario-free
# 18  SCENARIOS       — comma-separated scenario names (optional, for log messages)
#
# A single-scenario call passes one path in each list and behaves exactly as it
# did before scenarios existed.
#
# WHAT IS AND IS NOT HOISTED OUT OF THE SCENARIO LOOP.
# The genotype matrix load, the chr:pos index build, the candidate.loci vector,
# the present-climate tables and the raster template are scenario-invariant and
# are computed once. The LFMM2 fit is NOT hoisted: LEA::genetic.gap() takes
# env/new.env/pred.env together and fits and projects in a single call, so there
# is no fitted object to reuse across futures. (Stacking all scenarios into one
# call would fit once, but it also changes what gap$eigenvalues/$vectors are
# computed over, so it is not a drop-in substitution — it would need validating
# against the per-scenario numbers before being adopted.)
##############################################################################
LFMM_IMP_FULL  <- args[1]
VCFSNP         <- args[2]
REMOVED        <- args[3]
CANDIDATE_SNPS <- args[4]
ENV_SITE_PRES  <- args[5]
ENV_SITE_FUT   <- args[6]  %>% str_split(',') %>% unlist
ENV_ALL_PRES   <- args[7]
ENV_ALL_FUT    <- args[8]  %>% str_split(',') %>% unlist
PRES_RASTER    <- args[9]
SAMPLES        <- args[10]
PREDICTORS     <- args[11] %>% str_split(',') %>% unlist
K              <- args[12] %>% as.integer
SCALE          <- as.logical(args[13])
OUT_SITE_TSV   <- args[14] %>% str_split(',') %>% unlist
OUT_MAP_TSV    <- args[15] %>% str_split(',') %>% unlist
OUT_RASTER_TIF <- args[16] %>% str_split(',') %>% unlist
OUT_IMPORTANCE <- args[17]
SCENARIOS      <- if (length(args) >= 18 && nzchar(args[18])) {
    args[18] %>% str_split(',') %>% unlist
} else {
    paste0('scenario_', seq_along(ENV_SITE_FUT))
}

N_SCEN <- length(ENV_SITE_FUT)
stopifnot(length(ENV_ALL_FUT)   == N_SCEN,
          length(OUT_SITE_TSV)  == N_SCEN,
          length(OUT_MAP_TSV)   == N_SCEN,
          length(OUT_RASTER_TIF) == N_SCEN,
          length(SCENARIOS)     == N_SCEN)
##############################################################################

set.seed(42)
message('INFO: Geometric Genetic Offset (LEA::genetic.gap)')
message(paste0('INFO: K=', K, '  scale=', SCALE, '  scenarios=', N_SCEN))

# ── 1. Load full imputed genotype matrix ─────────────────────────────────────
message('INFO: Loading full imputed LFMM matrix')
lfmm_imp <- as.matrix(fread(LFMM_IMP_FULL))
message(paste0('INFO: LFMM matrix: ', nrow(lfmm_imp), ' samples x ', ncol(lfmm_imp), ' SNPs'))
# Keyed on the sentinel 9, not on `max > 2`, so a polyploid dataset (sNMF.ploidy > 2) is not
# falsely rejected. 9 is LFMM's reserved missing code and never a valid dosage in this format.
if (any(lfmm_imp == 9, na.rm = TRUE)) {
    stop("Geometric offset: genotype matrix contains the LFMM missing code 9 — this is an ",
         "UNIMPUTED LFMM file. LFMM2 would fit effect sizes against the missing-data code. ",
         "Requires an imputed matrix (*_imp*). Check the Snakemake input wiring.")
}

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

# ── 6. Load present climate (scenario-invariant) ─────────────────────────────
env_site_pres <- fread(ENV_SITE_PRES) %>% dplyr::select(!!PREDICTORS)
env_all_pres  <- fread(ENV_ALL_PRES)  %>% dplyr::select(ID, !!PREDICTORS)

n_site <- nrow(env_site_pres)
message(paste0('INFO: Site rows: ', n_site,
               '  Raster rows: ', nrow(env_all_pres)))

# Defense-in-depth: env_site_pres/env_site_fut and lfmm_imp are all supposed to already be
# row-aligned to the same climate-valid sample set (subset_lfmm_climate + the climate_site
# tables narrowed by download_climate_present/merge_climate_future). This should never fire --
# if it does, something upstream desynced the sample sets and feeding mismatched rows into
# LEA::genetic.gap() would silently bind the wrong sample's genotypes to the wrong climate
# values. Fail loud instead.
if (nrow(lfmm_imp) != n_site) {
    stop(paste0('FATAL: Sample-count mismatch between climate and genotype inputs -- ',
                'env_site_pres has ', n_site, ' rows, lfmm_imp has ', nrow(lfmm_imp),
                ' rows. These must match (both should already be subset to the same ',
                'climate-valid sample set). Re-run filter_climate_valid_samples/',
                'subset_lfmm_climate and check for a stale intermediate file.'))
}

env      <- as.matrix(env_site_pres)             # training environment = present at sites
clim_pres <- rast(PRES_RASTER)
samples <- fread(SAMPLES,
                 colClasses = c('site' = 'character', 'sample' = 'character',
                                'latitude' = 'numeric',  'longitude' = 'numeric'))
site_ids <- samples[, c('site', 'sample')] %>% unique()

# gap object of the FIRST scenario — the importance plot below is derived from
# the fit, which is scenario-free, so it is written once from this one.
gap_first <- NULL

# ── 7-10. Per-scenario: fit + project, then write site/raster/map ─────────────
for (i in seq_len(N_SCEN)) {
    message(paste0('INFO: [', i, '/', N_SCEN, '] scenario ', SCENARIOS[i]))

    env_site_fut <- fread(ENV_SITE_FUT[i]) %>% dplyr::select(!!PREDICTORS)
    env_all_fut  <- fread(ENV_ALL_FUT[i])  %>% dplyr::select(ID, !!PREDICTORS)

    if (nrow(env_site_fut) != n_site) {
        stop(paste0('FATAL: Sample-count mismatch -- env_site_pres has ', n_site,
                    ' rows but ', ENV_SITE_FUT[i], ' has ', nrow(env_site_fut), '.'))
    }

    # Drop raster cells with any NA across present+future predictors
    ok <- complete.cases(env_all_pres[, -1]) & complete.cases(env_all_fut[, -1])
    env_all_pres_ok <- env_all_pres[ok, ]
    env_all_fut_ok  <- env_all_fut[ok, ]
    message(paste0('INFO: Non-NA raster cells: ', sum(ok)))

    # Stack: sites first (top n_site rows), then landscape cells
    new.env  <- as.matrix(rbind(env_site_pres, env_all_pres_ok[, -1]))
    pred.env <- as.matrix(rbind(env_site_fut,  env_all_fut_ok[,  -1]))

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
    if (is.null(gap_first)) gap_first <- gap

    message(paste0('INFO: genetic.gap complete. Offset length: ', length(gap$offset)))

    # Slice: first n_site = per-site, remainder = landscape
    offset_site      <- gap$offset[seq_len(n_site)]
    offset_landscape <- gap$offset[(n_site + 1):length(gap$offset)]

    message(paste0('INFO: Site offset range: ',
                   round(min(offset_site, na.rm = TRUE), 4), ' – ',
                   round(max(offset_site, na.rm = TRUE), 4)))

    # ── Write site TSV (byte-compatible with GF output) ──────────────────────
    go_site <- site_ids %>%
        cbind(genetic_offset = offset_site) %>%
        dplyr::arrange(desc(genetic_offset))

    fwrite(go_site, OUT_SITE_TSV[i], sep = '\t')
    message(paste0('INFO: Saved site GO values: ', OUT_SITE_TSV[i]))

    # ── Write raster (mirror gradient_forest_offset.R) ───────────────────────
    # Re-derived from clim_pres each iteration: terra SpatRasters wrap external
    # pointers, so reusing a hoisted template risks mutating it.
    rast_offset <- clim_pres[[PREDICTORS[1]]]

    # Assign offset to the valid cell IDs
    rast_offset[] <- NA
    rast_offset[env_all_pres_ok$ID] <- offset_landscape

    writeRaster(rast_offset,
                filename  = OUT_RASTER_TIF[i],
                overwrite = TRUE,
                gdal      = c('INTERLEAVE=BAND', 'COMPRESS=LZW'))
    message(paste0('INFO: Saved offset raster: ', OUT_RASTER_TIF[i]))

    # ── Write map TSV (mirror gradient_forest_offset.R) ──────────────────────
    GO_matrix <- matrix(values(rast_offset), ncol = ncol(rast_offset), byrow = TRUE)
    GO_matrix %>%
        as.data.table() %>%
        fwrite(OUT_MAP_TSV[i], sep = '\t', col.names = FALSE)
    message(paste0('INFO: Saved map GO values: ', OUT_MAP_TSV[i]))
}

gap <- gap_first

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
