library(dplyr)
library(data.table)
library(gradientForest)
library(stringr)
library(qs)
args = commandArgs(trailingOnly=TRUE)
##############
LFMM = args[1]
SIGSNPS = args[2]       # Selected_SNPs.tsv with SNPID column
VCFSNP = args[3]        # .vcfsnp file
REMOVED = args[4]       # .removed file
SAMPLES = args[5]
CLIM_PRESENT_SITE = args[6]
PREDICTORS_SELECTED = args[7] %>% str_split(',') %>% unlist
NTREE = args[8] %>% as.numeric
COR_THRESHOLD = args[9] %>% as.numeric
SPATIAL_CORRECTION = args[10]         # "with" or "without"
MODEL_TYPE = args[11]   # "adaptive" or "random"
OUTPUT = args[12]
# Appended after OUTPUT (positions 13-14), matching rda.R's "append after
# outputs, never insert mid-list" convention — inserting earlier would
# silently shift MODEL_TYPE/OUTPUT. "NULL" (literal string) when
# SPATIAL_CORRECTION == 'without' (workflow/rules/maladaptation.smk's
# nospatial branch never reads these).
DBMEM_VECTORS = args[13]   # PreGEA/tables/spatial/dbmem_vectors.tsv, or "NULL"
DBMEM_SELECTED = args[14]  # PreGEA/tables/varpart/dbmem_selected.tsv, or "NULL"
##############

set.seed(42)
message(paste0('INFO: Building ', MODEL_TYPE, ' Gradient Forest model'))

# Load LFMM genotype matrix
lfmm_dt <- fread(LFMM)
message(paste0('INFO: LFMM matrix: ', nrow(lfmm_dt), ' samples x ', ncol(lfmm_dt), ' SNPs'))

# Load SNP IDs and removed SNPs
sigsnps <- fread(SIGSNPS, header = TRUE)$SNPID
vcfsnp <- fread(VCFSNP, header = FALSE) %>%
  dplyr::select(V1, V2) %>%
  setNames(c('chr', 'pos')) %>%
  dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>%
  .$chrpos

removed_dt <- fread(REMOVED, header = FALSE)
if (!is.null(removed_dt) && nrow(removed_dt) > 0) {
  removed <- removed_dt %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>%
    .$chrpos
  mask_adaptive <- (!vcfsnp %in% removed) & (vcfsnp %in% sigsnps)
} else {
  mask_adaptive <- vcfsnp %in% sigsnps
}

n_adaptive <- sum(mask_adaptive)
message(paste0('INFO: Adaptive SNPs: ', n_adaptive))

# Select SNPs based on model type
if (MODEL_TYPE == 'adaptive') {
  snp_subset <- lfmm_dt[, ..mask_adaptive]
  message(paste0('INFO: Using ', ncol(snp_subset), ' adaptive SNPs'))
} else if (MODEL_TYPE == 'random') {
  # [CHANGED 2026-08-10] Size-matched to the adaptive set. Was
  # min(max(n_adaptive, 300), ncol(lfmm_dt)) -- a 300-SNP floor, so a 185-SNP curated set was
  # "controlled" by a 300-SNP random model, 62% larger than the thing it controls for. The
  # neutral model only means anything as a like-for-like comparison, and nothing recorded the
  # size difference downstream. See docs/pipeline_improvement_requests.md request M.
  n_random <- min(n_adaptive, ncol(lfmm_dt))
  if (n_random < 50) {
    message(paste0('WARNING: random model has only ', n_random, ' SNPs (matched to the adaptive ',
                   'set). Gradient Forest turnover functions are unstable at this size; the ',
                   'previous 300-SNP floor existed to avoid this, and was removed because it ',
                   'broke the size match. Treat this neutral model as indicative only.'))
  }
  # [CHANGED 2026-08-10] Seed derived from the run label. Previously the global set.seed(42) at
  # the top of this script was the only seed, and the sampling pool is the full matrix regardless
  # of {run_label} -- so EVERY random_model.qs in a project was the identical draw. Verified on
  # Trifolium86FinalPC3: the three run_labels' response column sets were identical(), 300/300
  # overlap. The pipeline reported three neutral models and had one. The spatial tag is stripped
  # so that {label}_nospatial and {label}_spatial share a draw, keeping the spatial-correction
  # comparison like-for-like.
  run_label <- sub('_(no)?spatial$', '', basename(dirname(OUTPUT)))
  set.seed(42L + sum(utf8ToInt(run_label)))
  random_cols <- sample(names(lfmm_dt), n_random)
  snp_subset <- lfmm_dt[, ..random_cols]
  message(paste0('INFO: Using ', ncol(snp_subset), ' random SNPs (size-matched to the ',
                 n_adaptive, '-SNP adaptive set; run_label "', run_label, '", seed ',
                 42L + sum(utf8ToInt(run_label)), ')'))
} else {
  stop(paste0('ERROR: Unknown MODEL_TYPE: ', MODEL_TYPE))
}

# Guard: these columns become gradientForest() RESPONSE variables, and gradientForest() has no
# na.action. An unimputed LFMM file encodes a missing call as the literal 9, which would enter the
# regression forest as a genotype 3.5x further from the reference homozygote than the real
# alternate homozygote is -- so the forest learns which samples were CALLED, not which carry the
# allele. Checked on the subset actually passed to the model. Same guard as scripts/rda_offset.R:106.
# Keyed on the sentinel 9 rather than on `max > 2`, so a legitimately polyploid dataset
# (sNMF.ploidy > 2, dosages 3,4,...) is not falsely rejected. 9 is LFMM's reserved missing code,
# so it is never a valid dosage in this format regardless of ploidy.
subset_mat <- as.matrix(snp_subset)
if (any(subset_mat == 9, na.rm = TRUE)) {
  frac9 <- mean(subset_mat == 9, na.rm = TRUE)
  stop(sprintf(paste0('Gradient Forest: genotype matrix contains the LFMM missing code 9 in ',
                      '%.1f%% of cells - this is an UNIMPUTED LFMM file. gradientForest() has no ',
                      'na.action and would read 9 as a genotype, so the forest learns which ',
                      'samples were CALLED rather than which carry the allele. Requires an ',
                      'imputed matrix (*_imp*). Check the Snakemake input wiring.'), 100 * frac9))
}
rm(subset_mat)

# Compute maxLevel for gradient forest
maxLevel <- log2(0.368 * nrow(lfmm_dt) / 2)

# Load samples (also the sample-ID key for the dbMEM join below)
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# Load climate predictors
env.bio <- fread(CLIM_PRESENT_SITE) %>% dplyr::select(!!PREDICTORS_SELECTED)

# Build input matrix and predictor list. Spatial covariates are the climate
# module's forward-selected dbMEM vectors (adespatial::dbmem() + ordiR2step,
# climate.smk Block 2 — scripts/pregea_dbmem.R + scripts/pregea_varpart.R,
# moved out of PreGEA in the module split), NOT the old ad hoc
# vegan::pcnm(dist(raw lon/lat)) + "keep half the positive eigenvalues"
# heuristic this replaces. "Keep all positive MEMs" was considered and
# rejected: docs/rda_research.md's A20/C4 rows always pair dbMEM WITH forward
# selection, never recommend the raw set alone as a GF predictor block.
if (SPATIAL_CORRECTION == 'with') {
  dbmem_dt <- fread(DBMEM_VECTORS, colClasses = c(sample = "character"))
  sel_dt   <- fread(DBMEM_SELECTED)
  selected_mems <- sel_dt[selected == TRUE, mem]

  mem_cols <- grep("^MEM\\d+$", names(dbmem_dt), value = TRUE)
  if (length(mem_cols) == 0 || length(selected_mems) == 0) {
    stop("HARD STOP: dbMEM produced no usable spatial vectors for this project ",
         "(see PreGEA/tables/spatial/dbmem_diagnostics.tsv and ",
         "PreGEA/tables/varpart/dbmem_selected.tsv for why — e.g. too few ",
         "unique sample coordinates, or ordiR2step selected 0 MEMs). A ",
         "'spatial' Gradient Forest model cannot be built on this dataset. ",
         "Set Maladaptation.methods.gradient_forest.spatial_correction to ",
         "'without', or investigate why dbMEM/forward-selection degenerated.")
  }

  # dbmem_vectors.tsv is keyed by sample (preGEA's own row order); this
  # script's LFMM/env.bio matrices are positionally aligned to `samples`
  # (SAMPLES arg) — match() explicitly here rather than assuming row order,
  # since dbMEM comes from a different source file than the other two.
  dbmem_dt <- dbmem_dt[match(samples$sample, sample)]
  stopifnot(!anyNA(dbmem_dt$sample))  # every GF sample must resolve to a dbMEM row

  mem_scores <- as.data.frame(dbmem_dt[, ..selected_mems])
  colnames(mem_scores) <- selected_mems
  message(paste0('INFO: Using ', ncol(mem_scores), ' forward-selected dbMEM axes: ',
                 paste(selected_mems, collapse = ', ')))

  input_matrix <- cbind(env.bio, mem_scores, snp_subset)
  predictor_vars <- c(colnames(env.bio), colnames(mem_scores))
} else {
  input_matrix <- cbind(env.bio, snp_subset)
  predictor_vars <- colnames(env.bio)
}

message(paste0('INFO: Input matrix: ', nrow(input_matrix), ' samples x ', ncol(input_matrix), ' variables'))
message(paste0('INFO: Predictor variables: ', paste(predictor_vars, collapse = ', ')))
message(paste0('INFO: Response variables (SNPs): ', ncol(snp_subset)))

# Run Gradient Forest
gf <- gradientForest(input_matrix,
                     predictor.vars = predictor_vars,
                     response.vars = colnames(snp_subset),
                     ntree = NTREE,
                     maxLevel = maxLevel,
                     trace = TRUE,
                     corr.threshold = COR_THRESHOLD)

# Save model
qsave(gf, OUTPUT)
message(paste0('INFO: Saved ', MODEL_TYPE, ' GF model to: ', OUTPUT))
