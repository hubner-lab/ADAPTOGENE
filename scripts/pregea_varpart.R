#!/usr/bin/env Rscript
# pregea_varpart.R — shared climate/structure/geography variance-partitioning
# panel (preGEA Block 4b). Ports rda_varpart_lasky.Rmd's forward-selection,
# Px, varpart, and anova.cca blocks (see that file for the reference
# implementation this was verified against).
#
# CORRECTIONS applied on port (docs/rda_research.md C.0/C5):
#   - Px: the Rmd's CODE already divides by tot.chi (sigma^2, total SNP
#     variance) correctly — only a COMMENT in the Rmd was wrong ("/ r2").
#     Fixed here from the start; see compute_Px() below.
#   - ordiR2step gains Pin/R2permutations/R2scope explicitly (the Rmd omitted
#     the first two) — Blanchet's double stopping rule (A17/C4).
#   - A17: the FULL forward-selection path is written (cumulative adjusted R2
#     per step + the full-model ceiling), never just the final variable set.
#
# RESPONSE MATRIX defaults to genomic PCs (LEA .projections truncated at a
# cumulative-variance target) — a documented SCALABILITY ADAPTATION, not
# literal Lasky et al. 2012 methodology (which used raw SNPs). The
# partitioning FRAMEWORK (varpart/dbMEM/Px) is what is being ported, not the
# exact response matrix — see docs/rda_research.md C.0.
#
# STRUCTURE TABLE is the sNMF Q-matrix (C1..CK, LAST column dropped for
# compositional singularity) — explicitly NOT the LEA PCs used as the
# response, which would make X_struct = Y tautological.

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)      # load_pca_covariates() (emmax_core.R) uses %>% internally
    library(vegan)
    library(ggplot2)
    library(ggrepel)
})

source("/pipeline/scripts/R/utils/theme_adaptogene.R")
source("/pipeline/scripts/R/utils/emmax_core.R")   # load_pca_covariates()

args <- commandArgs(trailingOnly = TRUE)
################################################################################
PROJECTIONS      <- args[1]
EIGENVALUES      <- args[2]
DBMEM_VECTORS    <- args[3]
CLUSTERS         <- args[4]                  # "NULL" if structure_table=none
CLIMATE          <- args[5]
CLIMATE_VALID    <- args[6]
SAMPLES_ORDER    <- args[7]
PREDICTORS       <- args[8]
RESPONSE         <- args[9]                  # "pcs" | "snps"
RESPONSE_VAR     <- as.numeric(args[10])
RESPONSE_MAX_PCS <- as.integer(args[11])
RESPONSE_MIN_PCS <- as.integer(args[12])
LFMM_PRUNED      <- args[13]                 # "NULL" when RESPONSE="pcs"
ORDIR2STEP_PIN   <- as.numeric(args[14])
R2_PERMUTATIONS  <- as.numeric(args[15])     # ordiR2step's internal R2permutations (fixed constant, Climate.Varpart layer)
VARPART_PERMS    <- as.numeric(args[16])     # the one surviving Varpart permutations knob (varpart anova test)
CONFOUND_FLAG    <- toupper(args[17]) == "TRUE"
SEED             <- as.numeric(args[18])
OUT_SELECTION_TSV<- args[19]
OUT_SELECTED_TSV <- args[20]
OUT_FRACTIONS_TSV<- args[21]
OUT_ANOVA_TSV    <- args[22]
OUT_PX_TSV       <- args[23]
PLOT_DIR         <- args[24]
INTER_DIR        <- args[25]
################################################################################
# NOTE: SEL_PERMUTATIONS (old arg 16, "ordir2step_permutations" config key) was
# dropped here — it was assigned but never referenced (ordiR2step's own
# R2permutations control comes from R2_PERMUTATIONS above, arg 15). Removed
# together with the config key (module split cleanup) rather than left dead.

for (d in c(dirname(OUT_SELECTION_TSV), PLOT_DIR, INTER_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

predictor_names <- strsplit(PREDICTORS, ",")[[1]] |> trimws()
climate_valid_ids <- fread(CLIMATE_VALID, header = FALSE, colClasses = "character")$V2
sample_order_all <- readLines(SAMPLES_ORDER)

################################################################################
# 1. Response matrix Y — genomic PCs (default) or raw pruned SNPs
################################################################################
if (RESPONSE == "snps") {
    Y <- fread(LFMM_PRUNED, sep = " ", header = FALSE) |> as.matrix()
    response_desc <- sprintf("raw pruned SNPs (m=%d)", ncol(Y))
} else {
    eig <- as.numeric(readLines(EIGENVALUES))
    cum_var <- cumsum(eig / sum(eig))
    k_target <- which(cum_var >= RESPONSE_VAR)[1]
    if (is.na(k_target)) k_target <- length(eig)
    k_use <- max(RESPONSE_MIN_PCS, min(k_target, RESPONSE_MAX_PCS, length(eig)))
    Y <- load_pca_covariates(PROJECTIONS, k_use, SAMPLES_ORDER, climate_valid_ids)
    if (is.null(Y)) stop("preGEA varpart: could not build genomic-PC response matrix (k_use=", k_use, ")")
    response_desc <- sprintf("genomic PCs (k=%d, %.1f%% cumulative variance target)", k_use, RESPONSE_VAR * 100)
}
message("INFO: Response matrix Y: ", nrow(Y), " x ", ncol(Y), " [", response_desc, "]")
n_samples <- nrow(Y)

################################################################################
# 2. Climate predictors X_clim (climate-valid, same row order as Y)
################################################################################
climate_dt <- fread(CLIMATE, sep = "\t", header = TRUE)
X_clim <- as.data.frame(climate_dt[, ..predictor_names])
if (nrow(X_clim) != n_samples) stop("preGEA varpart: climate/response row mismatch (", nrow(X_clim), " vs ", n_samples, ")")
rownames(X_clim) <- climate_dt$sample %||% seq_len(n_samples)

################################################################################
# 3. Structure X_struct — sNMF Q-matrix, last column dropped (compositional)
################################################################################
X_struct <- NULL
if (!identical(CLUSTERS, "NULL") && file.exists(CLUSTERS)) {
    clusters_dt <- fread(CLUSTERS, colClasses = c(sample = "character"))
    q_cols <- grep("^C\\d+$", names(clusters_dt), value = TRUE)
    if (length(q_cols) >= 2) {
        q_cols_use <- q_cols[-length(q_cols)]   # drop last column (compositional singularity)
        clusters_dt <- clusters_dt[match(climate_valid_ids, sample)]
        X_struct <- as.data.frame(clusters_dt[, ..q_cols_use])
        message("INFO: Structure table: ", length(q_cols_use), "/", length(q_cols), " Q-matrix columns (last dropped)")
    }
}

################################################################################
# 4. Geography X_geo — dbMEM vectors, forward-selected (A17: full path kept)
################################################################################
dbmem_dt <- fread(DBMEM_VECTORS, colClasses = c(sample = "character"))
mem_cols <- grep("^MEM\\d+$", names(dbmem_dt), value = TRUE)
dbmem_status <- if (length(mem_cols) == 0) "no_spatial_vectors" else "ok"

X_geo_sel <- NULL
selection_dt <- data.table(step = integer(), variable = character(),
                           r2_adj_cumulative = numeric(), r2_adj_full_ceiling = numeric())
selected_dt <- data.table(mem = character(), selected = logical())

if (dbmem_status == "ok") {
    dbmem_dt <- dbmem_dt[match(climate_valid_ids, sample)]
    X_geo_full <- as.data.frame(dbmem_dt[, ..mem_cols])
    full_fit <- tryCatch(rda(as.formula(paste("Y ~", paste(mem_cols, collapse = " + "))), data = X_geo_full),
                         error = function(e) NULL)
    full_ceiling <- if (!is.null(full_fit)) suppressWarnings(RsquareAdj(full_fit)$adj.r.squared) else NA_real_

    selected_names <- character(0)
    if (!is.null(full_fit) && !is.na(full_ceiling) && full_ceiling > 0) {
        null_fit <- rda(Y ~ 1, data = X_geo_full)
        set.seed(SEED)
        step_res <- tryCatch(
            ordiR2step(null_fit, scope = formula(full_fit), Pin = ORDIR2STEP_PIN,
                      R2scope = TRUE, R2permutations = R2_PERMUTATIONS, trace = FALSE),
            error = function(e) { message("WARNING: ordiR2step (dbMEM) failed: ", e$message); NULL })
        if (!is.null(step_res)) selected_names <- attr(terms(step_res), "term.labels")
    }
    if (length(selected_names) == 0) {
        message("WARNING: forward selection chose 0 MEMs; falling back to all ", length(mem_cols), " MEM(s).")
        selected_names <- mem_cols
    }

    cum_terms <- character(0)
    for (i in seq_along(selected_names)) {
        cum_terms <- c(cum_terms, selected_names[i])
        step_fit <- tryCatch(rda(as.formula(paste("Y ~", paste(cum_terms, collapse = " + "))),
                                 data = X_geo_full), error = function(e) NULL)
        if (is.null(step_fit)) next
        selection_dt <- rbind(selection_dt, data.table(
            step = i, variable = selected_names[i],
            r2_adj_cumulative = suppressWarnings(RsquareAdj(step_fit)$adj.r.squared),
            r2_adj_full_ceiling = full_ceiling))
    }
    selected_dt <- data.table(mem = mem_cols, selected = mem_cols %in% selected_names)
    X_geo_sel <- X_geo_full[, selected_names, drop = FALSE]
    message("INFO: dbMEM forward selection: ", length(selected_names), "/", length(mem_cols), " MEM(s) selected")
} else {
    message("INFO: dbMEM status=", dbmem_status, " — geography fraction unavailable, running climate-vs-structure only")
}
fwrite(selection_dt, OUT_SELECTION_TSV, sep = "\t", quote = FALSE)
fwrite(selected_dt,  OUT_SELECTED_TSV,  sep = "\t", quote = FALSE)

################################################################################
# 5. Lasky Px per climate variable — compute_Px(), verbatim port (denominator
#    is tot.chi / (tot.chi - pCCA$tot.chi for partial models), NOT r2 — the
#    Rmd's comment was wrong, its code was always correct).
################################################################################
compute_Px <- function(rda_obj, X) {
    eig <- rda_obj$CCA$eig
    if (is.null(eig) || length(eig) == 0) return(NULL)
    denom <- if (!is.null(rda_obj$pCCA)) rda_obj$tot.chi - rda_obj$pCCA$tot.chi else rda_obj$tot.chi
    if (is.null(denom) || denom <= 0) return(NULL)
    lc <- scores(rda_obj, display = "lc", choices = seq_along(eig), scaling = 0)
    if (is.vector(lc)) lc <- matrix(lc, ncol = 1)
    Xs <- scale(as.matrix(X))
    R  <- suppressWarnings(cor(Xs, lc))
    Px <- as.numeric(abs(R) %*% eig) / denom
    data.table(variable = colnames(Xs), Px = Px, Px_pct = 100 * Px, rank = rank(-Px, ties.method = "min"))
}

rda_B <- tryCatch(rda(as.formula(paste("Y ~", paste(predictor_names, collapse = " + "))), data = X_clim),
                  error = function(e) NULL)
px_rows <- list()
if (!is.null(rda_B)) {
    px_b <- compute_Px(rda_B, X_clim)
    if (!is.null(px_b)) px_rows[[length(px_rows) + 1]] <- px_b[, model := "unconditional"]
}
if (!is.null(X_geo_sel) && ncol(X_geo_sel) > 0) {
    combined_df <- cbind(X_clim, X_geo_sel)
    mem_names <- paste(colnames(X_geo_sel), collapse = " + ")
    rda_C <- tryCatch(
        rda(as.formula(paste("Y ~", paste(predictor_names, collapse = " + "),
                             "+ Condition(", mem_names, ")")), data = combined_df),
        error = function(e) { message("WARNING: partial RDA (Px model C) failed: ", e$message); NULL })
    if (!is.null(rda_C)) {
        px_c <- compute_Px(rda_C, X_clim)
        if (!is.null(px_c)) px_rows[[length(px_rows) + 1]] <- px_c[, model := "partial_geo"]
    }
}
px_dt <- if (length(px_rows) > 0) rbindlist(px_rows) else data.table(
    variable = character(), Px = numeric(), Px_pct = numeric(), rank = integer(), model = character())
fwrite(px_dt, OUT_PX_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote Px table (", nrow(px_dt), " rows)")

################################################################################
# 6. Variance partitioning — 3-table (when struct+geo both available) and the
#    2-table climate-vs-geo Venn (the ONLY place the confounding flag is
#    computed — with 3 tables "shared" is ambiguous, d+e+f+g).
################################################################################
frac_rows <- list()
anova_rows <- list()

# vegan::varpart()'s indfract table names its adj-R2 column "Adj.R.squared" for a
# 2-table Venn but "Adj.R.square" (no trailing "d") for 3+ tables — confirmed on
# this Docker's vegan 2.6-8. data.frame[idx, "wrong_name"] returns NULL silently
# (no error), so hardcoding either spelling makes the OTHER table count go to NA
# with zero warning. Resolve the real column name at each call site instead.
adjr2_col <- function(df) {
    hit <- grep("^Adj\\.R\\.square", colnames(df), value = TRUE)
    if (!length(hit)) stop("varpart indfract table has no Adj.R.square(d) column — vegan API changed?")
    hit[1]
}

add_testable_fraction <- function(label, rda_obj, testable) {
    p <- NA_real_
    if (testable && !is.null(rda_obj)) {
        set.seed(SEED)
        av <- tryCatch(anova.cca(rda_obj, permutations = VARPART_PERMS), error = function(e) NULL)
        if (!is.null(av)) p <- as.data.frame(av)[["Pr(>F)"]][1]
    }
    anova_rows[[length(anova_rows) + 1]] <<- data.table(fraction = label, testable = testable, p_value = p)
}

have_struct <- !is.null(X_struct) && ncol(X_struct) > 0
have_geo    <- !is.null(X_geo_sel) && ncol(X_geo_sel) > 0

# 2-table: climate vs geo — the confounding-flag partition (A.2 point 4)
if (have_geo) {
    vp2 <- tryCatch(varpart(Y, X_clim, X_geo_sel), error = function(e) { message("WARNING: 2-table varpart failed: ", e$message); NULL })
    if (!is.null(vp2)) {
        fr <- vp2$part$indfract
        get_fr <- function(pat) { idx <- grep(pat, rownames(fr), fixed = TRUE); if (length(idx)) fr[idx, adjr2_col(fr)] else NA_real_ }
        clim_unique <- get_fr("[a]"); shared_2 <- get_fr("[b]"); geo_unique <- get_fr("[c]")
        frac_rows[[length(frac_rows) + 1]] <- data.table(
            partition = "2-table", fraction = c("climate_unique", "shared", "geo_unique"),
            adj_r2 = c(clim_unique, shared_2, geo_unique), fraction_type = c("unique", "shared", "unique"))
        rda_clim_u <- tryCatch(rda(Y, X_clim, X_geo_sel), error = function(e) NULL)
        rda_geo_u  <- tryCatch(rda(Y, X_geo_sel, X_clim), error = function(e) NULL)
        add_testable_fraction("2table_climate_unique", rda_clim_u, TRUE)
        add_testable_fraction("2table_geo_unique", rda_geo_u, TRUE)

        if (CONFOUND_FLAG && is.finite(shared_2) && is.finite(clim_unique) && is.finite(geo_unique)) {
            max_unique <- max(clim_unique, geo_unique, na.rm = TRUE)
            confounded <- shared_2 > max_unique
            frac_rows[[length(frac_rows) + 1]] <- data.table(
                partition = "2-table", fraction = "confounding_flag",
                adj_r2 = as.numeric(confounded), fraction_type = "diagnostic")
            if (confounded) message("WARNING: climate-geography CONFOUNDED — shared (", round(shared_2, 4),
                                    ") exceeds the largest unique fraction (", round(max_unique, 4), ")")
        }
    }
}

# 3-table: climate / structure / geography, when both auxiliary tables exist
status <- if (!have_geo) "no_spatial_vectors" else if (!have_struct) "no_structure_table" else "ok"
if (have_struct && have_geo) {
    vp3 <- tryCatch(varpart(Y, X_clim, X_struct, X_geo_sel), error = function(e) { message("WARNING: 3-table varpart failed: ", e$message); NULL })
    if (!is.null(vp3)) {
        fr3 <- vp3$part$indfract
        letters8 <- c("a", "b", "c", "d", "e", "f", "g")
        for (lt in letters8) {
            idx <- grep(paste0("\\[", lt, "\\]"), rownames(fr3), fixed = FALSE)
            if (length(idx)) frac_rows[[length(frac_rows) + 1]] <- data.table(
                partition = "3-table", fraction = paste0("[", lt, "]"),
                adj_r2 = fr3[idx, adjr2_col(fr3)],
                fraction_type = if (lt %in% c("a", "b", "c")) "unique" else "shared")
        }
        rda_clim3 <- tryCatch(rda(Y, X_clim, cbind(X_struct, X_geo_sel)), error = function(e) NULL)
        rda_str3  <- tryCatch(rda(Y, X_struct, cbind(X_clim, X_geo_sel)), error = function(e) NULL)
        rda_geo3  <- tryCatch(rda(Y, X_geo_sel, cbind(X_clim, X_struct)), error = function(e) NULL)
        add_testable_fraction("3table_climate_unique_a", rda_clim3, TRUE)
        add_testable_fraction("3table_structure_unique_b", rda_str3, TRUE)
        add_testable_fraction("3table_geo_unique_c", rda_geo3, TRUE)
    }
} else if (have_struct && !have_geo) {
    # Degenerate-input fallback (dbmem status != ok): climate-vs-structure 2-table instead.
    vp2s <- tryCatch(varpart(Y, X_clim, X_struct), error = function(e) NULL)
    if (!is.null(vp2s)) {
        frs <- vp2s$part$indfract
        get_fr <- function(pat) { idx <- grep(pat, rownames(frs), fixed = TRUE); if (length(idx)) frs[idx, adjr2_col(frs)] else NA_real_ }
        frac_rows[[length(frac_rows) + 1]] <- data.table(
            partition = "2-table-fallback", fraction = c("climate_unique", "shared", "structure_unique"),
            adj_r2 = c(get_fr("[a]"), get_fr("[b]"), get_fr("[c]")), fraction_type = c("unique", "shared", "unique"))
        rda_clim_u <- tryCatch(rda(Y, X_clim, X_struct), error = function(e) NULL)
        rda_str_u  <- tryCatch(rda(Y, X_struct, X_clim), error = function(e) NULL)
        add_testable_fraction("2table_fallback_climate_unique", rda_clim_u, TRUE)
        add_testable_fraction("2table_fallback_structure_unique", rda_str_u, TRUE)
    }
}

fractions_dt <- if (length(frac_rows) > 0) rbindlist(frac_rows, fill = TRUE) else data.table(
    partition = character(), fraction = character(), adj_r2 = numeric(), fraction_type = character())
fractions_dt[, status := status]
anova_dt <- if (length(anova_rows) > 0) rbindlist(anova_rows) else data.table(fraction = character(), testable = logical(), p_value = numeric())

fwrite(fractions_dt, OUT_FRACTIONS_TSV, sep = "\t", quote = FALSE)
fwrite(anova_dt,      OUT_ANOVA_TSV,      sep = "\t", quote = FALSE)
message("INFO: Wrote varpart fractions (", nrow(fractions_dt), " rows, status=", status, ") and anova (",
       nrow(anova_dt), " rows)")

################################################################################
# 7. Plots — all low-N (variables/fractions/steps), plain geoms are correct
#    (Rule 6 does not apply). No on-plot commentary (Rule 8) — the
#    "shared can be negative / untestable" caveat lives in the Shiny note.
################################################################################
# save_both()/empty_plot() now thin wrappers over the shared
# adapt_save_both()/adapt_empty_plot() helpers (theme_adaptogene.R) — hoisted
# there so pregea_rda_setup.R's duplicated ggsave() pairs can reuse them too.
save_both <- function(name, g, w = 7, h = 5) adapt_save_both(file.path(PLOT_DIR, name), g, w, h)
empty_plot <- adapt_empty_plot

if (nrow(selection_dt) > 0) {
    g_path <- ggplot(selection_dt, aes(x = step, y = r2_adj_cumulative)) +
        geom_hline(aes(yintercept = r2_adj_full_ceiling), color = ADAPT_THRESHOLD, linetype = "dashed") +
        geom_line(color = ADAPT_NEUTRAL) + geom_point(color = ADAPT_NEUTRAL, size = 2) +
        ggrepel::geom_text_repel(aes(label = variable), size = 3, color = ADAPT_COL$fg) +
        labs(x = "Forward-selection step", y = "Cumulative adjusted R2",
            title = "dbMEM forward-selection path (A17)") +
        theme_adaptogene()
    save_both("dbmem_selection_path", g_path)
} else save_both("dbmem_selection_path", empty_plot("No dbMEM selection path (geography unavailable)"))

venn_rows <- fractions_dt[partition == "2-table" & fraction_type != "diagnostic"]
if (nrow(venn_rows) > 0) {
    g_venn <- ggplot(venn_rows, aes(x = fraction, y = adj_r2, fill = fraction_type)) +
        geom_col() +
        scale_fill_manual(values = c(unique = ADAPT_RETAINED, shared = ADAPT_THRESHOLD), guide = "none") +
        labs(x = NULL, y = "Adjusted R2", title = "Climate vs geography variance partition (2-table)") +
        theme_adaptogene()
    save_both("varpart_venn", g_venn)
} else save_both("varpart_venn", empty_plot("2-table climate/geography partition unavailable"))

if (nrow(fractions_dt) > 0 && any(fractions_dt$partition == "3-table")) {
    g_frac <- ggplot(fractions_dt[partition == "3-table"], aes(x = fraction, y = adj_r2, fill = fraction_type)) +
        geom_col() +
        scale_fill_manual(values = c(unique = ADAPT_RETAINED, shared = ADAPT_THRESHOLD), guide = "none") +
        labs(x = NULL, y = "Adjusted R2", title = "Climate / structure / geography 3-table partition") +
        theme_adaptogene()
    save_both("varpart_fractions_bar", g_frac)
} else save_both("varpart_fractions_bar", empty_plot("3-table partition unavailable (needs structure + geography)"))

if (nrow(px_dt) > 0) {
    g_px <- ggplot(px_dt, aes(x = reorder(variable, Px), y = Px_pct, fill = model)) +
        geom_col(position = "dodge") + coord_flip() +
        scale_fill_adaptogene(name = NULL) +
        labs(x = NULL, y = "Px (% of available variance)", title = "Lasky Px per climate variable") +
        theme_adaptogene()
    save_both("px_barplot", g_px)
} else save_both("px_barplot", empty_plot("No Px values computed"))

message("INFO: preGEA varpart complete.")
