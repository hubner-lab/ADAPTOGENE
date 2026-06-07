# Tell data.table that this package uses its syntax (avoids cedta() errors)
.datatable.aware <- TRUE

#' Keep only canonical merged per-method pvalue/wza files from a glob result.
#'
#' GWAS DROP-mode (Path B) writes per-trait files ({trait}_pvalues_K{k}.tsv)
#' directly into the same methods/{METHOD}/ directory as the merged
#' {METHOD}_pvalues_K{k}.tsv. The canonical merged file always starts with the
#' directory (method) name followed by the regime token ("_pvalues_K" or "_wza_K").
#' Per-trait files never start with the method name, so this filter excludes them.
#'
#' @param files Character vector of glob results (already filtered for sig_snps/wza/block).
#' @param regime "snp" (default) or "wza".
#' @noRd
keep_merged_method_files <- function(files, regime = "snp") {
    if (length(files) == 0) return(files)
    token <- if (regime == "wza") "_wza_K" else "_pvalues_K"
    keep <- vapply(files, function(f) {
        parts <- strsplit(f, .Platform$file.sep, fixed = TRUE)[[1]]
        mi <- which(parts == "methods")
        if (length(mi) == 0) return(FALSE)
        startsWith(basename(f), paste0(parts[mi + 1L], token))
    }, logical(1))
    files[keep]
}

#' Check if a file path exists and is non-empty
#' @noRd
file_ok <- function(path) {
    !is.null(path) && length(path) == 1 && !is.na(path) &&
        nzchar(path) && file.exists(path) && file.size(path) > 0
}

#' Safely read a JSON file, returning NULL on error
#' @noRd
safe_read_json <- function(path) {
    if (!file_ok(path)) return(NULL)
    tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) NULL)
}

#' Resolve K_best: from config, then fall back to max available K
#' @noRd
resolve_k_best <- function(project, config = NULL) {
    k_cfg <- if (!is.null(config)) config_k_best(config) else NA_integer_
    if (!is.na(k_cfg) && k_cfg > 0L) return(k_cfg)
    ks <- find_k_values(project)
    if (length(ks) == 0) return(NA_integer_)
    # Flag the fallback: config k_best is absent/invalid so we guess max available K.
    # This is the condition that can pin a stale K (e.g. K6 while outputs are K5) — keep
    # the fallback (needed pre-sNMF) but make it visible in logs rather than silent.
    k_fallback <- max(ks)
    warning(sprintf(
        "resolve_k_best: config sNMF.k_best absent/invalid for project '%s'; falling back to max available K = %s",
        project, k_fallback))
    k_fallback
}

#' Resolve adjust string for a method from ASSOC_CONFIGS entries
#' Returns "bonf_0.05" style string or NULL
#' @noRd
resolve_adjust <- function(config, method, module = "GEA") {
    # Try module-specific configs first (e.g. GWAS.configs)
    configs <- config_get(config, module, "configs", default = list())
    # Fall back to GEA.configs (GWAS inherits from GEA)
    if (length(configs) == 0) configs <- config_assoc_configs(config)
    for (cfg in configs) {
        if (identical(cfg$method, method)) {
            return(paste0(cfg$adjust, "_", cfg$threshold))
        }
    }
    NULL
}

#' Create a reactive project data bundle for use across modules
#' @noRd
make_project_data <- function(project, pipeline_path = get_pipeline_path()) {
    config <- read_project_config(project, pipeline_path)
    k_best <- resolve_k_best(project, config)
    list(
        name    = project,
        path    = project_base(project, pipeline_path),
        config  = config,
        k_best  = k_best,
        status  = check_module_status(project)
    )
}

#' Map Shiny tab name to Snakemake mode string
#' Returns NULL for haplotype (handled by two separate runner instances)
#' @noRd
tab_to_mode <- function(tab) {
    switch(tab,
        processing    = "processing",
        prestructure  = "prestructure",
        structure     = "structure",
        gea           = "gea",
        gwas          = "gwas",
        gea_x_gwas    = "gea_x_gwas",
        maladaptation = "maladaptation",
        NULL
    )
}

#' Safe DT rendering with "no data" message fallback
#' @noRd
safe_datatable <- function(dt, ...) {
    if (is.null(dt) || nrow(dt) == 0) {
        dt <- data.frame(Message = "No data available")
    }
    DT::datatable(dt, ...)
}

#' Format p-value columns in a DT with signif formatting
#' @noRd
dt_format_pvals <- function(tbl, cols = c("pvalue", "p_adjust")) {
    cols_present <- intersect(cols, names(tbl))
    if (length(cols_present) == 0) return(tbl)
    DT::formatSignif(tbl, columns = cols_present, digits = 3)
}

#' Resolve snp_clumping_distance config value to a numeric bp for the Shiny filter bar.
#'
#' The pipeline config accepts "auto_per_chromosome", "auto_genome_wide", or an integer.
#' The filter bar numericInput always needs an integer. When auto_* is set, we read
#' the genome-wide r2_02_bp from the LD decay TSV (group "All"). Falls back to 1,000,000.
#' @noRd
resolve_ui_snp_clumping_distance <- function(config, module, project_name,
                                              fallback = 1000000L) {
    raw <- config_get(config, module, "snp_clumping_distance",
                      default = config_get(config, module, "region_distance",
                                           default = "auto_per_chromosome"))
    val <- suppressWarnings(as.integer(raw))
    if (!is.na(val) && val >= 1000L) return(val)
    # auto mode: read LD decay TSV
    tsv <- ld_decay_table_path(project_name)
    if (file_ok(tsv)) {
        dt <- tryCatch(
            data.table::fread(tsv, sep = "\t", header = TRUE),
            error = function(e) NULL
        )
        required <- c("scope", "group", "r2_02_bp")
        if (!is.null(dt) && nrow(dt) > 0 && all(required %in% names(dt))) {
            mask <- dt$scope == "genome_wide" & dt$group == "All"
            mask[is.na(mask)] <- FALSE
            gw <- dt$r2_02_bp[mask]
            if (length(gw) > 0 && !is.na(gw[1]) && gw[1] > 0)
                return(as.integer(gw[1]))
        }
    }
    fallback
}

#' Legacy alias — kept for any callers not yet updated
#' @noRd
resolve_ui_region_distance <- function(config, module, project_name,
                                       fallback = 1000000L) {
    resolve_ui_snp_clumping_distance(config, module, project_name, fallback)
}
