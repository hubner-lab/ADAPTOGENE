# Tell data.table that this package uses its syntax (avoids cedta() errors)
.datatable.aware <- TRUE

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
    if (length(ks) > 0) max(ks) else NA_integer_
}

#' Resolve adjust string for a method from ASSOC_CONFIGS entries
#' Returns "bonf_0.05" style string or NULL
#' @noRd
resolve_adjust <- function(config, method, module = "association") {
    # Try module-specific configs first (e.g. phenotype_association.configs)
    configs <- config_get(config, module, "configs", default = list())
    # Fall back to association.configs (phenotype_association inherits from association)
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
