#' Read a project config YAML file
#' @param project project name
#' @param pipeline_path pipeline root path
#' @return named list of config values, or empty list if not found
#' @noRd
read_project_config <- function(project, pipeline_path = get_pipeline_path()) {
    config_file <- file.path(pipeline_path, paste0("config_", project, ".yaml"))
    if (!file.exists(config_file)) {
        config_file <- file.path(pipeline_path, "config.yaml")
    }
    if (!file.exists(config_file)) return(list())
    tryCatch(yaml::read_yaml(config_file), error = function(e) list())
}

#' Safely extract nested config values with a default
#' @param config config list from read_project_config()
#' @param ... key path (character scalars)
#' @param default value to return if key not found
#' @noRd
config_get <- function(config, ..., default = NULL) {
    keys <- list(...)
    val <- config
    for (key in keys) {
        if (is.null(val) || !key %in% names(val)) return(default)
        val <- val[[key]]
    }
    val %||% default
}

#' Get K_best from config
#' @noRd
config_k_best <- function(config) {
    config_get(config, "sNMF", "k_best", default = NA_integer_)
}

#' Get association configs list from config
#' @noRd
config_assoc_configs <- function(config) {
    config_get(config, "GEA", "configs", default = list())
}

#' Get predictors (climate variables) from config
#' @noRd
config_predictors <- function(config) {
    preds <- config_get(config, "Climate", "predictors", default = character(0))
    if (is.character(preds)) preds else character(0)
}

#' Get top_regions setting from config
#' @noRd
config_top_regions <- function(config, module = "GEA") {
    config_get(config, module, "top_regions", default = 10L)
}

#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0) x else y
