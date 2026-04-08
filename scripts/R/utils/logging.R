# logging.R — lightweight logging helpers
# Usage: source("/pipeline/scripts/R/utils/logging.R")

log_info  <- function(...) message("INFO: ", ...)
log_warn  <- function(...) message("WARNING: ", ...)
log_error <- function(...) message("ERROR: ", ...)
