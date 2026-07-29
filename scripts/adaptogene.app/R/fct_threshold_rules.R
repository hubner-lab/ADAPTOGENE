#' Per-method significance rule helpers.
#'
#' A "rule" is list(type=, value=). The MASTER rule is the threshold bar's
#' scalar (type, value). An OVERRIDE map is METHOD -> rule for the subset of
#' methods that deviate; a method absent from the map follows master. An empty
#' map is byte-identical to the pre-rework single-global behaviour.
#' @noRd

#' Deterministic, order-independent cache key for an override map.
#' MUST be sorted: an unsorted list key produces cache misses that look like a
#' perf regression (bindCache in mod_gea.R, digest() in fct_data_loading.R).
#' @noRd
threshold_overrides_key <- function(overrides) {
    if (length(overrides) == 0) return("")
    ov <- overrides[order(names(overrides))]
    paste(vapply(names(ov), function(m) sprintf(
        "%s=%s:%s", m, ov[[m]]$type,
        format(as.numeric(ov[[m]]$value), digits = 15, scientific = TRUE)),
        character(1)), collapse = "|")
}

#' Resolve the effective rule for one method.
#' @noRd
effective_rule_for <- function(method, overrides, master_type, master_value) {
    r <- overrides[[method]]
    if (is.null(r)) return(list(type = master_type, value = as.numeric(master_value),
                                source = "master"))
    list(type = r$type %||% master_type,
         value = as.numeric(r$value %||% master_value),
         source = "override")
}

#' Normalise a map read back from region_params.json.
#' read_json(simplifyVector = FALSE) yields nested lists; drop malformed entries
#' and coerce value to numeric. Returns list() for absent/legacy files.
#' @noRd
normalize_threshold_overrides <- function(raw) {
    if (is.null(raw) || length(raw) == 0 || is.null(names(raw))) return(list())
    out <- list()
    for (m in names(raw)) {
        e  <- raw[[m]]
        ty <- tryCatch(as.character(e$type), error = function(...) NA_character_)
        if (length(ty) != 1) ty <- NA_character_
        # e$value on an entry missing "value" returns NULL (no error) — coerce
        # NULL/absent to a scalar NA_real_ before is.na(), or as.numeric(NULL)
        # (a zero-length vector) breaks the length-1 `if (... || is.na(va))`
        # condition below instead of taking the "skip this entry" branch.
        va_raw <- tryCatch(e$value, error = function(...) NULL)
        va <- if (is.null(va_raw)) NA_real_ else suppressWarnings(as.numeric(va_raw))
        if (length(va) != 1) va <- NA_real_
        if (is.na(ty) || !nzchar(ty) || is.na(va)) next
        if (!threshold_value_valid_for_type(ty, va)) next   # fct_combine.R:227
        out[[m]] <- list(type = ty, value = va)
    }
    out
}
