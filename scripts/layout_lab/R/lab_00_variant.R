# lab_00_variant.R — variant selection + theme override.
# Sourced first (files load alphabetically), so later lab files can call lab_variant().

.lab_variant_store <- new.env(parent = emptyenv())
.lab_variant_store$v <- "orig"
.lab_variant_store$chrome <- "bottom"

#' Which layout variant is this page load rendering?
#' @return one of "a", "b", "c", "orig"
lab_variant <- function() .lab_variant_store$v

#' Which chrome is this page load rendering?
#' "bottom" (default) = single merged bar at the foot of the page.
#' "top"              = the app's original navbar + module bar, for comparison.
lab_chrome <- function() .lab_variant_store$chrome

#' Lab UI entry: read ?v= from the request, then delegate to the REAL app_ui().
#'
#' app_ui() is a function(request), so Shiny hands us the HTTP request and the
#' variant can change on a browser reload -- no container restart.
lab_app_ui <- function(request) {
    q <- tryCatch(shiny::parseQueryString(request$QUERY_STRING), error = function(e) list())
    v <- tolower(as.character(q[["v"]] %||% "orig"))
    if (!length(v) || is.na(v) || !v %in% c("a", "b", "c", "orig")) v <- "orig"
    .lab_variant_store$v <- v

    ch <- tolower(as.character(q[["chrome"]] %||% "bottom"))
    if (!length(ch) || is.na(ch) || !ch %in% c("bottom", "top")) ch <- "bottom"
    .lab_variant_store$chrome <- ch

    message("[lab] variant: ", v, " | chrome: ", ch)
    app_ui(request)
}

# ── Theme override: inject lab.scss AFTER the app's custom.scss ───────────────
# Absolute mounted paths, so this never depends on how app_sys() resolves.
app_theme <- function() {
    bslib::bs_theme(
        version = 5,
        preset  = "shiny",
        bg = "#FAFBFC", fg = "#2D3748",
        primary = "#1B7A6E", secondary = "#4A5568",
        success = "#38A169", info = "#3182CE",
        warning = "#D69E2E", danger = "#E53E3E",
        base_font    = bslib::font_google("Source Sans 3"),
        heading_font = bslib::font_google("Inter", wght = c(500, 600, 700)),
        code_font    = bslib::font_google("JetBrains Mono"),
        font_scale   = 0.95,
        "navbar-bg"          = "#1A2332",
        "navbar-dark-color"  = "#E2E8F0",
        "card-border-radius" = "0.5rem",
        "card-border-color"  = "#E2E8F0",
        "border-radius"      = "0.375rem",
        "accordion-border-color" = "#E2E8F0"
    ) |>
        bslib::bs_add_rules(sass::sass_file(
            "/pipeline/scripts/adaptogene.app/inst/app/www/custom.scss")) |>
        bslib::bs_add_rules(sass::sass_file(
            "/pipeline/scripts/layout_lab/www/lab.scss"))
}
