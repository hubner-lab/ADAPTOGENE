#' RDA Diagnostics & Candidates panel
#'
#' Second accordion panel inside mod_gea.R's per-method accordion, shown only
#' when "RDA" is a configured GEA method. mod_gea.R itself stays fully
#' method-agnostic (it derives methods/traits/colors from discovered files);
#' this module is a separate, RDA-specific detail panel — the same pattern as
#' mod_maladaptation.R's per-method comparison tabs.
#'
#' @param id module namespace id
#' @noRd
mod_rda_details_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("panel"))
}

#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @param methods reactive character vector of configured GEA methods (from mod_gea.R)
#' @noRd
mod_rda_details_server <- function(id, project_data, methods) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        has_rda <- shiny::reactive(isTRUE("RDA" %in% methods()))

        diagnostics <- shiny::reactive({
            shiny::req(has_rda())
            load_rda_diagnostics(project_data()$name)
        })
        anova_dt <- shiny::reactive({
            shiny::req(has_rda())
            load_rda_anova(project_data()$name)
        })
        candidates_dt <- shiny::reactive({
            shiny::req(has_rda())
            load_rda_candidates(project_data()$name)
        })

        d <- function(key, default = NA) {
            v <- diagnostics()[[key]]
            if (is.null(v) || identical(v, "")) default else v
        }

        # ── Value-box summary row ───────────────────────────────────────────
        output$value_boxes <- shiny::renderUI({
            shiny::req(has_rda())
            fit_mode <- d("fit_mode", "unknown")
            bslib::layout_column_wrap(
                width = 1 / 5, heights_equal = "row",
                bslib::value_box(
                    title = "Fit mode", value = fit_mode,
                    theme = if (identical(fit_mode, "pruned")) "danger" else "success",
                    showcase = bsicons::bs_icon("cpu")
                ),
                bslib::value_box(
                    title = "Markers fitted", value = d("n_markers_fitted", "?"),
                    theme = "info", showcase = bsicons::bs_icon("database-fill")
                ),
                bslib::value_box(
                    title = "RDA axes (K)",
                    value = sprintf("%s / %s", d("rda_axes", "?"), d("rda_axes_max", "?")),
                    theme = "info", showcase = bsicons::bs_icon("bounding-box")
                ),
                bslib::value_box(
                    title = "Adj. R2", value = {
                        v <- suppressWarnings(as.numeric(d("adj_r_squared", NA)))
                        if (is.na(v)) "?" else sprintf("%.4f", v)
                    },
                    theme = "info", showcase = bsicons::bs_icon("graph-up")
                ),
                bslib::value_box(
                    title = "GIF lambda", value = {
                        v <- suppressWarnings(as.numeric(d("gif_lambda", NA)))
                        if (is.na(v)) "?" else sprintf("%.3f", v)
                    },
                    theme = if (!is.na(suppressWarnings(as.numeric(d("gif_lambda", NA)))) &&
                               abs(as.numeric(d("gif_lambda")) - 1) > 0.2) "warning" else "success",
                    showcase = bsicons::bs_icon("rulers")
                )
            )
        })

        # ── Warning card: pruned fallback and/or above the validated envelope ─
        output$warning_card <- shiny::renderUI({
            shiny::req(has_rda())
            fit_mode <- d("fit_mode", "full")
            envelope <- d("marker_envelope_status", "")
            msgs <- character(0)
            if (identical(fit_mode, "pruned")) msgs <- c(msgs, trimws(d("fallback_warning", "")))
            if (identical(envelope, "ABOVE_VALIDATED")) msgs <- c(msgs, d("marker_envelope_note", ""))
            msgs <- msgs[nzchar(msgs)]
            if (length(msgs) == 0) return(NULL)
            bslib::card(
                class = "border-danger",
                bslib::card_header(class = "bg-danger-subtle",
                                   bsicons::bs_icon("exclamation-triangle-fill"), " RDA run notice"),
                bslib::card_body(
                    lapply(msgs, function(m) htmltools::tags$pre(class = "small mb-2", m))
                )
            )
        })

        # ── ANOVA table (degenerate rows highlighted) ────────────────────────
        output$anova_table <- DT::renderDT({
            shiny::req(has_rda())
            dt <- anova_dt()
            tbl <- safe_datatable(dt, options = list(pageLength = 15, scrollX = TRUE))
            if ("status" %in% names(dt)) {
                tbl <- DT::formatStyle(tbl, "status",
                    backgroundColor = DT::styleEqual(c("degenerate", "ok"),
                                                     c("rgba(229, 62, 62, 0.12)", "")))
            }
            tbl
        })

        # ── Plot cards ────────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data(); shiny::req(has_rda())
            mod_image_card_server(
                "screeplot",
                path  = shiny::reactive(rda_plot_path(pd$name, MOD_GEA, "screeplot")),
                title = shiny::reactive("RDA Screeplot"),
                dl_name = shiny::reactive("rda_screeplot"),
                suggestion = shiny::reactive("Run GEA with RDA in GEA.configs"),
                note = shiny::reactive(help_note("rda_screeplot", results = sprintf(
                    "K=%s (max %s), condition_pcs=%s", d("rda_axes", "?"), d("rda_axes_max", "?"),
                    d("condition_pcs", "?"))))
            )
            mod_image_card_server(
                "pval_hist",
                path  = shiny::reactive(rda_plot_path(pd$name, MOD_GEA, "pvalue_histogram")),
                title = shiny::reactive("RDA p-value Histogram"),
                dl_name = shiny::reactive("rda_pvalue_histogram"),
                suggestion = shiny::reactive("Run GEA with RDA in GEA.configs"),
                note = shiny::reactive(help_note("rda_pval_hist", results = sprintf(
                    "GIF lambda=%s, qvalue method=%s", d("gif_lambda", "?"), d("qvalue_method", "?"))))
            )
            mod_image_card_server(
                "biplot",
                path  = shiny::reactive(rda_plot_path(pd$name, MOD_GEA, "loadings_biplot")),
                title = shiny::reactive("RDA Loadings Biplot"),
                dl_name = shiny::reactive("rda_loadings_biplot"),
                suggestion = shiny::reactive("Run GEA with RDA in GEA.configs"),
                note = shiny::reactive(help_note("rda_biplot", results = sprintf(
                    "%s candidate(s) at bonf 0.05 / q 0.1", d("n_candidates_total", "?"))))
            )
        })

        # ── Candidates table (server-side paged — can be large at WGS density) ─
        output$candidates_table <- DT::renderDT({
            shiny::req(has_rda())
            safe_datatable(
                candidates_dt(),
                filter = "top", server = TRUE,
                options = list(pageLength = 25, scrollX = TRUE,
                               order = list(list(which(names(candidates_dt()) == "rank") - 1, "asc")))
            )
        }, server = TRUE)

        # ── Panel shell — hidden entirely when RDA isn't configured ──────────
        output$panel <- shiny::renderUI({
            if (!has_rda()) return(NULL)
            htmltools::tagList(
                shiny::uiOutput(ns("value_boxes")),
                shiny::uiOutput(ns("warning_card")),
                bslib::layout_column_wrap(
                    width = 1 / 3,
                    mod_image_card_ui(ns("screeplot")),
                    mod_image_card_ui(ns("pval_hist")),
                    mod_image_card_ui(ns("biplot"))
                ),
                bslib::card(
                    bslib::card_header("ANOVA (full / by-axis / by-margin)"),
                    DT::DTOutput(ns("anova_table"))
                ),
                bslib::card(
                    bslib::card_header("Candidates (LD-unclumped — see Region Explorer below for clumped regions)"),
                    DT::DTOutput(ns("candidates_table"))
                )
            )
        })
    })
}
