# lab_70_climate_a.R — Environmental (climate) in the approved Variant A pattern.
#
# The user asked for BOTH the correlation heatmap and the variance-partitioning
# Venn in the left focus column: this module has two halves ("what do my
# predictors look like" and "how much do climate/structure/geography overlap"),
# and each half has exactly one headline figure.
#
# Right column carries the supporting plots: predictor distributions first,
# variance partitioning / spatial structure last.
#
# Additive server (lab_climate_extras) supplies the KPI strip and the
# fractions-bar card; the real mod_climate_server() runs untouched.

mod_climate_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "climate",

        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── LEFT: both headline figures, in focus ─────────────────────────
            htmltools::div(
                class = "lab-hero-col",

                shiny::uiOutput(ns("climate_invariant_warning")),
                lab_hero(mod_image_card_ui(ns("climate_heatmap"))),

                htmltools::div(
                    class = "lab-varpart-badges",
                    shiny::uiOutput(ns("variance_explained_badge")),
                    shiny::uiOutput(ns("dbmem_skip_warning")),
                    shiny::uiOutput(ns("confounding_badge"))
                ),
                lab_hero(mod_image_card_ui(ns("varpart_venn")))
            ),

            # ── RIGHT: plenty of plots ────────────────────────────────────────
            htmltools::div(
                class = "lab-multiples-col",

                lab_section_header("Predictor distributions", icon = "graph-up"),
                # density_plot_present is a multi-panel grid (5250x3600, 18.9 MP).
                # A full-width row keeps the per-predictor panels legible; a
                # standard tile height caps the width too and collapses it.
                htmltools::div(
                    class = "lab-density",
                    lab_thumb_grid(
                        cols = 1,
                        lab_thumb(mod_image_card_ui(ns("climate_density")))
                    )
                ),

                lab_section_header("Variance partitioning & spatial structure",
                                   icon = "pie-chart"),
                # Three ~1.40 plots -> three columns, matching their aspect ratio.
                lab_thumb_grid(
                    cols = 3,
                    lab_thumb(mod_image_card_ui(ns("px_barplot"))),
                    lab_thumb(mod_image_card_ui(ns("dbmem_screeplot"))),
                    lab_thumb(mod_image_card_ui(ns("dbmem_selection_path")))
                )
            )
        ),

        # Below the fold: the three varpart/dbMEM tables. Collapsed -> suspended.
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE, multiple = TRUE,
            bslib::accordion_panel("Variance partition", value = "vp",
                DT::DTOutput(ns("varpart_table"))),
            bslib::accordion_panel("Px per variable", value = "px",
                DT::DTOutput(ns("px_table"))),
            bslib::accordion_panel("dbMEM diagnostics", value = "dbmem",
                DT::DTOutput(ns("dbmem_diagnostics_table")))
        )
    )
}

#' Additive server: the KPI strip + the fractions-bar card the app never wired.
lab_climate_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {

        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "climate" & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))

            unexp     <- num(sv("varpart_unexplained_R2adj"))
            unexp_txt <- if (is.na(unexp)) "—" else paste0(round(unexp, 1), "%")
            # Nearly all variance unexplained is the interesting/bad case.
            unexp_theme <- if (is.na(unexp)) "secondary"
                           else if (unexp >= 90) "danger"
                           else if (unexp >= 75) "warning" else "success"

            # Climate UNIQUE variance is the headline number of this module:
            # what climate explains once structure and geography are removed.
            # No yes/no boxes here -- the KPI strip carries numbers, and the
            # confounding call is already made by confounding_badge below.
            clim      <- num(sv("varpart_climate_R2adj"))
            clim_txt  <- if (is.na(clim)) "—" else paste0(round(clim, 2), "%")
            clim_theme <- if (is.na(clim)) "secondary"
                          else if (clim >= 5) "success"
                          else if (clim >= 1) "warning" else "danger"

            n_inv       <- num(sv("n_invariant_predictors"))
            n_inv_theme <- if (!is.na(n_inv) && n_inv > 0) "warning" else "success"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("Climate unique (R²adj)", clim_txt, clim_theme, "thermometer-half"),
                lab_value_box("Unexplained", unexp_txt, unexp_theme, "pie-chart"),
                lab_value_box("Invariant predictors",
                              as.character(sv("n_invariant_predictors", "—")),
                              n_inv_theme, "exclamation-triangle"),
                lab_value_box("dbMEM vectors",
                              as.character(sv("dbMEM_n_mem_positive", "—")),
                              "info", "diagram-2")
            )
        })
    })
}
