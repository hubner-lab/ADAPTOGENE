# lab_80_traits_a.R — Phenotypic (traits) in the approved Variant A pattern.
#
# Anchor (left, 5/12): the trait correlogram / pairs figure -- the "how do these
# traits relate to each other (and to climate)" figure this module exists for --
# with its two controls above it and the trait summary table docked beneath,
# mirroring PreStructure's PCA + Q-matrix.
# Right (7/12): distributions across the top, trait map below.
#
# The two switches stay STATIC (outside renderUI) exactly as the real UI has
# them -- mod_traits.R:32-34 notes their state must survive re-renders, and
# shinyjs::hidden() toggles the climate one when no joint correlogram exists.
#
# Additive server (lab_traits_extras) supplies only the KPI strip.

mod_traits_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "traits",

        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── LEFT: the anchor ──────────────────────────────────────────────
            htmltools::div(
                class = "lab-hero-col",

                shiny::uiOutput(ns("trait_invariant_warning")),

                htmltools::div(
                    class = "control-bar lab-trait-bar d-flex align-items-center gap-3",
                    shiny::selectInput(
                        ns("figure"), "Figure",
                        choices  = c("Correlogram" = "correlogram",
                                     "Pairs plot"  = "pairs"),
                        selected = "correlogram"
                    ),
                    shinyjs::hidden(htmltools::div(
                        id = ns("climate_toggle_wrap"),
                        bslib::input_switch(ns("with_climate"),
                                            "Include climate factors", value = FALSE)
                    ))
                ),

                lab_hero(mod_image_card_ui(ns("trait_figure"))),

                bslib::card(
                    class = "lab-dock-card lab-dock-traitsum",
                    bslib::card_header(bsicons::bs_icon("table"), " Trait summary"),
                    bslib::card_body(DT::DTOutput(ns("trait_summary_table")))
                )
            ),

            # ── RIGHT: plenty of plots ────────────────────────────────────────
            htmltools::div(
                class = "lab-multiples-col",

                lab_section_header("Distributions", icon = "graph-up"),
                # density_plot_phenotypes is 3.5:1 -- a full-width row is the
                # one shape that fits it without stranding it.
                lab_thumb_grid(
                    cols = 1,
                    lab_thumb(htmltools::div(class = "lab-density-strip",
                                             mod_image_card_ui(ns("phenotype_density"))))
                ),

                lab_section_header("Trait map", icon = "geo-alt"),
                htmltools::div(
                    class = "control-bar lab-phenomap-bar d-flex align-items-center gap-3",
                    shiny::uiOutput(ns("phenomap_trait_selector")),
                    bslib::input_switch(ns("phenomap_points"), "Points", value = FALSE)
                ),
                htmltools::div(
                    class = "lab-traitmap-row",
                    lab_thumb_grid(
                        cols = 1,
                        lab_thumb(bslib::card(
                            class = "lab-phenomap-card",
                            bslib::card_header(bsicons::bs_icon("geo-alt"), " Trait map"),
                            bslib::card_body(
                                class = "p-2 text-center",
                                htmltools::div(class = "piemap-container",
                                               shiny::uiOutput(ns("phenomap_content")))
                            )
                        ))
                    )
                )
            )
        ),

        # Invariant traits below the fold. With 0 invariant traits the DT renders
        # an empty "No data available" box -- a message box where the KPI strip
        # already carries the number. Collapsed here, so it is reachable when a
        # project DOES have invariant traits and silent when it does not.
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE,
            bslib::accordion_panel(
                value = "invariant",
                title = htmltools::tagList(bsicons::bs_icon("exclamation-triangle"),
                                           " Invariant traits"),
                DT::DTOutput(ns("trait_invariant_table"))
            )
        )
    )
}

#' Additive server: the KPI strip only.
lab_traits_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "traits" & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))

            # Worst missingness across traits: not in pipeline_summary.tsv, so
            # derive it from trait_summary.tsv, which the app already loads.
            miss_txt <- "—"; miss_theme <- "secondary"
            p <- trait_table_path(pd$name, "trait_summary")
            if (file_ok(p)) {
                ts <- tryCatch(as.data.frame(data.table::fread(p)),
                               error = function(e) NULL)
                if (!is.null(ts) && "pct_missing" %in% names(ts) && nrow(ts) > 0) {
                    mx <- max(num(ts$pct_missing), na.rm = TRUE)
                    if (is.finite(mx)) {
                        miss_txt   <- paste0(round(mx, 1), "%")
                        miss_theme <- if (mx >= 20) "danger"
                                      else if (mx >= 5) "warning" else "success"
                    }
                }
            }

            n_inv       <- num(sv("n_invariant_traits"))
            n_inv_theme <- if (!is.na(n_inv) && n_inv > 0) "warning" else "success"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("Traits", as.character(sv("n_traits", "—")),
                              "primary", "rulers"),
                lab_value_box("With missing data",
                              as.character(sv("n_traits_with_missing", "—")),
                              "info", "question-octagon"),
                lab_value_box("Worst missingness", miss_txt, miss_theme, "percent"),
                lab_value_box("Invariant traits",
                              as.character(sv("n_invariant_traits", "—")),
                              n_inv_theme, "exclamation-triangle")
            )
        })
    })
}
