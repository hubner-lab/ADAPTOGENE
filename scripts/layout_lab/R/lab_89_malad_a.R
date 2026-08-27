# lab_89_malad_a.R — Maladaptation in the Variant A pattern.
#
# Anchor (left): the genetic-offset piemap -- the module's deliverable -- with
# the model selector and its badges above it, and the site-offset table docked
# beneath.
# Right: the two Gradient-Forest importance plots, then the cross-model
# comparison in a collapsed accordion (it is opt-in and expensive).
#
# mod_maladaptation.R has the longest UI in the app (137 lines) and many
# renderUI'd badges; all of their ids are preserved verbatim.

mod_maladaptation_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "malad",

        lab_kpi_row(ns, alert_id = NULL),

        # Model selector + the badge cluster the server renders
        htmltools::div(
            class = "control-bar lab-malad-modelbar d-flex align-items-center gap-2 flex-wrap",
            shiny::uiOutput(ns("model_selector")),
            shiny::uiOutput(ns("snp_count_badge")),
            shiny::uiOutput(ns("method_params_badge")),
            shiny::uiOutput(ns("climate_scenario_badge")),
            shiny::uiOutput(ns("sensitivity_badge")),
            shiny::uiOutput(ns("set_details_btn")),
            shiny::uiOutput(ns("set_delete_btn"))
        ),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── LEFT: the offset map + its table ──────────────────────────────
            htmltools::div(
                class = "lab-hero-col",
                htmltools::div(
                    class = "control-bar lab-piemap-bar",
                    bslib::layout_columns(
                        col_widths = c(6, 6),
                        shiny::uiOutput(ns("piemap_variant_ui")),
                        shiny::uiOutput(ns("zoom_selector"))
                    ),
                    htmltools::div(
                        class = "d-flex align-items-center gap-2 mt-1",
                        bslib::input_switch(ns("points"), "Points (clear map)",
                                            value = FALSE)
                    )
                ),
                lab_hero(htmltools::div(class = "piemap-container",
                                        mod_image_card_ui(ns("offset_piemap")))),
                bslib::card(
                    class = "lab-dock-card lab-dock-offset",
                    bslib::card_header(
                        class = "d-flex justify-content-between align-items-center",
                        htmltools::span(bsicons::bs_icon("table"), " Site offsets"),
                        shiny::downloadButton(ns("dl_site"), "CSV",
                                              class = "btn-sm btn-outline-secondary")
                    ),
                    bslib::card_body(
                        shiny::uiOutput(ns("offset_summary")),
                        DT::DTOutput(ns("site_table"))
                    )
                )
            ),

            # ── RIGHT: importance plots ───────────────────────────────────────
            htmltools::div(
                class = "lab-multiples-col",
                lab_section_header("Predictor importance", icon = "bar-chart-line"),
                lab_thumb_grid(
                    cols = 2,
                    lab_thumb(mod_image_card_ui(ns("overall_importance"))),
                    lab_thumb(mod_image_card_ui(ns("cumulative_importance")))
                )
            )
        ),

        # Cross-model comparison: opt-in and expensive -> collapsed.
        bslib::accordion(
            id = ns("malad_sections"),
            class = "lab-below-fold",
            open = FALSE, multiple = TRUE,
            bslib::accordion_panel(
                value = "compare",
                title = htmltools::tagList(bsicons::bs_icon("shuffle"),
                                           " Cross-model comparison"),
                compare_honesty_banner(),
                bslib::layout_columns(
                    col_widths = c(6, 6),
                    htmltools::div(class = "control-bar lab-cmp-a",
                                   shiny::uiOutput(ns("compare_model_a_ui"))),
                    htmltools::div(class = "control-bar lab-cmp-b",
                                   shiny::uiOutput(ns("compare_model_b_ui")))
                ),
                shiny::uiOutput(ns("compare_spatial_note")),
                htmltools::div(
                    class = "d-flex align-items-center gap-2 mb-2",
                    shiny::actionButton(ns("run_compare"), "Compare",
                                        class = "btn-sm btn-primary"),
                    shiny::uiOutput(ns("compare_cache_badge"))
                ),
                shiny::uiOutput(ns("compare_results_ui"))
            )
        )
    )
}

#' Additive server: the KPI strip only.
lab_malad_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = "—") {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "maladaptation" & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }
            # Which offset methods actually produced output.
            meths <- character(0)
            mp <- mod_path(pd$name, MOD_MALAD, "tables")
            if (dir.exists(mp)) meths <- list.dirs(mp, full.names = FALSE, recursive = FALSE)
            meths <- meths[nzchar(meths)]

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("Offset methods",
                              if (length(meths)) as.character(length(meths)) else "—",
                              "primary", "layers"),
                lab_value_box("SNP sets", sv("n_snp_sets"), "info", "collection"),
                lab_value_box("Sites", sv("n_sites"), "info", "geo-alt"),
                lab_value_box("Scenario", sv("climate_scenario"), "success", "thermometer-half")
            )
        })
    })
}
