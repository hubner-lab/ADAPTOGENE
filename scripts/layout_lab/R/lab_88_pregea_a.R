# lab_88_pregea_a.R — PreGEA in the Variant A pattern.
#
# This module exists to RECOMMEND hyperparameters, so the KPI strip carries the
# four recommendations themselves (LFMM K, EMMAX #PCs, RDA Condition() PCs, RDA
# axes) and the recommendations table is the left anchor -- it is the module's
# summary table in the most literal sense.
#
# Right column keeps the app's own navset (LFMM / EMMAX / RDA): three tabs of
# diagnostics, where the inactive tabs stay in the DOM but suspend, so 13 figures
# never load at once.
#
# Tile shapes follow the figures: lambda/hits are 1.56 landscape, the LFMM
# histogram/QQ grids are 0.91, and the EMMAX grids are 0.45 -- nearly 1:2.2
# portrait -- so they get their own taller row rather than a shared height.

.lab_pregea_note_table <- function(ns, note_id, table_id, title) {
    htmltools::div(
        htmltools::div(
            class = "d-flex align-items-center gap-2 mt-2 mb-1",
            htmltools::span(class = "lab-pregea-tbl-title", title),
            shiny::uiOutput(ns(note_id), inline = TRUE)
        ),
        DT::DTOutput(ns(table_id))
    )
}

mod_pregea_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "pregea",

        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── LEFT: the recommendations, this module's whole output ─────────
            htmltools::div(
                class = "lab-hero-col",
                lab_section_header("Recommended hyperparameters", icon = "sliders"),
                bslib::card(
                    class = "lab-dock-card lab-dock-recs",
                    bslib::card_header(bsicons::bs_icon("table"),
                                       " Recommendations"),
                    bslib::card_body(DT::DTOutput(ns("recommendations_table")))
                )
            ),

            # ── RIGHT: per-method ladders (inactive tabs suspend) ─────────────
            htmltools::div(
                class = "lab-multiples-col",
                bslib::navset_card_tab(
                    id = ns("pregea_tabs"),

                    bslib::nav_panel(
                        "LFMM K ladder", value = "lfmm",
                        lab_thumb_grid(
                            cols = 3,
                            lab_thumb(mod_image_card_ui(ns("lfmm_screeplot"))),
                            lab_thumb(mod_image_card_ui(ns("lfmm_lambda"))),
                            lab_thumb(mod_image_card_ui(ns("lfmm_hits")))
                        ),
                        htmltools::div(
                            class = "lab-pregea-grids",
                            lab_thumb_grid(
                                cols = 2,
                                lab_thumb(mod_image_card_ui(ns("lfmm_hist"))),
                                lab_thumb(mod_image_card_ui(ns("lfmm_qq")))
                            )
                        ),
                        .lab_pregea_note_table(ns, "lfmm_table_note", "lfmm_table",
                                               "LFMM ladder")
                    ),

                    bslib::nav_panel(
                        "EMMAX / GAPIT #PC ladder", value = "emmax",
                        lab_thumb_grid(
                            cols = 2,
                            lab_thumb(mod_image_card_ui(ns("emmax_lambda"))),
                            lab_thumb(mod_image_card_ui(ns("emmax_hits")))
                        ),
                        # 1720x3820 -- these are nearly 1:2.2 portrait strips
                        htmltools::div(
                            class = "lab-pregea-tallgrids",
                            lab_thumb_grid(
                                cols = 2,
                                lab_thumb(mod_image_card_ui(ns("emmax_hist"))),
                                lab_thumb(mod_image_card_ui(ns("emmax_qq")))
                            )
                        ),
                        .lab_pregea_note_table(ns, "emmax_table_note", "emmax_table",
                                               "EMMAX ladder")
                    ),

                    bslib::nav_panel(
                        "RDA setup", value = "rda",
                        htmltools::div(
                            class = "lab-pregea-wide",
                            lab_thumb_grid(
                                cols = 1,
                                lab_thumb(mod_image_card_ui(ns("rda_comparison")))
                            )
                        ),
                        htmltools::div(
                            class = "control-bar lab-k-bar",
                            shiny::uiOutput(ns("rda_model_selector"))
                        ),
                        lab_thumb_grid(
                            cols = 2,
                            lab_thumb(mod_image_card_ui(ns("rda_biplot"))),
                            lab_thumb(mod_image_card_ui(ns("rda_screeplot")))
                        ),
                        bslib::accordion(
                            open = FALSE, multiple = TRUE, class = "mt-2",
                            bslib::accordion_panel("Axis ANOVA", value = "axis",
                                DT::DTOutput(ns("rda_model_axis_table"))),
                            bslib::accordion_panel("Condition() ladder", value = "ladder",
                                DT::DTOutput(ns("rda_ladder_table"))),
                            bslib::accordion_panel("Predictor collinearity", value = "collin",
                                DT::DTOutput(ns("rda_collinearity_table")))
                        )
                    )
                )
            )
        ),

        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE,
            bslib::accordion_panel(
                value = "transfer",
                title = htmltools::tagList(bsicons::bs_icon("shield-check"),
                                           " Transfer guard (full-set validation)"),
                mod_image_card_ui(ns("transfer_guard_plot")),
                DT::DTOutput(ns("transfer_guard_table"))
            )
        )
    )
}

#' Additive server: the recommendations as a KPI strip.
lab_pregea_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = "—") {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "pregea" & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }
            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("LFMM K", sv("LFMM_K_recommended"),
                              "primary", "diagram-3"),
                lab_value_box("EMMAX #PCs", sv("EMMAX_n_pcs_recommended"),
                              "info", "sliders"),
                lab_value_box("RDA Condition() PCs", sv("RDA_condition_pcs_recommended"),
                              "info", "funnel"),
                lab_value_box("RDA axes", sv("RDA_axes_recommended"),
                              "success", "compass")
            )
        })
    })
}
