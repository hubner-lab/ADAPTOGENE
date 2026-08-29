# lab_89_malad_a.R — Maladaptation in the Variant A pattern.
#
# Three equal columns under the model bar (2026-08-27):
#   1. the genetic-offset piemap -- the module's deliverable -- with its map
#      controls inline above it;
#   2. the two importance figures, stacked;
#   3. the per-site offset table.
# Then the cross-model comparison in a collapsed accordion (opt-in, expensive).
#
# Equal thirds rather than a weighted split: all three are outputs of the SAME
# model, read together (map = where, importance = which predictor, table = the
# numbers). The previous 5/7 split docked the table under the map, which made it
# compete with the map for height and pushed the importance plots into a column
# of their own for no reason. All three columns share one height so the row ends
# flush -- see --lab-malad-row in lab.scss.
#
# mod_maladaptation.R has the longest UI in the app (137 lines) and many
# renderUI'd badges; all of their ids are preserved verbatim.
#
# There is NO KPI strip (dropped 2026-08-27). It carried four value boxes and
# three of them could never show a number: `n_snp_sets`, `n_sites` and
# `climate_scenario` are keys the lab invented, and nothing in the pipeline ever
# writes them to pipeline_summary.tsv -- measured on SIMDATA, they rendered
# "SNP sets — / Sites — / Scenario —". The fourth, "Offset methods 3", counted
# subdirectories of Maladaptation/tables/: a property of how the PIPELINE was
# configured, not of the analysis the user is looking at. Meanwhile the model
# bar directly below already carried six badges of real per-model fact (adaptive
# SNP count, ntree, cor.thr, spatial, random model, future scenario), so the
# strip was 69px of blanks sitting on top of the answer.
#
# What replaced it: `offset_summary` -- the genetic-offset range/mean/n-sites
# for the SELECTED model, which mod_maladaptation.R already renders as a
# filter_note badge -- moves OUT of the site-table card and onto the map strip,
# right of the Points switch. It is the one number that changes with the model
# selector, and it describes the map it now sits above.

mod_maladaptation_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "malad",

        # Model selector + the badge cluster the server renders. offset_summary
        # joins them here (it used to sit above the site table); it is the live
        # result for the selected model, so it reads as the bar's headline.
        htmltools::div(
            class = "control-bar lab-malad-modelbar d-flex align-items-center gap-2 flex-wrap",
            shiny::uiOutput(ns("model_selector")),
            shiny::uiOutput(ns("snp_count_badge")),
            shiny::uiOutput(ns("method_params_badge")),
            shiny::uiOutput(ns("climate_scenario_badge")),
            shiny::uiOutput(ns("sensitivity_badge")),
            htmltools::div(class = "ms-auto d-flex align-items-center gap-2",
                           shiny::uiOutput(ns("set_details_btn")),
                           shiny::uiOutput(ns("set_delete_btn")))
        ),

        bslib::layout_columns(
            col_widths = c(4, 4, 4),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── COLUMN 1: the offset map ──────────────────────────────────────
            htmltools::div(
                class = "lab-malad-col lab-malad-col-map",
                # One inline strip, not a control-bar box. The old version wrapped
                # the two selectors in a layout_columns and put the Points switch
                # on its own row beneath: on a project with neither a piemap
                # variant nor a zoom map (both renderUI's return NULL when there
                # is only one choice) that was a 64px bordered box containing two
                # zero-height grid columns and a 23px switch. Flex + gap means the
                # row is exactly as tall as whatever actually rendered.
                htmltools::div(
                    class = "lab-mapbar lab-malad-mapbar",
                    shiny::uiOutput(ns("piemap_variant_ui")),
                    shiny::uiOutput(ns("zoom_selector")),
                    bslib::input_switch(ns("points"), "Points", value = FALSE),
                    # The offset range/mean/n-sites badge sits with the map
                    # rather than in the model bar: it describes what this map
                    # is showing, and the map strip is where the eye already is
                    # when reading it.
                    htmltools::div(class = "lab-malad-offsetnote ms-auto",
                                   shiny::uiOutput(ns("offset_summary")))
                ),
                lab_hero(htmltools::div(class = "piemap-container",
                                        mod_image_card_ui(ns("offset_piemap"))))
            ),

            # ── COLUMN 2: predictor importance, stacked ───────────────────────
            htmltools::div(
                class = "lab-malad-col lab-malad-col-imp",
                mod_image_card_ui(ns("overall_importance")),
                mod_image_card_ui(ns("cumulative_importance"))
            ),

            # ── COLUMN 3: per-site offsets ────────────────────────────────────
            htmltools::div(
                class = "lab-malad-col lab-malad-col-tbl",
                bslib::card(
                    class = "lab-dock-card lab-dock-offset",
                    bslib::card_header(
                        class = "d-flex justify-content-between align-items-center",
                        htmltools::span(bsicons::bs_icon("geo-alt"), " Site offsets"),
                        shiny::downloadButton(ns("dl_site"), "CSV",
                                              class = "btn-sm btn-outline-secondary")
                    ),
                    bslib::card_body(DT::DTOutput(ns("site_table")))
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

#' Additive server for the Maladaptation lab layout.
#'
#' Nothing left to add: the KPI strip this used to feed is gone (see the header
#' note), and every remaining output is rendered by the real
#' mod_maladaptation_server(). Kept as a no-op so lab_90_dispatch.R's uniform
#' `<mod>_extras(id, project_data)` wiring does not need a special case.
lab_malad_extras <- function(id, project_data) {
    invisible(NULL)
}
