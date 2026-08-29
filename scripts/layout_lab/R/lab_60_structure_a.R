# lab_60_structure_a.R — Structure in the approved Variant A pattern.
#
# Anchor (left, 5/12): the ancestry PIEMAP -- the figure this module exists for --
# with the four controls that drive it on one inline strip directly above, since
# none of them affect anything in the right column.
# Right (7/12): plenty of plots.
#
# Like PreStructure, mod_structure_server() has no summary_boxes output, so the KPI
# strip comes from an additive server (lab_structure_extras). Nothing in
# scripts/adaptogene.app/ is modified.

mod_structure_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "structure",

        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # ── LEFT: controls + the piemap anchor ────────────────────────────
            htmltools::div(
                class = "lab-hero-col",

                # All four controls drive ONLY the piemap, so they belong in this
                # column rather than spanning the page -- and they sit on ONE
                # inline strip (.lab-mapbar), not in a bordered box.
                #
                # The boxed version measured 116px: a layout_columns of three
                # 52px selector cells, plus the Points switch on its own 23px
                # row, plus padding -- more chrome than the four controls need,
                # directly above the figure they drive. Shared with
                # Maladaptation's map strip; see .lab-mapbar in lab.scss.
                htmltools::div(
                    class = "lab-mapbar",
                    shiny::uiOutput(ns("bio_selector")),
                    shiny::uiOutput(ns("metric_selector")),
                    shiny::uiOutput(ns("zoom_selector")),
                    bslib::input_switch(ns("points"), "Points", value = FALSE)
                ),

                lab_hero(mod_piemap_viewer_ui(ns("piemap")))
            ),

            # ── RIGHT: plenty of plots ────────────────────────────────────────
            htmltools::div(
                class = "lab-multiples-col",

                # pop_stats_info is a filter_note badge ("3 pops" -> N
                # populations at K, samples per population). On its own line it
                # cost a full row to say one thing about the section beneath it,
                # so it rides the divider instead -- same aside slot as
                # PreStructure's K selector.
                lab_section_header("Population structure", icon = "diagram-3",
                                   aside = shiny::uiOutput(ns("pop_stats_info"))),
                # pca_structure_k is 1:1; mantel/amova are absent unless
                # Pop.calc_stats ran, so they render as placeholders -- kept
                # visible (not hidden) so the gap is legible.
                lab_thumb_grid(
                    cols = 3,
                    lab_thumb(mod_image_card_ui(ns("pca_structure_k"))),
                    lab_thumb(mod_image_card_ui(ns("mantel"))),
                    lab_thumb(mod_image_card_ui(ns("amova")))
                ),

                lab_section_header("Linkage disequilibrium", icon = "bar-chart-line"),
                shiny::uiOutput(ns("ld_decay_skipped_alert")),
                # Two ~1.4:1 curves at one column each, and the half-distance
                # table takes the remaining two columns -- the row fills exactly.
                lab_thumb_grid(
                    cols = 4,
                    lab_thumb(mod_image_card_ui(ns("ld_decay"))),
                    lab_thumb(mod_image_card_ui(ns("ld_decay_chr"))),
                    lab_thumb(
                        bslib::card(
                            class = "lab-dock-card lab-dock-ld",
                            bslib::card_header(bsicons::bs_icon("table"),
                                               " LD half-distances"),
                            bslib::card_body(DT::DTOutput(ns("ld_decay_table")))
                        ),
                        span = 2
                    )
                )
            )
        ),

        # Below the fold: the pop-stats table accordion. Collapsed -> suspended.
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE,
            bslib::accordion_panel(
                value = "pop_tables",
                title = htmltools::tagList(bsicons::bs_icon("table"),
                                           " Population statistics tables"),
                shiny::uiOutput(ns("tables_content"))
            )
        )
    )
}

#' Additive server: supplies ONLY the KPI strip, which the real server lacks.
lab_structure_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))

            sv <- function(step, metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == step & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))
            bp <- function(x) {
                v <- num(x)
                if (is.na(v)) "—"
                else if (v >= 1e6) paste0(round(v / 1e6, 2), " Mb")
                else if (v >= 1e3) paste0(round(v / 1e3, 1), " kb")
                else paste0(v, " bp")
            }

            k_best <- sv("structure", "K_best")
            if (is.na(k_best)) k_best <- if (!is.na(pd$k_best)) pd$k_best else "—"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("K best", as.character(k_best), "primary", "diagram-3"),
                lab_value_box("Climate predictors",
                              as.character(sv("structure", "n_climate_variables", "—")),
                              "info", "thermometer-half"),
                lab_value_box("LD half-decay",
                              bp(sv("structure", "ld_decay_half_decay_genome_wide")),
                              "success", "bar-chart-line"),
                lab_value_box("LD r² = 0.2",
                              bp(sv("structure", "ld_decay_r2_02_genome_wide")),
                              "info", "rulers")
            )
        })
    })
}
