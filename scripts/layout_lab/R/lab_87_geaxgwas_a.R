# lab_87_geaxgwas_a.R — GEAxGWAS in the Variant A pattern.
#
# The Miami plot is the hero (GEA above the zero line, GWAS below), so it takes
# the full width like the Manhattans in GEA/GWAS. The two filter panels move
# SIDE BY SIDE above it instead of sandwiching it: sandwiching is semantically
# cute (above/below matches the plot) but costs two full-width blocks of vertical
# space, and the panels are labelled anyway.
#
# There are no `gea_x_gwas` rows in pipeline_summary.tsv -- overlaps are computed
# interactively in Shiny, not pipeline-side (CLAUDE.md). So the KPI strip
# compares the two SIDES from their own summaries, which is what this module is
# actually about.
#
# The original's alert.alert-info explaining the fill buttons is dropped: static
# instructional prose, which rule 8 puts in a note badge, not a page banner.

mod_gea_x_gwas_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "geaxgwas",

        lab_kpi_row(ns, alert_id = NULL),

        bslib::layout_columns(
            col_widths = c(6, 6),
            fill = FALSE, fillable = FALSE,
            gap = "0.6rem",

            bslib::card(
                class = "lab-xg-filter lab-xg-gea",
                bslib::card_header(
                    class = "d-flex justify-content-between align-items-center",
                    htmltools::span(bsicons::bs_icon("arrow-up"),
                                    " GEA filter (above the zero line)"),
                    shiny::actionButton(ns("fill_gea"), "Fill from GEA",
                                        class = "btn-sm btn-outline-secondary")
                ),
                bslib::card_body(
                    shiny::uiOutput(ns("gea_threshold_bar")),
                    shiny::uiOutput(ns("gea_filter_bar"))
                )
            ),

            bslib::card(
                class = "lab-xg-filter lab-xg-gwas",
                bslib::card_header(
                    class = "d-flex justify-content-between align-items-center",
                    htmltools::span(bsicons::bs_icon("arrow-down"),
                                    " GWAS filter (below the zero line)"),
                    shiny::actionButton(ns("fill_gwas"), "Fill from GWAS",
                                        class = "btn-sm btn-outline-secondary")
                ),
                bslib::card_body(
                    shiny::uiOutput(ns("gwas_threshold_bar")),
                    shiny::uiOutput(ns("gwas_filter_bar"))
                )
            )
        ),

        # ── HERO: the Miami plot ──────────────────────────────────────────────
        htmltools::div(
            class = "lab-gea-hero",
            mod_manhattan_overlay_ui(ns("miami"), height = "100%")
        ),

        htmltools::div(
            class = "control-bar lab-xg-bounds d-flex align-items-center gap-3",
            htmltools::span(class = "lab-xg-bounds-label", "Overlap bounds"),
            shiny::radioButtons(
                ns("overlap_bounds"), label = NULL,
                choices = c("Union" = "union", "Intersection" = "intersection",
                            "GEA only" = "gea", "GWAS only" = "gwas"),
                selected = "union", inline = TRUE
            )
        ),

        mod_region_explorer_ui(ns("region_explorer")),

        bslib::accordion(
            id = ns("pairwise_accordion"),
            class = "lab-below-fold",
            open = FALSE,
            bslib::accordion_panel(
                value = "pairwise",
                title = htmltools::tagList(bsicons::bs_icon("grid-3x3"),
                                           " Pairwise trait overlap"),
                mod_pairwise_overlap_ui(ns("pairwise"))
            )
        )
    )
}

#' Additive server: a KPI strip comparing the two sides.
lab_geaxgwas_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(step, metric, default = "—") {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == step & s$metric == metric, ]
                if (nrow(row) == 0) default else as.character(row$value[1])
            }
            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("GEA SNPs",  sv("gea",  "selected_snps_total"),
                              "primary", "arrow-up"),
                lab_value_box("GEA regions", sv("gea", "regions_combined"),
                              "info", "geo-alt"),
                lab_value_box("GWAS SNPs", sv("gwas", "selected_snps_total"),
                              "success", "arrow-down"),
                lab_value_box("GWAS regions", sv("gwas", "regions_combined"),
                              "info", "geo-alt")
            )
        })
    })
}
