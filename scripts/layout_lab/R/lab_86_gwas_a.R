# lab_86_gwas_a.R — GWAS in the Variant A pattern.
#
# Same shape as GEA (full-width Manhattan hero, controls extracted above it,
# Region Explorer untouched at the bottom), with the two GWAS-only differences:
#   * a phenotype map with its own trait/points controls -> the left anchor of
#     the split below the hero,
#   * no Save-SNP-set button and no RDA panel.
#
# The original's companion "click a SNP..." help card is dropped: it is static
# instructional text, which CLAUDE.md rule 8 puts in a note badge, not a card
# occupying a third of a row.

mod_gwas_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "gwas",

        lab_kpi_row(ns, alert_id = NULL),

        shiny::uiOutput(ns("config_badges")),
        shiny::uiOutput(ns("pheno_missing_alert")),

        htmltools::div(
            class = "control-bar lab-gea-controls",
            shiny::uiOutput(ns("threshold_bar")),
            shiny::uiOutput(ns("wza_collapse_note")),
            shiny::uiOutput(ns("filter_bar"))
        ),

        htmltools::div(
            class = "lab-gea-hero",
            mod_manhattan_overlay_ui(ns("combined_manhattan"), height = "100%")
        ),

        shiny::uiOutput(ns("no_sig_snps_warning")),

        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",

            # LEFT: phenotype map + the two controls that drive only it
            htmltools::div(
                class = "lab-hero-col",
                lab_section_header("Phenotype map", icon = "geo-alt"),
                htmltools::div(
                    class = "control-bar lab-phenomap-bar d-flex align-items-center gap-3",
                    shiny::uiOutput(ns("trait_selector")),
                    bslib::input_switch(ns("points"), "Points", value = FALSE)
                ),
                bslib::card(
                    class = "lab-phenomap-card",
                    bslib::card_header(bsicons::bs_icon("geo-alt"), " Phenotype map"),
                    bslib::card_body(
                        class = "p-2 text-center",
                        htmltools::div(class = "piemap-container",
                                       shiny::uiOutput(ns("phenomap_content")))
                    )
                )
            ),

            # RIGHT: per-method detail
            htmltools::div(
                class = "lab-multiples-col",
                lab_section_header("Per-method detail", icon = "layers"),
                htmltools::div(
                    class = "lab-gea-methodbar d-flex align-items-center gap-3 flex-wrap",
                    shiny::uiOutput(ns("method_tabs_ui")),
                    shiny::uiOutput(ns("per_method_trait_ui"))
                ),
                bslib::layout_columns(
                    col_widths = c(4, 8),
                    fill = FALSE, fillable = FALSE,
                    gap = "0.5rem",
                    htmltools::div(class = "lab-gea-qq",
                                   mod_image_card_ui(ns("qq_plot"))),
                    htmltools::div(class = "lab-gea-methodman",
                                   mod_manhattan_overlay_ui(ns("method_manhattan"),
                                                            height = "100%"))
                )
            )
        ),

        mod_region_explorer_ui(ns("region_explorer"))
    )
}

#' Additive server: the KPI strip only.
lab_gwas_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "gwas" & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))

            meth <- if (nrow(s) == 0) character(0) else
                s$metric[s$step == "gwas" & grepl("^sig_snps_", s$metric)]
            n_snps     <- num(sv("selected_snps_total"))
            snps_theme <- if (is.na(n_snps)) "secondary"
                          else if (n_snps == 0) "danger" else "primary"

            bslib::layout_column_wrap(
                width = 1 / 4, fill = FALSE,
                lab_value_box("Significant SNPs",
                              as.character(sv("selected_snps_total", "—")),
                              snps_theme, "bullseye"),
                lab_value_box("Regions", as.character(sv("regions_combined", "—")),
                              "info", "geo-alt"),
                lab_value_box("Genes", as.character(sv("genes_found", "—")),
                              "success", "diagram-2"),
                lab_value_box("Methods",
                              if (length(meth) > 0) as.character(length(meth)) else "—",
                              "info", "layers")
            )
        })
    })
}
