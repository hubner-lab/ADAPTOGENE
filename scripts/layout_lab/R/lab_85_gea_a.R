# lab_85_gea_a.R — GEA in the Variant A pattern, adapted for an interactive module.
#
# The pattern bends here, deliberately:
#   * The combined Manhattan is a genome-wide x-axis. It gets the FULL width as
#     the hero rather than a 5/12 anchor column -- compressing it by 40% would
#     cost real resolution on the one plot that needs it.
#   * The threshold bar, WZA note, trait x method matrix and Save button move OUT
#     of the Manhattan card into a control bar above it. Today they are injected
#     between the card header and the plot (mod_gea.R filter_ui slot), which is
#     the one control convention that differs from every other module AND eats
#     the plot's height.
#   * Region Explorer is left exactly as-is, full width at the bottom. Its server
#     is 1342 lines; there is nothing to gain from re-housing it.
#
# Additive server (lab_gea_extras) supplies only the KPI strip.

mod_gea_ui_a <- function(id) {
    ns <- shiny::NS(id)

    lab_root("a", module = "gea",

        lab_kpi_row(ns, alert_id = NULL),

        shiny::uiOutput(ns("config_badges")),

        # ── Controls, above the plot they drive ───────────────────────────────
        htmltools::div(
            class = "control-bar lab-gea-controls",
            shiny::uiOutput(ns("threshold_bar")),
            shiny::uiOutput(ns("wza_collapse_note")),
            shiny::uiOutput(ns("filter_bar")),
            htmltools::div(
                class = "d-flex justify-content-end mt-1",
                shiny::actionButton(
                    ns("save_snp_set"), "Save SNP set for maladaptation",
                    class = "btn-sm btn-success"
                )
            )
        ),

        # ── HERO: combined Manhattan, full width, no injected filter_ui ───────
        htmltools::div(
            class = "lab-gea-hero",
            mod_manhattan_overlay_ui(ns("combined_manhattan"), height = "100%")
        ),

        shiny::uiOutput(ns("no_sig_snps_warning")),

        # ── Per-method detail ─────────────────────────────────────────────────
        lab_section_header("Per-method detail", icon = "layers"),
        htmltools::div(
            class = "lab-gea-methodbar d-flex align-items-center gap-3 flex-wrap",
            shiny::uiOutput(ns("method_tabs_ui")),
            shiny::uiOutput(ns("per_method_trait_ui"))
        ),
        bslib::layout_columns(
            col_widths = c(5, 7),
            fill = FALSE, fillable = FALSE,
            gap = "0.75rem",
            # QQ is square-ish -> the narrow column; the per-method Manhattan is
            # wide -> the broad one. Inverted vs the display modules on purpose.
            htmltools::div(class = "lab-gea-qq", mod_image_card_ui(ns("qq_plot"))),
            htmltools::div(class = "lab-gea-methodman",
                           mod_manhattan_overlay_ui(ns("method_manhattan"),
                                                    height = "100%"))
        ),

        # ── RDA diagnostics: collapsed (suspended until opened) ───────────────
        bslib::accordion(
            class = "lab-below-fold",
            open = FALSE,
            bslib::accordion_panel(
                value = "rda",
                title = htmltools::tagList(bsicons::bs_icon("compass"),
                                           " RDA diagnostics & candidates"),
                mod_rda_details_ui(ns("rda_details"))
            )
        ),

        # ── Region Explorer: unchanged, full width, bottom ────────────────────
        mod_region_explorer_ui(ns("region_explorer"))
    )
}

#' Additive server: the KPI strip only.
lab_gea_extras <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        output$summary_boxes <- shiny::renderUI({
            pd <- project_data()
            s  <- as.data.frame(load_pipeline_summary(pd$name))
            sv <- function(metric, default = NA) {
                if (nrow(s) == 0) return(default)
                row <- s[s$step == "gea" & s$metric == metric, ]
                if (nrow(row) == 0) default else row$value[1]
            }
            num <- function(x) suppressWarnings(as.numeric(x))

            # Methods actually run = the sig_snps_<METHOD> rows the summary wrote.
            meth <- if (nrow(s) == 0) character(0) else
                s$metric[s$step == "gea" & grepl("^sig_snps_", s$metric)]
            n_meth <- length(meth)

            n_snps      <- num(sv("selected_snps_total"))
            snps_theme  <- if (is.na(n_snps)) "secondary"
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
                              if (n_meth > 0) as.character(n_meth) else "—",
                              "info", "layers")
            )
        })
    })
}
