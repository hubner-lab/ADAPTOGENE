#' Climate tab UI
#'
#' Predictor characterization: correlation/density plots, invariant-predictor
#' detection, and spatial structure (dbMEM + variance partitioning). Split out
#' of the old preGEA module so each module answers one question — Climate
#' answers "what do my predictors look like and how much do they overlap with
#' geography/structure", PreGEA (mod_pregea.R) answers "how many K / #PCs".
#'
#' Producers: correlation heatmap / density plots / invariant-predictor table
#' come from structure.smk (unchanged, moved display only in an earlier pass);
#' dbMEM + varpart come from climate.smk (mode=climate, this tab's own Run
#' button) — see workflow/rules/climate.smk.
#'
#' @param id module namespace id
#' @noRd
mod_climate_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        bslib::card(
            bslib::card_header("Predictors"),
            shiny::uiOutput(ns("climate_invariant_warning")),
            mod_image_card_ui(ns("climate_heatmap")),
            bslib::layout_column_wrap(
                width = 1 / 2,
                mod_image_card_ui(ns("climate_density")),
                mod_image_card_ui(ns("phenotype_density"))
            )
        ),

        bslib::card(
            bslib::card_header("Variance Partitioning"),
            bslib::layout_column_wrap(
                width = 1 / 3,
                mod_image_card_ui(ns("dbmem_screeplot")),
                mod_image_card_ui(ns("dbmem_selection_path")),
                mod_image_card_ui(ns("varpart_venn"))
            ),
            mod_image_card_ui(ns("px_barplot")),
            htmltools::h6("Tables", class = "mt-3"),
            bslib::accordion(
                open = FALSE, multiple = TRUE,
                bslib::accordion_panel("Variance-partition fractions",
                    DT::DTOutput(ns("varpart_fractions_table"))),
                bslib::accordion_panel("Testable-fraction significance",
                    DT::DTOutput(ns("varpart_anova_table"))),
                bslib::accordion_panel("Px per variable",
                    DT::DTOutput(ns("px_table"))),
                bslib::accordion_panel("dbMEM diagnostics",
                    DT::DTOutput(ns("dbmem_diagnostics_table")))
            )
        )
    )
}

#' Climate tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_climate_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        render_tsv_table <- function(path) {
            dt <- if (!file_ok(path)) data.table::data.table() else
                tryCatch(data.table::fread(path, sep = "\t", header = TRUE),
                         error = function(e) data.table::data.table())
            safe_datatable(as.data.frame(dt),
                           options = list(dom = "t", scrollX = TRUE, pageLength = -1),
                           rownames = FALSE)
        }

        # ── Predictors ───────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Run mode=climate.")

            # Collinearity-drop badge, next to the standard help_note badge —
            # the old rda_predictor_collinearity.png is retired (redundant
            # with this heatmap); which predictors the |r| pre-screen dropped
            # (PreGEA.Advanced.collinearity_r) is still on disk in
            # rda_predictor_collinearity.tsv and shown here instead.
            collin_note <- shiny::reactive({
                path <- pregea_table_path(pd$name, "rda", "rda_predictor_collinearity")
                if (!file_ok(path)) return(help_note("climate_heatmap"))
                dt <- tryCatch(data.table::fread(path, sep = "\t", header = TRUE),
                               error = function(e) data.table::data.table())
                dropped <- if (nrow(dt) > 0 && "action" %in% names(dt)) dt[dt$action == "dropped", ] else dt[0]
                extra_badge <- if (nrow(dropped) > 0) {
                    filter_note(
                        paste0(nrow(dropped), " dropped (PreGEA)"),
                        htmltools::p(htmltools::strong(nrow(dropped), "predictor(s) dropped for collinearity"),
                            " before the PreGEA RDA fit (PreGEA.Advanced.collinearity_r): ",
                            paste(dropped$predictor, collapse = ", "), ". ",
                            "Never auto-applied to this heatmap — informational only."),
                        class = "bg-warning text-dark"
                    )
                } else NULL
                htmltools::tagList(help_note("climate_heatmap"), extra_badge)
            })

            mod_image_card_server("climate_heatmap",
                path       = shiny::reactive(climate_heatmap_path(pd$name)),
                title      = shiny::reactive("Climate Correlation Heatmap"),
                dl_name    = shiny::reactive("climate_heatmap"),
                suggestion = suggestion,
                note       = collin_note
            )
            mod_image_card_server("climate_density",
                path       = shiny::reactive(climate_density_path(pd$name)),
                title      = shiny::reactive("Climate Density (Present)"),
                dl_name    = shiny::reactive("climate_density_present"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_density"))
            )
            mod_image_card_server("phenotype_density",
                path       = shiny::reactive(phenotype_density_path(pd$name)),
                title      = shiny::reactive("Phenotype Density"),
                dl_name    = shiny::reactive("phenotype_density"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("phenotype_density"))
            )

            output$climate_invariant_warning <- shiny::renderUI({
                p <- climate_invariant_path(pd$name)
                if (!file.exists(p)) return(NULL)
                dt <- tryCatch(
                    data.table::fread(p, sep = "\t", header = TRUE),
                    error = function(e) data.table::data.table()
                )
                if (nrow(dt) == 0 || !"predictor" %in% names(dt)) return(NULL)
                preds <- dt$predictor
                items <- lapply(preds, function(pr) {
                    reason <- if ("reason" %in% names(dt)) {
                        r <- dt[dt$predictor == pr, ]$reason
                        if (length(r) > 0) paste0(" — ", r[1]) else ""
                    } else ""
                    htmltools::tags$li(htmltools::tags$code(pr), reason)
                })
                htmltools::div(
                    class = "d-flex justify-content-end mb-2",
                    filter_note(
                        paste0(length(preds), " excluded"),
                        htmltools::p(
                            htmltools::tags$strong(length(preds),
                                if (length(preds) == 1) "climate predictor" else "climate predictors",
                                "excluded due to zero variance at sample sites"),
                            " — not used in GEA or Gradient Forest analyses. ",
                            "Consider removing from ", htmltools::tags$code("Climate.predictors"), ".",
                            htmltools::tags$ul(class = "mb-0 mt-1", items)
                        ),
                        class = "bg-warning text-dark"
                    )
                )
            })
        })

        # ── Variance partitioning ────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Run mode=climate.")

            mod_image_card_server("dbmem_screeplot",
                path       = shiny::reactive(climate_plot_path(pd$name, "spatial", "dbmem_screeplot")),
                title      = shiny::reactive("dbMEM Screeplot"),
                dl_name    = shiny::reactive("dbmem_screeplot"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_dbmem_screeplot"))
            )
            mod_image_card_server("dbmem_selection_path",
                path       = shiny::reactive(climate_plot_path(pd$name, "varpart", "dbmem_selection_path")),
                title      = shiny::reactive("dbMEM Forward-selection Path"),
                dl_name    = shiny::reactive("dbmem_selection_path"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_dbmem_selection_path"))
            )
            mod_image_card_server("varpart_venn",
                path       = shiny::reactive(climate_plot_path(pd$name, "varpart", "varpart_venn")),
                title      = shiny::reactive("Variance-partition Venn"),
                dl_name    = shiny::reactive("varpart_venn"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_varpart_venn"))
            )
            mod_image_card_server("px_barplot",
                path       = shiny::reactive(climate_plot_path(pd$name, "varpart", "px_barplot")),
                title      = shiny::reactive("Px per Predictor (Lasky et al. 2012)"),
                dl_name    = shiny::reactive("px_barplot"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("climate_px_barplot"))
            )
        })
        output$varpart_fractions_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(climate_table_path(pd$name, "varpart", "varpart_fractions"))
        })
        output$varpart_anova_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(climate_table_path(pd$name, "varpart", "varpart_anova"))
        })
        output$px_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(climate_table_path(pd$name, "varpart", "px_per_variable"))
        })
        output$dbmem_diagnostics_table <- DT::renderDataTable({
            pd <- project_data()
            render_tsv_table(climate_table_path(pd$name, "spatial", "dbmem_diagnostics"))
        })
    })
}
