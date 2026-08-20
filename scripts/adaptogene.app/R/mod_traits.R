#' Phenotypic Factors tab UI
#'
#' The phenotypic counterpart of mod_climate.R: "what do my traits look like,
#' how do they covary, and does any of them carry no signal at all".
#'
#' Exists in BOTH regimes. mode=climate raises when Climate.enabled is false,
#' so in gwas_only — the regime where traits are the only factors a project
#' has — this is the only factor-characterization tab available. Everything
#' here comes from mode=traits (workflow/rules/traits.smk), which needs no
#' coordinates and no climate; the single climate-dependent product is the
#' joint traits x climate correlogram, and the switch that selects it hides
#' itself when that file is absent.
#'
#' Phenomaps are NOT produced here — they are mode=gwas outputs
#' (gwas.smk phenomap_{trait}) surfaced again on this tab, which is where a
#' user reasoning about traits expects to find them.
#'
#' @param id module namespace id
#' @noRd
mod_traits_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        bslib::card(
            bslib::card_header("Factors"),
            shiny::uiOutput(ns("trait_invariant_warning")),
            htmltools::div(
                class = "control-bar d-flex align-items-end gap-3 flex-wrap mb-2",
                shiny::selectInput(ns("figure"), "Figure",
                                   choices = c("Correlogram" = "corr",
                                               "Pairs plot"  = "pairs"),
                                   selected = "corr", width = "200px"),
                # Static (not inside renderUI) so its state survives re-renders;
                # shinyjs hides it when the project has no traits x climate file.
                shinyjs::hidden(htmltools::div(
                    id = ns("climate_toggle_wrap"),
                    bslib::input_switch(ns("with_climate"), "Include climate factors",
                                        value = FALSE)
                ))
            ),
            mod_image_card_ui(ns("trait_figure"))
        ),

        bslib::layout_column_wrap(
            width = 1 / 2,
            bslib::card(
                bslib::card_header("Distributions"),
                mod_image_card_ui(ns("phenotype_density"))
            ),
            bslib::card(
                bslib::card_header("Trait map"),
                htmltools::div(
                    class = "control-bar d-flex align-items-end gap-3 flex-wrap mb-2",
                    shiny::uiOutput(ns("phenomap_trait_selector")),
                    bslib::input_switch(ns("phenomap_points"), "Points", value = FALSE)
                ),
                htmltools::div(class = "piemap-container",
                               shiny::uiOutput(ns("phenomap_content")))
            )
        ),

        bslib::card(
            bslib::card_header("Tables"),
            bslib::accordion(
                open = FALSE, multiple = TRUE,
                bslib::accordion_panel("Trait summary",
                    DT::DTOutput(ns("trait_summary_table"))),
                bslib::accordion_panel("Invariant traits",
                    DT::DTOutput(ns("trait_invariant_table")))
            )
        )
    )
}

#' Phenotypic Factors tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_traits_server <- function(id, project_data) {
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

        trait_names <- shiny::reactive({
            p <- trait_table_path(project_data()$name, "trait_summary")
            if (!file_ok(p)) return(character(0))
            dt <- tryCatch(data.table::fread(p, sep = "\t", header = TRUE),
                           error = function(e) data.table::data.table())
            if (nrow(dt) == 0 || !"trait" %in% names(dt)) character(0) else as.character(dt$trait)
        })

        has_climate_corr <- shiny::reactive(
            file_ok(trait_corr_path(project_data()$name, with_climate = TRUE))
        )

        # ── Figure card ──────────────────────────────────────────────────────
        # The traits-only correlogram always exists once mode=traits has run;
        # the joint one only in a project that also has climate. Rather than
        # offering a switch that resolves to a missing file, hide it.
        shiny::observe({
            if (has_climate_corr()) {
                shinyjs::show("climate_toggle_wrap")
            } else {
                shinyjs::hide("climate_toggle_wrap")
                bslib::toggle_switch("with_climate", value = FALSE)
            }
        })

        figure_is_pairs <- shiny::reactive(identical(input$figure, "pairs"))

        shiny::observe({
            pd <- project_data()
            suggestion <- shiny::reactive("Run mode=traits.")

            fig_path <- shiny::reactive({
                if (figure_is_pairs()) {
                    trait_plot_path(pd$name, "trait_pairs")
                } else {
                    trait_corr_path(pd$name,
                                    with_climate = isTRUE(input$with_climate) && has_climate_corr())
                }
            })

            fig_title <- shiny::reactive({
                if (figure_is_pairs()) "Trait Pairs Plot"
                else if (isTRUE(input$with_climate) && has_climate_corr())
                    "Trait x Climate Correlogram"
                else "Trait Correlogram"
            })

            fig_note <- shiny::reactive({
                if (figure_is_pairs()) {
                    n <- length(trait_names())
                    # The pipeline replaces the grid with a placeholder above
                    # Traits.pairs_max_factors; explain that here rather than
                    # baking the count onto the image (CLAUDE.md rule 8).
                    help_note("trait_pairs",
                              results = if (n > 0) paste0(n, " trait(s) in this project") else NULL)
                } else if (isTRUE(input$with_climate) && has_climate_corr()) {
                    help_note("trait_correlogram_climate")
                } else {
                    help_note("trait_correlogram")
                }
            })

            mod_image_card_server("trait_figure",
                path       = fig_path,
                title      = fig_title,
                dl_name    = shiny::reactive(if (figure_is_pairs()) "trait_pairs" else "trait_correlogram"),
                suggestion = suggestion,
                note       = fig_note
            )

            mod_image_card_server("phenotype_density",
                path       = shiny::reactive(phenotype_density_path(pd$name)),
                title      = shiny::reactive("Phenotype Density"),
                dl_name    = shiny::reactive("phenotype_density"),
                suggestion = suggestion,
                note       = shiny::reactive(help_note("phenotype_density"))
            )
        })

        # ── Invariant-trait badge ────────────────────────────────────────────
        # Same shape as mod_climate.R's climate_invariant_warning: an invariant
        # factor is silently useless downstream, so it gets a visible count.
        output$trait_invariant_warning <- shiny::renderUI({
            p <- trait_table_path(project_data()$name, "trait_invariant")
            if (!file_ok(p)) return(NULL)
            dt <- tryCatch(data.table::fread(p, sep = "\t", header = TRUE),
                           error = function(e) data.table::data.table())
            if (nrow(dt) == 0 || !"predictor" %in% names(dt)) return(NULL)
            traits <- dt$predictor
            items <- lapply(traits, function(tr) {
                reason <- if ("reason" %in% names(dt)) {
                    r <- dt[dt$predictor == tr, ]$reason
                    if (length(r) > 0) paste0(" — ", r[1]) else ""
                } else ""
                htmltools::tags$li(htmltools::tags$code(tr), reason)
            })
            htmltools::div(
                class = "d-flex justify-content-end mb-2",
                filter_note(
                    paste0(length(traits), " invariant"),
                    htmltools::p(
                        htmltools::tags$strong(length(traits),
                            if (length(traits) == 1) "trait" else "traits",
                            "with no variance across samples"),
                        " — constant or all-missing. They are dropped from the ",
                        "correlogram and the pairs plot, and a GWAS on one cannot ",
                        "produce a signal.",
                        htmltools::tags$ul(class = "mb-0 mt-1", items)
                    ),
                    class = "bg-warning text-dark"
                )
            )
        })

        # ── Phenomap (mode=gwas product, shown here too) ─────────────────────
        output$phenomap_trait_selector <- shiny::renderUI({
            tr <- trait_names()
            if (length(tr) == 0) return(NULL)
            shiny::selectInput(ns("phenomap_trait"), "Trait",
                               choices = tr, selected = tr[1], width = "200px")
        })

        selected_phenomap_trait <- shiny::reactive({
            tr <- trait_names()
            if (length(tr) == 0) return(NULL)
            input$phenomap_trait %||% tr[1]
        })

        output$phenomap_content <- shiny::renderUI({
            tr <- selected_phenomap_trait()
            if (is.null(tr)) return(plot_placeholder("No traits found",
                                                     "Run mode=traits."))
            pd   <- project_data()
            path <- pheno_piemap_path(pd$name, tr, points = isTRUE(input$phenomap_points))
            if (file_ok(path)) {
                shiny::imageOutput(ns("phenomap_img"), height = "auto", width = "100%")
            } else {
                plot_placeholder("Trait map not available",
                                 "Run mode=gwas — phenotype piemaps need coordinates.")
            }
        })

        output$phenomap_img <- shiny::renderImage({
            tr <- shiny::req(selected_phenomap_trait())
            pd <- project_data()
            p  <- pheno_piemap_path(pd$name, tr, points = isTRUE(input$phenomap_points))
            shiny::validate(shiny::need(file_ok(p), "Trait map not found"))
            list(src = p, contentType = "image/png", width = "100%",
                 alt = paste("Phenomap", tr))
        }, deleteFile = FALSE)

        # ── Tables ───────────────────────────────────────────────────────────
        output$trait_summary_table <- DT::renderDataTable({
            render_tsv_table(trait_table_path(project_data()$name, "trait_summary"))
        })
        output$trait_invariant_table <- DT::renderDataTable({
            render_tsv_table(trait_table_path(project_data()$name, "trait_invariant"))
        })
    })
}
