#' Maladaptation tab UI
#'
#' Gradient Forest results: importance plots, genetic offset piemaps,
#' zoom maps, and site-level offset table.
#'
#' @param id module namespace id
#' @noRd
mod_maladaptation_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # GF Run selector at top (controls everything)
        htmltools::div(
            class = "control-bar",
            shiny::uiOutput(ns("suffix_selector"))
        ),

        # Importance plots
        bslib::layout_column_wrap(
            width = 1 / 2,
            mod_image_card_ui(ns("overall_importance")),
            mod_image_card_ui(ns("cumulative_importance"))
        ),

        shiny::hr(),

        # Piemap Type + Zoom right before piemap
        htmltools::div(
            class = "control-bar",
            bslib::layout_columns(
                col_widths = c(6, 6),
                shiny::uiOutput(ns("piemap_variant_ui")),
                shiny::uiOutput(ns("zoom_selector"))
            )
        ),

        # Genetic offset piemap (constrained to 2/3 width; full-screen for detail)
        bslib::layout_columns(
            col_widths = c(8, 4),
            htmltools::div(
                class = "piemap-container",
                mod_image_card_ui(ns("offset_piemap"))
            ),
            bslib::card(
                bslib::card_body(
                    class = "text-muted small",
                    htmltools::p(bsicons::bs_icon("info-circle"), " Genetic offset estimates local population vulnerability to projected climate change."),
                    htmltools::p(bsicons::bs_icon("arrows-fullscreen"), " Use the expand icon for a full-resolution view of the map."),
                    htmltools::p(bsicons::bs_icon("table"), " Site-level offset values are in the table below.")
                )
            )
        ),

        # Site table in accordion
        bslib::accordion(
            id       = ns("malad_sections"),
            open     = FALSE,
            multiple = TRUE,

            bslib::accordion_panel(
                "Site Offset Table",
                value = "site_table",
                icon  = bsicons::bs_icon("table"),
                bslib::card_body(
                    shiny::downloadButton(ns("dl_site"), "Download CSV",
                                          class = "btn-sm btn-outline-secondary mb-2"),
                    DT::DTOutput(ns("site_table"))
                )
            )
        )
    )
}

#' Maladaptation tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_maladaptation_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        gf_suffixes <- shiny::reactive(find_gf_suffixes(project_data()$name))

        output$suffix_selector <- shiny::renderUI({
            suf <- gf_suffixes()
            if (length(suf) == 0)
                return(shiny::p("No GF runs found.", class = "text-muted small"))
            shiny::selectInput(ns("suffix"), "GF Run", choices = suf, selected = suf[1])
        })

        selected_suffix <- shiny::reactive({
            s <- input$suffix
            if (is.null(s)) {
                suf <- gf_suffixes()
                if (length(suf) == 0) return(NULL)
                suf[1]
            } else s
        })

        # ── Piemap variant selector (only shows variants that exist) ───────────
        output$piemap_variant_ui <- shiny::renderUI({
            suf <- shiny::req(selected_suffix())
            pd  <- project_data()
            choices <- c("Genetic Offset" = "base")
            if (file_ok(gf_offset_piemap_path(pd$name, suf, "tajima_d")))
                choices <- c(choices, c("Offset \u00d7 Tajima's D" = "tajima_d"))
            if (file_ok(gf_offset_piemap_path(pd$name, suf, "pi_diversity")))
                choices <- c(choices, c("Offset \u00d7 Pi Diversity" = "pi_diversity"))
            if (length(choices) == 1) return(NULL)  # only base — no selector needed
            shiny::selectInput(ns("piemap_variant"), "Piemap Type",
                               choices = choices, selected = "base")
        })

        # ── Zoom selector (only shows if zoom maps exist) ──────────────────────
        output$zoom_selector <- shiny::renderUI({
            suf   <- shiny::req(selected_suffix())
            pd    <- project_data()
            zooms <- find_gf_zooms(pd$name, suf)
            if (length(zooms) == 0) return(NULL)
            shiny::selectInput(ns("zoom"), "Zoom Region",
                choices  = c("Full view" = "none", setNames(zooms, zooms)),
                selected = "none"
            )
        })

        selected_variant <- shiny::reactive(input$piemap_variant %||% "base")
        selected_zoom    <- shiny::reactive(input$zoom %||% "none")

        # ── Importance images ──────────────────────────────────────────────────
        shiny::observe({
            suf <- shiny::req(selected_suffix())
            pd  <- project_data()

            mod_image_card_server("overall_importance",
                path    = shiny::reactive(gf_importance_path(pd$name, suf, "overall")),
                title   = shiny::reactive("Overall Variable Importance"),
                dl_name = shiny::reactive(paste0("overall_importance_", suf))
            )
            mod_image_card_server("cumulative_importance",
                path    = shiny::reactive(gf_importance_path(pd$name, suf, "cumulative")),
                title   = shiny::reactive("Cumulative Importance"),
                dl_name = shiny::reactive(paste0("cumulative_importance_", suf))
            )
        })

        # ── Single selector-driven offset piemap ───────────────────────────────
        piemap_path <- shiny::reactive({
            suf     <- shiny::req(selected_suffix())
            pd      <- project_data()
            variant <- selected_variant()
            zoom    <- selected_zoom()
            if (zoom != "none") {
                mod_path(pd$name, MOD_MALAD, "plots", suf, "zoom",
                         paste0(zoom, ".png"))
            } else {
                gf_offset_piemap_path(pd$name, suf, variant)
            }
        })

        piemap_title <- shiny::reactive({
            variant <- selected_variant()
            zoom    <- selected_zoom()
            label   <- c(
                base         = "Genetic Offset",
                tajima_d     = "Offset \u00d7 Tajima's D",
                pi_diversity = "Offset \u00d7 Pi Diversity"
            )[variant]
            if (zoom != "none") paste0(label, " \u2014 Zoom: ", zoom) else label
        })

        mod_image_card_server("offset_piemap",
            path        = piemap_path,
            title       = piemap_title,
            dl_name     = shiny::reactive(paste0("offset_piemap_", selected_variant(),
                                                  "_", selected_suffix() %||% "gf")),
            placeholder = shiny::reactive("Run mode=maladaptation to generate Gradient Forest results")
        )

        # ── Site table ─────────────────────────────────────────────────────────
        site_data <- shiny::reactive({
            suf <- selected_suffix()
            if (is.null(suf)) return(data.table::data.table())
            load_cached(paste0("gf_site_", project_data()$name, "_", suf), function() {
                p <- gf_site_table_path(project_data()$name, suf)
                if (!file_ok(p)) return(data.table::data.table())
                dt <- data.table::fread(p, colClasses = c("site" = "character",
                                                            "sample" = "character"))
                # Join pop stats (TajimaD, PI) if available
                tajima_p <- mod_path(project_data()$name, MOD_SK,
                                     "tables", "pop_stats", "tajima_d_by_pop.tsv")
                pi_p     <- mod_path(project_data()$name, MOD_SK,
                                     "tables", "pop_stats", "pi_diversity_by_pop.tsv")
                if (file_ok(tajima_p)) {
                    taj <- data.table::fread(tajima_p, colClasses = c("site" = "character"))
                    dt  <- taj[dt, on = "site"]
                }
                if (file_ok(pi_p)) {
                    pi_dt <- data.table::fread(pi_p, colClasses = c("site" = "character"))
                    dt    <- pi_dt[dt, on = "site"]
                }
                dt
            })
        })

        output$site_table <- DT::renderDataTable({
            dt <- site_data()
            safe_datatable(
                as.data.frame(dt),
                extensions = "Buttons",
                options    = list(
                    dom       = "Bfrtip",
                    buttons   = list("csv"),
                    scrollX   = TRUE,
                    pageLength = 20
                ),
                rownames = FALSE
            )
        })

        output$dl_site <- shiny::downloadHandler(
            filename = function() {
                paste0("genetic_offset_site_", project_data()$name,
                       "_", selected_suffix(), ".csv")
            },
            content = function(file) {
                dt <- site_data()
                utils::write.csv(as.data.frame(dt), file, row.names = FALSE)
            }
        )
    })
}
