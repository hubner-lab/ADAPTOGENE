#' Maladaptation tab UI
#'
#' Gradient Forest results: importance plots, genetic offset piemaps,
#' zoom maps, and site-level offset table.
#'
#' @param id module namespace id
#' @noRd
mod_maladaptation_ui <- function(id) {
    ns <- shiny::NS(id)
    bslib::page_sidebar(
        sidebar = bslib::sidebar(
            title = "Controls",
            shiny::uiOutput(ns("suffix_selector"))
        ),

        # Importance plots
        bslib::layout_column_wrap(
            width = 1 / 2,
            mod_image_card_ui(ns("overall_importance")),
            mod_image_card_ui(ns("cumulative_importance"))
        ),

        shiny::hr(),

        # Genetic Offset piemaps
        shiny::h5("Genetic Offset Piemaps", class = "text-muted mt-2"),
        bslib::layout_column_wrap(
            width = 1 / 3,
            mod_image_card_ui(ns("piemap_base")),
            mod_image_card_ui(ns("piemap_tajima")),
            mod_image_card_ui(ns("piemap_pi"))
        ),

        # Zoom maps + site table in accordion
        bslib::accordion(
            id       = ns("malad_sections"),
            open     = FALSE,
            multiple = TRUE,

            bslib::accordion_panel(
                "Zoom Maps",
                value = "zoom",
                icon  = bsicons::bs_icon("zoom-in"),
                shiny::uiOutput(ns("zoom_content"))
            ),

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
            shiny::selectInput(ns("suffix"), "GF run", choices = suf, selected = suf[1])
        })

        selected_suffix <- shiny::reactive({
            s <- input$suffix
            if (is.null(s)) {
                suf <- gf_suffixes()
                if (length(suf) == 0) return(NULL)
                suf[1]
            } else s
        })

        # ── Importance + piemap images ─────────────────────────────────────────
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
            mod_image_card_server("piemap_base",
                path    = shiny::reactive(gf_offset_piemap_path(pd$name, suf, "base")),
                title   = shiny::reactive("Genetic Offset (Base)"),
                dl_name = shiny::reactive(paste0("offset_piemap_base_", suf))
            )
            mod_image_card_server("piemap_tajima",
                path    = shiny::reactive(gf_offset_piemap_path(pd$name, suf, "tajima_d")),
                title   = shiny::reactive("Offset \u00d7 Tajima's D"),
                dl_name = shiny::reactive(paste0("offset_piemap_tajima_", suf))
            )
            mod_image_card_server("piemap_pi",
                path    = shiny::reactive(gf_offset_piemap_path(pd$name, suf, "pi_diversity")),
                title   = shiny::reactive("Offset \u00d7 Pi Diversity"),
                dl_name = shiny::reactive(paste0("offset_piemap_pi_", suf))
            )
        })

        # ── Zoom maps ──────────────────────────────────────────────────────────
        output$zoom_content <- shiny::renderUI({
            suf  <- shiny::req(selected_suffix())
            pd   <- project_data()
            zooms <- find_gf_zooms(pd$name, suf)

            if (length(zooms) == 0) {
                return(htmltools::p("No zoom maps available.", class = "text-muted small"))
            }

            zoom_cards <- lapply(seq_along(zooms), function(i) {
                z_id <- paste0("zoom_", i)
                path <- mod_path(pd$name, MOD_MALAD, "plots", suf, "zoom",
                                  paste0(zooms[i], ".png"))
                mod_image_card_server(z_id,
                    path    = shiny::reactive(path),
                    title   = shiny::reactive(paste0("Zoom: ", zooms[i])),
                    dl_name = shiny::reactive(paste0("offset_zoom_", zooms[i], "_", suf))
                )
                mod_image_card_ui(ns(z_id))
            })

            do.call(bslib::layout_column_wrap, c(list(width = 1 / 2), zoom_cards))
        })

        # ── Site table ─────────────────────────────────────────────────────────
        site_data <- shiny::reactive({
            suf <- selected_suffix()
            if (is.null(suf)) return(data.table::data.table())
            load_cached(paste0("gf_site_", project_data()$name, "_", suf), function() {
                p <- gf_site_table_path(project_data()$name, suf)
                if (!file_ok(p)) return(data.table::data.table())
                data.table::fread(p)
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
