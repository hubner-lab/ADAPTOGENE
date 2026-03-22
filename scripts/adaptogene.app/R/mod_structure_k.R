#' Structure K tab UI
#'
#' Piemaps (bio variable + metric + zoom), climate plots, population stats,
#' LD decay, and summary tables.
#'
#' @param id module namespace id
#' @noRd
mod_structure_k_ui <- function(id) {
    ns <- shiny::NS(id)
    bslib::page_sidebar(
        sidebar = bslib::sidebar(
            title = "Controls",
            shiny::uiOutput(ns("bio_selector")),
            shiny::uiOutput(ns("metric_selector")),
            shiny::uiOutput(ns("zoom_selector"))
        ),

        # Piemap (main)
        mod_piemap_viewer_ui(ns("piemap")),

        # Tajima D + Pi piemaps (optional, in accordion)
        bslib::accordion(
            id       = ns("sk_sections"),
            open     = FALSE,
            multiple = TRUE,

            bslib::accordion_panel(
                "Climate Plots",
                value = "climate",
                icon  = bsicons::bs_icon("thermometer-half"),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("climate_heatmap")),
                    mod_image_card_ui(ns("climate_density"))
                )
            ),

            bslib::accordion_panel(
                "Population Statistics",
                value = "pop_stats",
                icon  = bsicons::bs_icon("graph-up"),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("mantel")),
                    mod_image_card_ui(ns("amova"))
                ),
                bslib::layout_column_wrap(
                    width = 1 / 2,
                    mod_image_card_ui(ns("ld_decay")),
                    mod_image_card_ui(ns("ld_decay_chr"))
                )
            ),

            bslib::accordion_panel(
                "Tables",
                value  = "tables",
                icon   = bsicons::bs_icon("table"),
                shiny::uiOutput(ns("tables_content"))
            )
        )
    )
}

#' Structure K tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_structure_k_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        bio_values <- shiny::reactive(find_bio_values(project_data()$name))

        # ── Sidebar selectors ──────────────────────────────────────────────────
        output$bio_selector <- shiny::renderUI({
            bios <- bio_values()
            if (length(bios) == 0) return(shiny::p("No piemaps found.", class = "text-muted small"))
            choices <- setNames(bios, paste0("bio", bios))
            shiny::selectInput(ns("bio"), "Bio variable", choices = choices, selected = bios[1])
        })

        output$metric_selector <- shiny::renderUI({
            shiny::selectInput(ns("metric"), "Scale by metric",
                choices = c(
                    "None"           = "none",
                    "Tajima's D"     = "tajima_d",
                    "Pi Diversity"   = "pi_diversity"
                ),
                selected = "none"
            )
        })

        output$zoom_selector <- shiny::renderUI({
            # Zoom directories are subdirectories of piemap/zoom/
            pd       <- project_data()
            zoom_dir <- mod_path(pd$name, MOD_SK, "plots", "piemap", "zoom")
            if (!dir.exists(zoom_dir)) return(NULL)
            tags <- list.dirs(zoom_dir, full.names = FALSE, recursive = FALSE)
            tags <- tags[nzchar(tags)]
            if (length(tags) == 0) return(NULL)
            shiny::selectInput(ns("zoom"), "Zoom region",
                choices  = c("Global" = "none", setNames(tags, tags)),
                selected = "none"
            )
        })

        selected_bio <- shiny::reactive({
            b <- input$bio
            if (is.null(b)) {
                bios <- bio_values()
                if (length(bios) == 0) return(1L)
                bios[1]
            } else as.integer(b)
        })

        selected_metric <- shiny::reactive(input$metric %||% "none")
        selected_zoom   <- shiny::reactive(input$zoom   %||% "none")

        # ── Piemap viewer ──────────────────────────────────────────────────────
        mod_piemap_viewer_server("piemap",
            project_data = project_data,
            bio          = selected_bio,
            metric       = selected_metric,
            zoom         = selected_zoom
        )

        # ── Climate images ─────────────────────────────────────────────────────
        shiny::observe({
            pd <- project_data()
            mod_image_card_server("climate_heatmap",
                path    = shiny::reactive(climate_heatmap_path(pd$name)),
                title   = shiny::reactive("Climate Correlation Heatmap"),
                dl_name = shiny::reactive("climate_heatmap")
            )
            mod_image_card_server("climate_density",
                path    = shiny::reactive(climate_density_path(pd$name)),
                title   = shiny::reactive("Climate Density (Present)"),
                dl_name = shiny::reactive("climate_density_present")
            )
        })

        # ── Population stats images ────────────────────────────────────────────
        shiny::observe({
            pd       <- project_data()
            stats_dir <- mod_path(pd$name, MOD_SK, "plots", "pop_stats")

            mantel_p <- file.path(stats_dir, "mantel_test.png")
            amova_p  <- file.path(stats_dir, "amova.png")

            mod_image_card_server("mantel",
                path    = shiny::reactive(if (file_ok(mantel_p)) mantel_p else NULL),
                title   = shiny::reactive("Mantel Test"),
                dl_name = shiny::reactive("mantel_test")
            )
            mod_image_card_server("amova",
                path    = shiny::reactive(if (file_ok(amova_p)) amova_p else NULL),
                title   = shiny::reactive("AMOVA"),
                dl_name = shiny::reactive("amova")
            )
            mod_image_card_server("ld_decay",
                path    = shiny::reactive(ld_decay_path(pd$name, per_chr = FALSE)),
                title   = shiny::reactive("LD Decay (Genome-wide)"),
                dl_name = shiny::reactive("ld_decay")
            )
            mod_image_card_server("ld_decay_chr",
                path    = shiny::reactive(ld_decay_path(pd$name, per_chr = TRUE)),
                title   = shiny::reactive("LD Decay (Per Chromosome)"),
                dl_name = shiny::reactive("ld_decay_per_chr")
            )
        })

        # ── Tables ─────────────────────────────────────────────────────────────
        output$tables_content <- shiny::renderUI({
            pd        <- project_data()
            tables_dir <- mod_path(pd$name, MOD_SK, "tables", "pop_stats")
            if (!dir.exists(tables_dir)) {
                return(htmltools::p("No population statistics tables available.",
                                    class = "text-muted small"))
            }
            files <- list.files(tables_dir, pattern = "\\.tsv$", full.names = TRUE)
            if (length(files) == 0) {
                return(htmltools::p("No tables available.", class = "text-muted small"))
            }

            # Build accordion panels for each table file
            panels <- lapply(files, function(f) {
                tbl_name <- tools::file_path_sans_ext(basename(f))
                tbl_id   <- gsub("[^A-Za-z0-9]", "_", tbl_name)
                bslib::accordion_panel(
                    gsub("_", " ", tbl_name),
                    DT::DTOutput(ns(paste0("tbl_", tbl_id)))
                )
            })

            # Render each table server-side
            lapply(files, function(f) {
                tbl_name <- tools::file_path_sans_ext(basename(f))
                tbl_id   <- gsub("[^A-Za-z0-9]", "_", tbl_name)
                output[[paste0("tbl_", tbl_id)]] <- DT::renderDataTable({
                    dt <- data.table::fread(f)
                    safe_datatable(
                        as.data.frame(dt),
                        extensions = "Buttons",
                        options    = list(
                            dom       = "Bfrtip",
                            buttons   = list("csv"),
                            scrollX   = TRUE,
                            pageLength = 15
                        ),
                        rownames = FALSE
                    )
                })
            })

            do.call(bslib::accordion, c(panels, list(open = FALSE, multiple = TRUE)))
        })
    })
}
