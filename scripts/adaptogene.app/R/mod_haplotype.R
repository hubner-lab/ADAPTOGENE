#' Haplotype Analysis tab UI
#'
#' Scan browser: tag selector, region selector, status summary, clustree plots,
#' per-region visualization (custom source only), and data tables.
#'
#' @param id module namespace id
#' @noRd
mod_haplotype_ui <- function(id) {
    ns <- shiny::NS(id)
    htmltools::tagList(
        # Tag selector inline above status + clustree
        htmltools::div(
            class = "control-bar",
            shiny::uiOutput(ns("tag_selector"))
        ),

        # Status summary value boxes
        shiny::uiOutput(ns("status_boxes")),

        shiny::br(),

        # Clustree region selector + plots
        htmltools::div(
            class = "control-bar",
            shiny::uiOutput(ns("clustree_region_selector"))
        ),
        bslib::layout_column_wrap(
            width = 1 / 2,
            mod_image_card_ui(ns("clustree_mg")),
            mod_image_card_ui(ns("clustree_hap"))
        ),

        # Per-region viz + tables in accordions
        bslib::accordion(
            id       = ns("hap_sections"),
            open     = FALSE,
            multiple = TRUE,

            bslib::accordion_panel(
                "Tables",
                value = "tables",
                icon  = bsicons::bs_icon("table"),
                shiny::uiOutput(ns("tables_content"))
            )
        )
    )
}

#' Haplotype Analysis tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_haplotype_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ── Tag discovery ──────────────────────────────────────────────────────
        hap_tags <- shiny::reactive(find_haplotype_tags(project_data()$name))

        output$tag_selector <- shiny::renderUI({
            tags <- hap_tags()
            if (length(tags) == 0) {
                return(shiny::p("No haplotype analyses found.", class = "text-muted small"))
            }
            labels <- setNames(tags, vapply(tags, format_hap_tag, character(1)))
            shiny::selectInput(ns("tag"), "Analysis tag", choices = labels, selected = tags[1])
        })

        selected_tag <- shiny::reactive(input$tag %||% hap_tags()[1])

        # ── Clustree region discovery ──────────────────────────────────────────
        clustree_regions <- shiny::reactive({
            tag <- selected_tag()
            if (is.null(tag)) return(character(0))
            find_haplotype_scan_regions(project_data()$name, tag)
        })

        output$clustree_region_selector <- shiny::renderUI({
            rids <- clustree_regions()
            if (length(rids) == 0) return(NULL)
            if (length(rids) == 1) return(NULL)  # no selector needed for single region
            choices <- setNames(rids, vapply(rids, format_region_id, character(1)))
            shiny::selectInput(ns("clustree_region"), "Region (clustree)",
                               choices = choices, selected = rids[1])
        })

        selected_clustree_region <- shiny::reactive({
            cr <- clustree_regions()
            input$clustree_region %||% (if (length(cr) > 0) cr[1] else NULL)
        })

        # ── Status summary ─────────────────────────────────────────────────────
        scan_status <- shiny::reactive({
            tag <- selected_tag()
            if (is.null(tag)) return(data.table::data.table())
            pd <- project_data()
            load_cached(paste0("hap_status_", pd$name, "_", tag), function() {
                load_haplotype_status(pd$name, tag)
            })
        })

        output$status_boxes <- shiny::renderUI({
            st  <- scan_status()
            tag <- selected_tag()
            if (nrow(st) == 0 || is.null(tag)) return(NULL)

            summary_row <- st[st$region_id == "SUMMARY", ]
            if (nrow(summary_row) == 0) return(NULL)

            # Summary row format: chr="total", start=n_total, end=n_success,
            # snp_count=n_skipped, status=n_failed
            n_total   <- as.integer(summary_row$start[1])
            n_success <- as.integer(summary_row$end[1])
            n_skip    <- as.integer(summary_row$snp_count[1])
            n_fail    <- as.integer(summary_row$status[1])

            status_theme <- if (n_fail > 0 && n_success == 0) "danger"
                            else if (n_fail > 0 || n_skip > 0) "warning"
                            else "success"

            bslib::layout_column_wrap(
                width = "180px",
                fill  = FALSE,
                bslib::value_box(
                    title    = "Total regions",
                    value    = n_total,
                    theme    = "info",
                    showcase = bsicons::bs_icon("geo")
                ),
                bslib::value_box(
                    title    = "Succeeded",
                    value    = n_success,
                    theme    = "success",
                    showcase = bsicons::bs_icon("check-circle")
                ),
                bslib::value_box(
                    title    = "Skipped",
                    value    = n_skip,
                    theme    = "warning",
                    showcase = bsicons::bs_icon("skip-forward")
                ),
                bslib::value_box(
                    title    = "Failed",
                    value    = n_fail,
                    theme    = status_theme,
                    showcase = bsicons::bs_icon("x-circle")
                )
            )
        })

        # ── Clustree plots (per selected clustree region) ──────────────────────
        shiny::observe({
            tag  <- shiny::req(selected_tag())
            crid <- shiny::req(selected_clustree_region())
            pd   <- project_data()

            mod_image_card_server("clustree_mg",
                path    = shiny::reactive(hap_clustree_path(pd$name, tag, crid, "MG")),
                title   = shiny::reactive(paste0("Clustree: Marker Groups (",
                                                  format_region_id(crid), ")")),
                dl_name = shiny::reactive(paste0("clustree_MG_", tag, "_", crid))
            )
            mod_image_card_server("clustree_hap",
                path    = shiny::reactive(hap_clustree_path(pd$name, tag, crid, "hap")),
                title   = shiny::reactive(paste0("Clustree: Haplotypes (",
                                                  format_region_id(crid), ")")),
                dl_name = shiny::reactive(paste0("clustree_hap_", tag, "_", crid))
            )
        })

        # ── Tables ─────────────────────────────────────────────────────────────
        output$tables_content <- shiny::renderUI({
            tag <- selected_tag()
            if (is.null(tag)) return(NULL)
            pd <- project_data()

            panels <- list(
                bslib::accordion_panel(
                    "Selected Regions",
                    DT::DTOutput(ns("tbl_selected"))
                ),
                bslib::accordion_panel(
                    "Scan Status Detail",
                    DT::DTOutput(ns("tbl_status"))
                )
            )

            # Add assignment/frequency tables when haplotype viz has been run
            if (file_ok(hap_assignments_path(pd$name, tag))) {
                panels <- c(panels, list(
                    bslib::accordion_panel(
                        "Haplotype Assignments",
                        DT::DTOutput(ns("tbl_assignments"))
                    ),
                    bslib::accordion_panel(
                        "Haplotype Frequencies",
                        DT::DTOutput(ns("tbl_frequencies"))
                    )
                ))
            }

            do.call(bslib::accordion, c(panels, list(open = FALSE, multiple = TRUE)))
        })

        output$tbl_selected <- DT::renderDataTable({
            tag <- shiny::req(selected_tag())
            pd  <- project_data()
            p   <- hap_selected_regions_path(pd$name, tag)
            dt  <- if (file_ok(p)) data.table::fread(p) else data.table::data.table()
            safe_datatable(
                as.data.frame(dt),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                   scrollX = TRUE, pageLength = 15),
                rownames   = FALSE
            )
        })

        output$tbl_status <- DT::renderDataTable({
            st <- scan_status()
            if (nrow(st) > 0 && "region_id" %in% names(st)) {
                st <- st[st$region_id != "SUMMARY", ]
            }
            safe_datatable(
                as.data.frame(st),
                options  = list(scrollX = TRUE, pageLength = 20),
                rownames = FALSE
            )
        })

        output$tbl_assignments <- DT::renderDataTable({
            tag <- shiny::req(selected_tag())
            pd  <- project_data()
            p   <- hap_assignments_path(pd$name, tag)
            dt  <- if (file_ok(p)) data.table::fread(p) else data.table::data.table()
            safe_datatable(
                as.data.frame(dt),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                   scrollX = TRUE, pageLength = 20),
                rownames   = FALSE
            )
        })

        output$tbl_frequencies <- DT::renderDataTable({
            tag <- shiny::req(selected_tag())
            pd  <- project_data()
            p   <- hap_frequencies_path(pd$name, tag)
            dt  <- if (file_ok(p)) data.table::fread(p) else data.table::data.table()
            safe_datatable(
                as.data.frame(dt),
                extensions = "Buttons",
                options    = list(dom = "Bfrtip", buttons = list("csv"),
                                   scrollX = TRUE, pageLength = 20),
                rownames   = FALSE
            )
        })
    })
}
