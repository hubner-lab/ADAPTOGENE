#' Overlapping Regions tab UI
#'
#' Miami Manhattan overlay, overlap pairs table, and region detail panel.
#'
#' @param id module namespace id
#' @noRd
mod_overlapping_ui <- function(id) {
    ns <- shiny::NS(id)
    bslib::page_sidebar(
        sidebar = bslib::sidebar(
            title = "Controls",
            shiny::uiOutput(ns("region_selector")),
            shiny::actionButton(ns("clear_region"), "Clear selection",
                                class = "btn-sm btn-outline-secondary w-100 mt-1")
        ),

        # Miami Manhattan (always visible)
        mod_manhattan_overlay_ui(ns("miami"), height = "550px"),

        # Overlap pairs table
        bslib::card(
            full_screen = TRUE,
            bslib::card_header(
                class = "d-flex justify-content-between align-items-center",
                htmltools::span(
                    bsicons::bs_icon("table"),
                    " Overlap Summary"
                ),
                shiny::downloadButton(ns("dl_overlap"), "CSV",
                                      class = "btn-sm btn-outline-secondary")
            ),
            bslib::card_body(
                DT::DTOutput(ns("overlap_table"))
            )
        ),

        # Region detail panel (conditional)
        shiny::uiOutput(ns("region_detail_ui"))
    )
}

#' Overlapping Regions tab server
#'
#' @param id module namespace id
#' @param project_data reactive project data bundle
#' @noRd
mod_overlapping_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        module <- MOD_OVERLAP

        # ── Data loading ───────────────────────────────────────────────────────
        regions_data <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("regions_combined_", pd$name, "_", module), function() {
                load_regions_combined(pd$name, module)
            })
        })

        genes_data <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("genes_", pd$name, "_", module), function() {
                load_genes(pd$name, module)
            })
        })

        overlap_data <- shiny::reactive({
            pd <- project_data()
            load_cached(paste0("overlap_summary_", pd$name), function() {
                load_overlap_summary(pd$name)
            })
        })

        # ── Region selection state ─────────────────────────────────────────────
        selected_region <- shiny::reactiveVal(NULL)

        shiny::observeEvent(input$clear_region, {
            selected_region(NULL)
            shiny::updateSelectInput(session, "region_id", selected = "")
        })

        shiny::observeEvent(input$region_id, {
            rid <- input$region_id
            selected_region(if (nzchar(rid %||% "")) rid else NULL)
        }, ignoreInit = TRUE)

        # ── Miami Manhattan ────────────────────────────────────────────────────
        miami_click <- mod_manhattan_overlay_server("miami",
            project_data         = project_data,
            module               = module,
            is_miami             = TRUE,
            combined             = TRUE,
            title_label          = shiny::reactive("Miami Plot (GEA \u2191 | GWAS \u2193)"),
            regions              = regions_data,
            show_regions_control = TRUE
        )

        shiny::observeEvent(miami_click(), {
            rid <- miami_click()
            if (!is.null(rid) && nzchar(rid)) {
                selected_region(rid)
                shiny::updateSelectInput(session, "region_id", selected = rid)
            }
        }, ignoreNULL = TRUE)

        # ── Region dropdown ────────────────────────────────────────────────────
        output$region_selector <- shiny::renderUI({
            rg <- regions_data()
            if (nrow(rg) == 0) {
                return(htmltools::p("No overlapping regions found.", class = "text-muted small"))
            }
            genes   <- genes_data()
            labels  <- build_region_labels(rg, genes)
            current <- selected_region()
            shiny::selectInput(ns("region_id"), "Select region",
                choices  = c(setNames("", ""), labels),
                selected = current %||% ""
            )
        })

        # ── Overlap table ──────────────────────────────────────────────────────
        output$overlap_table <- DT::renderDataTable({
            ov <- overlap_data()
            # Filter out SUMMARY rows for main table display
            if (nrow(ov) > 0 && "region_id" %in% names(ov)) {
                ov <- ov[ov$region_id != "SUMMARY", ]
            }
            safe_datatable(
                as.data.frame(ov),
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

        output$dl_overlap <- shiny::downloadHandler(
            filename = function() paste0("overlap_summary_", project_data()$name, ".csv"),
            content  = function(file) {
                utils::write.csv(as.data.frame(overlap_data()), file, row.names = FALSE)
            }
        )

        # ── Haplotype tag ──────────────────────────────────────────────────────
        hap_tag <- shiny::reactive({
            pd   <- project_data()
            tags <- find_haplotype_tags(pd$name)
            matched <- Filter(function(tag) {
                parts <- strsplit(tag, "_")[[1]]
                if (length(parts) < 2) return(FALSE)
                paste(parts[-1], collapse = "_") == "overlapping"
            }, tags)
            if (length(matched) > 0) matched[[1]] else NULL
        })

        # ── Region detail panel ────────────────────────────────────────────────
        output$region_detail_ui <- shiny::renderUI({
            rid <- selected_region()
            if (is.null(rid)) return(NULL)
            mod_region_detail_ui(ns("region_detail"))
        })

        region_trait <- shiny::reactive({
            rid <- selected_region()
            rg  <- regions_data()
            if (is.null(rid) || nrow(rg) == 0) return(NULL)
            row <- rg[rg$region_id == rid, ]
            if (nrow(row) == 0 || !"trait" %in% names(row)) return(NULL)
            row$trait[1]
        })

        mod_region_detail_server("region_detail",
            project_data = project_data,
            region_id    = selected_region,
            module       = module,
            trait        = region_trait,
            genes_data   = genes_data,
            hap_tag      = hap_tag
        )
    })
}
